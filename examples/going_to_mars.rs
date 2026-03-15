use poliastrs::bodies::{EARTH, SUN, Body};
use poliastrs::twobody::orbit::Orbit;
use poliastrs::core::elements::ClassicalElements;
use poliastrs::core::angles::{m_to_e, e_to_nu};
use poliastrs::iod::lambert::izzo;
use poliastrs::frames::Plane;
use nalgebra::Vector3;
use chrono::{TimeZone, Utc, DateTime};
use plotters::prelude::*;

// Constants
const AU_KM: f64 = 149_597_870.7;

fn mean_to_true(mean: f64, ecc: f64) -> f64 {
    let e_anom = m_to_e(mean, ecc);
    e_to_nu(e_anom, ecc)
}

fn get_j2000_orbit(
    attractor: Body,
    a_au: f64,
    ecc: f64,
    inc_deg: f64,
    raan_deg: f64,
    argp_deg: f64,
    mean_anom_deg: f64,
) -> Orbit {
    let a_km = a_au * AU_KM;
    let mean_rad = mean_anom_deg.to_radians();
    let nu_rad = mean_to_true(mean_rad, ecc);
    
    Orbit::from_classical_at(
        attractor,
        ClassicalElements {
            p_km: a_km * (1.0 - ecc * ecc),
            ecc,
            inc_rad: inc_deg.to_radians(),
            raan_rad: raan_deg.to_radians(),
            argp_rad: argp_deg.to_radians(),
            nu_rad,
        },
        0.0, // J2000
        Plane::EarthEcliptic,
    ).unwrap()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Going to Mars with Rust and poliastrs");

    // 1. Initial Data
    let date_launch = Utc.with_ymd_and_hms(2011, 11, 26, 15, 2, 0).unwrap();
    let date_arrival = Utc.with_ymd_and_hms(2012, 8, 6, 5, 17, 0).unwrap();
    
    let j2000 = Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap();
    let t_launch = (date_launch - j2000).num_milliseconds() as f64 / 1000.0;
    let t_arrival = (date_arrival - j2000).num_milliseconds() as f64 / 1000.0;
    let tof = t_arrival - t_launch;
    
    println!("Time of Flight: {:.2} days", tof / 86400.0);

    // 2. Ephemerides (Approximate J2000)
    // Earth J2000: a=1.00000011 AU, e=0.01671022, i=0.00005 deg, L=100.46435 deg, w_bar=102.94719 deg, raan=-11.26064 deg
    // M = L - w_bar = -2.48284 deg
    // w = w_bar - raan = 114.20783 deg
    let earth_j2000 = get_j2000_orbit(
        SUN, 
        1.00000011, 
        0.01671022, 
        0.00005, 
        -11.26064, 
        114.20783, 
        -2.48284
    );
    
    // Mars J2000: a=1.523679 AU, e=0.093405, i=1.8497 deg, raan=49.558 deg, w=286.502 deg, M=19.412 deg
    // Using values from previous thought
    let mars_j2000 = get_j2000_orbit(
        SUN,
        1.523679,
        0.093405,
        1.8497,
        49.558,
        286.502,
        19.412
    );
    
    // Propagate to dates
    let earth_at_launch = earth_j2000.propagate_seconds(t_launch).unwrap();
    let mars_at_arrival = mars_j2000.propagate_seconds(t_arrival).unwrap();
    
    let r0 = Vector3::new(earth_at_launch.state.r_km.x, earth_at_launch.state.r_km.y, earth_at_launch.state.r_km.z);
    let rf = Vector3::new(mars_at_arrival.state.r_km.x, mars_at_arrival.state.r_km.y, mars_at_arrival.state.r_km.z);
    
    // 3. Solve Lambert
    println!("Solving Lambert Transfer...");
    let (v0, _vf) = izzo(SUN.mu_km3_s2, r0, rf, tof, 0, true, false, 35, 1e-8).map_err(|e| format!("{:?}", e))?;
    
    let v_earth = Vector3::new(earth_at_launch.state.v_km_s.x, earth_at_launch.state.v_km_s.y, earth_at_launch.state.v_km_s.z);
    let dv_launch = (v0 - v_earth).norm();
    
    println!("Launch C3: {:.4} km^2/s^2", dv_launch.powi(2));
    println!("Launch Delta-V (Hyperbolic Excess): {:.4} km/s", dv_launch);

    // 4. Plotting in 3D
    let root = BitMapBackend::new("going_to_mars.png", (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("MSL Mission: Earth to Mars", ("sans-serif", 30))
        .margin(20)
        .build_cartesian_3d(
            -2.0 * AU_KM..2.0 * AU_KM, 
            -2.0 * AU_KM..2.0 * AU_KM, 
            -0.5 * AU_KM..0.5 * AU_KM
        )?;

    chart.configure_axes().draw()?;

    // Sun
    chart.draw_series(PointSeries::of_element(
        vec![(0.0, 0.0, 0.0)],
        5,
        &YELLOW,
        &|c, s, st| {
            return EmptyElement::at(c) + Circle::new((0,0), s, st.filled());
        },
    ))?;
    
    // Earth Orbit (at launch)
    let earth_points: Vec<(f64, f64, f64)> = (0..360).map(|i| {
        let nu = (i as f64).to_radians();
        let state = earth_j2000.propagate_to_anomaly(nu).unwrap();
        (state.state.r_km.x, state.state.r_km.y, state.state.r_km.z)
    }).collect();
    chart.draw_series(LineSeries::new(earth_points, &BLUE))?
        .label("Earth Orbit")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
        
    // Mars Orbit (at arrival)
    let mars_points: Vec<(f64, f64, f64)> = (0..360).map(|i| {
        let nu = (i as f64).to_radians();
        let state = mars_j2000.propagate_to_anomaly(nu).unwrap();
        (state.state.r_km.x, state.state.r_km.y, state.state.r_km.z)
    }).collect();
    chart.draw_series(LineSeries::new(mars_points, &RED))?
        .label("Mars Orbit")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Transfer Trajectory
    let transfer_orbit = Orbit::from_vectors_at(
        SUN,
        [r0.x, r0.y, r0.z],
        [v0.x, v0.y, v0.z],
        t_launch,
        Plane::EarthEcliptic
    );
    
    let transfer_points: Vec<(f64, f64, f64)> = (0..100).map(|i| {
        let t = i as f64 * tof / 100.0;
        let state = transfer_orbit.propagate_seconds(t_launch + t).unwrap();
        (state.state.r_km.x, state.state.r_km.y, state.state.r_km.z)
    }).collect();
    chart.draw_series(LineSeries::new(transfer_points, &GREEN))?
        .label("Transfer")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    // Markers for Launch and Arrival
    chart.draw_series(PointSeries::of_element(
        vec![(r0.x, r0.y, r0.z)],
        5,
        &BLUE,
        &|c, s, st| {
            return EmptyElement::at(c) + Circle::new((0,0), s, st.filled());
        },
    ))?; // Launch point
    
    chart.draw_series(PointSeries::of_element(
        vec![(rf.x, rf.y, rf.z)],
        5,
        &RED,
        &|c, s, st| {
            return EmptyElement::at(c) + Circle::new((0,0), s, st.filled());
        },
    ))?; // Arrival point

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    println!("Generated going_to_mars.png");
    Ok(())
}
