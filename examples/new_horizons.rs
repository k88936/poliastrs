use poliastrs::bodies::{EARTH, SUN, Body};
use poliastrs::twobody::orbit::Orbit;
use poliastrs::core::elements::ClassicalElements;
use poliastrs::core::angles::{m_to_e, e_to_nu};
use poliastrs::iod::lambert::izzo;
use poliastrs::frames::Plane;
use nalgebra::Vector3;
use chrono::{TimeZone, Utc, DateTime};
use plotters::prelude::*;

fn mean_to_true(mean: f64, ecc: f64) -> f64 {
    let e_anom = m_to_e(mean, ecc);
    e_to_nu(e_anom, ecc)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Exploring the New Horizons launch with Rust and poliastrs");

    // 1. Parking Orbit
    // r_p = 165 km altitude, r_a = 215 km altitude
    let r_p_alt: f64 = 165.0;
    let r_a_alt: f64 = 215.0;
    let r_earth = EARTH.equatorial_radius_km; // km
    let r_p = r_earth + r_p_alt;
    let r_a = r_earth + r_a_alt;
    
    let a_parking = (r_p + r_a) / 2.0;
    let ecc_parking = 1.0 - r_p / a_parking;
    
    println!("Parking Orbit: a = {:.2} km, e = {:.6}", a_parking, ecc_parking);
    
    let parking = Orbit::from_classical(
        EARTH,
        ClassicalElements {
            p_km: a_parking * (1.0 - ecc_parking.powi(2)),
            ecc: ecc_parking,
            inc_rad: 0.0,
            raan_rad: 0.0,
            argp_rad: 0.0,
            nu_rad: 0.0,
        }
    );
    
    println!("Parking Orbit Velocity at Perigee: {:.4} km/s", parking.state.v_km_s.norm());
    
    // 2. Hyperbolic Exit (Option A: C3 from design)
    let c3_design = 157.6561; // km^2/s^2
    let a_exit = -(EARTH.mu_km3_s2 / c3_design);
    let ecc_exit = 1.0 - r_p / a_exit;
    
    println!("Hyperbolic Exit (Design): a = {:.2} km, e = {:.6}", a_exit, ecc_exit);
    
    let exit_orbit = Orbit::from_classical(
        EARTH,
        ClassicalElements {
            p_km: a_exit * (1.0 - ecc_exit.powi(2)),
            ecc: ecc_exit,
            inc_rad: 0.0,
            raan_rad: 0.0,
            argp_rad: 0.0,
            nu_rad: 0.0,
        }
    );
    
    let v_exit_perigee = exit_orbit.state.v_km_s.norm();
    println!("Hyperbolic Exit Velocity at Perigee: {:.4} km/s", v_exit_perigee);
    
    let v_estimated = 16.2;
    let error = (v_exit_perigee - v_estimated).abs() / v_estimated * 100.0;
    println!("Relative Error vs Estimate (16.2 km/s): {:.2}%", error);
    
    // 3. Option B: Lambert from Earth to Jupiter
    let date_launch = Utc.with_ymd_and_hms(2006, 1, 19, 19, 0, 0).unwrap();
    let date_flyby = Utc.with_ymd_and_hms(2007, 2, 28, 5, 43, 40).unwrap();
    
    let j2000 = Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap();
    let t_launch = (date_launch - j2000).num_milliseconds() as f64 / 1000.0;
    let t_flyby = (date_flyby - j2000).num_milliseconds() as f64 / 1000.0;
    let tof = t_flyby - t_launch;
    
    // Ephemerides (Approximate J2000)
    let au_km = 149_597_870.7;
    
    // Earth
    let earth_a = 1.00000011 * au_km;
    let earth_e = 0.01671022;
    let earth_i = 0.00005_f64.to_radians();
    let earth_raan = (-11.26064_f64).to_radians();
    let earth_w = 114.20783_f64.to_radians();
    let earth_m = (-2.48284_f64).to_radians(); // at J2000
    
    let earth_j2000 = Orbit::from_classical_at(
        SUN,
        ClassicalElements {
            p_km: earth_a * (1.0 - earth_e * earth_e),
            ecc: earth_e,
            inc_rad: earth_i,
            raan_rad: earth_raan,
            argp_rad: earth_w,
            nu_rad: mean_to_true(earth_m, earth_e),
        },
        0.0,
        Plane::EarthEcliptic
    ).unwrap();
    
    // Jupiter
    let jup_a = 5.20336301 * au_km;
    let jup_e = 0.04839266;
    let jup_i = 1.30530_f64.to_radians();
    let jup_raan = 100.55615_f64.to_radians();
    let jup_w = (-85.8023_f64).to_radians();
    let jup_m = 19.65053_f64.to_radians();
    
    let jup_j2000 = Orbit::from_classical_at(
        SUN,
        ClassicalElements {
            p_km: jup_a * (1.0 - jup_e * jup_e),
            ecc: jup_e,
            inc_rad: jup_i,
            raan_rad: jup_raan,
            argp_rad: jup_w,
            nu_rad: mean_to_true(jup_m, jup_e),
        },
        0.0,
        Plane::EarthEcliptic
    ).unwrap();
    
    let earth_state = earth_j2000.propagate_seconds(t_launch).unwrap();
    let jup_state = jup_j2000.propagate_seconds(t_flyby).unwrap();
    
    let r0 = Vector3::new(earth_state.state.r_km.x, earth_state.state.r_km.y, earth_state.state.r_km.z);
    let rf = Vector3::new(jup_state.state.r_km.x, jup_state.state.r_km.y, jup_state.state.r_km.z);
    
    println!("Solving Lambert...");
    // prograde = true for typical interplanetary transfer
    let (v0, _vf) = izzo(SUN.mu_km3_s2, r0, rf, tof, 0, true, false, 35, 1e-8).map_err(|e| format!("{:?}", e))?;
    
    let v_earth = Vector3::new(earth_state.state.v_km_s.x, earth_state.state.v_km_s.y, earth_state.state.v_km_s.z);
    
    let v_inf_vec = v0 - v_earth;
    let c3_lambert = v_inf_vec.norm_squared();
    
    println!("C3 (Lambert): {:.4} km^2/s^2", c3_lambert);
    println!("Relative Error vs Design ({:.4}): {:.2}%", c3_design, (c3_lambert - c3_design).abs() / c3_design * 100.0);

    // Plotting
    let root = BitMapBackend::new("new_horizons.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("New Horizons Transfer", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-10.0 * au_km..10.0 * au_km, -10.0 * au_km..10.0 * au_km)?;

    chart.configure_mesh().draw()?;
    
    // Plot Sun
    chart.draw_series(PointSeries::of_element(
        vec![(0.0, 0.0)],
        5,
        &YELLOW,
        &|c, s, st| {
            return EmptyElement::at(c) + Circle::new((0,0), s, st.filled());
        },
    ))?;

    // Plot Earth Orbit
    let earth_points: Vec<(f64, f64)> = (0..360).map(|i| {
        let nu = (i as f64).to_radians();
        let state = earth_j2000.propagate_to_anomaly(nu);
        (state.unwrap().state.r_km.x, state.unwrap().state.r_km.y)
    }).collect();
    chart.draw_series(LineSeries::new(earth_points, &BLUE))?
        .label("Earth")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    // Plot Jupiter Orbit
    let jup_points: Vec<(f64, f64)> = (0..360).map(|i| {
        let nu = (i as f64).to_radians();
        let state = jup_j2000.propagate_to_anomaly(nu);
        (state.unwrap().state.r_km.x, state.unwrap().state.r_km.y)
    }).collect();
    chart.draw_series(LineSeries::new(jup_points, &RED))?
        .label("Jupiter")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Plot Transfer
    let transfer_orbit = Orbit::from_vectors_at(
        SUN, 
        [earth_state.state.r_km.x, earth_state.state.r_km.y, earth_state.state.r_km.z], 
        [v0.x, v0.y, v0.z], 
        t_launch,
        Plane::EarthEcliptic 
    );
    let transfer_points: Vec<(f64, f64)> = (0..100).map(|i| {
        let t = i as f64 * tof / 100.0;
        let state = transfer_orbit.propagate_seconds(t_launch + t).unwrap();
        (state.state.r_km.x, state.state.r_km.y)
    }).collect();
    chart.draw_series(LineSeries::new(transfer_points, &GREEN))?
        .label("New Horizons")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));
        
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    println!("Generated new_horizons.png");

    Ok(())
}
