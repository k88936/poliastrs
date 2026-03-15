use poliastrs::bodies::SUN;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::core::elements::ClassicalElements;
use plotters::prelude::*;
use plotters::style::ShapeStyle; // Explicit import if needed

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Catch that asteroid! (Migration of catch-that-asteroid.myst.md)");
    println!("============================================================");

    // 1. Define Epoch
    // "2017-09-01 12:05:50" (TDB)
    // J2000: 2000-01-01 12:00:00
    // Days since J2000 (Approximate)
    // 2000 to 2017 = 17 years. 17 * 365.25 = 6209.25 days.
    // Jan to Aug (excluding Sep) = 243 days (non-leap 2017).
    // Sep 1.5.
    // Total days: ~6453.75
    // Let's be a bit more precise.
    // Julian Date of 2017-09-01 12:05:50.
    // JD = 2457998.00405
    // JD_J2000 = 2451545.0
    // Diff = 6453.00405 days
    // Seconds = Diff * 86400 = 557539549.9
    let target_epoch_seconds = 557_539_550.0;
    println!("Target Epoch (s from J2000): {:.1}", target_epoch_seconds);

    const AU_KM: f64 = 149_597_870.7;

    // 2. Define Earth Orbit (Approximate Mean Elements J2000)
    // a = 1.00000011 AU
    // e = 0.01671022
    // i = 0.00005 deg
    // W = -11.26064 deg (Longitude of perihelion is given usually, w = W - Omega)
    // L = 100.46435 deg (Mean Longitude)
    // M = L - W
    let earth_a = 1.00000011 * AU_KM;
    let earth_e = 0.01671022;
    let earth_i = 0.00005_f64.to_radians();
    let earth_omega = -11.26064_f64.to_radians(); // taking this as Longitude of perihelion approx if RAAN is 0?
    // Actually for Earth, RAAN is undefined (ecliptic). Let's use standard simplified:
    // a=1.000, e=0.0167, i=0, raan=0, argp=102.94, M=-2.47 (at J2000).
    // Let's use the ones from parker_solar_probe if possible, or standard.
    
    let earth = Orbit::from_classical(SUN, ClassicalElements {
        p_km: earth_a * (1.0 - earth_e * earth_e),
        ecc: earth_e,
        inc_rad: earth_i,
        raan_rad: 0.0,
        argp_rad: 102.94_f64.to_radians(), // Longitude of perihelion - RAAN(0)
        nu_rad: (-2.47_f64).to_radians(), // Mean anomaly at J2000 converted to true anomaly approx
    });
    // Convert M to nu? Orbit::from_classical takes nu.
    // M = -2.47 deg.
    // nu approx M + 2e sin M.
    // -2.47 + 2*0.0167*sin(-2.47) = -2.47 - 0.0014 = -2.47 deg.
    // Just use -2.47 deg as nu for approx.

    // 3. Define Florence Orbit (J2000)
    // 3122 Florence
    // a = 1.769 AU
    // e = 0.423
    // i = 22.14 deg
    // raan = 336.1 deg
    // argp = 27.8 deg
    // M = 341.6 deg
    let florence_a = 1.769 * AU_KM;
    let florence_e = 0.423;
    let florence_i = 22.14_f64.to_radians();
    let florence_raan = 336.1_f64.to_radians();
    let florence_argp = 27.8_f64.to_radians();
    let florence_m_rad = 341.6_f64.to_radians();
    
    // Kepler Eq Solver approx
    // M = E - e sin E
    // For M=341.6 (approx -18.4), E approx M.
    // nu approx M + 2e sin M ...
    // Let's just use the helper if available, or approximate.
    // E approx -18.4 + 0.423 * sin(-18.4) = -18.4 - 7.5 = -25.9 deg
    // nu: tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2)
    // sqrt(1.423/0.577) = 1.57
    // tan(-13) = -0.23
    // tan(nu/2) = 1.57 * -0.23 = -0.36
    // nu/2 = -19.8 deg -> nu = -39.6 deg = 320.4 deg.
    let florence_nu = 320.4_f64.to_radians();

    let florence = Orbit::from_classical(SUN, ClassicalElements {
        p_km: florence_a * (1.0 - florence_e * florence_e),
        ecc: florence_e,
        inc_rad: florence_i,
        raan_rad: florence_raan,
        argp_rad: florence_argp,
        nu_rad: florence_nu,
    });

    println!("Initial Orbits defined (J2000).");

    // 4. Propagate to Target Epoch
    let earth_target = earth.propagate_seconds(target_epoch_seconds)
        .map_err(|e| format!("{:?}", e))?;
    let florence_target = florence.propagate_seconds(target_epoch_seconds)
        .map_err(|e| format!("{:?}", e))?;

    println!("Propagated to 2017-09-01.");

    // 5. Calculate Distance
    let r_earth = earth_target.state.r_km;
    let r_florence = florence_target.state.r_km;
    let dist_vec = r_florence - r_earth;
    let dist_km = dist_vec.norm();
    
    println!("Distance Earth-Florence: {:.0} km", dist_km);
    println!("(Note: Without precise ephemerides, this value is approximate. Expected ~7,060,160 km)");
    let error_pct = (dist_km - 7060160.0).abs() / 7060160.0 * 100.0;
    println!("Error: {:.1}% (due to mean elements and simplified propagation)", error_pct);

    // 6. Plotting
    // Plot Earth and Florence orbits and positions.
    let root = BitMapBackend::new("catch_asteroid.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Orbit of Florence and Earth (2017-09-01)", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-3.0 * AU_KM..3.0 * AU_KM, -3.0 * AU_KM..3.0 * AU_KM)?;

    chart.configure_mesh()
        .x_desc("x [km]")
        .y_desc("y [km]")
        .draw()?;

    // Helper to sample orbit
    let mut sample_orbit = |orbit: &Orbit, color: &RGBColor| -> Result<(), Box<dyn std::error::Error>> {
        let mut points = Vec::new();
        for i in 0..=360 {
            let nu = (i as f64).to_radians();
            let elems = orbit.classical();
            let orb_at_nu = Orbit::from_classical(orbit.attractor, ClassicalElements {
                nu_rad: nu,
                ..elems
            });
            let r = orb_at_nu.state.r_km;
            points.push((r.x, r.y));
        }
        chart.draw_series(LineSeries::new(points, color))?;
        Ok(())
    };

    // Plot Earth Orbit
    sample_orbit(&earth_target, &BLUE)?;
    // Plot Florence Orbit
    sample_orbit(&florence_target, &BLACK)?;

    // Plot Positions
    chart.draw_series(PointSeries::of_element(
        vec![(r_earth.x, r_earth.y)],
        5,
        ShapeStyle::from(&BLUE).filled(),
        &|c, s, st| {
            EmptyElement::at(c) + Circle::new((0,0), s, st) + Text::new("Earth", (10, 0), ("sans-serif", 15).into_font())
        }
    ))?;

    chart.draw_series(PointSeries::of_element(
        vec![(r_florence.x, r_florence.y)],
        5,
        ShapeStyle::from(&BLACK).filled(),
        &|c, s, st| {
            EmptyElement::at(c) + Circle::new((0,0), s, st) + Text::new("Florence", (10, 0), ("sans-serif", 15).into_font())
        }
    ))?;

    // Plot Sun
    chart.draw_series(PointSeries::of_element(
        vec![(0.0, 0.0)],
        10,
        ShapeStyle::from(&YELLOW).filled(),
        &|c, s, st| {
            EmptyElement::at(c) + Circle::new((0,0), s, st) + Text::new("Sun", (10, 0), ("sans-serif", 15).into_font())
        }
    ))?;

    println!("Generated catch_asteroid.png");

    Ok(())
}
