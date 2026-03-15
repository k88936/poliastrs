use poliastrs::bodies::{EARTH, SUN};
use poliastrs::twobody::orbit::Orbit;
use poliastrs::core::elements::ClassicalElements;
use poliastrs::twobody::maneuver::Maneuver;
use poliastrs::frames::Plane;

fn main() {
    println!("Poliastro Rust Quickstart");
    println!("=========================");

    // 1. Define Orbit from vectors
    println!("\n## 1. Defining orbit from vectors");
    // Data from Curtis, example 4.3
    let r_km = [-6045.0, -3490.0, 2500.0];
    let v_km_s = [-3.457, 6.618, 2.533];
    
    let orb = Orbit::from_vectors(EARTH, r_km, v_km_s);
    
    // Print orbit properties similar to Python output
    // Python: 7283 x 10293 km x 153.2 deg (GCRS) orbit around Earth (♁) at epoch J2000.000 (TT)
    let rp = orb.r_p_km();
    let ra = orb.r_a_km().unwrap_or(f64::INFINITY); // Handle hyperbolic/parabolic if needed
    let inc_deg = orb.inc_rad().to_degrees();
    
    println!("Orbit created from vectors:");
    println!("  Perigee: {:.0} km", rp);
    if ra.is_finite() {
        println!("  Apogee:  {:.0} km", ra);
    } else {
        println!("  Apogee:  Infinite");
    }
    println!("  Inclination: {:.1} deg", inc_deg);
    println!("  Frame: {:?}", orb.plane);
    println!("  Attractor: {}", orb.attractor.name);

    let (r, v) = orb.rv();
    println!("  Position vector: {:?}", r);
    println!("  Velocity vector: {:?}", v);

    // 2. From classical orbital elements
    println!("\n## 2. Defining orbit from classical elements");
    // Data for Mars at J2000 from JPL HORIZONS
    const AU_KM: f64 = 149_597_870.7;
    let a_km = 1.523679 * AU_KM;
    let ecc = 0.093315;
    let inc_rad = 1.85_f64.to_radians();
    let raan_rad = 49.562_f64.to_radians();
    let argp_rad = 286.537_f64.to_radians();
    let nu_rad = 23.33_f64.to_radians();

    // Calculate p from a and ecc
    let p_km = a_km * (1.0 - ecc * ecc);

    let coe = ClassicalElements {
        p_km,
        ecc,
        inc_rad,
        raan_rad,
        argp_rad,
        nu_rad,
    };

    let orb_mars = Orbit::from_classical(SUN, coe);
    println!("Mars Orbit created from classical elements:");
    println!("  Semimajor axis: {:.3} AU ({:.3} km)", orb_mars.a_km() / AU_KM, orb_mars.a_km());
    println!("  Eccentricity: {:.6}", orb_mars.ecc());
    println!("  Inclination: {:.2} deg", orb_mars.inc_rad().to_degrees());
    if let Some(period) = orb_mars.period_seconds() {
        println!("  Period: {:.2} days", period / 86400.0);
    } else {
        println!("  Period: Infinite");
    }

    // 3. Propagation
    println!("\n## 3. Propagation");
    // Example ISS orbit
    // ISS approx: 400km altitude, circular
    let r_iss = EARTH.mean_radius_km + 400.0;
    
    let iss_orbit = Orbit::from_classical(EARTH, ClassicalElements {
        p_km: r_iss, // circular, p = a
        ecc: 0.0,
        inc_rad: 51.6_f64.to_radians(),
        raan_rad: 0.0,
        argp_rad: 0.0,
        nu_rad: 0.0 // Start at periapsis (doesn't matter for circular)
    });
    
    println!("Initial ISS Orbit:");
    println!("  Epoch: {:.3} s (TDB relative to J2000)", iss_orbit.epoch_tdb_seconds);
    println!("  True Anomaly: {:.2} deg", iss_orbit.nu_rad().to_degrees());
    
    let dt_seconds = 30.0 * 60.0; // 30 minutes
    match iss_orbit.propagate_seconds(dt_seconds) {
        Ok(iss_30m) => {
             println!("ISS Orbit after 30 min propagation:");
             println!("  Epoch: {:.3} s", iss_30m.epoch_tdb_seconds);
             println!("  True Anomaly: {:.2} deg", iss_30m.nu_rad().to_degrees());
             
             // Check if mean motion matches Python example: 3.887 deg/min
             // n = sqrt(mu/a^3)
             let n_rad_s = (EARTH.mu_km3_s2 / iss_orbit.a_km().powi(3)).sqrt();
             let n_deg_min = n_rad_s.to_degrees() * 60.0;
             println!("  Mean motion: {:.3} deg/min", n_deg_min);
        },
        Err(e) => println!("Propagation error: {:?}", e),
    }

    // 4. Maneuvers (Hohmann)
    println!("\n## 4. Maneuvers (Hohmann Transfer)");
    // Initial orbit: 700 km altitude circular
    let r_i = EARTH.mean_radius_km + 700.0;
    let orb_i = Orbit::from_classical(EARTH, ClassicalElements {
        p_km: r_i,
        ecc: 0.0,
        inc_rad: 0.0,
        raan_rad: 0.0,
        argp_rad: 0.0,
        nu_rad: 0.0
    });
    
    println!("Initial Orbit Altitude: 700 km");
    
    // Target altitude: 36000 km (GEO)
    let alt_f = 36000.0;
    let r_f = EARTH.mean_radius_km + alt_f;
    println!("Target Orbit Altitude: {:.0} km", alt_f);
    
    // Hohmann transfer
    let hoh = Maneuver::hohmann(&orb_i, r_f);
    
    println!("Hohmann Transfer calculated:");
    println!("  Total Cost (Delta V): {:.3} km/s", hoh.get_total_cost());
    println!("  Transfer Time: {:.2} s ({:.2} hours)", hoh.get_total_time(), hoh.get_total_time() / 3600.0);
    
    let (dv1_t, dv1_v) = hoh.impulses[0];
    let (dv2_t, dv2_v) = hoh.impulses[1];
    
    println!("  Impulse 1 at t={:.0} s: |dv| = {:.3} km/s", dv1_t, dv1_v.norm());
    println!("  Impulse 2 at t={:.0} s: |dv| = {:.3} km/s", dv2_t, dv2_v.norm());

    // Apply maneuver
    let orb_f = orb_i.apply_maneuver(&hoh);
    println!("Final Orbit after maneuver:");
    println!("  Semimajor axis: {:.3} km", orb_f.a_km());
    println!("  Eccentricity: {:.6}", orb_f.ecc());

    // 5. Lambert Problem
    println!("\n## 5. Lambert Problem (Example 5.2 from Curtis)");
    // r0 = [5000, 10000, 2100] km
    // rf = [-14600, 2500, 7000] km
    // dt = 1 hour
    
    let r0_km = [5000.0, 10000.0, 2100.0];
    let rf_km = [-14600.0, 2500.0, 7000.0];
    let tof_seconds = 3600.0; // 1 hour

    // Create orbits at t0 and tf (velocity doesn't matter for Lambert input, only position and time)
    let orb0 = Orbit::from_vectors_at(EARTH, r0_km, [0.0, 0.0, 0.0], 0.0, Plane::EarthEquator);
    let orbf = Orbit::from_vectors_at(EARTH, rf_km, [0.0, 0.0, 0.0], tof_seconds, Plane::EarthEquator);

    println!("Solving Lambert problem for dt = {:.0} s", tof_seconds);
    match Maneuver::lambert(&orb0, &orbf) {
        Ok(man_lambert) => {
            let (_dv_a_t, dv_a_v) = man_lambert.impulses[0];
            let (_dv_b_t, dv_b_v) = man_lambert.impulses[1];
            
            // Expected results from test_iod.py test_curtis52:
            // va = [-5.9925, 1.9254, 3.2456] km/s
            // vb = [-3.3125, -4.1966, -0.38529] km/s
            
            println!("  Velocity at departure (va): {:.4?} km/s", dv_a_v);
            // Since target velocity was 0, dv_b = v_final - v_transfer_arrival = 0 - v_transfer_arrival = -v_transfer_arrival.
            // So v_transfer_arrival = -dv_b.
            // But we want vb (velocity at arrival on transfer orbit).
            // Maneuver impulses are (v_transfer_dep - v_initial) and (v_final - v_transfer_arr).
            // Here v_initial = 0, v_final = 0.
            // So dv_a = v_transfer_dep - 0 = va.
            // dv_b = 0 - v_transfer_arr = -vb.
            // So vb = -dv_b.
            
            println!("  Velocity at arrival (vb):   {:.4?} km/s", -dv_b_v);
        },
        Err(e) => println!("Lambert error: {}", e),
    }

    println!("\nQuickstart completed successfully!");
}
