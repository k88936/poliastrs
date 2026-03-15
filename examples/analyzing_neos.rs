use poliastrs::bodies::SUN;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::core::elements::ClassicalElements;
use poliastrs::frames::Plane;

fn main() {
    println!("Analyzing Near-Earth Objects (NEOs)");
    println!("===================================");
    println!("Note: Orbit.from_sbdb() is not yet implemented in Rust.");
    println!("This example uses hardcoded orbital elements for demonstration.\n");

    const AU_KM: f64 = 149_597_870.7;

    // 1. Define Eros (433)
    // Approximate elements (Epoch J2000)
    let eros_a = 1.4582 * AU_KM;
    let eros_e = 0.2226;
    let eros_i = 10.82_f64.to_radians();
    let eros_omega = 304.32_f64.to_radians();
    let eros_w = 178.82_f64.to_radians();
    let eros_nu = 0.0_f64.to_radians(); // Arbitrary true anomaly

    let eros = Orbit::from_classical(SUN, ClassicalElements {
        p_km: eros_a * (1.0 - eros_e * eros_e),
        ecc: eros_e,
        inc_rad: eros_i,
        raan_rad: eros_omega,
        argp_rad: eros_w,
        nu_rad: eros_nu,
    });

    println!("Defined Orbit for Eros:");
    println!("  Semimajor Axis: {:.3} AU", eros.a_km() / AU_KM);
    println!("  Eccentricity:   {:.3}", eros.ecc());
    println!("  Inclination:    {:.1} deg", eros.inc_rad().to_degrees());
    println!("");

    // 2. Define Ganymed (1036) - Not the moon of Jupiter!
    // Largest NEO.
    let ganymed_a = 2.6628 * AU_KM;
    let ganymed_e = 0.5337;
    let ganymed_i = 26.69_f64.to_radians();
    let ganymed_omega = 130.68_f64.to_radians();
    let ganymed_w = 132.43_f64.to_radians();
    let ganymed_nu = 0.0;

    let ganymed = Orbit::from_classical(SUN, ClassicalElements {
        p_km: ganymed_a * (1.0 - ganymed_e * ganymed_e),
        ecc: ganymed_e,
        inc_rad: ganymed_i,
        raan_rad: ganymed_omega,
        argp_rad: ganymed_w,
        nu_rad: ganymed_nu,
    });

    println!("Defined Orbit for Ganymed (1036):");
    println!("  Semimajor Axis: {:.3} AU", ganymed.a_km() / AU_KM);
    println!("  Eccentricity:   {:.3}", ganymed.ecc());
    println!("  Inclination:    {:.1} deg", ganymed.inc_rad().to_degrees());
    println!("");

    // 3. Define Amor (1221)
    let amor_a = 1.9186 * AU_KM;
    let amor_e = 0.4346;
    let amor_i = 11.87_f64.to_radians();
    let amor_omega = 171.32_f64.to_radians();
    let amor_w = 26.31_f64.to_radians();
    let amor_nu = 0.0;

    let amor = Orbit::from_classical(SUN, ClassicalElements {
        p_km: amor_a * (1.0 - amor_e * amor_e),
        ecc: amor_e,
        inc_rad: amor_i,
        raan_rad: amor_omega,
        argp_rad: amor_w,
        nu_rad: amor_nu,
    });

    println!("Defined Orbit for Amor (1221):");
    println!("  Semimajor Axis: {:.3} AU", amor.a_km() / AU_KM);
    println!("  Eccentricity:   {:.3}", amor.ecc());
    println!("  Inclination:    {:.1} deg", amor.inc_rad().to_degrees());
    println!("");

    // 4. Propagation example
    // Propagate Eros to a new epoch (e.g., +30 days)
    let dt_seconds = 30.0 * 86400.0;
    println!("Propagating Eros by 30 days...");
    match eros.propagate_seconds(dt_seconds) {
        Ok(eros_future) => {
            println!("  Initial Epoch: {:.1} s (TDB)", eros.epoch_tdb_seconds);
            println!("  Final Epoch:   {:.1} s (TDB)", eros_future.epoch_tdb_seconds);
            println!("  Final True Anomaly: {:.1} deg", eros_future.nu_rad().to_degrees());
        },
        Err(e) => println!("Propagation failed: {:?}", e),
    }
}
