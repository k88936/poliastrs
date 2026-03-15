use poliastrs::plotting::porkchop::PorkchopPlotter;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::bodies::SUN;
use poliastrs::core::elements::ClassicalElements;
use poliastrs::frames::Plane;

#[test]
fn test_porkchop_generation() {
    let earth_a = 149_597_870.7; // km
    let earth_ecc = 0.0167086;
    let earth_p = earth_a * (1.0 - earth_ecc * earth_ecc);
    
    let mars_a = 227_939_200.0; // km
    let mars_ecc = 0.0934025; // approx
    let mars_p = mars_a * (1.0 - mars_ecc * mars_ecc);
    
    let earth_orbit = Orbit::from_classical_at(
        SUN,
        ClassicalElements {
            p_km: earth_p,
            ecc: earth_ecc,
            inc_rad: 0.0,
            raan_rad: 0.0,
            argp_rad: 0.0,
            nu_rad: 0.0,
        },
        0.0,
        Plane::EarthEcliptic
    ).unwrap();
    
    let mars_orbit = Orbit::from_classical_at(
        SUN,
        ClassicalElements {
            p_km: mars_p,
            ecc: mars_ecc,
            inc_rad: 1.85_f64.to_radians(),
            raan_rad: 0.0,
            argp_rad: 0.0,
            nu_rad: 1.0, // Some initial separation
        },
        0.0,
        Plane::EarthEcliptic
    ).unwrap();
    
    // Launch span: 0 to 30 days
    // Arrival span: 200 to 400 days
    let launch_span = (0.0, 30.0 * 86400.0);
    let arrival_span = (200.0 * 86400.0, 400.0 * 86400.0);
    
    let plotter = PorkchopPlotter::new(earth_orbit, mars_orbit, launch_span, arrival_span);
    
    let path = std::path::Path::new("test_porkchop.png");
    let res = plotter.plot(path);
    if let Err(e) = &res {
        println!("Error generating porkchop plot: {:?}", e);
    }
    assert!(res.is_ok());
    
    // if path.exists() {
    //     let _ = std::fs::remove_file(path);
    // }
}
