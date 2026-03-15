use poliastrs::plotting::gabbard::GabbardPlotter;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::bodies::EARTH;
use poliastrs::core::elements::ClassicalElements;
use poliastrs::frames::Plane;

fn create_orbit(a: f64, ecc: f64, inc_deg: f64, raan_deg: f64, argp_deg: f64, nu_deg: f64) -> Orbit {
    let p = a * (1.0 - ecc * ecc);
    let inc = inc_deg.to_radians();
    let raan = raan_deg.to_radians();
    let argp = argp_deg.to_radians();
    let nu = nu_deg.to_radians();
    
    Orbit::from_classical_at(
        EARTH,
        ClassicalElements {
            p_km: p,
            ecc,
            inc_rad: inc,
            raan_rad: raan,
            argp_rad: argp,
            nu_rad: nu,
        },
        0.0,
        Plane::EarthEquator
    ).unwrap()
}

#[test]
fn test_gabbard_plot_generation() {
    let mut plotter = GabbardPlotter::new(false);
    
    // Create debris orbits (subset of Python test data)
    let orbits = vec![
        create_orbit(6828.193338988509, 0.0062140354170737815, 82.69387440482602, 37.33894561668519, 200.62393574484153, -117.55203086408737),
        create_orbit(6821.922877133498, 0.0023241234515223646, 82.65684766470754, 36.3401421924121, 125.29597430617513, -151.64963315597913),
        create_orbit(6836.825166360441, 0.004635589373624103, 82.69764910622918, 36.757861621556614, 44.219092511353594, 133.63349740950568),
    ];

    plotter.plot_orbits(&orbits, Some("COSMOS 1408 DEB"));
    
    let path = std::path::Path::new("test_gabbard_plot.png");
    let res = plotter.save(path);
    if let Err(e) = &res {
        println!("Error saving plot: {:?}", e);
    }
    assert!(res.is_ok());
    assert!(path.exists());
    
    // Clean up
    // let _ = std::fs::remove_file(path);
}

#[test]
fn test_gabbard_dark_mode() {
    let plotter = GabbardPlotter::new(true);
    assert!(plotter.dark_mode);
}
