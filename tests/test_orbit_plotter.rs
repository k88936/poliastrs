use poliastrs::plotting::orbit_plotter::OrbitPlotter;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::bodies::EARTH;

#[test]
fn test_orbit_plotter_basic() {
    let mut plotter = OrbitPlotter::new();
    
    // Create an orbit
    let orbit = Orbit::circular(EARTH, 7000.0).unwrap();
    
    plotter.plot(&orbit, Some("Test Orbit"));
    
    let path = std::path::Path::new("test_orbit_plot_2d.png");
    let res = plotter.save_2d(path);
    if let Err(e) = &res {
        println!("Error saving plot: {:?}", e);
    }
    assert!(res.is_ok());
    assert!(path.exists());
    
    // Clean up
    // let _ = std::fs::remove_file(path);
}

#[test]
fn test_orbit_plotter_3d_generation() {
    let mut plotter = OrbitPlotter::new();
    let orbit = Orbit::circular(EARTH, 7000.0).unwrap();
    plotter.plot(&orbit, Some("3D Orbit"));
    
    // We haven't implemented save_3d yet, but we can verify trajectories_3d is populated
    assert_eq!(plotter.trajectories_3d.len(), 1);
    assert!(!plotter.trajectories_3d[0].0.is_empty());
}

#[test]
fn test_orbit_plotter_set_attractor() {
    let mut plotter = OrbitPlotter::new();
    plotter.set_attractor(EARTH);
    assert_eq!(plotter.attractor, Some(EARTH));
}

#[test]
#[should_panic]
fn test_orbit_plotter_attractor_mismatch() {
    let mut plotter = OrbitPlotter::new();
    plotter.set_attractor(EARTH);
    
    // Try to set different attractor
    plotter.set_attractor(poliastrs::bodies::MARS);
}
