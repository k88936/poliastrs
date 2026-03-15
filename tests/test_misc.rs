use poliastrs::plotting::misc::plot_solar_system;

#[test]
fn test_plot_inner_solar_system() {
    let plotter = plot_solar_system(false, 0.0);
    
    // Verify 4 planets
    assert_eq!(plotter.trajectories.len(), 4);
    
    let path = std::path::Path::new("test_inner_solar_system.png");
    let res = plotter.save_2d(path);
    assert!(res.is_ok());
    // if path.exists() {
    //     let _ = std::fs::remove_file(path);
    // }
}

#[test]
fn test_plot_outer_solar_system() {
    let plotter = plot_solar_system(true, 0.0);
    
    // Verify 5 planets
    assert_eq!(plotter.trajectories.len(), 5);
    
    let path = std::path::Path::new("test_outer_solar_system.png");
    let res = plotter.save_2d(path);
    assert!(res.is_ok());
    // if path.exists() {
    //     let _ = std::fs::remove_file(path);
    // }
}
