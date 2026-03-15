use poliastrs::plotting::tisserand::{TisserandPlotter, TisserandKind};
use poliastrs::bodies::{EARTH, MARS, VENUS};

#[test]
fn test_tisserand_generation() {
    let mut plotter = TisserandPlotter::new(TisserandKind::Energy);
    
    plotter.plot(VENUS, (1.0, 14.0), 14, None);
    plotter.plot(EARTH, (1.0, 14.0), 14, None);
    plotter.plot(MARS, (1.0, 14.0), 14, None);
    
    let path = std::path::Path::new("test_tisserand_plot.png");
    let res = plotter.save(path);
    if let Err(e) = &res {
        println!("Error generating tisserand plot: {:?}", e);
    }
    assert!(res.is_ok());
    
    // if path.exists() {
    //     let _ = std::fs::remove_file(path);
    // }
}
