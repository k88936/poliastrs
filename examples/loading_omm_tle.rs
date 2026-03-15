use sgp4::{Elements, Constants};
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Loading OMM and TLE satellite data");

    // Hardcoded TLE for ISS (ZARYA)
    let tle_line1 = "1 25544U 98067A   26073.80519411  .00010373  00000+0  20105-3 0  9991";
    let tle_line2 = "2 25544  51.6329  47.0659 0006371 190.9724 169.1126 15.48281436557126";
    
    let elements = Elements::from_tle(
        Some("ISS (ZARYA)".to_string()),
        tle_line1.as_bytes(),
        tle_line2.as_bytes(),
    )?;

    println!("Loaded TLE for: {:?}", elements.object_name);

    // Propagate
    // Epoch of TLE is 2020 day 262.19325852
    // Let's propagate for 90 minutes from epoch.
    
    // sgp4 crate uses minutes since epoch for propagation?
    // Constants::from_elements(&elements) creates a predictor.
    let constants = Constants::from_elements(&elements)?;

    // Plotting
    let root = BitMapBackend::new("loading_omm_tle.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("ISS Orbit from TLE", ("sans-serif", 30))
        .margin(20)
        .build_cartesian_3d(-7000.0..7000.0, -7000.0..7000.0, -7000.0..7000.0)?;

    chart.configure_axes().draw()?;

    // Earth
    chart.draw_series(PointSeries::of_element(
        vec![(0.0, 0.0, 0.0)],
        5,
        &BLUE,
        &|c, s, st| {
            return EmptyElement::at(c) + Circle::new((0,0), s, st.filled());
        },
    ))?;

    // Propagate and collect points
    let mut points = Vec::new();
    for i in 0..100 {
        let t_min = (i as f64) * 90.0 / 100.0; // 0 to 90 minutes
        let prediction = constants.propagate(sgp4::MinutesSinceEpoch(t_min))?;
        let r = prediction.position; // [f64; 3] in km
        points.push((r[0], r[1], r[2]));
    }

    chart.draw_series(LineSeries::new(
        points,
        &RED,
    ))?.label("ISS").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    println!("Generated loading_omm_tle.png");
    
    Ok(())
}
