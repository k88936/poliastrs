use plotters::prelude::*;
use poliastrs::earth::atmosphere::COESA76;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let coesa76 = COESA76::new();

    // Generate data
    let mut data = Vec::new();
    let steps = 1000;
    let max_alt = 1000.0;
    
    for i in 0..=steps {
        let alt = (i as f64 / steps as f64) * max_alt;
        let (t, p, rho) = coesa76.properties(alt);
        data.push((alt, t, p, rho));
    }

    let root = BitMapBackend::new("atmospheric_models.png", (1200, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let areas = root.split_evenly((1, 3));

    // 1. Temperature
    let mut chart_t = ChartBuilder::on(&areas[0])
        .caption("Temperature vs Altitude", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0f64..2000.0f64, 0f64..1000.0f64)?; // T range 0-2000K, Alt 0-1000km

    chart_t.configure_mesh()
        .x_desc("Temperature [K]")
        .y_desc("Altitude [km]")
        .draw()?;

    chart_t.draw_series(LineSeries::new(
        data.iter().map(|(alt, t, _, _)| (*t, *alt)),
        &RED,
    ))?;

    // Add layer lines
    let layers = [
        (11.0, "Troposphere"),
        (47.0, "Stratosphere"),
        (86.0, "Mesosphere"),
        (1000.0, "Thermosphere"), // Just label
    ];
    
    // Annotations for layers (simplified)
    for &(h, name) in layers.iter() {
        if h < 1000.0 {
            chart_t.draw_series(std::iter::once(PathElement::new(
                vec![(0.0, h), (2000.0, h)],
                BLACK.mix(0.5).filled(),
            )))?;
        }
        
        let y_pos = if h == 1000.0 { 950.0 } else { h + 10.0 };
        chart_t.draw_series(std::iter::once(Text::new(
            name,
            (1500.0, y_pos),
            ("sans-serif", 15).into_font(),
        )))?;
    }


    // 2. Pressure (Log scale)
    // Pressure ranges from ~10^5 Pa to ~10^-11 Pa
    let min_p = data.last().unwrap().2.max(1e-12); // Avoid 0
    let max_p = data.first().unwrap().2;

    let mut chart_p = ChartBuilder::on(&areas[1])
        .caption("Pressure vs Altitude", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d((min_p..max_p).log_scale(), 0f64..1000.0f64)?;

    chart_p.configure_mesh()
        .x_desc("Pressure [Pa]")
        .y_desc("Altitude [km]")
        .draw()?;

    chart_p.draw_series(LineSeries::new(
        data.iter().map(|(alt, _, p, _)| (*p, *alt)),
        &BLUE,
    ))?;

    // 3. Density (Log scale)
    // Density ranges from ~1.2 to ~10^-13
    let min_rho = data.last().unwrap().3.max(1e-14);
    let max_rho = data.first().unwrap().3;

    let mut chart_rho = ChartBuilder::on(&areas[2])
        .caption("Density vs Altitude", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d((min_rho..max_rho).log_scale(), 0f64..1000.0f64)?;

    chart_rho.configure_mesh()
        .x_desc("Density [kg/m^3]")
        .y_desc("Altitude [km]")
        .draw()?;

    chart_rho.draw_series(LineSeries::new(
        data.iter().map(|(alt, _, _, rho)| (*rho, *alt)),
        &GREEN,
    ))?;
    
    println!("Generated atmospheric_models.png");

    Ok(())
}
