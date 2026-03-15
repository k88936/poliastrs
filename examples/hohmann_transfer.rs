use poliastrs::bodies::EARTH;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::core::elements::ClassicalElements;
use poliastrs::twobody::maneuver::Maneuver;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Studying Hohmann Transfers");
    
    // Initial orbit: circular, 800 km altitude
    let r_i_km = EARTH.mean_radius_km + 800.0;
    let orb_i = Orbit::from_classical(EARTH, ClassicalElements {
        p_km: r_i_km,
        ecc: 0.0,
        inc_rad: 0.0,
        raan_rad: 0.0,
        argp_rad: 0.0,
        nu_rad: 0.0,
    });
    
    // Initial velocity magnitude
    // For circular orbit v = sqrt(mu/r)
    let v_i = (EARTH.mu_km3_s2 / r_i_km).sqrt();
    
    println!("Initial Orbit Radius: {:.2} km", r_i_km);
    println!("Initial Velocity: {:.4} km/s", v_i);
    
    // Range of R = r_f / r_i from 1 to 100
    let n_points = 1000;
    let r_ratios: Vec<f64> = (0..n_points).map(|i| 1.0 + (99.0 * i as f64 / (n_points - 1) as f64)).collect();
    
    let mut dv_a_normalized: Vec<f64> = Vec::with_capacity(n_points);
    let mut dv_b_normalized: Vec<f64> = Vec::with_capacity(n_points);
    let mut total_cost_normalized: Vec<f64> = Vec::with_capacity(n_points);
    
    for &ratio in &r_ratios {
        let r_f = r_i_km * ratio;
        let man = Maneuver::hohmann(&orb_i, r_f);
        
        let dv_a = man.impulses[0].1.norm();
        let dv_b = man.impulses[1].1.norm();
        
        dv_a_normalized.push(dv_a / v_i);
        dv_b_normalized.push(dv_b / v_i);
        total_cost_normalized.push((dv_a + dv_b) / v_i);
    }
    
    // Plotting
    let root = BitMapBackend::new("hohmann_transfer.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("Hohmann Transfer Costs", ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(1.0f64..100.0f64, 0.0f64..0.7f64)?;
        
    chart.configure_mesh()
        .x_desc("R (r_f / r_i)")
        .y_desc("Delta V / v_i")
        .draw()?;
        
    // Plot First Impulse
    chart.draw_series(LineSeries::new(
        r_ratios.iter().zip(dv_a_normalized.iter()).map(|(&x, &y)| (x, y)),
        &BLUE,
    ))?
    .label("First Impulse")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    
    // Plot Second Impulse
    chart.draw_series(LineSeries::new(
        r_ratios.iter().zip(dv_b_normalized.iter()).map(|(&x, &y)| (x, y)),
        &RED,
    ))?
    .label("Second Impulse")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    
    // Plot Total Cost
    chart.draw_series(LineSeries::new(
        r_ratios.iter().zip(total_cost_normalized.iter()).map(|(&x, &y)| (x, y)),
        &GREEN,
    ))?
    .label("Total Cost")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));
    
    // Reference lines
    // sqrt(2) - 1 approx 0.414
    let limit_val = 2.0_f64.sqrt() - 1.0;
    chart.draw_series(LineSeries::new(
        vec![(1.0, limit_val), (100.0, limit_val)],
        &BLACK.mix(0.5),
    ))?;
    
    // 1 / sqrt(R) line?
    // The Python code plots: (1 / np.sqrt(r_f_vector / r_i)).value
    // But normalized by v_i?
    // Wait, the Python code:
    // ax.plot((r_f_vector / r_i).value, (1 / np.sqrt(r_f_vector / r_i)).value, "k--")
    // This is 1/sqrt(R).
    // Let's plot it too.
    chart.draw_series(LineSeries::new(
        r_ratios.iter().map(|&x| (x, 1.0 / x.sqrt())),
        &BLACK.mix(0.5),
    ))?;
    
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
        
    println!("Plot saved to hohmann_transfer.png");
    
    Ok(())
}
