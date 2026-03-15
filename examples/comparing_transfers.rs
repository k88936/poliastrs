use poliastrs::bodies::EARTH;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::twobody::maneuver::Maneuver;
use poliastrs::core::elements::ClassicalElements;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Comparing Hohmann and Bi-elliptic Transfers");
    
    // Initial Orbit: Circular 7000 km (approx LEO)
    let r_i_km = 7000.0;
    let orb_i = Orbit::from_classical(EARTH, ClassicalElements {
        p_km: r_i_km, // circular, p = a
        ecc: 0.0,
        inc_rad: 0.0,
        raan_rad: 0.0,
        argp_rad: 0.0,
        nu_rad: 0.0,
    });

    let v_i = orb_i.state.v_km_s.norm();

    // R range: 2 to 75
    let r_ratios: Vec<f64> = (0..100).map(|i| 2.0 + (75.0 - 2.0) * (i as f64 / 99.0)).collect();
    
    // Rstar values (ratio of r_b / r_i)
    let r_stars = vec![15.58, 40.0, 60.0, 100.0, 200.0, 10000.0]; // 10000 as approx for inf
    let colors = vec![&RED, &BLUE, &GREEN, &MAGENTA, &CYAN, &BLACK];
    let labels = vec!["15.58", "40", "60", "100", "200", "inf"];

    // Collect Data
    let mut hohmann_costs = Vec::new();
    let mut bielliptic_costs = vec![Vec::new(); r_stars.len()];

    for &r in &r_ratios {
        let r_f = r * r_i_km;
        
        // Hohmann
        let man_h = Maneuver::hohmann(&orb_i, r_f);
        hohmann_costs.push((r, man_h.get_total_cost() / v_i));

        // Bi-elliptic
        for (j, &r_star) in r_stars.iter().enumerate() {
            let r_b = r_star * r_i_km;
            // Check if r_b > r_f. Bi-elliptic implies going out further.
            // If r_b < r_f, it's just a less efficient transfer? 
            // The formula works regardless, but typically r_b > max(r_i, r_f).
            // Example sets Rstar fixed.
            
            let man_b = Maneuver::bielliptic(&orb_i, r_b, r_f);
            bielliptic_costs[j].push((r, man_b.get_total_cost() / v_i));
        }
    }

    // Plotting
    let root = BitMapBackend::new("comparing_transfers.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Hohmann vs Bi-elliptic Transfers", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(2.0f64..75.0f64, 0.35f64..0.60f64)?;

    chart.configure_mesh()
        .x_desc("R (r_f / r_i)")
        .y_desc("Delta-V / v_i")
        .draw()?;

    // Plot Hohmann (Thick Red line in example, but let's use standard color)
    // Actually example uses 'l' color for bielliptic too?
    // "for jj in range(len(Rstar)): ax.plot(..., color=l.get_color())"
    // Wait, example plots Hohmann and then plots bielliptic with same color?
    // No, `(l,) = ax.plot(R, hohmann_data, lw=2)` returns line `l`.
    // Then `color=l.get_color()` reuses the color.
    // So all lines are same color? That's confusing.
    // Ah, usually they are different to distinguish Rstar.
    // But the code says `color=l.get_color()`.
    // Maybe they rely on line style? Or maybe it's just one color for all?
    // Let's use different colors for clarity.

    chart.draw_series(LineSeries::new(
        hohmann_costs.clone(),
        BLACK.stroke_width(2),
    ))?.label("Hohmann").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK.stroke_width(2)));

    for (j, costs) in bielliptic_costs.iter().enumerate() {
        let color = colors[j % colors.len()];
        chart.draw_series(LineSeries::new(
            costs.clone(),
            color,
        ))?.label(format!("R*={}", labels[j]))
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color));
    }
    
    // Vertical lines
    
    // Draw vertical line at 11.94
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(11.94, 0.35), (11.94, 0.6)],
        BLACK.mix(0.5),
    )))?;
    
    // Max Hohmann cost R
    // Finding max of hohmann_costs
    let max_hohmann = hohmann_costs.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap();
    let r_max = max_hohmann.0;
    
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(r_max, 0.35), (r_max, 0.6)],
        BLACK.mix(0.5),
    )))?;

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    println!("Generated comparing_transfers.png");

    Ok(())
}
