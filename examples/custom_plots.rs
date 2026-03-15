use poliastrs::bodies::{EARTH, MARS, JUPITER, SUN, Body};
use poliastrs::twobody::orbit::Orbit;
use poliastrs::core::elements::ClassicalElements;
use plotters::prelude::*;
use std::f64::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Customising Static Orbit Plots");
    
    // 1. Define Orbits (Approximate Mean Elements J2000)
    let earth = get_mean_orbit(EARTH, 1.000, 0.0167, 0.0, 0.0, 102.9, 100.5);
    let mars = get_mean_orbit(MARS, 1.524, 0.0934, 1.85, 49.5, 286.5, 19.4);
    let jupiter = get_mean_orbit(JUPITER, 5.203, 0.0484, 1.30, 100.5, 274.3, 20.0);

    // 2. Plotting
    let root = BitMapBackend::new("custom_plots.png", (800, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Customized Orbit Plots", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-6.0f64..6.0f64, -6.0f64..6.0f64)?; // AU units

    chart.configure_mesh()
        .x_desc("x [AU]")
        .y_desc("y [AU]")
        .draw()?;

    const AU_KM: f64 = 149_597_870.7;

    // Helper to sample orbit points in AU
    let sample_orbit = |orbit: &Orbit| -> Vec<(f64, f64)> {
        let mut points = Vec::new();
        for i in 0..=360 {
            let nu = (i as f64).to_radians();
            let elems = orbit.classical();
            let orb_at_nu = Orbit::from_classical(orbit.attractor, ClassicalElements {
                nu_rad: nu,
                ..elems
            });
            let r = orb_at_nu.state.r_km;
            points.push((r.x / AU_KM, r.y / AU_KM));
        }
        points
    };

    // Plot Earth: Solid Line, Hexagon Marker
    let earth_points = sample_orbit(&earth);
    let earth_pos = earth_points[0]; // Just pick first point as "current" position
    
    chart.draw_series(LineSeries::new(
        earth_points.clone(),
        BLUE.stroke_width(2),
    ))?.label("Earth").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE.stroke_width(2)));

    // Hexagon Marker
    chart.draw_series(std::iter::once(earth_pos).map(|(x, y)| {
        let mut points = Vec::new();
        let s = 10.0;
        for i in 0..6 {
            let angle = (i as f64) * PI / 3.0;
            points.push((
                (s * angle.cos()) as i32,
                (s * angle.sin()) as i32,
            ));
        }
        EmptyElement::at((x, y)) + Polygon::new(points, BLUE.filled())
    }))?;


    // Plot Mars: Dashed Line (Simulated), Diamond Marker
    let mars_points = sample_orbit(&mars);
    let mars_pos = mars_points[90]; // 90 deg position
    
    // Simulated dashed line (using dots)
    chart.draw_series(
        mars_points.iter().step_by(5).map(|(x, y)| Circle::new((*x, *y), 1, RED.filled()))
    )?.label("Mars").legend(|(x, y)| Circle::new((x + 10, y), 2, RED.filled()));

    // Diamond Marker
    chart.draw_series(std::iter::once(mars_pos).map(|(x, y)| {
        let s = 10;
        let points = vec![
            (0, -s),
            (s, 0),
            (0, s),
            (-s, 0),
        ];
        EmptyElement::at((x, y)) + Polygon::new(points, RED.filled())
    }))?;

    // Plot Jupiter: No Line (just marker), Star Marker
    let jupiter_points = sample_orbit(&jupiter);
    let jupiter_pos = jupiter_points[180];
    
    // Star Marker
    chart.draw_series(std::iter::once(jupiter_pos).map(|(x, y)| {
        let mut points = Vec::new();
        let s = 15.0;
        let outer_r = s;
        let inner_r = s * 0.4;
        for i in 0..10 {
            let angle = (i as f64) * PI / 5.0 - PI / 2.0;
            let r = if i % 2 == 0 { outer_r } else { inner_r };
            points.push((
                (r * angle.cos()) as i32,
                (r * angle.sin()) as i32,
            ));
        }
        EmptyElement::at((x, y)) + Polygon::new(points, GREEN.filled())
    }))?.label("Jupiter").legend(|(x, y)| {
             let mut points = Vec::new();
            let s = 5.0;
            let outer_r = s;
            let inner_r = s * 0.4;
            for i in 0..10 {
                let angle = (i as f64) * PI / 5.0 - PI / 2.0;
                let r = if i % 2 == 0 { outer_r } else { inner_r };
                points.push((
                    x + (r * angle.cos()) as i32,
                    y + (r * angle.sin()) as i32,
                ));
            }
            Polygon::new(points, GREEN.filled())
    });

    // Sun
    chart.draw_series(std::iter::once((0.0, 0.0)).map(|(x, y)| {
        Circle::new((x, y), 10, YELLOW.filled())
    }))?;

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    println!("Generated custom_plots.png");
    Ok(())
}

fn get_mean_orbit(_body: Body, a_au: f64, ecc: f64, inc_deg: f64, raan_deg: f64, argp_deg: f64, nu_deg: f64) -> Orbit {
    const AU_KM: f64 = 149_597_870.7;
    let a_km = a_au * AU_KM;
    Orbit::from_classical(SUN, ClassicalElements {
        p_km: a_km * (1.0 - ecc * ecc),
        ecc,
        inc_rad: inc_deg.to_radians(),
        raan_rad: raan_deg.to_radians(),
        argp_rad: argp_deg.to_radians(),
        nu_rad: nu_deg.to_radians(),
    })
}
