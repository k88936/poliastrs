use poliastrs::bodies::EARTH;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::spheroid_location::SpheroidLocation;
use nalgebra::{Vector3, Rotation3};
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ISS Orbit (Approximate)
    // Altitude ~400km, Inclination ~51.6 deg.
    let h = 400.0;
    let r = EARTH.equatorial_radius_km + h;
    let v_circ = (EARTH.mu_km3_s2 / r).sqrt();
    
    // Create orbit using vectors (easier)
    // Position at x-axis, Velocity inclined
    let inc = 51.6_f64.to_radians();
    let r0 = Vector3::new(r, 0.0, 0.0);
    let v0 = Vector3::new(0.0, v_circ * inc.cos(), v_circ * inc.sin());
    
    let iss = Orbit::from_vectors(EARTH, [r0.x, r0.y, r0.z], [v0.x, v0.y, v0.z]);
    
    // Simulation parameters
    let duration = 3.0 * 3600.0; // 3 hours
    let steps = 300;
    let dt = duration / steps as f64;
    
    let mut groundtrack = Vec::new();
    
    // Earth rotation
    let w_earth = EARTH.angular_velocity_rad_s();
    // Initial GMST (assume 0 for simplicity)
    let gmst0 = 0.0;
    
    for i in 0..=steps {
        let t = i as f64 * dt;
        let orbit_t = iss.propagate_seconds(t).unwrap();
        let r_eci = orbit_t.state.r_km; // Vector3
        
        // Convert to ECEF
        // Angle increases with time.
        // Rotation of ECEF frame relative to ECI is theta.
        // r_ecef = Rot(-theta) * r_eci
        let theta = gmst0 + w_earth * t;
        let rot = Rotation3::from_axis_angle(&Vector3::z_axis(), -theta);
        let r_ecef = rot * r_eci;
        
        // Convert to Lat/Lon
        let (lon_rad, lat_rad, _h) = SpheroidLocation::cartesian_to_ellipsoidal(
            EARTH, r_ecef.x, r_ecef.y, r_ecef.z
        );
        
        // Convert to degrees
        let lat_deg = lat_rad.to_degrees();
        let lon_deg = lon_rad.to_degrees();
        
        groundtrack.push((lon_deg, lat_deg));
    }
    
    // Plot
    let root = BitMapBackend::new("groundtrack.png", (1024, 512)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption("ISS Ground Track", ("sans-serif", 30))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-180.0..180.0, -90.0..90.0)?;
        
    chart.configure_mesh()
        .x_desc("Longitude")
        .y_desc("Latitude")
        .draw()?;
        
    // Points
    chart.draw_series(PointSeries::of_element(
        groundtrack,
        2,
        &RED,
        &|c, s, st| {
            return EmptyElement::at(c)    + Circle::new((0,0),s,st.filled());
        },
    ))?;
    
    println!("Generated groundtrack.png");
    
    Ok(())
}
