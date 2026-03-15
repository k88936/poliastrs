use poliastrs::bodies::EARTH;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::twobody::propagation::cowell::propagate;
use nalgebra::Vector3;

fn main() {
    println!("Cowell Propagation with Perturbations");
    
    // Initial Orbit from Python Quickstart
    let r0 = [-2384.46, 5729.01, 3050.46];
    let v0 = [-7.36138, -2.98997, 1.64354];
    
    let initial = Orbit::from_vectors(EARTH, r0, v0);
    println!("Initial Orbit:");
    println!("  Perigee: {:.0} km", initial.r_p_km());
    println!("  Apogee:  {:.0} km", initial.r_a_km().unwrap_or(0.0));
    println!("  Inclination: {:.2} deg", initial.inc_rad().to_degrees());
    
    let r0_vec = Vector3::new(r0[0], r0[1], r0[2]);
    let v0_vec = Vector3::new(v0[0], v0[1], v0[2]);
    
    let tof = 3.0 * 86400.0; // 3 days
    
    // Define perturbing acceleration: constant acceleration aligned with velocity
    // a = 1e-5 * v / |v|
    // Note: 1e-5 km/s^2 is actually quite large thrust! 
    // In Python example: 1e-5 * v_vec / norm_v. Units are km/s^2 if v is km/s.
    let ad = |_t: f64, _r: &Vector3<f64>, v: &Vector3<f64>| -> Vector3<f64> {
        let v_norm = v.norm();
        if v_norm > 0.0 {
            1e-5 * v / v_norm
        } else {
            Vector3::zeros()
        }
    };
    
    println!("Propagating for {:.0} seconds (3 days)...", tof);
    
    match propagate(EARTH.mu_km3_s2, r0_vec, v0_vec, tof, ad) {
        Ok((r_final, v_final)) => {
            println!("Propagation successful");
            // Create final orbit
            let final_orbit = Orbit::from_vectors_at(
                EARTH, 
                [r_final.x, r_final.y, r_final.z], 
                [v_final.x, v_final.y, v_final.z], 
                initial.epoch_tdb_seconds + tof, 
                initial.plane
            );
            
            println!("Final Orbit after 3 days:");
            let rp = final_orbit.r_p_km();
            let ra = final_orbit.r_a_km().unwrap_or(0.0);
            
            println!("  Perigee Radius: {:.0} km", rp);
            println!("  Apogee Radius:  {:.0} km", ra);
            println!("  Inclination: {:.1} deg", final_orbit.inc_rad().to_degrees());
            
            // Python expected: 18255 x 21848 km x 28.0 deg
            // Let's see if we match.
        },
        Err(e) => println!("Propagation failed: {}", e),
    }
}
