use poliastrs::bodies::{EARTH, Body};
use poliastrs::twobody::orbit::Orbit;
use poliastrs::twobody::propagation::cowell::{propagate_with_events, Event};
use poliastrs::earth::atmosphere::COESA76;
use poliastrs::core::elements::ClassicalElements;
use nalgebra::Vector3;
use std::f64::consts::PI;

struct AltitudeCrossEvent {
    thresh_alt_km: f64,
    body_radius_km: f64,
    terminal: bool,
}

impl Event for AltitudeCrossEvent {
    fn evaluate(&self, _t: f64, r: &Vector3<f64>, _v: &Vector3<f64>) -> f64 {
        let r_norm = r.norm();
        let alt = r_norm - self.body_radius_km;
        alt - self.thresh_alt_km
    }

    fn is_terminal(&self) -> bool {
        self.terminal
    }
}

struct NodeCrossEvent {
    terminal: bool,
}

impl Event for NodeCrossEvent {
    fn evaluate(&self, _t: f64, r: &Vector3<f64>, _v: &Vector3<f64>) -> f64 {
        r.z
    }

    fn is_terminal(&self) -> bool {
        self.terminal
    }
}

fn main() {
    println!("Detecting Events with Cowell Propagator");
    
    // 1. Altitude Cross Event
    println!("\n1. Altitude Cross Event (Atmospheric Drag)");
    
    let orbit = Orbit::circular(EARTH, 150.0).unwrap(); // 150 km initial altitude
    let coesa = COESA76::new();
    
    // Parameters for drag
    let cd = 2.2;
    let area_m2 = PI / 4.0; 
    let mass_kg = 100.0;
    let area_km2 = area_m2 * 1e-6;
    let am_ratio = area_km2 / mass_kg;
    
    // Drag function closure
    let drag_perturbation = |t: f64, r: &Vector3<f64>, v: &Vector3<f64>| -> Vector3<f64> {
        let r_norm = r.norm();
        let alt_km = r_norm - EARTH.mean_radius_km;
        
        // Atmosphere limit
        if alt_km < 0.0 || alt_km > 1000.0 {
            return Vector3::zeros();
        }
        
        let rho_kg_m3 = coesa.density(alt_km);
        let rho_kg_km3 = rho_kg_m3 * 1e9;
        
        let w = Vector3::new(0.0, 0.0, EARTH.angular_velocity_rad_s());
        let v_atm = w.cross(r);
        let v_rel = v - v_atm;
        let v_rel_norm = v_rel.norm();
        
        if v_rel_norm < 1e-6 { return Vector3::zeros(); }
        
        -0.5 * rho_kg_km3 * v_rel_norm * cd * am_ratio * v_rel
    };
    
    let thresh_alt = 50.0; // Stop at 50 km
    let alt_event = AltitudeCrossEvent {
        thresh_alt_km: thresh_alt,
        body_radius_km: EARTH.mean_radius_km,
        terminal: true,
    };
    
    let events = vec![alt_event];
    
    // Propagate
    let res = propagate_with_events(
        EARTH.mu_km3_s2,
        orbit.state.r_km,
        orbit.state.v_km_s,
        24000.0, // Increased time to ensure decay happens (2400s might be too short for 150km->50km depending on drag)
        drag_perturbation,
        &events
    );
    
    match res {
        Ok((_r, _v, detected)) => {
             if let Some(ev) = detected.first() {
                 println!("Altitude threshold detected at t = {:.2} s. Event index: {}", ev.t, ev.event_index);
                 let alt = ev.r.norm() - EARTH.mean_radius_km;
                 println!("Altitude at detection: {:.4} km (Expected: {:.4} km)", alt, thresh_alt);
             } else {
                 println!("No event detected within duration.");
             }
        },
        Err(e) => println!("Error: {}", e),
    }

    // 2. Node Cross Event
    println!("\n2. Node Cross Event (No perturbation)");
    
    // Inclined orbit
    let coe = ClassicalElements {
        p_km: 7000.0,
        ecc: 0.001,
        inc_rad: 45.0_f64.to_radians(),
        raan_rad: 0.0,
        argp_rad: 0.0,
        nu_rad: 0.0,
    };
    let orbit2 = Orbit::from_classical(EARTH, coe);
    
    let node_event = NodeCrossEvent { terminal: true };
    let events2 = vec![node_event];
    
    let period = orbit2.period_seconds().unwrap();
    println!("Orbit period: {:.2} s", period);
    
    let res2 = propagate_with_events(
        EARTH.mu_km3_s2,
        orbit2.state.r_km,
        orbit2.state.v_km_s,
        period * 1.5, // > 1 orbit
        |_, _, _| Vector3::zeros(),
        &events2
    );
    
    match res2 {
        Ok((_r, _v, detected)) => {
             for ev in detected {
                 println!("Node crossed at t = {:.2} s. Z = {:.6} km. Event index: {}", ev.t, ev.r.z, ev.event_index);
             }
        },
        Err(e) => println!("Error: {}", e),
    }
    
    // 3. Multiple Events
    println!("\n3. Multiple Events (Altitude + Node)");
    // Use Box<dyn Event>
    let event1 = Box::new(AltitudeCrossEvent {
        thresh_alt_km: 145.0, // Higher threshold (150 -> 145)
        body_radius_km: EARTH.mean_radius_km,
        terminal: false, // Don't stop
    });
    let event2 = Box::new(NodeCrossEvent { terminal: false });
    
    let events3: Vec<Box<dyn Event>> = vec![event1, event2];
    
    // Reuse orbit and drag
    // Initial alt is 150. Orbit decays.
    let res3 = propagate_with_events(
        EARTH.mu_km3_s2,
        orbit.state.r_km, // Use first orbit
        orbit.state.v_km_s,
        6000.0, // 100 mins
        drag_perturbation, // Use drag
        &events3
    );
    
     match res3 {
        Ok((_r, _v, detected)) => {
             println!("Detected {} events:", detected.len());
             for ev in detected {
                 let type_str = if ev.event_index == 0 { "Altitude (145km)" } else { "Node" };
                 println!("Event: {} at t = {:.2} s", type_str, ev.t);
             }
        },
        Err(e) => println!("Error: {}", e),
    }
}
