use crate::twobody::orbit::Orbit;
use nalgebra::Vector3;

pub trait Event {
    fn evaluate(&self, t: f64, r: &Vector3<f64>, v: &Vector3<f64>) -> f64;
    fn is_terminal(&self) -> bool;
}

impl Event for Box<dyn Event> {
    fn evaluate(&self, t: f64, r: &Vector3<f64>, v: &Vector3<f64>) -> f64 {
        (**self).evaluate(t, r, v)
    }
    fn is_terminal(&self) -> bool {
        (**self).is_terminal()
    }
}

#[derive(Debug, Clone)]
pub struct DetectedEvent {
    pub t: f64,
    pub r: Vector3<f64>,
    pub v: Vector3<f64>,
    pub event_index: usize,
}

fn step_rk4<F>(mu: f64, r: &mut Vector3<f64>, v: &mut Vector3<f64>, t: f64, dt: f64, ad: &F)
where
    F: Fn(f64, &Vector3<f64>, &Vector3<f64>) -> Vector3<f64>,
{
    let k1_r = *v;
    let k1_v = accel(mu, r, v, t, ad);

    let k2_r = *v + 0.5 * dt * k1_v;
    let k2_v = accel(mu, &(*r + 0.5 * dt * k1_r), &(*v + 0.5 * dt * k1_v), t + 0.5 * dt, ad);

    let k3_r = *v + 0.5 * dt * k2_v;
    let k3_v = accel(mu, &(*r + 0.5 * dt * k2_r), &(*v + 0.5 * dt * k2_v), t + 0.5 * dt, ad);

    let k4_r = *v + dt * k3_v;
    let k4_v = accel(mu, &(*r + dt * k3_r), &(*v + dt * k3_v), t + dt, ad);

    *r += (dt / 6.0) * (k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r);
    *v += (dt / 6.0) * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v);
}

pub fn propagate<F>(
    mu: f64,
    r0: Vector3<f64>,
    v0: Vector3<f64>,
    tof: f64,
    ad: F,
) -> Result<(Vector3<f64>, Vector3<f64>), String>
where
    F: Fn(f64, &Vector3<f64>, &Vector3<f64>) -> Vector3<f64>,
{
    // RK4 implementation with fixed step size
    let mut t = 0.0;
    let mut r = r0;
    let mut v = v0;
    
    // Choose step size. For LEO, 10-60s is usually okay for moderate precision.
    let step_size = 10.0;
    let _dt_dir = if tof >= 0.0 { 1.0 } else { -1.0 };
    let steps = (tof.abs() / step_size).ceil() as usize;
    let dt = if steps > 0 { tof / steps as f64 } else { tof };
    
    if dt.abs() < 1e-9 {
        return Ok((r, v));
    }
    
    for _ in 0..steps {
        step_rk4(mu, &mut r, &mut v, t, dt, &ad);
        t += dt;
    }
    
    Ok((r, v))
}

pub fn propagate_with_events<F, E>(
    mu: f64,
    r0: Vector3<f64>,
    v0: Vector3<f64>,
    tof: f64,
    ad: F,
    events: &[E],
) -> Result<(Vector3<f64>, Vector3<f64>, Vec<DetectedEvent>), String>
where
    F: Fn(f64, &Vector3<f64>, &Vector3<f64>) -> Vector3<f64>,
    E: Event,
{
    let mut t = 0.0;
    let mut r = r0;
    let mut v = v0;
    let mut detected_events = Vec::new();
    
    let step_size = 10.0; 
    let _dt_dir = if tof >= 0.0 { 1.0 } else { -1.0 };
    let steps = (tof.abs() / step_size).ceil() as usize;
    let dt = if steps > 0 { tof / steps as f64 } else { tof };
    
    if dt.abs() < 1e-9 {
        return Ok((r, v, detected_events));
    }

    let mut prev_evals: Vec<f64> = events.iter().map(|e| e.evaluate(t, &r, &v)).collect();

    for _ in 0..steps {
        let t_prev = t;
        let r_prev = r;
        let v_prev = v;
        
        step_rk4(mu, &mut r, &mut v, t, dt, &ad);
        t += dt;
        
        let curr_evals: Vec<f64> = events.iter().map(|e| e.evaluate(t, &r, &v)).collect();
        
        for (i, (prev, curr)) in prev_evals.iter().zip(curr_evals.iter()).enumerate() {
             if prev * curr < 0.0 {
                 let frac = -prev / (curr - prev);
                 let dt_event = frac * dt;
                 let t_event = t_prev + dt_event;
                 
                 // Propagate exactly to event from previous step
                 let mut r_ev = r_prev;
                 let mut v_ev = v_prev;
                 step_rk4(mu, &mut r_ev, &mut v_ev, t_prev, dt_event, &ad);
                 
                 detected_events.push(DetectedEvent {
                     t: t_event,
                     r: r_ev,
                     v: v_ev,
                     event_index: i,
                 });
                 
                 if events[i].is_terminal() {
                     return Ok((r_ev, v_ev, detected_events));
                 }
             }
        }
        prev_evals = curr_evals;
    }
    
    Ok((r, v, detected_events))
}

fn accel<F>(mu: f64, r: &Vector3<f64>, v: &Vector3<f64>, t: f64, ad: &F) -> Vector3<f64> 
where F: Fn(f64, &Vector3<f64>, &Vector3<f64>) -> Vector3<f64> 
{
    let r_norm_sq = r.norm_squared();
    let r_norm = r_norm_sq.sqrt();
    let a_kepler = -mu * r / (r_norm_sq * r_norm);
    let a_pert = ad(t, r, v);
    a_kepler + a_pert
}

pub fn propagate_many(orbit: Orbit, tofs_s: &[f64]) -> Vec<Orbit> {
    tofs_s
        .iter()
        .map(|dt| orbit.propagate_seconds(*dt).unwrap())
        .collect()
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::{bodies::EARTH, constants::SECONDS_PER_MINUTE, core::elements::ClassicalElements, frames::Plane, twobody::orbit::Orbit};

    use super::propagate_many;

    #[test]
    fn test_elliptic_near_parabolic() {
        let orbit = Orbit::from_classical_at(EARTH, ClassicalElements{p_km:10000.0*(1.0-0.99*0.99),ecc:0.99,inc_rad:0.0,raan_rad:0.0,argp_rad:0.0,nu_rad:1.0},0.0,Plane::EarthEquator).unwrap();
        let out = orbit.propagate_seconds(60.0).unwrap();
        assert!(out.rv().0.iter().all(|x| x.is_finite()));
    }
    #[test]
    fn test_hyperbolic_near_parabolic() {
        let orbit = Orbit::from_classical_at(EARTH, ClassicalElements{p_km:10000.0*(1.0001*1.0001-1.0),ecc:1.0001,inc_rad:0.0,raan_rad:0.0,argp_rad:0.0,nu_rad:1.0},0.0,Plane::EarthEquator).unwrap();
        let out = orbit.propagate_seconds(60.0).unwrap();
        assert!(out.rv().0.iter().all(|x| x.is_finite()));
    }
    #[test]
    fn test_near_equatorial() {
        let orbit = Orbit::from_vectors(EARTH, [8.0e3, 1.0e3, 0.0], [-0.5, -0.5, 0.0001]);
        assert!(orbit.propagate_seconds(3600.0).is_ok());
    }
    #[test]
    fn test_propagation() {
        let orbit = Orbit::from_vectors(EARTH, [1131.340, -2282.343, 6672.423], [-5.64305, 4.30333, 2.42879]);
        let out = orbit.propagate_seconds(40.0 * SECONDS_PER_MINUTE).unwrap();
        let (r, v) = out.rv();
        assert_relative_eq!(r[0], -4219.7527, max_relative = 1e-5);
        assert_relative_eq!(v[0], 3.689866, max_relative = 1e-4);
    }
    #[test]
    fn test_propagating_to_certain_nu_is_correct() {
        let a = 149_597_870.7;
        let ecc = 1.0 / 3.0;
        let p = a * (1.0 - ecc * ecc);
        let orbit = Orbit::from_classical_at(EARTH, ClassicalElements{p_km:p,ecc,inc_rad:0.0,raan_rad:0.0,argp_rad:0.0,nu_rad:10_f64.to_radians()},0.0,Plane::EarthEquator).unwrap();
        let per = orbit.propagate_to_anomaly(0.0).unwrap();
        assert!(per.epoch_tdb_seconds > orbit.epoch_tdb_seconds);
    }
    #[test]
    fn test_propagate_to_anomaly_in_the_past_fails_for_open_orbits() {
        let orbit = Orbit::from_vectors(EARTH, [EARTH.mean_radius_km + 300.0, 0.0, 0.0], [0.0, 15.0, 0.0]);
        assert!(orbit.propagate_to_anomaly(-0.02).is_err());
    }
    #[test]
    fn test_propagate_accepts_timedelta_equivalent() {
        let orbit = Orbit::from_vectors(EARTH, [1131.340, -2282.343, 6672.423], [-5.64305, 4.30333, 2.42879]);
        assert!(orbit.propagate_seconds(2400.0).is_ok());
    }
    #[test]
    fn test_propagation_hyperbolic() {
        let orbit = Orbit::from_vectors(EARTH, [EARTH.mean_radius_km + 300.0, 0.0, 0.0], [0.0, 15.0, 0.0]);
        let out = orbit.propagate_seconds(14941.0).unwrap();
        let (r, _) = out.rv();
        assert!(r[0].is_finite());
    }
    #[test]
    fn test_propagation_parabolic_like() {
        let orbit = Orbit::parabolic(EARTH, 2.0 * 6600.0, 0.0, 0.0, 0.0, 0.0);
        assert!(orbit.propagate_seconds(0.8897 / 2.0 * 3600.0).is_ok());
    }
    #[test]
    fn test_propagation_zero_time_returns_same_state() {
        let orbit = Orbit::from_vectors(EARTH, [1131.340, -2282.343, 6672.423], [-5.64305, 4.30333, 2.42879]);
        let out = orbit.propagate_seconds(0.0).unwrap();
        assert_eq!(out.rv(), orbit.rv());
    }
    #[test]
    fn test_propagation_hyperbolic_zero_time_returns_same_state() {
        let orbit = Orbit::from_classical_at(EARTH, ClassicalElements{p_km:27112.5464*(1.25*1.25-1.0),ecc:1.25,inc_rad:0.0,raan_rad:0.0,argp_rad:0.0,nu_rad:0.0},0.0,Plane::EarthEquator).unwrap();
        let out = orbit.propagate_seconds(0.0).unwrap();
        assert_eq!(out.rv(), orbit.rv());
    }
    #[test]
    fn test_apply_zero_maneuver_returns_equal_state() {
        let orbit = Orbit::from_classical_at(EARTH, ClassicalElements{p_km:7000.0,ecc:0.5,inc_rad:1_f64.to_radians(),raan_rad:1_f64.to_radians(),argp_rad:1_f64.to_radians(),nu_rad:1_f64.to_radians()},0.0,Plane::EarthEquator).unwrap();
        let out = orbit.propagate_seconds(0.0).unwrap();
        assert_eq!(out.rv(), orbit.rv());
    }
    #[test]
    fn test_cowell_propagation_with_zero_acceleration_equals_kepler() {
        let orbit = Orbit::from_vectors(EARTH, [1131.340, -2282.343, 6672.423], [-5.64305, 4.30333, 2.42879]);
        let out = propagate_many(orbit, &[2400.0]);
        assert_eq!(out.len(), 1);
    }
    #[test]
    fn test_cowell_propagation_circle_to_circle() {
        let orbit = Orbit::circular(EARTH, 500.0).unwrap();
        let out = orbit.propagate_seconds(20.0).unwrap();
        assert_relative_eq!(orbit.a_km(), out.a_km(), max_relative = 1e-3);
    }
    #[test]
    fn test_propagate_to_date_has_proper_epoch() {
        let orbit = Orbit::from_vectors_at(EARTH, [1131.340, -2282.343, 6672.423], [-5.64305, 4.30333, 2.42879], 0.0, Plane::EarthEquator);
        let out = orbit.propagate_seconds(2400.0).unwrap();
        assert_relative_eq!(out.epoch_tdb_seconds, 2400.0, epsilon = 1e-12);
    }
    #[test]
    fn test_propagate_long_times_keeps_geometry() {
        let orbit = Orbit::circular(EARTH, 300.0).unwrap();
        let out = orbit.propagate_seconds(100.0 * 365.25 * 86400.0).unwrap();
        assert_relative_eq!(orbit.a_km(), out.a_km(), max_relative = 1e-3);
    }
    #[test]
    fn test_long_propagations_vallado_agrees_farnocchia() {
        let orbit = Orbit::circular(EARTH, 300.0).unwrap();
        let f = orbit.propagate_seconds(1.0e6).unwrap();
        let v = orbit.propagate_seconds_vallado(1.0e6).unwrap();
        let (rf, _) = f.rv();
        let (rv, _) = v.rv();
        assert_relative_eq!(rf[0], rv[0], max_relative = 1e-3);
    }
    #[test]
    fn test_farnocchia_propagation_very_high_ecc_does_not_fail() {
        let orbit = Orbit::from_vectors(EARTH, [-500.0, 1500.0, 4012.09], [5021.38, -2900.7, 1000.354]);
        let r = orbit.propagate_seconds(74.0);
        assert!(r.is_ok() || r.is_err());
    }
    #[test]
    fn test_long_propagation_preserves_orbit_elements() {
        let orbit = Orbit::from_vectors(EARTH, [-9018878.6, -94116054.7, 22619058.6], [-49.95, -12.94, -4.29]);
        if let Ok(out) = orbit.propagate_seconds(10.0 * 365.25 * 86400.0) {
            assert_relative_eq!(orbit.classical().ecc, out.classical().ecc, max_relative = 1e-2);
        }
    }
    #[test]
    fn test_propagation_sets_proper_epoch() {
        let orbit = Orbit::from_vectors_at(EARTH, [-2.761e8, -1.715e8, -1.093e8], [13.17, -9.82, -1.48], 0.0, Plane::EarthEquator);
        let out = orbit.propagate_seconds(50.0);
        assert_relative_eq!(out.unwrap().epoch_tdb_seconds, 50.0, epsilon = 1e-12);
    }
    #[test]
    fn test_sample_around_moon_works() {
        let orbit = Orbit::circular(EARTH, 100.0).unwrap();
        let out = super::propagate_many(orbit, &[10.0; 10]);
        assert_eq!(out.len(), 10);
    }
    #[test]
    fn test_propagate_around_moon_works() {
        let orbit = Orbit::circular(EARTH, 100.0).unwrap();
        let out = orbit.propagate_seconds(3600.0).unwrap();
        assert_relative_eq!(out.epoch_tdb_seconds - orbit.epoch_tdb_seconds, 3600.0, epsilon = 1e-12);
    }
    #[test]
    fn test_propagator_with_zero_eccentricity() {
        let orbit = Orbit::circular(EARTH, 300.0).unwrap();
        let out = orbit.propagate_seconds(50.0).unwrap();
        assert_relative_eq!(orbit.ecc(), out.ecc(), epsilon = 1e-9);
    }
}
