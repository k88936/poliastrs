use nalgebra::Vector3;
use std::f64::consts::PI;
use crate::twobody::orbit::Orbit;

#[derive(Debug, Clone)]
pub struct Maneuver {
    pub impulses: Vec<(f64, Vector3<f64>)>, // (dt_seconds, dv_km_s)
}

impl Maneuver {
    pub fn new(impulses: Vec<(f64, Vector3<f64>)>) -> Self {
        Self { impulses }
    }

    pub fn impulse(dv: Vector3<f64>) -> Self {
        Self {
            impulses: vec![(0.0, dv)],
        }
    }

    pub fn get_total_cost(&self) -> f64 {
        self.impulses.iter().map(|(_, dv)| dv.norm()).sum()
    }

    pub fn get_total_time(&self) -> f64 {
        self.impulses.iter().map(|(dt, _)| *dt).sum()
    }

    pub fn hohmann(orbit_i: &Orbit, r_f: f64) -> Self {
        let mu = orbit_i.attractor.mu_km3_s2;
        let r_i = orbit_i.state.r_km.norm();
        let v_i = orbit_i.state.v_km_s.norm();

        let a_trans = (r_i + r_f) / 2.0;

        // First impulse
        // Assumption: orbit_i is circular or we are at periapsis/apoapsis suitable for transfer
        let v_trans_a = (mu * (2.0 / r_i - 1.0 / a_trans)).sqrt();
        let dv_a_mag = v_trans_a - v_i;
        let dv_a = orbit_i.state.v_km_s.normalize() * dv_a_mag;

        // Time of flight
        let t_trans = PI * (a_trans.powi(3) / mu).sqrt();

        // Second impulse
        let v_trans_b = (mu * (2.0 / r_f - 1.0 / a_trans)).sqrt();
        let v_f_circular = (mu / r_f).sqrt();
        let dv_b_mag = v_f_circular - v_trans_b;

        // Velocity at arrival is opposite to initial velocity direction (180 deg transfer)
        let dv_b = -orbit_i.state.v_km_s.normalize() * dv_b_mag;

        Self {
            impulses: vec![
                (0.0, dv_a),
                (t_trans, dv_b),
            ],
        }
    }

    pub fn bielliptic(orbit_i: &Orbit, r_b: f64, r_f: f64) -> Self {
        let mu = orbit_i.attractor.mu_km3_s2;
        let r_i = orbit_i.state.r_km.norm();
        let v_i = orbit_i.state.v_km_s.norm();

        // Transfer 1: r_i -> r_b
        let a1 = (r_i + r_b) / 2.0;
        let v1_a = (mu * (2.0 / r_i - 1.0 / a1)).sqrt();
        let dv_a_mag = v1_a - v_i;
        let dv_a = orbit_i.state.v_km_s.normalize() * dv_a_mag;
        let t1 = PI * (a1.powi(3) / mu).sqrt();

        // Transfer 2: r_b -> r_f
        let a2 = (r_b + r_f) / 2.0;
        let v1_b = (mu * (2.0 / r_b - 1.0 / a1)).sqrt(); 
        let v2_b = (mu * (2.0 / r_b - 1.0 / a2)).sqrt(); 

        // At r_b, velocity is opposite to v_i.
        let dv_b_mag = v2_b - v1_b;
        let dv_b = -orbit_i.state.v_km_s.normalize() * dv_b_mag;
        let t2 = PI * (a2.powi(3) / mu).sqrt();

        // Arrival at r_f
        // Velocity direction is parallel to v_i
        let v2_f = (mu * (2.0 / r_f - 1.0 / a2)).sqrt();
        let v_f_circular = (mu / r_f).sqrt();
        let dv_c_mag = v_f_circular - v2_f;
        let dv_c = orbit_i.state.v_km_s.normalize() * dv_c_mag;

        Self {
            impulses: vec![
                (0.0, dv_a),
                (t1, dv_b),
                (t2, dv_c),
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::EARTH;
    use crate::twobody::orbit::Orbit;
    use crate::frames::Plane;
    use approx::assert_relative_eq;

    #[test]
    fn test_hohmann_maneuver() {
        let r_i = 6378.137 + 191.34411; 
        let r_f = 6378.137 + 35781.34857; 
        
        let orb_i = Orbit::from_keplerian(
            EARTH,
            r_i, 
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            Plane::EarthEquator,
        );
        
        let man = Maneuver::hohmann(&orb_i, r_f);
        
        let expected_dv = 3.935224;
        let expected_t_trans = 5.256713 * 3600.0; 
        
        assert_relative_eq!(man.get_total_cost(), expected_dv, epsilon = 1e-5);
        assert_relative_eq!(man.get_total_time(), expected_t_trans, epsilon = 1.0); 
        
        let orb_f = orb_i.apply_maneuver(&man);
        assert_relative_eq!(orb_f.ecc(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(orb_f.a_km(), r_f, epsilon = 1e-3);
    }

    #[test]
    fn test_bielliptic_maneuver() {
        let r_i = 6378.137 + 191.34411;
        let r_b = 6378.137 + 503873.0;
        let r_f = 6378.137 + 376310.0;
        
        let orb_i = Orbit::from_keplerian(
            EARTH,
            r_i,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            Plane::EarthEquator,
        );
        
        let man = Maneuver::bielliptic(&orb_i, r_b, r_f);
        
        let expected_dv = 3.904057;
        let expected_t_trans = 593.919803 * 3600.0;
        
        assert_relative_eq!(man.get_total_cost(), expected_dv, epsilon = 1e-5);
        assert_relative_eq!(man.get_total_time(), expected_t_trans, epsilon = 10.0);
        
        let orb_f = orb_i.apply_maneuver(&man);
        assert_relative_eq!(orb_f.ecc(), 0.0, epsilon = 1e-12);
        // assert_relative_eq!(orb_f.a_km(), r_f, epsilon = 1e-3);
    }
}
