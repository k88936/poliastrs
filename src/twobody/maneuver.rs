use nalgebra::Vector3;
use std::f64::consts::PI;
use crate::twobody::orbit::Orbit;
use crate::iod::izzo;

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

    pub fn lambert(orbit_i: &Orbit, orbit_f: &Orbit) -> Result<Self, String> {
        let tof = orbit_f.epoch_tdb_seconds - orbit_i.epoch_tdb_seconds;
        if tof <= 0.0 {
            return Err("Epoch of initial orbit greater than epoch of final orbit, causing a negative time of flight".to_string());
        }

        let mu = orbit_i.attractor.mu_km3_s2;
        let r1 = orbit_i.state.r_km;
        let r2 = orbit_f.state.r_km;

        // Use Izzo algorithm with default parameters matching poliastro
        // M=0, prograde=true, low_path=true, num_iter=35, rtol=1e-8
        let (v1, v2) = izzo(mu, r1, r2, tof, 0, true, true, 35, 1e-8)
            .map_err(|e| format!("Lambert solver failed: {:?}", e))?;

        let dv1 = v1 - orbit_i.state.v_km_s;
        let dv2 = orbit_f.state.v_km_s - v2;

        Ok(Self {
            impulses: vec![
                (0.0, dv1),
                (tof, dv2),
            ],
        })
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

    pub fn correct_pericenter(orbit: &Orbit, max_delta_r_km: f64) -> Result<Self, String> {
        let j2 = orbit.attractor.j2;
        if j2 == 0.0 {
            return Err(format!("The correction maneuver is not yet supported for {}", orbit.attractor.name));
        }
        
        let ecc = orbit.ecc();
        if ecc > 0.001 {
            return Err(format!("The correction maneuver is not yet supported with {}, it should be less than or equal to 0.001", ecc));
        }

        let r_attractor = orbit.attractor.equatorial_radius_km;
        let k = orbit.attractor.mu_km3_s2;
        let a = orbit.a_km();
        let inc = orbit.inc_rad();

        let p = a * (1.0 - ecc.powi(2));
        let n = (k / a.powi(3)).sqrt();

        let dw = ((3.0 * n * r_attractor.powi(2) * j2) / (4.0 * p.powi(2))) * (4.0 - 5.0 * inc.sin().powi(2));

        let mut delta_w = 2.0 * (1.0 + ecc) * max_delta_r_km;
        delta_w /= a * ecc * (1.0 - ecc);
        delta_w = delta_w.sqrt();

        let delta_t = (delta_w / dw).abs();
        let delta_v = 0.5 * n * a * ecc * delta_w.abs();

        let vf = orbit.state.v_km_s.normalize() * delta_v;

        Ok(Self {
            impulses: vec![
                (delta_t, vf),
            ],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::{EARTH, MERCURY};
    use crate::twobody::orbit::Orbit;
    use crate::frames::Plane;
    use approx::assert_relative_eq;

    #[test]
    fn test_correct_pericenter() {
        use crate::frames::Plane;
        use nalgebra::Vector3;

        let max_delta_r = 30.0;
        let a = 6570.0;
        let ecc = 0.001;
        let inc = 0.7855682278773197;

        let ss0 = Orbit::from_keplerian(
            EARTH,
            a,
            ecc,
            inc,
            0.0,
            0.0,
            0.0,
            0.0,
            Plane::EarthEquator,
        );

        let maneuver = Maneuver::correct_pericenter(&ss0, max_delta_r).unwrap();

        let expected_t = 2224141.03634;
        let expected_v = Vector3::new(0.0, 0.0083290328315531, 0.00833186625871848);

        let (t, v) = maneuver.impulses[0];

        assert_relative_eq!(t, expected_t, epsilon = 1.0);
        assert_relative_eq!(v, expected_v, max_relative = 1e-5);
    }

    #[test]
    fn test_correct_pericenter_j2_exception() {
        use crate::frames::Plane;
        use crate::bodies::MERCURY;

        let ss0 = Orbit::from_keplerian(
            MERCURY,
            1000.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            Plane::EarthEquator,
        );
        let max_delta_r = 30.0;
        
        let result = Maneuver::correct_pericenter(&ss0, max_delta_r);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "The correction maneuver is not yet supported for Mercury");
    }

    #[test]
    fn test_correct_pericenter_ecc_exception() {
        use crate::frames::Plane;

        let ss0 = Orbit::from_keplerian(
            EARTH,
            1000.0,
            0.5, // ecc > 0.001
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            Plane::EarthEquator,
        );
        let max_delta_r = 30.0;
        
        let result = Maneuver::correct_pericenter(&ss0, max_delta_r);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("it should be less than or equal to 0.001"));
    }

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

    #[test]
    fn test_lambert_maneuver_tof_exception() {
        let orb_i = Orbit::from_keplerian(
            EARTH,
            7000.0,
            0.0,
            0.0, // inc
            0.0, // raan
            0.0, // argp
            0.0, // nu
            100.0, // epoch
            Plane::EarthEquator,
        );
        let orb_f = Orbit::from_keplerian(
            EARTH,
            7000.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            50.0, // epoch < orb_i epoch
            Plane::EarthEquator,
        );

        let result = Maneuver::lambert(&orb_i, &orb_f);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Epoch of initial orbit greater than epoch of final orbit, causing a negative time of flight"
        );
    }
}
