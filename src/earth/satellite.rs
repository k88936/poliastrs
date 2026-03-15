use crate::twobody::orbit::Orbit;
use crate::spacecraft::Spacecraft;
use crate::earth::enums::EarthGravity;
use crate::bodies::EARTH;
use crate::twobody::propagation::cowell;
use crate::earth::atmosphere::Atmosphere;
use crate::twobody::perturbations::j2_accel;
use nalgebra::Vector3;

const J2_EARTH: f64 = 1.08263e-3;
const R_EARTH_EQ_KM: f64 = 6378.137;

#[derive(Clone, Debug)]
pub struct EarthSatellite {
    pub orbit: Orbit,
    pub spacecraft: Option<Spacecraft>,
}

impl EarthSatellite {
    pub fn new(orbit: Orbit, spacecraft: Option<Spacecraft>) -> Result<Self, String> {
        if (orbit.attractor.mu_km3_s2 - EARTH.mu_km3_s2).abs() > 1e-6 {
             return Err("The attractor must be Earth".to_string());
        }
        Ok(Self { orbit, spacecraft })
    }
    
    pub fn propagate(&self, tof_s: f64, gravity: Option<EarthGravity>, atmosphere: Option<&dyn Atmosphere>) -> Result<Self, String> {
        let r0 = self.orbit.state.r_km;
        let v0 = self.orbit.state.v_km_s;
        let mu = self.orbit.attractor.mu_km3_s2;
        
        let j2_val = if let Some(EarthGravity::J2) = gravity {
            Some((J2_EARTH, R_EARTH_EQ_KM))
        } else {
            None
        };
        
        // Prepare drag params if spacecraft and atmosphere exist
        let drag_params = if let (Some(sc), Some(atm)) = (self.spacecraft, atmosphere) {
            Some((sc.ballistic_coefficient(), atm))
        } else {
            None
        };

        if j2_val.is_none() && drag_params.is_none() {
             // Keplerian propagation
             let new_orbit = self.orbit.propagate_seconds(tof_s).map_err(|e| format!("{:?}", e))?;
             return Ok(Self { orbit: new_orbit, spacecraft: self.spacecraft });
        }
        
        // Cowell propagation
        let ad = move |_t: f64, r: &Vector3<f64>, v: &Vector3<f64>| -> Vector3<f64> {
            let mut acc = Vector3::new(0.0, 0.0, 0.0);
            let r_arr = [r.x, r.y, r.z];
            
            if let Some((j2, r_eq)) = j2_val {
                let acc_j2 = j2_accel(r_arr, mu, j2, r_eq);
                acc += Vector3::new(acc_j2[0], acc_j2[1], acc_j2[2]);
            }
            
            if let Some((b_star, atm)) = drag_params {
                let r_norm = r.norm();
                let v_norm = v.norm();
                // Use R_EARTH_EQ_KM for geometric altitude
                let z = r_norm - R_EARTH_EQ_KM;
                
                if z > 0.0 && z <= 1000.0 {
                     // Clamp to safe range to avoid floating point issues
                     let safe_z = z.clamp(0.0, 1000.0);
                     let rho_kg_m3 = atm.density(safe_z);
                     let rho_kg_km3 = rho_kg_m3 * 1e9;
                     
                     // Drag acceleration: -1/2 * rho * v * B* * v_vec
                     let drag_acc = -0.5 * rho_kg_km3 * v_norm * b_star * v;
                     acc += drag_acc;
                }
            }
            
            acc
        };
        
        let (r_new, v_new) = cowell::propagate(mu, r0, v0, tof_s, ad)?;
                     
        let new_orb = Orbit::from_vectors_at(
           self.orbit.attractor, 
           [r_new.x, r_new.y, r_new.z], 
           [v_new.x, v_new.y, v_new.z], 
           self.orbit.epoch_tdb_seconds + tof_s, 
           self.orbit.plane
       );
        Ok(Self { orbit: new_orb, spacecraft: self.spacecraft })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::MARS;
    use crate::frames::Plane;
    use std::f64::consts::PI;
    use crate::earth::atmosphere::COESA76;
    
    #[test]
    fn test_earth_satellite_orbit() {
        let r = [3539.088, 5310.199, 3066.313];
        let v = [-6.497, 3.249, 1.875];
        let orb = Orbit::from_vectors_at(EARTH, r, v, 0.0, Plane::EarthEquator);
        
        let sc = Spacecraft::new(PI/4.0 * 1e-6, 2.2, 100.0);
        let sat = EarthSatellite::new(orb.clone(), Some(sc)).unwrap();
        
        assert_eq!(sat.orbit.attractor.mu_km3_s2, EARTH.mu_km3_s2);
    }
    
    #[test]
    fn test_orbit_attractor() {
        let r = [3539.088, 5310.199, 3066.313];
        let v = [-6.497, 3.249, 1.875];
        let orb = Orbit::from_vectors_at(MARS, r, v, 0.0, Plane::EarthEquator);
        
        let sc = Spacecraft::new(PI/4.0 * 1e-6, 2.2, 100.0);
        let res = EarthSatellite::new(orb, Some(sc));
        assert!(res.is_err());
        assert_eq!(res.unwrap_err(), "The attractor must be Earth");
    }
    
    #[test]
    fn test_propagate_instance() {
        // Setup initial orbit
        let r = [7000.0, 0.0, 0.0];
        let v = [0.0, 7.5, 0.0];
        let orb = Orbit::from_vectors_at(EARTH, r, v, 0.0, Plane::EarthEquator);
        let sc = Spacecraft::new(1e-6, 2.2, 100.0);
        let sat = EarthSatellite::new(orb, Some(sc)).unwrap();
        
        let tof = 60.0;
        
        // Propagate J2
        let res_j2 = sat.propagate(tof, Some(EarthGravity::J2), None);
        assert!(res_j2.is_ok());
        
        // Propagate Kepler
        let res_kep = sat.propagate(tof, None, None);
        assert!(res_kep.is_ok());

        // Propagate J2 + Atmosphere
        let coesa76 = COESA76::new();
        let res_atm = sat.propagate(tof, Some(EarthGravity::J2), Some(&coesa76 as &dyn crate::earth::atmosphere::Atmosphere));
        assert!(res_atm.is_ok());
    }
}
