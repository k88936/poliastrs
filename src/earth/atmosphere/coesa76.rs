use crate::earth::atmosphere::util::{check_altitude, get_index};
use crate::earth::atmosphere::Atmosphere;

pub const R_AIR: f64 = 287.053; // J / kg / K
pub const GAMMA: f64 = 1.4;
pub const BETA: f64 = 1.458e-6; // kg / s / m / K^0.5
pub const S: f64 = 110.4; // K
pub const R0: f64 = 6356.766; // km
pub const T_INF: f64 = 1000.0; // K

// Data tables
const ZB_LEVELS: [f64; 13] = [0.0, 11.019, 20.063, 32.162, 47.350, 51.413, 71.802, 86.0, 91.0, 110.0, 120.0, 500.0, 1000.0];
const HB_LEVELS: [f64; 13] = [0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852, 89.716, 108.129, 117.777, 463.340, 864.071];
const TB_LEVELS: [f64; 13] = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.95, 187.36, 254.93, 397.91, 2019.69, 7351.15];
const LB_LEVELS: [f64; 13] = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0, 0.0, 12.0, 0.0, 0.0, 0.0];
const PB_LEVELS: [f64; 13] = [1013.25, 226.32, 54.748, 8.6801, 1.109, 0.66938, 0.039564, 0.0037338, 0.0015381, 7.1042e-05, 2.5382e-05, 3.0236e-09, 7.5138e-11]; // mbar

const Z_COEFF: [f64; 11] = [86.0, 91.0, 100.0, 110.0, 120.0, 150.0, 200.0, 300.0, 500.0, 750.0, 1000.0];

const P_COEFF: [[f64; 5]; 11] = [
    [0.0, 2.159582e-06, -0.0004836957, -0.1425192, 13.4753],
    [0.0, 3.304895e-05, -0.00906273, 0.6516698, -11.03037],
    [0.0, 6.693926e-05, -0.01945388, 1.71908, -47.7503],
    [0.0, -6.539316e-05, 0.02485568, -3.22362, 135.9355],
    [2.283506e-07, -0.0001343221, 0.02999016, -3.055446, 113.5764],
    [1.209434e-08, -9.692458e-06, 0.003002041, -0.4523015, 19.19151],
    [8.113942e-10, -9.822568e-07, 0.0004687616, -0.123171, 3.067409],
    [9.814674e-11, -1.654439e-07, 0.0001148115, -0.05431334, -2.011365],
    [-7.835161e-11, 1.964589e-07, -0.0001657213, 0.04305869, -14.77132],
    [2.813255e-11, -1.120689e-07, 0.0001695568, -0.1188941, 14.56718],
    [2.813255e-11, -1.120689e-07, 0.0001695568, -0.1188941, 14.56718], // Duplicate for upper bound check safety
];

const RHO_COEFF: [[f64; 5]; 11] = [
    [0.0, -3.322622e-06, 0.000911146, -0.2609971, 5.944694],
    [0.0, 2.873405e-05, -0.008492037, 0.6541179, -23.6201],
    [-1.240774e-05, 0.005162063, -0.8048342, 55.55996, -1443.338],
    [0.0, -8.854164e-05, 0.03373254, -4.390837, 176.5294],
    [3.661771e-07, -0.0002154344, 0.04809214, -4.884744, 172.3597],
    [1.906032e-08, -1.527799e-05, 0.004724294, -0.699234, 20.50921],
    [1.199282e-09, -1.451051e-06, 0.0006910474, -0.173622, -5.321644],
    [1.140564e-10, -2.130756e-07, 0.0001570762, -0.07029296, -12.89844],
    [8.105631e-12, -2.358417e-09, -2.63511e-06, -0.01562608, -20.02246],
    [-3.701195e-12, -8.608611e-09, 5.118829e-05, -0.06600998, -6.137674],
    [-3.701195e-12, -8.608611e-09, 5.118829e-05, -0.06600998, -6.137674],
];

#[derive(Clone, Copy, Debug)]
pub struct COESA76;

impl COESA76 {
    pub fn new() -> Self {
        Self
    }

    fn get_coefficients_above_86(z: f64, coeffs: &[[f64; 5]]) -> [f64; 5] {
        // Find index in Z_COEFF
        // Since we are above 86km, we can assume z >= 86.0
        // But get_index handles general case.
        // Z_COEFF = [86, 91, 100, ...]
        // If z=90, index=0 (86).
        // If z=86, index=0 (86).
        // If z=1000, index=10?
        // Let's use get_index logic.
        
        let mut i = 0;
        for (idx, &val) in Z_COEFF.iter().enumerate() {
            if val > z {
                if idx == 0 { i = 0; break; }
                i = idx - 1;
                break;
            }
            i = idx;
        }
        
        // Handle upper bound
        if i >= coeffs.len() {
            i = coeffs.len() - 1;
        }

        coeffs[i]
    }

    pub fn temperature(&self, alt: f64) -> f64 {
        let (z, h) = check_altitude(alt, R0, true, 0.0, 1000.0).unwrap();
        
        let i = get_index(z, &ZB_LEVELS);
        let tb = TB_LEVELS[i];
        let lb = LB_LEVELS[i];
        let hb = HB_LEVELS[i];

        if z < ZB_LEVELS[7] {
            // Below 86 km
            let tm = tb + lb * (h - hb);
            tm
        } else if z >= ZB_LEVELS[7] && z < ZB_LEVELS[8] {
            // [86, 91)
            186.87
        } else if z >= ZB_LEVELS[8] && z < ZB_LEVELS[9] {
            // [91, 110)
            let tc = 263.1905;
            let a = -76.3232;
            let a_coeff = -19.9429;
            tc + a * (1.0 - ((z - ZB_LEVELS[8]) / a_coeff).powi(2)).sqrt()
        } else if z >= ZB_LEVELS[9] && z < ZB_LEVELS[10] {
            // [110, 120)
            240.0 + lb * (z - ZB_LEVELS[9])
        } else {
             // >= 120
             let t10 = 360.0;
             let gamma_val = LB_LEVELS[9] / (T_INF - t10);
             let epsilon = (z - ZB_LEVELS[10]) * (R0 + ZB_LEVELS[10]) / (R0 + z);
             T_INF - (T_INF - t10) * (-gamma_val * epsilon).exp()
        }
    }

    pub fn pressure(&self, alt: f64) -> f64 {
        let (z, h) = check_altitude(alt, R0, true, 0.0, 1000.0).unwrap();

        let i = get_index(z, &ZB_LEVELS);
        let tb = TB_LEVELS[i];
        let lb = LB_LEVELS[i];
        let hb = HB_LEVELS[i];
        let pb = PB_LEVELS[i];

        if z < 86.0 {
            let p_mbar = if lb == 0.0 {
                pb * (-34.1632 * (h - hb) / tb).exp()
            } else {
                let t = self.temperature(z);
                pb * (tb / t).powf(34.1632 / lb)
            };
            p_mbar * 100.0
        } else {
             let coeffs = Self::get_coefficients_above_86(z, &P_COEFF);
             let a = coeffs[0];
             let b = coeffs[1];
             let c = coeffs[2];
             let d = coeffs[3];
             let e = coeffs[4];
             (a * z.powi(4) + b * z.powi(3) + c * z.powi(2) + d * z + e).exp()
        }
    }

    pub fn density(&self, alt: f64) -> f64 {
        let (z, _h) = check_altitude(alt, R0, true, 0.0, 1000.0).unwrap();

        if z <= 86.0 {
             let t = self.temperature(z);
             let p = self.pressure(z);
             p / (R_AIR * t)
        } else {
             let coeffs = Self::get_coefficients_above_86(z, &RHO_COEFF);
             let a = coeffs[0];
             let b = coeffs[1];
             let c = coeffs[2];
             let d = coeffs[3];
             let e = coeffs[4];
             (a * z.powi(4) + b * z.powi(3) + c * z.powi(2) + d * z + e).exp()
        }
    }
    
    pub fn properties(&self, alt: f64) -> (f64, f64, f64) {
        (self.temperature(alt), self.pressure(alt), self.density(alt))
    }

    pub fn sound_speed(&self, alt: f64) -> Result<f64, String> {
        let (z, _) = check_altitude(alt, R0, true, 0.0, 1000.0).unwrap();
        if z > 86.0 {
            return Err("Speed of sound in COESA76 has just been implemented up to 86km.".to_string());
        }
        let t = self.temperature(z);
        Ok((GAMMA * R_AIR * t).sqrt())
    }

    pub fn viscosity(&self, alt: f64) -> Result<f64, String> {
        let (z, _) = check_altitude(alt, R0, true, 0.0, 1000.0).unwrap();
        if z > 86.0 {
            return Err("Dynamic Viscosity in COESA76 has just been implemented up to 86km.".to_string());
        }
        let t = self.temperature(z);
        Ok(BETA * t.powf(1.5) / (t + S))
    }

    pub fn thermal_conductivity(&self, alt: f64) -> Result<f64, String> {
        let (z, _) = check_altitude(alt, R0, true, 0.0, 1000.0).unwrap();
        if z > 86.0 {
            return Err("Thermal conductivity in COESA76 has just been implemented up to 86km.".to_string());
        }
        let t = self.temperature(z);
        Ok((2.64638e-3 * t.powf(1.5)) / (t + 245.4 * 10.0_f64.powf(-12.0 / t)))
    }
}

impl Atmosphere for COESA76 {
    fn density(&self, alt_km: f64) -> f64 {
        self.density(alt_km)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_properties_coesa76_0_5km() {
        let coesa = COESA76::new();
        let (t, p, rho) = coesa.properties(0.5);
        // T=284.90 K, p=9.5461e2 mbar = 95461 Pa, rho=1.1673 kg/m3
        assert_relative_eq!(t, 284.90, max_relative = 1e-4);
        assert_relative_eq!(p, 95461.0, max_relative = 1e-3);
        assert_relative_eq!(rho, 1.1673, max_relative = 1e-3);
    }
    
    #[test]
    fn test_properties_coesa76_86km() {
        let coesa = COESA76::new();
        let (t, p, rho) = coesa.properties(86.0);
        // T=186.87 K, p=3.7338e-3 mbar = 0.37338 Pa, rho=6.958e-6 kg/m3
        assert_relative_eq!(t, 186.87, max_relative = 1e-4);
        assert_relative_eq!(p, 0.37338, max_relative = 1e-4);
        assert_relative_eq!(rho, 6.958e-6, max_relative = 1e-3);
    }

    #[test]
    fn test_sound_speed_0_5km() {
        let coesa = COESA76::new();
        let cs = coesa.sound_speed(0.5).unwrap();
        // 338.37 m/s
        assert_relative_eq!(cs, 338.37, epsilon = 1e-2);
    }

    #[test]
    #[should_panic(expected = "Geometric altitude must be in range")]
    fn test_outside_altitude_range() {
        let coesa = COESA76::new();
        let _ = coesa.temperature(1001.0);
    }
    
    #[test]
    fn test_sound_speed_over_86km() {
        let coesa = COESA76::new();
        let res = coesa.sound_speed(87.0);
        assert!(res.is_err());
    }
}
