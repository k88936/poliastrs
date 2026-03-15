use nalgebra::Vector3;

use crate::bodies::Body;
use crate::core::elements::{ClassicalElements, EquinoctialElements, coe2mee};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CartesianState {
    pub r_km: Vector3<f64>,
    pub v_km_s: Vector3<f64>,
}

impl CartesianState {
    pub fn new(r_km: Vector3<f64>, v_km_s: Vector3<f64>) -> Self {
        Self { r_km, v_km_s }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClassicalState {
    pub attractor: Body,
    pub elements: ClassicalElements,
}

impl ClassicalState {
    pub fn new(attractor: Body, elements: ClassicalElements) -> Self {
        Self {
            attractor,
            elements,
        }
    }

    pub fn n(&self) -> f64 {
        let a = self.elements.p_km / (1.0 - self.elements.ecc * self.elements.ecc);
        (self.attractor.mu_km3_s2 / (a * a * a)).sqrt()
    }

    pub fn to_equinoctial(&self) -> Result<EquinoctialElements, &'static str> {
        coe2mee(self.elements)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RVState {
    pub attractor: Body,
    pub cartesian: CartesianState,
}

impl RVState {
    pub fn new(attractor: Body, cartesian: CartesianState) -> Self {
        Self {
            attractor,
            cartesian,
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    use crate::bodies::{EARTH, SUN};
    use crate::core::elements::ClassicalElements;

    use super::{CartesianState, ClassicalState, RVState};

    #[test]
    fn state_has_attractor_given_in_constructor() {
        let ss = ClassicalState::new(
            SUN,
            ClassicalElements {
                p_km: 149_597_870.7,
                ecc: 0.5,
                inc_rad: 1_f64.to_radians(),
                raan_rad: 1_f64.to_radians(),
                argp_rad: 1_f64.to_radians(),
                nu_rad: 1_f64.to_radians(),
            },
        );
        assert_eq!(ss.attractor, SUN);
    }

    #[test]
    fn classical_state_has_elements_given_in_constructor() {
        let p = 227_939_200.0 * (1.0 - 0.093315_f64.powi(2));
        let ss = ClassicalState::new(
            SUN,
            ClassicalElements {
                p_km: p,
                ecc: 0.093315,
                inc_rad: 1.85_f64.to_radians(),
                raan_rad: 49.562_f64.to_radians(),
                argp_rad: 286.537_f64.to_radians(),
                nu_rad: 23.33_f64.to_radians(),
            },
        );
        assert_relative_eq!(ss.elements.p_km, p, max_relative = 1e-12);
        assert_relative_eq!(ss.elements.ecc, 0.093315, max_relative = 1e-12);
    }

    #[test]
    fn rv_state_has_rv_given_in_constructor() {
        let ss = RVState::new(
            SUN,
            CartesianState::new(
                Vector3::new(149_597_870.7, 0.0, 0.0),
                Vector3::new(0.0, 0.001, 0.0),
            ),
        );
        assert_eq!(ss.cartesian.r_km.x, 149_597_870.7);
        assert_eq!(ss.cartesian.v_km_s.y, 0.001);
    }

    #[test]
    fn mean_motion() {
        let ss = ClassicalState::new(
            EARTH,
            ClassicalElements {
                p_km: 42_164.1696,
                ecc: 0.0,
                inc_rad: 1.85_f64.to_radians(),
                raan_rad: 50_f64.to_radians(),
                argp_rad: 200_f64.to_radians(),
                nu_rad: 20_f64.to_radians(),
            },
        );
        let expected = 2.0 * std::f64::consts::PI / 86164.090518;
        assert_relative_eq!(ss.n(), expected, max_relative = 1e-8);
    }

    #[test]
    fn to_equinoctial_retrograde_equatorial_raises_error() {
        let ss = ClassicalState::new(
            SUN,
            ClassicalElements {
                p_km: 10000.0 * (1.0 - 0.3 * 0.3),
                ecc: 0.3,
                inc_rad: std::f64::consts::PI,
                raan_rad: 49.562_f64.to_radians(),
                argp_rad: 286.537_f64.to_radians(),
                nu_rad: 23.33_f64.to_radians(),
            },
        );
        let err = ss.to_equinoctial().unwrap_err();
        assert_eq!(
            err,
            "Cannot compute modified equinoctial set for 180 degrees orbit inclination due to `h` and `k` singularity."
        );
    }
}
