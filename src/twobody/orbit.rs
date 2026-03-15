use nalgebra::Vector3;

use crate::bodies::Body;
use crate::twobody::propagation::farnocchia::{PropagationError, propagate_two_body};
use crate::twobody::states::CartesianState;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orbit {
    pub attractor: Body,
    pub state: CartesianState,
}

impl Orbit {
    pub fn from_vectors(attractor: Body, r_km: [f64; 3], v_km_s: [f64; 3]) -> Self {
        Self {
            attractor,
            state: CartesianState::new(
                Vector3::new(r_km[0], r_km[1], r_km[2]),
                Vector3::new(v_km_s[0], v_km_s[1], v_km_s[2]),
            ),
        }
    }

    pub fn rv(self) -> ([f64; 3], [f64; 3]) {
        (
            [self.state.r_km.x, self.state.r_km.y, self.state.r_km.z],
            [
                self.state.v_km_s.x,
                self.state.v_km_s.y,
                self.state.v_km_s.z,
            ],
        )
    }

    pub fn propagate_seconds(self, dt_seconds: f64) -> Result<Self, PropagationError> {
        let propagated = propagate_two_body(self.attractor.mu_km3_s2, &self.state, dt_seconds)?;
        Ok(Self {
            attractor: self.attractor,
            state: propagated,
        })
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::bodies::EARTH;
    use crate::constants::SECONDS_PER_MINUTE;

    use super::Orbit;

    #[test]
    fn propagates_vallado_example_2_4() {
        let orbit = Orbit::from_vectors(
            EARTH,
            [1131.340, -2282.343, 6672.423],
            [-5.64305, 4.30333, 2.42879],
        );
        let propagated = orbit.propagate_seconds(40.0 * SECONDS_PER_MINUTE).unwrap();
        let (r, v) = propagated.rv();

        let expected_r = [-4219.7527, 4363.0292, -3958.7666];
        let expected_v = [3.689866, -1.916735, -6.112511];

        for i in 0..3 {
            assert_relative_eq!(r[i], expected_r[i], max_relative = 1e-5);
            assert_relative_eq!(v[i], expected_v[i], max_relative = 1e-4);
        }
    }

    #[test]
    fn zero_time_returns_same_state() {
        let orbit = Orbit::from_vectors(
            EARTH,
            [1131.340, -2282.343, 6672.423],
            [-5.64305, 4.30333, 2.42879],
        );
        let propagated = orbit.propagate_seconds(0.0).unwrap();
        let (r, v) = propagated.rv();
        let (r0, v0) = orbit.rv();
        assert_eq!(r, r0);
        assert_eq!(v, v0);
    }

    #[test]
    fn hyperbolic_case_matches_reference_norms() {
        let orbit = Orbit::from_vectors(EARTH, [6678.1363, 0.0, 0.0], [0.0, 15.0, 0.0]);
        let propagated = orbit.propagate_seconds(14941.0).unwrap();
        let (r, v) = propagated.rv();
        let r_norm = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
        let v_norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert_relative_eq!(r_norm, 163180.0, max_relative = 1e-4);
        assert_relative_eq!(v_norm, 10.51, max_relative = 1e-3);
    }
}
