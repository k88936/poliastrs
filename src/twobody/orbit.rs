use nalgebra::Vector3;

use crate::bodies::Body;
use crate::core::angles::{e_to_m, nu_to_e, wrap_to_pi};
use crate::core::elements::{ClassicalElements, coe2rv, rv2coe};
use crate::twobody::propagation::farnocchia::{PropagationError, propagate_two_body};
use crate::twobody::propagation::vallado::{ValladoError, propagate_vallado};
use crate::twobody::states::CartesianState;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orbit {
    pub attractor: Body,
    pub state: CartesianState,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OrbitError {
    Propagation(PropagationError),
    Vallado(ValladoError),
    UnreachableAnomaly,
}

impl From<PropagationError> for OrbitError {
    fn from(value: PropagationError) -> Self {
        Self::Propagation(value)
    }
}

impl From<ValladoError> for OrbitError {
    fn from(value: ValladoError) -> Self {
        Self::Vallado(value)
    }
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

    pub fn from_classical(attractor: Body, coe: ClassicalElements) -> Self {
        Self {
            attractor,
            state: coe2rv(attractor.mu_km3_s2, coe),
        }
    }

    pub fn classical(self) -> ClassicalElements {
        rv2coe(self.attractor.mu_km3_s2, &self.state)
    }

    pub fn propagate_to_anomaly(self, target_nu_rad: f64) -> Result<Self, OrbitError> {
        let coe = self.classical();
        let mu = self.attractor.mu_km3_s2;
        if coe.ecc >= 1.0 {
            return Err(OrbitError::UnreachableAnomaly);
        }

        let a = coe.p_km / (1.0 - coe.ecc * coe.ecc);
        let n = (mu / (a * a * a)).sqrt();
        let e_now = nu_to_e(wrap_to_pi(coe.nu_rad), coe.ecc);
        let e_target = nu_to_e(wrap_to_pi(target_nu_rad), coe.ecc);
        let m_now = e_to_m(e_now, coe.ecc);
        let m_target = e_to_m(e_target, coe.ecc);
        let mut delta_m = (m_target - m_now).rem_euclid(std::f64::consts::TAU);
        if delta_m.abs() < 1e-14 {
            delta_m = std::f64::consts::TAU;
        }
        self.propagate_seconds(delta_m / n).map_err(Into::into)
    }

    pub fn propagate_seconds_vallado(self, dt_seconds: f64) -> Result<Self, OrbitError> {
        let propagated = propagate_vallado(self.attractor.mu_km3_s2, &self.state, dt_seconds)?;
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
    use crate::core::angles::{e_to_m, f_to_m, m_to_e, m_to_f, nu_to_e, nu_to_f};
    use crate::core::elements::ClassicalElements;

    use super::{Orbit, OrbitError};

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

    #[test]
    fn from_classical_roundtrip_preserves_elements() {
        let coe = ClassicalElements {
            p_km: 7000.0,
            ecc: 0.1,
            inc_rad: 0.3,
            raan_rad: 0.4,
            argp_rad: 0.5,
            nu_rad: 1.0,
        };
        let orbit = Orbit::from_classical(EARTH, coe);
        let out = orbit.classical();
        assert_relative_eq!(out.p_km, coe.p_km, max_relative = 1e-9);
        assert_relative_eq!(out.ecc, coe.ecc, max_relative = 1e-9);
        assert_relative_eq!(out.inc_rad, coe.inc_rad, max_relative = 1e-9);
    }

    #[test]
    fn propagate_to_anomaly_matches_periapsis_apoapsis_radius() {
        let a = 149_597_870.7;
        let ecc = 1.0 / 3.0;
        let p = a * (1.0 - ecc * ecc);
        let orbit = Orbit::from_classical(
            EARTH,
            ClassicalElements {
                p_km: p,
                ecc,
                inc_rad: 0.0,
                raan_rad: 0.0,
                argp_rad: 0.0,
                nu_rad: 10f64.to_radians(),
            },
        );
        let per = orbit.propagate_to_anomaly(0.0).unwrap();
        let apo = orbit.propagate_to_anomaly(std::f64::consts::PI).unwrap();
        let (r_per, _) = per.rv();
        let (r_apo, _) = apo.rv();
        let rp = (r_per[0] * r_per[0] + r_per[1] * r_per[1] + r_per[2] * r_per[2]).sqrt();
        let ra = (r_apo[0] * r_apo[0] + r_apo[1] * r_apo[1] + r_apo[2] * r_apo[2]).sqrt();
        assert_relative_eq!(rp, a * (1.0 - ecc), max_relative = 1e-8);
        assert_relative_eq!(ra, a * (1.0 + ecc), max_relative = 1e-8);
    }

    #[test]
    fn propagate_to_anomaly_open_orbit_returns_error() {
        let orbit = Orbit::from_vectors(EARTH, [6678.1363, 0.0, 0.0], [0.0, 15.0, 0.0]);
        let err = orbit.propagate_to_anomaly(0.1).unwrap_err();
        assert_eq!(err, OrbitError::UnreachableAnomaly);
    }

    #[test]
    fn circular_orbit_keeps_geometry_short_step() {
        let orbit = Orbit::from_vectors(EARTH, [6678.1363, 0.0, 0.0], [0.0, 7.7257602321, 0.0]);
        let res = orbit.propagate_seconds(50.0).unwrap();
        let in_coe = orbit.classical();
        let out_coe = res.classical();
        assert_relative_eq!(in_coe.p_km, out_coe.p_km, max_relative = 1e-8);
        assert_relative_eq!(out_coe.ecc, 0.0, epsilon = 1e-6);
        assert_relative_eq!(in_coe.inc_rad, out_coe.inc_rad, epsilon = 1e-10);
    }

    #[test]
    fn anomaly_conversion_roundtrip_elliptic_hyperbolic() {
        let ecc = 0.3;
        let nu = 1.0;
        let e = nu_to_e(nu, ecc);
        let m = e_to_m(e, ecc);
        let e2 = m_to_e(m, ecc);
        assert_relative_eq!(e, e2, epsilon = 1e-11);

        let ecc_h = 1.2;
        let nu_h = 1.0;
        let f = nu_to_f(nu_h, ecc_h);
        let m_h = f_to_m(f, ecc_h);
        let f2 = m_to_f(m_h, ecc_h);
        assert_relative_eq!(f, f2, epsilon = 1e-11);
    }

    #[test]
    fn vallado_agrees_with_farnocchia_for_reference_case() {
        let orbit = Orbit::from_vectors(
            EARTH,
            [1131.340, -2282.343, 6672.423],
            [-5.64305, 4.30333, 2.42879],
        );
        let f = orbit.propagate_seconds(2400.0).unwrap();
        let v = orbit.propagate_seconds_vallado(2400.0).unwrap();
        let (rf, vf) = f.rv();
        let (rv, vv) = v.rv();
        for i in 0..3 {
            assert_relative_eq!(rf[i], rv[i], max_relative = 5e-5);
            assert_relative_eq!(vf[i], vv[i], max_relative = 5e-5);
        }
    }
}
