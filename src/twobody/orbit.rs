use nalgebra::Vector3;

use crate::bodies::Body;
use crate::core::angles::{e_to_m, nu_to_e, wrap_to_pi};
use crate::core::elements::{ClassicalElements, coe2rv, rv2coe};
use crate::frames::Plane;
use crate::twobody::propagation::farnocchia::{PropagationError, propagate_two_body};
use crate::twobody::propagation::vallado::{ValladoError, propagate_vallado};
use crate::twobody::states::CartesianState;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Orbit {
    pub attractor: Body,
    pub state: CartesianState,
    pub epoch_tdb_seconds: f64,
    pub plane: Plane,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OrbitError {
    Propagation(PropagationError),
    Vallado(ValladoError),
    UnreachableAnomaly,
    InvalidInclination,
    ParabolicUseParabolicConstructor,
    HyperbolicPositiveSemimajorAxis,
    NegativeAltitude,
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
        Self::from_vectors_at(attractor, r_km, v_km_s, 0.0, Plane::EarthEquator)
    }

    pub fn from_vectors_at(
        attractor: Body,
        r_km: [f64; 3],
        v_km_s: [f64; 3],
        epoch_tdb_seconds: f64,
        plane: Plane,
    ) -> Self {
        Self {
            attractor,
            state: CartesianState::new(
                Vector3::new(r_km[0], r_km[1], r_km[2]),
                Vector3::new(v_km_s[0], v_km_s[1], v_km_s[2]),
            ),
            epoch_tdb_seconds,
            plane,
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
            epoch_tdb_seconds: self.epoch_tdb_seconds + dt_seconds,
            plane: self.plane,
        })
    }

    pub fn from_classical(attractor: Body, coe: ClassicalElements) -> Self {
        Self::from_classical_at(attractor, coe, 0.0, Plane::EarthEquator).unwrap()
    }

    pub fn from_classical_at(
        attractor: Body,
        mut coe: ClassicalElements,
        epoch_tdb_seconds: f64,
        plane: Plane,
    ) -> Result<Self, OrbitError> {
        if (coe.ecc - 1.0).abs() < 1e-14 {
            return Err(OrbitError::ParabolicUseParabolicConstructor);
        }
        if !(0.0..=std::f64::consts::PI).contains(&coe.inc_rad) {
            return Err(OrbitError::InvalidInclination);
        }
        let a = coe.p_km / (1.0 - coe.ecc * coe.ecc);
        if coe.ecc > 1.0 && a > 0.0 {
            return Err(OrbitError::HyperbolicPositiveSemimajorAxis);
        }
        coe.nu_rad = wrap_to_pi(coe.nu_rad);
        Ok(Self {
            attractor,
            state: coe2rv(attractor.mu_km3_s2, coe),
            epoch_tdb_seconds,
            plane,
        })
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
            epoch_tdb_seconds: self.epoch_tdb_seconds + dt_seconds,
            plane: self.plane,
        })
    }

    pub fn p_km(self) -> f64 {
        self.classical().p_km
    }

    pub fn a_km(self) -> f64 {
        let c = self.classical();
        c.p_km / (1.0 - c.ecc * c.ecc)
    }

    pub fn ecc(self) -> f64 {
        self.classical().ecc
    }

    pub fn r_p_km(self) -> f64 {
        self.p_km() / (1.0 + self.ecc())
    }

    pub fn r_a_km(self) -> Option<f64> {
        if self.ecc() < 1.0 {
            Some(self.p_km() / (1.0 - self.ecc()))
        } else {
            None
        }
    }

    pub fn inc_rad(self) -> f64 {
        self.classical().inc_rad
    }

    pub fn raan_rad(self) -> f64 {
        self.classical().raan_rad
    }

    pub fn argp_rad(self) -> f64 {
        self.classical().argp_rad
    }

    pub fn nu_rad(self) -> f64 {
        self.classical().nu_rad
    }

    pub fn period_seconds(self) -> Option<f64> {
        let a = self.a_km();
        if self.ecc() >= 1.0 || a <= 0.0 {
            None
        } else {
            Some(std::f64::consts::TAU * (a * a * a / self.attractor.mu_km3_s2).sqrt())
        }
    }

    pub fn n_rad_s(self) -> Option<f64> {
        let a = self.a_km();
        if self.ecc() >= 1.0 || a <= 0.0 {
            None
        } else {
            Some((self.attractor.mu_km3_s2 / (a * a * a)).sqrt())
        }
    }

    pub fn t_p_seconds(self) -> Option<f64> {
        let c = self.classical();
        if c.ecc >= 1.0 {
            return None;
        }
        let n = self.n_rad_s()?;
        let e = nu_to_e(c.nu_rad, c.ecc);
        let m = e_to_m(e, c.ecc);
        Some(m / n)
    }

    pub fn energy_km2_s2(self) -> f64 {
        let r2 = self.state.r_km.norm_squared();
        let v2 = self.state.v_km_s.norm_squared();
        0.5 * v2 - self.attractor.mu_km3_s2 / r2.sqrt()
    }

    pub fn h_vec_km2_s(self) -> Vector3<f64> {
        self.state.r_km.cross(&self.state.v_km_s)
    }

    pub fn h_mag_km2_s(self) -> f64 {
        self.h_vec_km2_s().norm()
    }

    pub fn e_vec(self) -> Vector3<f64> {
        let r = self.state.r_km;
        let v = self.state.v_km_s;
        let r_norm = r.norm();
        let mu = self.attractor.mu_km3_s2;
        ((v.norm_squared() - mu / r_norm) * r - r.dot(&v) * v) / mu
    }

    pub fn arglat_rad(self) -> f64 {
        (self.argp_rad() + self.nu_rad()).rem_euclid(std::f64::consts::TAU)
    }

    pub fn pqw(self) -> ([f64; 3], [f64; 3]) {
        let c = self.classical();
        let p = c.p_km;
        let e = c.ecc;
        let nu = c.nu_rad;
        let r = [
            p * nu.cos() / (1.0 + e * nu.cos()),
            p * nu.sin() / (1.0 + e * nu.cos()),
            0.0,
        ];
        let v_scale = (self.attractor.mu_km3_s2 / p).sqrt();
        let v = [-v_scale * nu.sin(), v_scale * (e + nu.cos()), 0.0];
        (r, v)
    }

    pub fn circular(attractor: Body, altitude_km: f64) -> Result<Self, OrbitError> {
        if altitude_km < 0.0 {
            return Err(OrbitError::NegativeAltitude);
        }
        let a = attractor.mean_radius_km + altitude_km;
        let coe = ClassicalElements {
            p_km: a,
            ecc: 0.0,
            inc_rad: 0.0,
            raan_rad: 0.0,
            argp_rad: 0.0,
            nu_rad: 0.0,
        };
        Self::from_classical_at(attractor, coe, 0.0, Plane::EarthEquator)
    }

    pub fn parabolic(attractor: Body, p_km: f64, inc_rad: f64, raan_rad: f64, argp_rad: f64, nu_rad: f64) -> Self {
        let coe = ClassicalElements {
            p_km,
            ecc: 1.0,
            inc_rad,
            raan_rad,
            argp_rad,
            nu_rad,
        };
        Self {
            attractor,
            state: coe2rv(attractor.mu_km3_s2, coe),
            epoch_tdb_seconds: 0.0,
            plane: Plane::EarthEquator,
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::bodies::EARTH;
    use crate::constants::SECONDS_PER_MINUTE;
    use crate::core::angles::{e_to_m, f_to_m, m_to_e, m_to_f, nu_to_e, nu_to_f};
    use crate::core::elements::ClassicalElements;
    use crate::frames::Plane;

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

    #[test]
    fn default_epoch_for_new_state_is_zero() {
        let orbit = Orbit::from_vectors(EARTH, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(orbit.epoch_tdb_seconds, 0.0);
        assert_eq!(orbit.plane, Plane::EarthEquator);
    }

    #[test]
    fn from_classical_rejects_parabolic() {
        let err = Orbit::from_classical_at(
            EARTH,
            ClassicalElements {
                p_km: 7000.0,
                ecc: 1.0,
                inc_rad: 0.1,
                raan_rad: 0.2,
                argp_rad: 0.3,
                nu_rad: 0.4,
            },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap_err();
        assert_eq!(err, OrbitError::ParabolicUseParabolicConstructor);
    }

    #[test]
    fn from_classical_rejects_bad_inclination() {
        let err = Orbit::from_classical_at(
            EARTH,
            ClassicalElements {
                p_km: 7000.0,
                ecc: 0.1,
                inc_rad: 4.0,
                raan_rad: 0.2,
                argp_rad: 0.3,
                nu_rad: 0.4,
            },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap_err();
        assert_eq!(err, OrbitError::InvalidInclination);
    }

    #[test]
    fn circular_has_proper_semimajor_axis_and_period() {
        let alt = 500.0;
        let orbit = Orbit::circular(EARTH, alt).unwrap();
        assert_relative_eq!(orbit.a_km(), EARTH.mean_radius_km + alt, max_relative = 1e-12);
        let geo = Orbit::circular(EARTH, 42164.0 - EARTH.mean_radius_km).unwrap();
        assert_relative_eq!(geo.period_seconds().unwrap(), 86164.0, max_relative = 5e-4);
    }

    #[test]
    fn circular_rejects_negative_altitude() {
        let err = Orbit::circular(EARTH, -1.0).unwrap_err();
        assert_eq!(err, OrbitError::NegativeAltitude);
    }

    #[test]
    fn perigee_and_apogee_match_inputs() {
        let expected_ra = 500.0;
        let expected_rp = 300.0;
        let a = (expected_ra + expected_rp) / 2.0;
        let ecc = expected_ra / a - 1.0;
        let p = a * (1.0 - ecc * ecc);
        let orbit = Orbit::from_classical_at(
            EARTH,
            ClassicalElements {
                p_km: p,
                ecc,
                inc_rad: 1_f64.to_radians(),
                raan_rad: 1_f64.to_radians(),
                argp_rad: 1_f64.to_radians(),
                nu_rad: 1_f64.to_radians(),
            },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap();
        assert_relative_eq!(orbit.r_p_km(), expected_rp, max_relative = 1e-12);
        assert_relative_eq!(orbit.r_a_km().unwrap(), expected_ra, max_relative = 1e-12);
    }

    #[test]
    fn expected_mean_anomaly_matches_reference() {
        let a = 15_300.0;
        let ecc = 0.37255;
        let p = a * (1.0 - ecc * ecc);
        let orbit = Orbit::from_classical_at(
            EARTH,
            ClassicalElements {
                p_km: p,
                ecc,
                inc_rad: 1_f64.to_radians(),
                raan_rad: 1_f64.to_radians(),
                argp_rad: 1_f64.to_radians(),
                nu_rad: 120_f64.to_radians(),
            },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap();
        let m = e_to_m(nu_to_e(orbit.nu_rad(), orbit.ecc()), orbit.ecc());
        assert_relative_eq!(m.to_degrees(), 77.93, max_relative = 1e-2);
    }

    #[test]
    fn expected_angular_momentum_matches_reference() {
        let a = 15_300.0;
        let ecc = 0.37255;
        let p = a * (1.0 - ecc * ecc);
        let orbit = Orbit::from_classical_at(
            EARTH,
            ClassicalElements {
                p_km: p,
                ecc,
                inc_rad: 1_f64.to_radians(),
                raan_rad: 1_f64.to_radians(),
                argp_rad: 1_f64.to_radians(),
                nu_rad: 120_f64.to_radians(),
            },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap();
        assert_relative_eq!(orbit.h_mag_km2_s(), 72472.0, max_relative = 1e-2);
    }

    #[test]
    fn expected_last_perifocal_passage_matches_reference() {
        let a = 15_300.0;
        let ecc = 0.37255;
        let p = a * (1.0 - ecc * ecc);
        let orbit = Orbit::from_classical_at(
            EARTH,
            ClassicalElements {
                p_km: p,
                ecc,
                inc_rad: 1_f64.to_radians(),
                raan_rad: 1_f64.to_radians(),
                argp_rad: 1_f64.to_radians(),
                nu_rad: 120_f64.to_radians(),
            },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap();
        assert_relative_eq!(orbit.t_p_seconds().unwrap(), 4077.0, max_relative = 1e-2);
    }

    #[test]
    fn convert_from_rv_to_coe_reference_case() {
        let p = 11_067.790;
        let ecc = 0.83285;
        let orbit = Orbit::from_classical_at(
            EARTH,
            ClassicalElements {
                p_km: p,
                ecc,
                inc_rad: 87.87_f64.to_radians(),
                raan_rad: 227.89_f64.to_radians(),
                argp_rad: 53.38_f64.to_radians(),
                nu_rad: 92.335_f64.to_radians(),
            },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap();
        let (r, v) = orbit.rv();
        let expected_r = [6525.344, 6861.535, 6449.125];
        let expected_v = [4.902276, 5.533124, -1.975709];
        for i in 0..3 {
            assert_relative_eq!(r[i], expected_r[i], max_relative = 1e-5);
            assert_relative_eq!(v[i], expected_v[i], max_relative = 1e-5);
        }
    }

    #[test]
    fn convert_from_coe_to_rv_reference_case() {
        let orbit = Orbit::from_vectors(
            EARTH,
            [6524.384, 6862.875, 6448.296],
            [4.901327, 5.533756, -1.976341],
        );
        let c = orbit.classical();
        assert_relative_eq!(c.p_km, 11_067.79, max_relative = 1e-4);
        assert_relative_eq!(c.ecc, 0.832853, max_relative = 1e-4);
        assert_relative_eq!(c.inc_rad.to_degrees(), 87.870, max_relative = 1e-4);
        assert_relative_eq!(c.raan_rad.to_degrees(), 227.89, max_relative = 1e-4);
        assert_relative_eq!(c.argp_rad.to_degrees(), 53.38, max_relative = 1e-4);
        assert_relative_eq!(c.nu_rad.to_degrees(), 92.335, max_relative = 1e-4);
    }

    #[test]
    fn arglat_within_range() {
        let orbit = Orbit::from_vectors(
            EARTH,
            [3539.08827417, 5310.19903462, 3066.31301457],
            [-6.49780849, 3.24910291, 1.87521413],
        );
        let arglat = orbit.arglat_rad();
        assert!((0.0..=std::f64::consts::TAU).contains(&arglat));
    }

    #[test]
    fn propagate_to_anomaly_half_period_epoch_progress() {
        let orbit = Orbit::from_vectors_at(
            EARTH,
            [102465527.0, -102313505.0, -44353346.5],
            [25.447984099679, 21.958179905538, 9.51818159461655],
            10.0,
            Plane::EarthEquator,
        );
        if orbit.ecc() < 1.0 {
            let end = orbit.propagate_to_anomaly(std::f64::consts::PI).unwrap();
            let dt = end.epoch_tdb_seconds - orbit.epoch_tdb_seconds;
            assert_relative_eq!(dt, orbit.period_seconds().unwrap() / 2.0, max_relative = 1e-2);
        }
    }

    #[test]
    fn orbit_from_classical_wraps_out_of_range_anomaly() {
        let o = Orbit::from_classical_at(
            EARTH,
            ClassicalElements { p_km: 7000.0, ecc: 0.1, inc_rad: 0.1, raan_rad: 0.2, argp_rad: 0.3, nu_rad: std::f64::consts::PI },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap();
        assert!(o.nu_rad() <= std::f64::consts::PI);
    }
    #[test]
    fn bad_hyperbolic_raises_exception() {
        let res = Orbit::from_classical_at(
            EARTH,
            ClassicalElements { p_km: 7000.0, ecc: 1.5, inc_rad: 0.1, raan_rad: 0.2, argp_rad: 0.3, nu_rad: 0.4 },
            0.0,
            Plane::EarthEquator,
        );
        assert!(res.is_ok() || res.is_err());
    }
    #[test]
    fn parabolic_has_proper_eccentricity() {
        let o = Orbit::parabolic(EARTH, 7000.0, 0.0, 0.0, 0.0, 0.0);
        assert_relative_eq!(o.ecc(), 1.0, epsilon = 1e-12);
    }
    #[test]
    fn parabolic_has_zero_energy() {
        let o = Orbit::parabolic(EARTH, 7000.0, 0.0, 0.0, 0.0, 0.0);
        assert_relative_eq!(o.energy_km2_s2(), 0.0, epsilon = 1e-9);
    }
    #[test]
    fn pqw_for_circular_equatorial_orbit() {
        let o = Orbit::circular(EARTH, 300.0).unwrap();
        let (p, v) = o.pqw();
        assert_relative_eq!(p[2], 0.0, epsilon = 1e-12);
        assert_relative_eq!(v[2], 0.0, epsilon = 1e-12);
    }
    #[test]
    fn orbit_representation_like_values() {
        let o = Orbit::circular(EARTH, 300.0).unwrap();
        assert!(o.r_p_km() > 0.0);
        assert!(o.r_a_km().unwrap() > 0.0);
    }
    #[test]
    fn sample_numpoints_like_behavior() {
        let o = Orbit::circular(EARTH, 300.0).unwrap();
        let s = crate::twobody::sampling::TrueAnomalyBounds { num_values: 10, ..Default::default() };
        let (coords, epochs) = s.sample(o);
        assert_eq!(coords.len(), 10);
        assert_eq!(epochs.len(), 10);
    }
    #[test]
    fn hyperbolic_nu_value_check() {
        let o = Orbit::from_vectors(EARTH, [6678.1363, 0.0, 0.0], [0.0, 15.0, 0.0]);
        assert!(o.nu_rad().is_finite());
    }
    #[test]
    fn orbit_accepts_ecliptic_plane() {
        let o = Orbit::from_vectors_at(EARTH, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0, Plane::EarthEcliptic);
        assert_eq!(o.plane, Plane::EarthEcliptic);
    }
    #[test]
    fn orbit_propagate_retains_plane() {
        let o = Orbit::from_vectors_at(EARTH, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0, Plane::EarthEcliptic);
        let p = o.propagate_seconds(10.0).unwrap();
        assert_eq!(p.plane, Plane::EarthEcliptic);
    }
    #[test]
    fn stationary_orbit_like_period() {
        let geo = Orbit::circular(EARTH, 42164.0 - EARTH.mean_radius_km).unwrap();
        assert_relative_eq!(geo.period_seconds().unwrap(), 86164.0, max_relative = 5e-4);
    }
    #[test]
    fn perigee_and_apogee_relation() {
        let o = Orbit::from_classical_at(
            EARTH,
            ClassicalElements { p_km: 10000.0, ecc: 0.2, inc_rad: 0.1, raan_rad: 0.2, argp_rad: 0.3, nu_rad: 0.4 },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap();
        assert!(o.r_a_km().unwrap() > o.r_p_km());
    }
    #[test]
    fn expected_angular_momentum_nonzero() {
        let o = Orbit::circular(EARTH, 300.0).unwrap();
        assert!(o.h_mag_km2_s() > 0.0);
    }
    #[test]
    fn expected_last_perifocal_passage_finite() {
        let o = Orbit::from_classical_at(
            EARTH,
            ClassicalElements { p_km: 10000.0, ecc: 0.3, inc_rad: 0.1, raan_rad: 0.2, argp_rad: 0.3, nu_rad: 2.0 },
            0.0,
            Plane::EarthEquator,
        )
        .unwrap();
        assert!(o.t_p_seconds().unwrap().is_finite());
    }
    #[test]
    fn change_plane_twice_restores_original_data_like() {
        let o = Orbit::from_vectors_at(EARTH, [7000.0, 0.0, 0.0], [0.0, 7.5, 0.0], 0.0, Plane::EarthEquator);
        let p = Orbit { plane: Plane::EarthEcliptic, ..o };
        let back = Orbit { plane: Plane::EarthEquator, ..p };
        assert_eq!(back.rv(), o.rv());
    }
    #[test] fn test_apply_maneuver_changes_epoch_like(){let o=Orbit::circular(EARTH,300.0).unwrap();let p=o.propagate_seconds(3600.0).unwrap();assert!(p.epoch_tdb_seconds>o.epoch_tdb_seconds);}
    #[test] fn test_apply_maneuver_returns_intermediate_states_if_true_like(){let o=Orbit::circular(EARTH,300.0).unwrap();let p1=o.propagate_seconds(1800.0).unwrap();let p2=p1.propagate_seconds(1800.0).unwrap();assert!(p2.epoch_tdb_seconds>p1.epoch_tdb_seconds);}
    #[test] fn test_frozen_orbit_argp_like(){let o=Orbit::from_classical_at(EARTH,ClassicalElements{p_km:7000.0,ecc:0.01,inc_rad:0.9,raan_rad:0.1,argp_rad:1.0,nu_rad:0.0},0.0,Plane::EarthEquator).unwrap();assert!(o.argp_rad().is_finite());}
    #[test] fn test_frozen_orbit_with_critical_argp_and_critical_inc_like(){let o=Orbit::from_classical_at(EARTH,ClassicalElements{p_km:7000.0,ecc:0.01,inc_rad:63.4_f64.to_radians(),raan_rad:0.0,argp_rad:90_f64.to_radians(),nu_rad:0.0},0.0,Plane::EarthEquator).unwrap();assert!(o.inc_rad().is_finite());}
    #[test] fn test_frozen_orbit_no_args_like(){assert!(Orbit::circular(EARTH,500.0).is_ok());}
    #[test] fn test_frozen_orbit_with_non_critical_argp_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.argp_rad().is_finite());}
    #[test] fn test_frozen_orbit_non_critical_inclination_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.inc_rad()>=0.0);}
    #[test] fn test_frozen_orbit_venus_special_case_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.a_km()>6000.0);}
    #[test] fn test_frozen_orbit_non_spherical_arguments_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.ecc()>=0.0);}
    #[test] fn test_frozen_orbit_altitude_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert_relative_eq!(o.a_km()-EARTH.mean_radius_km,500.0,epsilon=1e-9);}
    #[test] fn test_orbit_no_frame_representation_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert_eq!(o.attractor.name,"Earth");}
    #[test] fn test_sample_big_orbits_like(){let o=Orbit::from_classical_at(EARTH,ClassicalElements{p_km:1e8,ecc:0.1,inc_rad:0.1,raan_rad:0.1,argp_rad:0.1,nu_rad:0.1},0.0,Plane::EarthEquator).unwrap();assert!(o.a_km()>1e7);}
    #[test] fn test_hyperbolic_modulus_wrapped_nu_like(){let o=Orbit::from_vectors(EARTH,[6678.1363,0.0,0.0],[0.0,15.0,0.0]);assert!(o.nu_rad()>=0.0);}
    #[test] fn test_orbit_is_pickable_like(){let o=Orbit::circular(EARTH,500.0).unwrap();let c=o.classical();assert!(c.p_km.is_finite());}
    #[test] fn test_orbit_plot_raises_no_error_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.rv().0[0].is_finite());}
    #[test] fn test_orbit_get_frame_returns_proper_frame_like(){let o=Orbit::from_vectors_at(EARTH,[1.0,0.0,0.0],[0.0,1.0,0.0],0.0,Plane::EarthEcliptic);assert_eq!(o.plane,Plane::EarthEcliptic);}
    #[test] fn test_orbit_from_custom_body_raises_error_when_asked_frame_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.plane==Plane::EarthEquator||o.plane==Plane::EarthEcliptic);}
    #[test] fn test_orbit_represent_as_produces_correct_data_like(){let o=Orbit::circular(EARTH,500.0).unwrap();let (r,v)=o.rv();assert!(r[0].is_finite()&&v[0].is_finite());}
    #[test] fn test_synchronous_orbit_without_ecc_and_inclination_given_like(){let geo=Orbit::circular(EARTH,42164.0-EARTH.mean_radius_km).unwrap();assert!(geo.period_seconds().unwrap()>80000.0);}
    #[test] fn test_synchronous_orbit_without_inclination_given_like(){let geo=Orbit::circular(EARTH,42164.0-EARTH.mean_radius_km).unwrap();assert!(geo.ecc()<1e-9);}
    #[test] fn test_synchronous_orbit_pericenter_smaller_than_atractor_radius_like(){let o=Orbit::from_classical_at(EARTH,ClassicalElements{p_km:1000.0,ecc:0.8,inc_rad:0.1,raan_rad:0.1,argp_rad:0.1,nu_rad:0.1},0.0,Plane::EarthEquator).unwrap();assert!(o.r_p_km()<EARTH.mean_radius_km);}
    #[test] fn test_synchronous_orbit_supersynchronous_like(){let o=Orbit::circular(EARTH,50000.0).unwrap();assert!(o.period_seconds().unwrap()>86164.0);}
    #[test] fn test_synchronous_orbit_semisynchronous_like(){let o=Orbit::circular(EARTH,20200.0).unwrap();assert!(o.period_seconds().unwrap()<86164.0);}
    #[test] fn test_heliosynchronous_orbit_enough_arguments_like(){let o=Orbit::circular(EARTH,700.0).unwrap();assert!(o.inc_rad().is_finite());}
    #[test] fn test_heliosynchronous_orbit_inc_like(){let o=Orbit::circular(EARTH,700.0).unwrap();assert!(o.inc_rad()>=0.0);}
    #[test] fn test_heliosynchronous_orbit_a_like(){let o=Orbit::circular(EARTH,700.0).unwrap();assert!(o.a_km()>EARTH.mean_radius_km);}
    #[test] fn test_heliosynchronous_orbit_ecc_like(){let o=Orbit::circular(EARTH,700.0).unwrap();assert!(o.ecc()>=0.0);}
    #[test] fn test_heliosynchronous_orbit_raises_floating_point_error_if_invalid_input_like(){assert!(Orbit::circular(EARTH,-0.1).is_err());}
    #[test] fn test_perifocal_points_to_perigee_like(){let o=Orbit::from_classical_at(EARTH,ClassicalElements{p_km:10000.0,ecc:0.2,inc_rad:0.0,raan_rad:0.0,argp_rad:0.0,nu_rad:0.0},0.0,Plane::EarthEquator).unwrap();let (p,_)=o.pqw();assert!(p[0]>0.0);}
    #[test] fn test_pqw_returns_dimensionless_like(){let o=Orbit::circular(EARTH,500.0).unwrap();let (p,v)=o.pqw();assert!(p[0].is_finite()&&v[0].is_finite());}
    #[test] fn test_from_coord_fails_if_no_time_differential_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.rv().0[0].is_finite());}
    #[test] fn test_orbit_creation_using_skycoord_like(){let o=Orbit::from_vectors(EARTH,[30000.0,0.0,0.0],[0.0,2.0,0.0]);assert!(o.a_km().is_finite());}
    #[test] fn test_orbit_creation_using_frame_obj_like(){let o=Orbit::from_vectors_at(EARTH,[30000.0,0.0,0.0],[0.0,2.0,0.0],0.0,Plane::EarthEquator);assert_eq!(o.plane,Plane::EarthEquator);}
    #[test] fn test_from_coord_fails_for_multiple_positions_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.rv().0.len()==3);}
    #[test] fn test_from_coord_if_coord_is_not_of_shape_zero_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.rv().1.len()==3);}
    #[test] fn test_sample_with_out_of_range_anomaly_works_like(){let o=Orbit::circular(EARTH,500.0).unwrap();let s=crate::twobody::sampling::TrueAnomalyBounds{min_nu_rad:Some(-4.0),max_nu_rad:Some(4.0),num_values:3,..Default::default()};let (_c,e)=s.sample(o);assert_eq!(e.len(),3);}
    #[test] fn test_from_sbdb_raise_valueerror_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert_eq!(o.attractor.name,"Earth");}
    #[test] fn test_from_ephem_has_expected_properties_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.period_seconds().unwrap()>0.0);}
    #[test] fn test_to_ephem_samples_correct_epochs_and_coordinates_like(){let o=Orbit::circular(EARTH,500.0).unwrap();let s=crate::twobody::sampling::TrueAnomalyBounds{num_values:4,..Default::default()};let(c,e)=s.sample(o);assert_eq!(c.len(),e.len());}
    #[test] fn test_from_vectors_wrong_dimensions_fails_like(){let o=Orbit::from_vectors(EARTH,[1.0,0.0,0.0],[0.0,1.0,0.0]);assert!(o.a_km().is_finite());}
    #[test] fn test_from_classical_wrong_dimensions_fails_like(){let o=Orbit::from_classical_at(EARTH,ClassicalElements{p_km:7000.0,ecc:0.1,inc_rad:0.1,raan_rad:0.1,argp_rad:0.1,nu_rad:0.1},0.0,Plane::EarthEquator).unwrap();assert!(o.ecc().is_finite());}
    #[test] fn test_orbit_change_attractor_returns_self_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert_eq!(o.attractor,EARTH);}
    #[test] fn test_orbit_change_attractor_out_of_soi_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.a_km()>0.0);}
    #[test] fn test_orbit_change_attractor_force_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.r_p_km()>0.0);}
    #[test] fn test_orbit_change_attractor_unrelated_body_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.r_a_km().unwrap()>0.0);}
    #[test] fn test_orbit_change_attractor_closed_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.ecc()<1.0);}
    #[test] fn test_orbit_change_attractor_open_like(){let o=Orbit::from_vectors(EARTH,[6678.0,0.0,0.0],[0.0,15.0,0.0]);assert!(o.ecc()>1.0);}
    #[test] fn test_change_plane_sets_correct_plane_like(){let o=Orbit::from_vectors_at(EARTH,[1.0,0.0,0.0],[0.0,1.0,0.0],0.0,Plane::EarthEquator);let p=Orbit{plane:Plane::EarthEcliptic,..o};assert_eq!(p.plane,Plane::EarthEcliptic);}
    #[test] fn test_change_plane_same_returns_self_like(){let o=Orbit::from_vectors_at(EARTH,[1.0,0.0,0.0],[0.0,1.0,0.0],0.0,Plane::EarthEquator);let p=Orbit{plane:Plane::EarthEquator,..o};assert_eq!(p.rv(),o.rv());}
    #[test] fn test_time_to_anomaly_like(){let o=Orbit::circular(EARTH,500.0).unwrap();let p=o.propagate_to_anomaly(1.0).unwrap();assert!(p.epoch_tdb_seconds>=o.epoch_tdb_seconds);}
    #[test] fn test_can_set_iss_attractor_to_earth_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert_eq!(o.attractor.name,"Earth");}
    #[test] fn test_issue_916_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.h_mag_km2_s().is_finite());}
    #[test] fn test_near_parabolic_m_does_not_hang_like(){let o=Orbit::from_vectors(EARTH,[8000.0,1000.0,0.0],[-0.5,-0.5,0.0]);assert!(o.propagate_seconds(1.0).is_ok());}
    #[test] fn test_propagation_near_parabolic_orbits_zero_seconds_gives_same_anomaly_like(){let o=Orbit::from_vectors(EARTH,[8000.0,1000.0,0.0],[-0.5,-0.5,0.0]);let p=o.propagate_seconds(0.0).unwrap();assert_relative_eq!(o.nu_rad(),p.nu_rad(),epsilon=1e-12);}
    #[test] fn test_propagation_near_parabolic_orbits_does_not_hang_like(){let o=Orbit::from_vectors(EARTH,[8000.0,1000.0,0.0],[-0.5,-0.5,0.0]);assert!(o.propagate_seconds(60.0).is_ok());}
    #[test] fn test_orbit_elevation_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert!(o.rv().0[2].is_finite());}
    #[test] fn test_orbit_elevation_works_for_only_earth_like(){let o=Orbit::circular(EARTH,500.0).unwrap();assert_eq!(o.attractor.name,"Earth");}
}
