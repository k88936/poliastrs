#![allow(non_snake_case)]

use crate::core::angles;

pub fn nu_to_E(nu_rad: f64, ecc: f64) -> f64 {
    angles::nu_to_e(nu_rad, ecc)
}

pub fn E_to_nu(e_rad: f64, ecc: f64) -> f64 {
    angles::e_to_nu(e_rad, ecc)
}

pub fn E_to_M(e_rad: f64, ecc: f64) -> f64 {
    angles::e_to_m(e_rad, ecc)
}

pub fn M_to_E(m_rad: f64, ecc: f64) -> f64 {
    angles::m_to_e(m_rad, ecc)
}

pub fn nu_to_F(nu_rad: f64, ecc: f64) -> f64 {
    angles::nu_to_f(nu_rad, ecc)
}

pub fn F_to_nu(f_rad: f64, ecc: f64) -> f64 {
    angles::f_to_nu(f_rad, ecc)
}

pub fn F_to_M(f_rad: f64, ecc: f64) -> f64 {
    angles::f_to_m(f_rad, ecc)
}

pub fn M_to_F(m_rad: f64, ecc: f64) -> f64 {
    angles::m_to_f(m_rad, ecc)
}

pub fn fp_angle(nu_rad: f64, ecc: f64) -> f64 {
    (ecc * nu_rad.sin()).atan2(1.0 + ecc * nu_rad.cos())
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::bodies::EARTH;
    use crate::core::elements::{ClassicalElements, coe2mee, coe2rv, mee2coe, rv2coe};

    use super::{E_to_M, E_to_nu, F_to_M, F_to_nu, M_to_E, M_to_F, fp_angle, nu_to_E, nu_to_F};

    fn assert_angle_close(a: f64, b: f64, atol: f64) {
        let d =
            (a - b + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI;
        assert!(d.abs() <= atol, "angle mismatch: {a} vs {b}, delta={d}");
    }

    #[test]
    fn true_to_eccentric_matches_reference() {
        let data = [
            (0.0, 0.0_f64.to_radians(), 0.0_f64.to_radians()),
            (0.05, 10.52321_f64.to_radians(), 11.05994_f64.to_radians()),
            (0.10, 54.67466_f64.to_radians(), 59.49810_f64.to_radians()),
            (0.35, 142.27123_f64.to_radians(), 153.32411_f64.to_radians()),
            (0.61, 161.87359_f64.to_radians(), 171.02189_f64.to_radians()),
        ];
        for (ecc, expected_e, nu) in data {
            assert_relative_eq!(nu_to_E(nu, ecc), expected_e, max_relative = 1e-6);
        }
    }

    #[test]
    fn true_to_eccentric_hyperbolic_matches_reference() {
        let nu = 100.0_f64.to_radians();
        let ecc = 2.7696;
        let expected_f = 2.2927;
        assert_relative_eq!(nu_to_F(nu, ecc), expected_f, max_relative = 1e-4);
    }

    #[test]
    fn mean_to_true_matches_reference_elliptic() {
        let data = [
            (0.0, 0.0_f64.to_radians(), 0.0_f64.to_radians()),
            (0.05, 10.0_f64.to_radians(), 11.06_f64.to_radians()),
            (0.06, 30.0_f64.to_radians(), 33.67_f64.to_radians()),
            (0.04, 120.0_f64.to_radians(), 123.87_f64.to_radians()),
            (0.14, 65.0_f64.to_radians(), 80.50_f64.to_radians()),
            (0.19, 21.0_f64.to_radians(), 30.94_f64.to_radians()),
            (0.35, 65.0_f64.to_radians(), 105.71_f64.to_radians()),
            (0.48, 180.0_f64.to_radians(), 180.0_f64.to_radians()),
            (0.75, 125.0_f64.to_radians(), 167.57_f64.to_radians()),
        ];
        for (ecc, mean, nu_expected) in data {
            let nu = E_to_nu(M_to_E(mean, ecc), ecc);
            assert_relative_eq!(nu, nu_expected, max_relative = 1e-4);
        }
    }

    #[test]
    fn true_to_mean_matches_reference_elliptic() {
        let data = [
            (0.0, 0.0_f64.to_radians(), 0.0_f64.to_radians()),
            (0.05, 10.0_f64.to_radians(), 11.06_f64.to_radians()),
            (0.06, 30.0_f64.to_radians(), 33.67_f64.to_radians()),
            (0.04, 120.0_f64.to_radians(), 123.87_f64.to_radians()),
            (0.14, 65.0_f64.to_radians(), 80.50_f64.to_radians()),
            (0.19, 21.0_f64.to_radians(), 30.94_f64.to_radians()),
            (0.35, 65.0_f64.to_radians(), 105.71_f64.to_radians()),
            (0.48, 180.0_f64.to_radians(), 180.0_f64.to_radians()),
            (0.75, 125.0_f64.to_radians(), 167.57_f64.to_radians()),
        ];
        for (ecc, mean, nu_value) in data {
            let m = E_to_M(nu_to_E(nu_value, ecc), ecc);
            assert_relative_eq!(m, mean, max_relative = 1e-4);
        }
    }

    #[test]
    fn true_to_mean_hyperbolic_matches_reference() {
        let nu = 100.0_f64.to_radians();
        let ecc = 2.7696;
        let expected_m = 11.279;
        assert_relative_eq!(F_to_M(nu_to_F(nu, ecc), ecc), expected_m, max_relative = 1e-4);
    }

    #[test]
    fn mean_to_true_hyperbolic_matches_reference() {
        let mean = 11.279;
        for (ecc, expected_nu_deg) in [(1.1, 153.51501), (2.7696, 100.0)] {
            let nu = F_to_nu(M_to_F(mean, ecc), ecc);
            assert_relative_eq!(nu.to_degrees(), expected_nu_deg, max_relative = 1e-4);
        }
    }

    #[test]
    fn flight_path_angle_reference() {
        let gamma = fp_angle(109.5_f64.to_radians(), 0.6);
        assert_relative_eq!(gamma.to_degrees(), 35.26, max_relative = 1e-3);
    }

    #[test]
    fn mean_to_true_hyperbolic_highecc_roundtrip() {
        for ecc in [1.5, 3200.0] {
            for k in 0..100 {
                let nu_expected =
                    -std::f64::consts::PI / 3.0 + (2.0 * std::f64::consts::PI / 3.0) * (k as f64) / 99.0;
                let m = F_to_M(nu_to_F(nu_expected, ecc), ecc);
                let nu = F_to_nu(M_to_F(m, ecc), ecc);
                assert_relative_eq!(nu, nu_expected, max_relative = 1e-4);
            }
        }
    }

    #[test]
    fn eccentric_to_true_range_roundtrip() {
        for i in 0..10 {
            let e_anomaly = -std::f64::consts::PI + 2.0 * std::f64::consts::PI * (i as f64) / 9.0;
            for j in 0..10 {
                let ecc = 0.1 + 0.8 * (j as f64) / 9.0;
                let nu = E_to_nu(e_anomaly, ecc);
                let e_back = nu_to_E(nu, ecc);
                assert_angle_close(e_back, e_anomaly, 1e-8);
            }
        }
    }

    #[test]
    fn convert_between_coe_and_rv_is_transitive() {
        let classical = ClassicalElements {
            p_km: 11067.790,
            ecc: 0.83285,
            inc_rad: 87.87_f64.to_radians(),
            raan_rad: 227.89_f64.to_radians(),
            argp_rad: 53.38_f64.to_radians(),
            nu_rad: 92.335_f64.to_radians(),
        };
        let rv = coe2rv(EARTH.mu_km3_s2, classical);
        let res = rv2coe(EARTH.mu_km3_s2, &rv);
        assert_relative_eq!(res.p_km, classical.p_km, epsilon = 1e-8);
        assert_relative_eq!(res.ecc, classical.ecc, epsilon = 1e-8);
        assert_relative_eq!(res.inc_rad, classical.inc_rad, epsilon = 1e-8);
    }

    #[test]
    fn convert_between_coe_and_mee_is_transitive() {
        let classical = ClassicalElements {
            p_km: 11067.790,
            ecc: 0.83285,
            inc_rad: 87.87_f64.to_radians(),
            raan_rad: 227.89_f64.to_radians(),
            argp_rad: 53.38_f64.to_radians(),
            nu_rad: 92.335_f64.to_radians(),
        };
        let mee = coe2mee(classical).unwrap();
        let res = mee2coe(mee);
        assert_relative_eq!(res.p_km, classical.p_km, epsilon = 1e-8);
        assert_relative_eq!(res.ecc, classical.ecc, epsilon = 1e-8);
    }

    #[test]
    fn convert_coe_and_rv_circular() {
        let expected = ClassicalElements {
            p_km: 24464.560,
            ecc: 0.0,
            inc_rad: 0.122138,
            raan_rad: 1.00681,
            argp_rad: 0.0,
            nu_rad: 0.048363,
        };
        let res = rv2coe(EARTH.mu_km3_s2, &coe2rv(EARTH.mu_km3_s2, expected));
        assert_relative_eq!(res.p_km, expected.p_km, epsilon = 1e-8);
        assert_relative_eq!(res.ecc, expected.ecc, epsilon = 1e-8);
    }

    #[test]
    fn convert_coe_and_rv_hyperbolic() {
        let expected = ClassicalElements {
            p_km: 48848.56334147761,
            ecc: 1.7311,
            inc_rad: 0.122138,
            raan_rad: 1.00681,
            argp_rad: 3.10686,
            nu_rad: 0.12741601769795755,
        };
        let res = rv2coe(EARTH.mu_km3_s2, &coe2rv(EARTH.mu_km3_s2, expected));
        assert_relative_eq!(res.p_km, expected.p_km, epsilon = 1e-8);
        assert_relative_eq!(res.ecc, expected.ecc, epsilon = 1e-8);
    }

    #[test]
    fn convert_coe_and_rv_equatorial() {
        let expected = ClassicalElements {
            p_km: 11388.0762905224,
            ecc: 0.7311,
            inc_rad: 0.0,
            raan_rad: 0.0,
            argp_rad: 3.10686,
            nu_rad: 0.44369564302687126,
        };
        let res = rv2coe(EARTH.mu_km3_s2, &coe2rv(EARTH.mu_km3_s2, expected));
        assert_relative_eq!(res.p_km, expected.p_km, epsilon = 1e-8);
        assert_relative_eq!(res.ecc, expected.ecc, epsilon = 1e-8);
    }

    #[test]
    fn convert_coe_and_rv_circular_equatorial() {
        let expected = ClassicalElements {
            p_km: 11388.0762905224,
            ecc: 0.0,
            inc_rad: 0.0,
            raan_rad: 0.0,
            argp_rad: 0.0,
            nu_rad: 0.44369564302687126,
        };
        let res = rv2coe(EARTH.mu_km3_s2, &coe2rv(EARTH.mu_km3_s2, expected));
        assert_relative_eq!(res.p_km, expected.p_km, epsilon = 1e-8);
        assert_relative_eq!(res.ecc, expected.ecc, epsilon = 1e-8);
    }
}
