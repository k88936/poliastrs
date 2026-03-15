pub fn stumpff_c2(psi: f64) -> f64 {
    if psi > 1e-12 {
        let sqrt_psi = psi.sqrt();
        (1.0 - sqrt_psi.cos()) / psi
    } else if psi < -1e-12 {
        let sqrt_minus_psi = (-psi).sqrt();
        (sqrt_minus_psi.cosh() - 1.0) / -psi
    } else {
        // Taylor series: 1/2! - psi/4! + psi^2/6! ...
        0.5 - psi / 24.0 + psi * psi / 720.0
    }
}

pub fn stumpff_c3(psi: f64) -> f64 {
    if psi > 1e-12 {
        let sqrt_psi = psi.sqrt();
        (sqrt_psi - sqrt_psi.sin()) / psi.powf(1.5)
    } else if psi < -1e-12 {
        let sqrt_minus_psi = (-psi).sqrt();
        (sqrt_minus_psi.sinh() - sqrt_minus_psi) / (-psi).powf(1.5)
    } else {
        // Taylor series: 1/3! - psi/5! + psi^2/7! ...
        1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use super::*;

    #[test]
    fn test_stumpff_functions_near_zero() {
        let psi = 0.5f64;
        let expected_c2 = (1.0 - psi.sqrt().cos()) / psi;
        let expected_c3 = (psi.sqrt() - psi.sqrt().sin()) / psi.powf(1.5);

        assert_relative_eq!(stumpff_c2(psi), expected_c2, epsilon = 1e-12);
        assert_relative_eq!(stumpff_c3(psi), expected_c3, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_functions_above_zero() {
        let psi = 3.0f64;
        let expected_c2 = (1.0 - psi.sqrt().cos()) / psi;
        let expected_c3 = (psi.sqrt() - psi.sqrt().sin()) / psi.powf(1.5);

        assert_relative_eq!(stumpff_c2(psi), expected_c2, epsilon = 1e-12);
        assert_relative_eq!(stumpff_c3(psi), expected_c3, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_functions_under_zero() {
        let psi = -3.0f64;
        let sqrt_neg_psi = (-psi).sqrt();
        let expected_c2 = (sqrt_neg_psi.cosh() - 1.0) / -psi;
        let expected_c3 = (sqrt_neg_psi.sinh() - sqrt_neg_psi) / (-psi).powf(1.5);

        assert_relative_eq!(stumpff_c2(psi), expected_c2, epsilon = 1e-12);
        assert_relative_eq!(stumpff_c3(psi), expected_c3, epsilon = 1e-12);
    }
}
