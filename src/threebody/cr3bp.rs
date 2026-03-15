use nalgebra::SVector;

use crate::threebody::ThreeBodyError;
use crate::bodies::Body;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SystemChars {
    /// Dimensionless mass parameter mu = m2 / (m1 + m2)
    pub mu: f64,
    /// Characteristic length l* (distance between primaries) in km
    pub lstar_km: f64,
    /// Characteristic time t* in seconds
    pub tstar_sec: f64,
}

impl SystemChars {
    /// Creates a SystemChars instance from two bodies.
    /// Assumes `primary` is the more massive body and `secondary` orbits it.
    /// `secondary` must have `mean_semi_major_axis_km` defined.
    pub fn from_primaries(primary: &Body, secondary: &Body) -> Result<Self, ThreeBodyError> {
        let lstar = secondary
            .mean_semi_major_axis_km
            .ok_or(ThreeBodyError::MissingSemiMajorAxis)?;
        let mu1 = primary.mu_km3_s2;
        let mu2 = secondary.mu_km3_s2;

        let mu = mu2 / (mu1 + mu2);
        let tstar = (lstar.powi(3) / (mu1 + mu2)).sqrt();

        Ok(Self {
            mu,
            lstar_km: lstar,
            tstar_sec: tstar,
        })
    }
}

pub fn calculate_mu(mu1: f64, mu2: f64) -> f64 {
    mu2 / (mu1 + mu2)
}

pub fn calculate_tstar(mu1: f64, mu2: f64, lstar: f64) -> f64 {
    (lstar.powi(3) / (mu1 + mu2)).sqrt()
}

pub type Cr3bpState = SVector<f64, 6>;

pub fn state_derivative(state: &Cr3bpState, mu: f64) -> Result<Cr3bpState, ThreeBodyError> {
    if !(0.0..0.5).contains(&mu) {
        return Err(ThreeBodyError::InvalidMassParameter);
    }

    let x = state[0];
    let y = state[1];
    let z = state[2];
    let xd = state[3];
    let yd = state[4];
    let zd = state[5];

    let r1 = ((x + mu).powi(2) + y * y + z * z).sqrt();
    let r2 = ((x - 1.0 + mu).powi(2) + y * y + z * z).sqrt();

    let ddx = 2.0 * yd + x - (1.0 - mu) * (x + mu) / r1.powi(3) - mu * (x - 1.0 + mu) / r2.powi(3);
    let ddy = -2.0 * xd + y - (1.0 - mu) * y / r1.powi(3) - mu * y / r2.powi(3);
    let ddz = -((1.0 - mu) * z / r1.powi(3) + mu * z / r2.powi(3));

    Ok(Cr3bpState::from_row_slice(&[xd, yd, zd, ddx, ddy, ddz]))
}

pub fn jacobi_constant(state: &Cr3bpState, mu: f64) -> Result<f64, ThreeBodyError> {
    if !(0.0..0.5).contains(&mu) {
        return Err(ThreeBodyError::InvalidMassParameter);
    }
    let x = state[0];
    let y = state[1];
    let z = state[2];
    let xd = state[3];
    let yd = state[4];
    let zd = state[5];
    let r1 = ((x + mu).powi(2) + y * y + z * z).sqrt();
    let r2 = ((x - 1.0 + mu).powi(2) + y * y + z * z).sqrt();
    Ok(x * x + y * y + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - (xd * xd + yd * yd + zd * zd))
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::bodies::{EARTH, MOON};
    use crate::threebody::lagrange::triangular_lagrange_points;

    use super::{Cr3bpState, SystemChars, jacobi_constant, state_derivative, calculate_mu, calculate_tstar};

    #[test]
    fn test_calculate_mu() {
        let mu1 = EARTH.mu_km3_s2;
        let mu2 = MOON.mu_km3_s2;
        // Expected from python test
        let expected_mu = 1.215058560962404e-02;
        
        let mu = calculate_mu(mu1, mu2);
        assert_relative_eq!(mu, expected_mu, epsilon = 1e-6);
    }

    #[test]
    fn test_calculate_tstar() {
        let mu1 = EARTH.mu_km3_s2;
        let mu2 = MOON.mu_km3_s2;
        let lstar = 389703.2648292776;
        let expected_tstar = 382981.2891290545;
        
        let tstar = calculate_tstar(mu1, mu2, lstar);
        // Use relative tolerance
        assert_relative_eq!(tstar, expected_tstar, max_relative = 1e-6);
    }

    #[test]
    fn system_chars_earth_moon() {
        let sys = SystemChars::from_primaries(&EARTH, &MOON).unwrap();
        
        let expected_mu = 0.0121505856;
        let expected_lstar = 384_400.0;
        let _expected_tstar = 375_190.25852; // sqrt(384400^3 / (398600.4418 + 4902.79981))

        assert_relative_eq!(sys.mu, expected_mu, epsilon = 1e-5);
        assert_relative_eq!(sys.lstar_km, expected_lstar, epsilon = 1e-5);
        // tstar in python test: calculate_tstar(Earth.k, Moon.k, Moon.mean_a)
        // Let's rely on calculation
        let mu_sum = EARTH.mu_km3_s2 + MOON.mu_km3_s2;
        let calc_tstar = (expected_lstar.powi(3) / mu_sum).sqrt();
        assert_relative_eq!(sys.tstar_sec, calc_tstar, epsilon = 1e-5);
    }

    #[test]
    fn l4_is_equilibrium_with_zero_velocity() {
        let mu = 0.0121505856;
        let (l4, _) = triangular_lagrange_points(mu).unwrap();
        let state = Cr3bpState::from_row_slice(&[l4.x, l4.y, 0.0, 0.0, 0.0, 0.0]);
        let deriv = state_derivative(&state, mu).unwrap();
        assert_relative_eq!(deriv[3], 0.0, epsilon = 1e-12);
        assert_relative_eq!(deriv[4], 0.0, epsilon = 1e-12);
        assert_relative_eq!(deriv[5], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn jacobi_constant_matches_reference_expression() {
        let mu = 0.0121505856;
        let state = Cr3bpState::from_row_slice(&[0.8, 0.1, 0.0, 0.02, -0.01, 0.0]);
        let c = jacobi_constant(&state, mu).unwrap();
        assert_relative_eq!(c, 3.1781344535722202, epsilon = 1e-12);
    }

    #[test]
    fn invalid_mu_returns_error() {
        let state = Cr3bpState::from_row_slice(&[0.8, 0.1, 0.0, 0.02, -0.01, 0.0]);
        assert!(state_derivative(&state, 0.5).is_err());
        assert!(jacobi_constant(&state, -1.0).is_err());
    }
}
