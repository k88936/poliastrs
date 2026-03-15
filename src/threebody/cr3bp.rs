use nalgebra::SVector;

use crate::threebody::ThreeBodyError;

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

    use crate::threebody::lagrange::triangular_lagrange_points;

    use super::{Cr3bpState, jacobi_constant, state_derivative};

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
