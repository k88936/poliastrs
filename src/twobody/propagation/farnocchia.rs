use nalgebra::Vector3;

use crate::core::stumpff::{stumpff_c2, stumpff_c3};
use crate::twobody::states::CartesianState;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PropagationError {
    DidNotConverge,
}

pub fn propagate_two_body(
    mu_km3_s2: f64,
    state: &CartesianState,
    dt_seconds: f64,
) -> Result<CartesianState, PropagationError> {
    if dt_seconds == 0.0 {
        return Ok(*state);
    }

    let r0 = state.r_km;
    let v0 = state.v_km_s;
    let r0_norm = r0.norm();
    let vr0 = r0.dot(&v0) / r0_norm;
    let alpha = 2.0 / r0_norm - v0.norm_squared() / mu_km3_s2;
    let sqrt_mu = mu_km3_s2.sqrt();

    // Handle multi-revolution for elliptic orbits to improve convergence
    let mut dt = dt_seconds;
    if alpha > 1e-10 {
        let a = 1.0 / alpha;
        let period = 2.0 * std::f64::consts::PI * (a * a * a / mu_km3_s2).sqrt();
        dt %= period;
        // Keep dt in [-T/2, T/2]
        if dt > period / 2.0 {
            dt -= period;
        } else if dt < -period / 2.0 {
            dt += period;
        }
    }

    let mut x = if alpha.abs() > 1e-10 {
        sqrt_mu * dt * alpha.abs()
    } else {
        sqrt_mu * dt / r0_norm
    };

    let max_iter = 200;
    let tol = 1e-10;

    for _ in 0..max_iter {
        let z = alpha * x * x;
        let c2 = stumpff_c2(z);
        let c3 = stumpff_c3(z);

        let f = r0_norm * vr0 / sqrt_mu * x * x * c2
            + (1.0 - alpha * r0_norm) * x * x * x * c3
            + r0_norm * x
            - sqrt_mu * dt; // Use modified dt

        let fp = r0_norm * vr0 / sqrt_mu * x * (1.0 - z * c3)
            + (1.0 - alpha * r0_norm) * x * x * c2
            + r0_norm;

        let dx = f / fp;
        x -= dx;
        if dx.abs() < tol {
            let z = alpha * x * x;
            let c2 = stumpff_c2(z);
            let c3 = stumpff_c3(z);

            let f_lagrange = 1.0 - (x * x / r0_norm) * c2;
            let g_lagrange = dt - (x * x * x / sqrt_mu) * c3; // Use modified dt
            let r: Vector3<f64> = f_lagrange * r0 + g_lagrange * v0;
            let r_norm = r.norm();

            let fdot = sqrt_mu / (r_norm * r0_norm) * (z * c3 - 1.0) * x;
            let gdot = 1.0 - (x * x / r_norm) * c2;
            let v = fdot * r0 + gdot * v0;

            return Ok(CartesianState::new(r, v));
        }
    }

    Err(PropagationError::DidNotConverge)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use crate::twobody::states::CartesianState;
    use nalgebra::Vector3;

    #[test]
    fn test_propagation_one_period_returns_same_state() {
        let mu = 398600.4418; // Earth
        let r = Vector3::new(-6045.0, -3490.0, 2500.0);
        let v = Vector3::new(-3.457, 6.618, 2.533);
        let state = CartesianState::new(r, v);
        
        // Calculate period
        let r_norm = r.norm();
        let v_norm_sq = v.norm_squared();
        let energy = v_norm_sq / 2.0 - mu / r_norm;
        let a = -mu / (2.0 * energy);
        let period = 2.0 * std::f64::consts::PI * (a.powi(3) / mu).sqrt();
        
        let new_state = propagate_two_body(mu, &state, period).unwrap();
        
        assert_relative_eq!(new_state.r_km, state.r_km, epsilon = 1e-6);
        assert_relative_eq!(new_state.v_km_s, state.v_km_s, epsilon = 1e-6);
    }
}
