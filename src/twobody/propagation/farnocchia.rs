use nalgebra::Vector3;

use crate::twobody::states::CartesianState;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PropagationError {
    DidNotConverge,
}

fn stumpff_c2(z: f64) -> f64 {
    if z > 1e-8 {
        let s = z.sqrt();
        (1.0 - s.cos()) / z
    } else if z < -1e-8 {
        let s = (-z).sqrt();
        (s.cosh() - 1.0) / (-z)
    } else {
        0.5
    }
}

fn stumpff_c3(z: f64) -> f64 {
    if z > 1e-8 {
        let s = z.sqrt();
        (s - s.sin()) / (s * s * s)
    } else if z < -1e-8 {
        let s = (-z).sqrt();
        (s.sinh() - s) / (s * s * s)
    } else {
        1.0 / 6.0
    }
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

    let mut x = if alpha.abs() > 1e-10 {
        sqrt_mu * dt_seconds * alpha.abs()
    } else {
        sqrt_mu * dt_seconds / r0_norm
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
            - sqrt_mu * dt_seconds;

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
            let g_lagrange = dt_seconds - (x * x * x / sqrt_mu) * c3;
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
