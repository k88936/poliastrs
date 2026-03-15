use nalgebra::Vector3;

use crate::threebody::ThreeBodyError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CollinearPoint {
    L1,
    L2,
    L3,
}

pub fn collinear_lagrange_x(mu: f64, point: CollinearPoint) -> Result<f64, ThreeBodyError> {
    if !(0.0..0.5).contains(&mu) {
        return Err(ThreeBodyError::InvalidMassParameter);
    }

    let mut x = match point {
        CollinearPoint::L1 => 1.0 - (mu / 3.0).cbrt() - mu,
        CollinearPoint::L2 => 1.0 + (mu / 3.0).cbrt() - mu,
        CollinearPoint::L3 => -1.0 - (5.0 / 12.0) * mu,
    };

    for _ in 0..100 {
        let f = collinear_equation(x, mu);
        let df = collinear_equation_derivative(x, mu);
        let dx = f / df;
        x -= dx;
        if dx.abs() < 1e-13 {
            return Ok(x);
        }
    }

    Err(ThreeBodyError::NonConvergentSolver)
}

pub fn triangular_lagrange_points(mu: f64) -> Result<(Vector3<f64>, Vector3<f64>), ThreeBodyError> {
    if !(0.0..0.5).contains(&mu) {
        return Err(ThreeBodyError::InvalidMassParameter);
    }
    let x = 0.5 - mu;
    let y = (3.0_f64).sqrt() / 2.0;
    Ok((Vector3::new(x, y, 0.0), Vector3::new(x, -y, 0.0)))
}

pub(crate) fn collinear_equation(x: f64, mu: f64) -> f64 {
    let t1 = (1.0 - mu) * (x + mu) / (x + mu).abs().powi(3);
    let t2 = mu * (x - 1.0 + mu) / (x - 1.0 + mu).abs().powi(3);
    x - t1 - t2
}

fn collinear_equation_derivative(x: f64, mu: f64) -> f64 {
    let d1 = (1.0 - mu) / (x + mu).abs().powi(3);
    let d2 = mu / (x - 1.0 + mu).abs().powi(3);
    1.0 + 2.0 * (d1 + d2)
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::{
        CollinearPoint, collinear_equation, collinear_lagrange_x, triangular_lagrange_points,
    };

    #[test]
    fn triangular_points_are_equidistant_to_primaries() {
        let mu = 0.0121505856;
        let (l4, l5) = triangular_lagrange_points(mu).unwrap();
        for point in [l4, l5] {
            let r1 = ((point.x + mu).powi(2) + point.y.powi(2)).sqrt();
            let r2 = ((point.x - 1.0 + mu).powi(2) + point.y.powi(2)).sqrt();
            assert_relative_eq!(r1, 1.0, epsilon = 1e-12);
            assert_relative_eq!(r2, 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn collinear_points_satisfy_equilibrium_equation() {
        let mu = 0.0121505856;
        for p in [CollinearPoint::L1, CollinearPoint::L2, CollinearPoint::L3] {
            let x = collinear_lagrange_x(mu, p).unwrap();
            assert_relative_eq!(collinear_equation(x, mu), 0.0, epsilon = 1e-11);
        }
    }
}
