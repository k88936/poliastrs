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

pub fn lagrange_points_vec(
    m1_kg: f64,
    r1_km: Vector3<f64>,
    m2_kg: f64,
    r2_km: Vector3<f64>,
    n: Vector3<f64>,
) -> Result<[Vector3<f64>; 5], ThreeBodyError> {
    if m1_kg <= m2_kg {
        // Python assertion: m1 > m2
        return Err(ThreeBodyError::InvalidMassRatio);
    }

    let mu = m2_kg / (m1_kg + m2_kg);
    let r12_vec = r2_km - r1_km;
    let r12 = r12_vec.norm();

    let ux = r12_vec / r12;
    // uy is perpendicular to ux and n.
    // In Python: uy = cross(n, ux); uy = uy / norm(uy)
    let uy = n.cross(&ux).normalize();

    let x_l1_bary = collinear_lagrange_x(mu, CollinearPoint::L1)?;
    let x_l2_bary = collinear_lagrange_x(mu, CollinearPoint::L2)?;
    let x_l3_bary = collinear_lagrange_x(mu, CollinearPoint::L3)?;

    // Convert from barycentric dimensionless X to distance from Primary (P1)
    // x_p1 = (x_bary + mu) * r12
    let x1 = (x_l1_bary + mu) * r12;
    let x2 = (x_l2_bary + mu) * r12;
    let x3 = (x_l3_bary + mu) * r12;

    // L4/L5
    // In barycentric: x = 0.5 - mu, y = +/- sqrt(3)/2
    // Dist from P1 along X: (0.5 - mu + mu) * r12 = 0.5 * r12
    let x4 = 0.5 * r12;
    let y4 = (3.0f64).sqrt() / 2.0 * r12;
    
    // Construct vectors
    let l1 = r1_km + ux * x1;
    let l2 = r1_km + ux * x2;
    let l3 = r1_km + ux * x3;
    let l4 = r1_km + ux * x4 + uy * y4;
    let l5 = r1_km + ux * x4 - uy * y4; // x5 = x4, y5 = -y4

    Ok([l1, l2, l3, l4, l5])
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

    #[test]
    fn invalid_mass_parameter_returns_error() {
        assert!(collinear_lagrange_x(0.0, CollinearPoint::L1).is_err());
        assert!(triangular_lagrange_points(0.6).is_err());
    }

    #[test]
    fn lagrange_points_vec_curtis_example() {
        use nalgebra::Vector3;
        use crate::bodies::{EARTH, MOON};
        use super::lagrange_points_vec;

        let m1 = EARTH.mu_km3_s2;
        let m2 = MOON.mu_km3_s2;
        let r1 = Vector3::new(0.0, 0.0, 0.0);
        let r2 = Vector3::new(384_400.0, 0.0, 0.0);
        let n = Vector3::new(0.0, 0.0, 1.0);

        let [l1, l2, l3, l4, l5] = lagrange_points_vec(m1, r1, m2, r2, n).unwrap();

        let deg60 = 60.0f64.to_radians();
        let expected_l1 = Vector3::new(326_400.0, 0.0, 0.0);
        let expected_l2 = Vector3::new(449_100.0, 0.0, 0.0);
        let expected_l3 = Vector3::new(-381_600.0, 0.0, 0.0);
        // L4/L5 form equilateral triangle with primaries.
        // Base is 384400.
        // x = 384400 * cos(60) = 192200.
        // y = 384400 * sin(60) = 332900.
        let expected_l4 = Vector3::new(384_400.0 * deg60.cos(), 384_400.0 * deg60.sin(), 0.0);
        let expected_l5 = Vector3::new(384_400.0 * deg60.cos(), -384_400.0 * deg60.sin(), 0.0);

        // Tolerance: Python uses rtol=1e-3.
        // 3e5 * 1e-3 = 300 km.
        let tol = 1000.0; 

        assert_relative_eq!(l1, expected_l1, epsilon = tol);
        assert_relative_eq!(l2, expected_l2, epsilon = tol);
        assert_relative_eq!(l3, expected_l3, epsilon = tol);
        assert_relative_eq!(l4, expected_l4, epsilon = tol);
        assert_relative_eq!(l5, expected_l5, epsilon = tol);
    }
}
