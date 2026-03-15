use nalgebra::Vector3;
use std::f64::consts::PI;
use crate::core::hyper::hyp2f1b;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LambertError {
    CollinearVectors,
    NoFeasibleSolution,
    ConvergenceFailed,
    DerivativeZero,
    TimeOfFlightNegative,
}

pub fn izzo(
    mu: f64,
    r1: Vector3<f64>,
    r2: Vector3<f64>,
    tof: f64,
    revolutions: i32,
    prograde: bool,
    low_path: bool,
    num_iter: usize,
    rtol: f64,
) -> Result<(Vector3<f64>, Vector3<f64>), LambertError> {
    if tof <= 0.0 {
        return Err(LambertError::TimeOfFlightNegative);
    }
    
    // Check collinearity
    let c_vec = r2 - r1;
    let cross_r1_r2 = r1.cross(&r2);
    if cross_r1_r2.norm() < 1e-12 {
        return Err(LambertError::CollinearVectors);
    }

    let r1_norm = r1.norm();
    let r2_norm = r2.norm();
    let c_norm = c_vec.norm();

    let s = (r1_norm + r2_norm + c_norm) * 0.5;

    let i_r1 = r1 / r1_norm;
    let i_r2 = r2 / r2_norm;
    let i_h = i_r1.cross(&i_r2).normalize();

    let mut lambda = (1.0 - (c_norm / s).min(1.0)).sqrt();

    let (i_t1, i_t2) = if i_h.z < 0.0 {
        lambda = -lambda;
        (i_r1.cross(&i_h), i_r2.cross(&i_h))
    } else {
        (i_h.cross(&i_r1), i_h.cross(&i_r2))
    };

    let (lambda, i_t1, i_t2) = if prograde {
        (lambda, i_t1, i_t2)
    } else {
        (-lambda, -i_t1, -i_t2)
    };

    let t_nd = (2.0 * mu / s.powi(3)).sqrt() * tof;

    let (x, y) = find_xy(lambda, t_nd, revolutions, num_iter, low_path, rtol)?;

    let gamma = (mu * s / 2.0).sqrt();
    let rho = (r1_norm - r2_norm) / c_norm;
    let sigma = (1.0 - rho * rho).sqrt();

    let (v_r1, v_r2, v_t1, v_t2) = reconstruct(x, y, r1_norm, r2_norm, lambda, gamma, rho, sigma);

    let v1 = i_r1 * v_r1 + i_t1 * v_t1;
    let v2 = i_r2 * v_r2 + i_t2 * v_t2;

    Ok((v1, v2))
}

fn reconstruct(
    x: f64,
    y: f64,
    r1: f64,
    r2: f64,
    lambda: f64,
    gamma: f64,
    rho: f64,
    sigma: f64,
) -> (f64, f64, f64, f64) {
    let v_r1 = gamma * ((lambda * y - x) - rho * (lambda * y + x)) / r1;
    let v_r2 = -gamma * ((lambda * y - x) + rho * (lambda * y + x)) / r2;
    let v_t1 = gamma * sigma * (y + lambda * x) / r1;
    let v_t2 = gamma * sigma * (y + lambda * x) / r2;
    (v_r1, v_r2, v_t1, v_t2)
}

fn find_xy(
    lambda: f64,
    t_nd: f64,
    m: i32,
    num_iter: usize,
    low_path: bool,
    rtol: f64,
) -> Result<(f64, f64), LambertError> {
    if lambda.abs() >= 1.0 {
        // Technically this should be covered by c_norm / s logic
    }
    
    let m_max = (t_nd / PI).floor() as i32;
    let t_00 = lambda.acos() + lambda * (1.0 - lambda * lambda).sqrt();

    let mut m_max_refined = m_max;
    if t_nd < t_00 + (m_max as f64) * PI && m_max > 0 {
        let (_, t_min) = compute_t_min(lambda, m_max, num_iter, rtol)?;
        if t_nd < t_min {
            m_max_refined -= 1;
        }
    }

    if m > m_max_refined {
        return Err(LambertError::NoFeasibleSolution);
    }

    let x_0 = initial_guess(t_nd, lambda, m, low_path);
    let x = householder(x_0, t_nd, lambda, m, rtol, num_iter)?;
    let y = compute_y(x, lambda);

    Ok((x, y))
}

fn compute_y(x: f64, lambda: f64) -> f64 {
    (1.0 - lambda * lambda * (1.0 - x * x)).sqrt()
}

fn compute_psi(x: f64, y: f64, lambda: f64) -> f64 {
    if x >= -1.0 && x < 1.0 {
        (x * y + lambda * (1.0 - x * x)).acos()
    } else if x > 1.0 {
        ((y - x * lambda) * (x * x - 1.0).sqrt()).asinh()
    } else {
        0.0
    }
}

fn tof_equation_y(x: f64, y: f64, t0: f64, lambda: f64, m: i32) -> f64 {
    let t_val = if m == 0 && x > 0.6f64.sqrt() && x < 1.4f64.sqrt() {
        let eta = y - lambda * x;
        let s_1 = (1.0 - lambda - x * eta) * 0.5;
        let q = 4.0 / 3.0 * hyp2f1b(s_1);
        (eta.powi(3) * q + 4.0 * lambda * eta) * 0.5
    } else {
        let psi = compute_psi(x, y, lambda);
        let num = psi + (m as f64) * PI;
        let den = (1.0 - x * x).abs().sqrt();
        (num / den - x + lambda * y) / (1.0 - x * x)
    };
    t_val - t0
}

fn tof_equation_p(x: f64, y: f64, t_val: f64, lambda: f64) -> f64 {
    (3.0 * t_val * x - 2.0 + 2.0 * lambda.powi(3) * x / y) / (1.0 - x * x)
}

fn tof_equation_p2(x: f64, y: f64, t_val: f64, dt: f64, lambda: f64) -> f64 {
    (3.0 * t_val + 5.0 * x * dt + 2.0 * (1.0 - lambda * lambda) * lambda.powi(3) / y.powi(3))
        / (1.0 - x * x)
}

fn tof_equation_p3(x: f64, y: f64, _t_val: f64, dt: f64, ddt: f64, lambda: f64) -> f64 {
    (7.0 * x * ddt + 8.0 * dt - 6.0 * (1.0 - lambda * lambda) * lambda.powi(5) * x / y.powi(5))
        / (1.0 - x * x)
}

fn householder(
    mut p0: f64,
    t0: f64,
    lambda: f64,
    m: i32,
    tol: f64,
    max_iter: usize,
) -> Result<f64, LambertError> {
    for _ in 0..max_iter {
        let y = compute_y(p0, lambda);
        let fval = tof_equation_y(p0, y, t0, lambda, m);
        let t_val = fval + t0;
        let fder = tof_equation_p(p0, y, t_val, lambda);
        let fder2 = tof_equation_p2(p0, y, t_val, fder, lambda);
        let fder3 = tof_equation_p3(p0, y, t_val, fder, fder2, lambda);

        // Householder step (quartic)
        // p = p0 - fval * (
        //     (fder**2 - fval * fder2 / 2)
        //     / (fder * (fder**2 - fval * fder2) + fder3 * fval**2 / 6)
        // )
        let num = fder.powi(2) - fval * fder2 / 2.0;
        let den = fder * (fder.powi(2) - fval * fder2) + fder3 * fval.powi(2) / 6.0;
        
        let p = p0 - fval * (num / den);

        if (p - p0).abs() < tol {
            return Ok(p);
        }
        p0 = p;
    }
    Err(LambertError::ConvergenceFailed)
}

fn halley(
    mut p0: f64,
    t0: f64,
    lambda: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64, LambertError> {
    for _ in 0..max_iter {
        let y = compute_y(p0, lambda);
        let fder = tof_equation_p(p0, y, t0, lambda);
        let fder2 = tof_equation_p2(p0, y, t0, fder, lambda);
        if fder2 == 0.0 {
            return Err(LambertError::DerivativeZero);
        }
        let fder3 = tof_equation_p3(p0, y, t0, fder, fder2, lambda);

        // Halley step (cubic)
        // p = p0 - 2 * fder * fder2 / (2 * fder2**2 - fder * fder3)
        let p = p0 - 2.0 * fder * fder2 / (2.0 * fder2.powi(2) - fder * fder3);

        if (p - p0).abs() < tol {
            return Ok(p);
        }
        p0 = p;
    }
    Err(LambertError::ConvergenceFailed)
}

fn compute_t_min(
    lambda: f64,
    m: i32,
    num_iter: usize,
    rtol: f64,
) -> Result<(f64, f64), LambertError> {
    if (lambda - 1.0).abs() < 1e-12 {
        let x_t_min = 0.0;
        let t_min = tof_equation_y(x_t_min, compute_y(x_t_min, lambda), 0.0, lambda, m);
        Ok((x_t_min, t_min))
    } else {
        if m == 0 {
            Ok((f64::INFINITY, 0.0))
        } else {
            let x_i = 0.1;
            let t_i = tof_equation_y(x_i, compute_y(x_i, lambda), 0.0, lambda, m);
            // Wait, t_i is calculated as fval (T - T0) where T0=0, so it's T.
            // Halley finds minimum of T(x). We search for T'(x) = 0?
            // No, the python code calls _halley(x_i, t_i, ll, rtol, numiter).
            // But `halley` minimizes?
            // `halley` in python finds a *minimum* of time of flight equation.
            // Wait, the halley implementation iterates to solve f(p) = 0?
            // The Python doc says "Find a minimum of time of flight equation".
            // But it iterates `p = p0 - 2 * fder * fder2 / ...`.
            // This is finding a root of f'(x)=0? Or f(x)=0?
            // Standard Halley is for root finding.
            // If we want minimum of T(x), we want root of T'(x).
            // In python `_halley`:
            // fder = _tof_equation_p(...)
            // fder2 = ...
            // Wait, if it's finding minimum, it should use derivatives of T.
            // The `halley` function uses `_tof_equation_p` (1st derivative), `p2` (2nd), `p3` (3rd).
            // So it finds root of `_tof_equation_p`? No, standard Halley finds root of f(x).
            // If f(x) = T(x), then it finds T(x)=0.
            // If we want minimum T, we need T'(x)=0.
            // Let's re-read Python code.
            // `_halley(p0, T0, ll, tol, maxiter)`
            // `y = _compute_y(p0, ll)`
            // `fder = _tof_equation_p(p0, y, T0, ll)`
            // It calculates fder.
            // Then fder2, fder3.
            // Then update p using fder, fder2, fder3.
            // This looks like finding root of `fder`?
            // If we are looking for root of `g(x) = T'(x)`, we need g, g', g''.
            // `fder` is T'. `fder2` is T''. `fder3` is T'''.
            // Yes, so it finds where T'(x) = 0. Correct.
            // So my `halley` function matches Python logic which uses derivatives starting from `tof_equation_p`.
            
            let x_t_min = halley(x_i, t_i, lambda, rtol, num_iter)?;
            let t_min = tof_equation_y(x_t_min, compute_y(x_t_min, lambda), 0.0, lambda, m);
            Ok((x_t_min, t_min))
        }
    }
}

fn initial_guess(t_nd: f64, lambda: f64, m: i32, low_path: bool) -> f64 {
    if m == 0 {
        let t_0 = lambda.acos() + lambda * (1.0 - lambda * lambda).sqrt() + (m as f64) * PI;
        let t_1 = 2.0 * (1.0 - lambda.powi(3)) / 3.0;
        
        if t_nd >= t_0 {
            (t_0 / t_nd).powf(2.0 / 3.0) - 1.0
        } else if t_nd < t_1 {
            5.0 / 2.0 * t_1 / t_nd * (t_1 - t_nd) / (1.0 - lambda.powi(5)) + 1.0
        } else {
            ((2.0f64.ln() * (t_nd / t_0).ln() / (t_1 / t_0).ln()).exp()) - 1.0
        }
    } else {
        let x_0l = (((m as f64 * PI + PI) / (8.0 * t_nd)).powf(2.0 / 3.0) - 1.0) /
                   (((m as f64 * PI + PI) / (8.0 * t_nd)).powf(2.0 / 3.0) + 1.0);
        let x_0r = (((8.0 * t_nd) / (m as f64 * PI)).powf(2.0 / 3.0) - 1.0) /
                   (((8.0 * t_nd) / (m as f64 * PI)).powf(2.0 / 3.0) + 1.0);

        if low_path {
            x_0l.max(x_0r)
        } else {
            x_0l.min(x_0r)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vallado75() {
        let mu = 398600.4418;
        let r0 = Vector3::new(15945.34, 0.0, 0.0);
        let r = Vector3::new(12214.83399, 10249.46731, 0.0);
        let tof = 76.0 * 60.0;

        let expected_va = Vector3::new(2.058925, 2.915956, 0.0);
        let expected_vb = Vector3::new(-3.451569, 0.910301, 0.0);

        // prograde = true
        let (va, vb) = izzo(mu, r0, r, tof, 0, true, false, 35, 1e-8).unwrap();

        assert_relative_eq!(va, expected_va, max_relative = 1e-5);
        assert_relative_eq!(vb, expected_vb, max_relative = 1e-4);
    }

    #[test]
    fn test_curtis52() {
        let mu = 398600.4418;
        let r0 = Vector3::new(5000.0, 10000.0, 2100.0);
        let r = Vector3::new(-14600.0, 2500.0, 7000.0);
        let tof = 3600.0;

        let expected_va = Vector3::new(-5.9925, 1.9254, 3.2456);
        let expected_vb = Vector3::new(-3.3125, -4.1966, -0.38529);

        // prograde = true
        let (va, vb) = izzo(mu, r0, r, tof, 0, true, false, 35, 1e-8).unwrap();

        assert_relative_eq!(va, expected_va, max_relative = 1e-4);
        assert_relative_eq!(vb, expected_vb, max_relative = 1e-4);
    }

    #[test]
    fn test_curtis53() {
        let mu = 398600.4418;
        let r0 = Vector3::new(273378.0, 0.0, 0.0);
        let r = Vector3::new(145820.0, 12758.0, 0.0);
        let tof = 13.5 * 3600.0;

        let expected_va = Vector3::new(-2.4356, 0.26741, 0.0);

        // prograde = true
        let (va, _) = izzo(mu, r0, r, tof, 0, true, false, 100, 1e-8).unwrap();
        assert_relative_eq!(va, expected_va, max_relative = 1e-4);
    }

    #[test]
    fn test_molniya_der_zero_full_revolution() {
        let mu = 398600.4418;
        let r0 = Vector3::new(22592.145603, -1599.915239, -19783.950506);
        let r = Vector3::new(1922.067697, 4054.157051, -8925.727465);
        let tof = 10.0 * 3600.0;

        let expected_va = Vector3::new(2.000652697, 0.387688615, -2.666947760);
        let expected_vb = Vector3::new(-3.79246619, -1.77707641, 6.856814395);

        // prograde = true
        let (va, vb) = izzo(mu, r0, r, tof, 0, true, false, 35, 1e-8).unwrap();

        assert_relative_eq!(va, expected_va, max_relative = 1e-5);
        assert_relative_eq!(vb, expected_vb, max_relative = 1e-5);
    }
}
