use nalgebra::Vector3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LambertError {
    MultiRevolutionUnsupported,
    PhaseAngle180Deg,
    MaximumIterationsReached,
}

pub fn lambert(
    mu_km3_s2: f64,
    r0_km: [f64; 3],
    r_km: [f64; 3],
    tof_seconds: f64,
    m: u32,
    prograde: bool,
    _lowpath: bool,
    numiter: usize,
    rtol: f64,
) -> Result<([f64; 3], [f64; 3]), LambertError> {
    if m > 0 {
        return Err(LambertError::MultiRevolutionUnsupported);
    }

    let r0 = Vector3::new(r0_km[0], r0_km[1], r0_km[2]);
    let rf = Vector3::new(r_km[0], r_km[1], r_km[2]);
    let norm_r0 = r0.norm();
    let norm_rf = rf.norm();
    let cos_dnu = r0.dot(&rf) / (norm_r0 * norm_rf);
    let tm = if prograde { 1.0 } else { -1.0 };
    let a = tm * (norm_r0 * norm_rf * (1.0 + cos_dnu)).sqrt();

    if a.abs() < 1e-15 {
        return Err(LambertError::PhaseAngle180Deg);
    }

    let mut psi = 0.0f64;
    let mut psi_low = -4.0 * std::f64::consts::PI.powi(2);
    let mut psi_up = 4.0 * std::f64::consts::PI.powi(2);

    let norm_sum = norm_r0 + norm_rf;
    let norm_prod = norm_r0 * norm_rf;

    for _ in 0..numiter {
        let c2_psi = c2(psi);
        let c3_psi = c3(psi);
        let mut y = norm_sum + a * (psi * c3_psi - 1.0) / c2_psi.sqrt();

        if a > 0.0 {
            while y < 0.0 {
                psi_low = psi;
                psi = 0.8 * (1.0 / c3(psi)) * (1.0 - norm_prod * c2(psi).sqrt() / a);
                let c2_new = c2(psi);
                let c3_new = c3(psi);
                y = norm_sum + a * (psi * c3_new - 1.0) / c2_new.sqrt();
            }
        }

        let xi = (y / c2(psi)).sqrt();
        let tof_new = (xi.powi(3) * c3(psi) + a * y.sqrt()) / mu_km3_s2.sqrt();

        if ((tof_new - tof_seconds) / tof_seconds).abs() < rtol {
            let f = 1.0 - y / norm_r0;
            let g = a * (y / mu_km3_s2).sqrt();
            let gdot = 1.0 - y / norm_rf;
            let v0 = (rf - f * r0) / g;
            let vf = (gdot * rf - r0) / g;
            return Ok(([v0.x, v0.y, v0.z], [vf.x, vf.y, vf.z]));
        }

        if tof_new <= tof_seconds {
            psi_low = psi;
        } else {
            psi_up = psi;
        }
        psi = 0.5 * (psi_up + psi_low);
    }

    Err(LambertError::MaximumIterationsReached)
}

fn c2(psi: f64) -> f64 {
    if psi > 1e-8 {
        (1.0 - psi.sqrt().cos()) / psi
    } else if psi < -1e-8 {
        (1.0 - (-psi).sqrt().cosh()) / psi
    } else {
        0.5
    }
}

fn c3(psi: f64) -> f64 {
    if psi > 1e-8 {
        (psi.sqrt() - psi.sqrt().sin()) / psi.powf(1.5)
    } else if psi < -1e-8 {
        (((-psi).sqrt()).sinh() - (-psi).sqrt()) / (-psi).powf(1.5)
    } else {
        1.0 / 6.0
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::bodies::EARTH;

    use super::lambert;

    #[test]
    fn lambert_vallado75_reference_case() {
        let r0 = [15945.34, 0.0, 0.0];
        let r = [12214.83399, 10249.46731, 0.0];
        let tof = 76.0 * 60.0;

        let (va, vb) = lambert(EARTH.mu_km3_s2, r0, r, tof, 0, true, true, 35, 1e-8).unwrap();
        let expected_va = [2.058925, 2.915956, 0.0];
        let expected_vb = [-3.451569, 0.910301, 0.0];

        for i in 0..3 {
            assert_relative_eq!(va[i], expected_va[i], max_relative = 1e-5);
            assert_relative_eq!(vb[i], expected_vb[i], max_relative = 1e-4);
        }
    }
}
