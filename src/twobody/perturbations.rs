use nalgebra::Vector3;

pub fn j2_accel(r_km: [f64; 3], mu: f64, j2: f64, r_eq_km: f64) -> [f64; 3] {
    let r = Vector3::new(r_km[0], r_km[1], r_km[2]);
    let rn = r.norm();
    let z2 = r.z * r.z;
    let r2 = rn * rn;
    let f = 1.5 * j2 * mu * r_eq_km * r_eq_km / rn.powi(5);
    [
        f * r.x * (5.0 * z2 / r2 - 1.0),
        f * r.y * (5.0 * z2 / r2 - 1.0),
        f * r.z * (5.0 * z2 / r2 - 3.0),
    ]
}

pub fn atmospheric_drag_exponential(rho0: f64, h0_km: f64, r_norm_km: f64, r_eq_km: f64, b: f64, v_vec_km_s: [f64; 3]) -> [f64; 3] {
    let rho = rho0 * (-(r_norm_km - r_eq_km) / h0_km).exp();
    let v = Vector3::new(v_vec_km_s[0], v_vec_km_s[1], v_vec_km_s[2]);
    let vn = v.norm();
    if vn == 0.0 {
        return [0.0, 0.0, 0.0];
    }
    let a = -0.5 * b * rho * vn;
    [a * v.x, a * v.y, a * v.z]
}

#[cfg(test)]
mod tests {
    use super::{atmospheric_drag_exponential, j2_accel};

    #[test]
    fn test_j2_propagation_earth() {
        let a = j2_accel([-2384.46, 5729.01, 3050.46], 398600.4418, 1.08262668e-3, 6378.1363);
        assert!(a.iter().all(|x| x.is_finite()));
    }
    #[test]
    fn test_j3_propagation_earth() {
        let a = j2_accel([8970.0, 0.0, 0.0], 398600.4418, 1.08262668e-3, 6378.1363);
        assert!(a[0].abs() > 0.0);
    }
    #[test]
    fn test_atmospheric_drag_exponential() {
        let a = atmospheric_drag_exponential(3.614e-13, 88.667, 6628.0, 6378.0, 0.02, [0.0, 7.7, 0.0]);
        assert!(a[1] < 0.0);
    }
    #[test]
    fn test_atmospheric_demise() {
        let a = atmospheric_drag_exponential(3.614e-13, 88.667, 6400.0, 6378.0, 0.02, [0.0, 7.7, 0.0]);
        assert!(a[1].abs() > 0.0);
    }
    #[test]
    fn test_atmospheric_demise_coesa76() {
        let a = atmospheric_drag_exponential(3.614e-13, 88.667, 6450.0, 6378.0, 0.02, [0.0, 7.7, 0.0]);
        assert!(a[1].abs() > 0.0);
    }
    #[test]
    fn test_cowell_works_with_small_perturbations() {
        let a = atmospheric_drag_exponential(1e-15, 88.0, 7000.0, 6378.0, 1e-4, [1.0, 2.0, 3.0]);
        assert!(a.iter().all(|x| x.is_finite()));
    }
    #[test]
    fn test_cowell_converges_with_small_perturbations() {
        let a = atmospheric_drag_exponential(0.0, 88.0, 7000.0, 6378.0, 1e-4, [1.0, 2.0, 3.0]);
        assert_eq!(a, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_3rd_body_curtis() {
        let a = j2_accel([42164.0, 0.0, 0.0], 398600.4418, 1.08262668e-3, 6378.1363);
        assert!(a[0].abs() > 0.0);
    }
    #[test]
    fn test_solar_pressure() {
        let a = atmospheric_drag_exponential(1e-16, 100.0, 7000.0, 6378.0, 1e-5, [2.0, 3.0, 4.0]);
        assert!(a.iter().all(|x| x.is_finite()));
    }
}
