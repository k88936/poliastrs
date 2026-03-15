#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HohmannResult {
    pub dv1_km_s: f64,
    pub dv2_km_s: f64,
    pub total_dv_km_s: f64,
}

pub fn hohmann_coplanar_circular(mu_km3_s2: f64, r1_km: f64, r2_km: f64) -> HohmannResult {
    let v1 = (mu_km3_s2 / r1_km).sqrt();
    let v2 = (mu_km3_s2 / r2_km).sqrt();
    let a_t = 0.5 * (r1_km + r2_km);
    let vt1 = (mu_km3_s2 * (2.0 / r1_km - 1.0 / a_t)).sqrt();
    let vt2 = (mu_km3_s2 * (2.0 / r2_km - 1.0 / a_t)).sqrt();
    let dv1 = (vt1 - v1).abs();
    let dv2 = (v2 - vt2).abs();
    HohmannResult {
        dv1_km_s: dv1,
        dv2_km_s: dv2,
        total_dv_km_s: dv1 + dv2,
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::bodies::EARTH;

    use super::hohmann_coplanar_circular;

    #[test]
    fn hohmann_lowers_total_dv_vs_single_impulse_speed_gap() {
        let r1 = 7000.0;
        let r2 = 12000.0;
        let res = hohmann_coplanar_circular(EARTH.mu_km3_s2, r1, r2);
        let v1 = (EARTH.mu_km3_s2 / r1).sqrt();
        let v2 = (EARTH.mu_km3_s2 / r2).sqrt();
        assert!(res.total_dv_km_s < (v1 - v2).abs() + 0.5);
    }

    #[test]
    fn hohmann_is_symmetric_for_swap_of_radii() {
        let up = hohmann_coplanar_circular(EARTH.mu_km3_s2, 7000.0, 12000.0);
        let down = hohmann_coplanar_circular(EARTH.mu_km3_s2, 12000.0, 7000.0);
        assert_relative_eq!(up.total_dv_km_s, down.total_dv_km_s, epsilon = 1e-12);
    }
}
