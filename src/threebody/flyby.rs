use nalgebra::Vector3;

use crate::threebody::ThreeBodyError;

pub fn compute_flyby(
    v_spacecraft_in_km_s: Vector3<f64>,
    v_body_km_s: Vector3<f64>,
    mu_body_km3_s2: f64,
    periapsis_radius_km: f64,
    b_plane_angle_rad: f64,
) -> Result<(Vector3<f64>, f64), ThreeBodyError> {
    if mu_body_km3_s2 <= 0.0 || periapsis_radius_km <= 0.0 {
        return Err(ThreeBodyError::InvalidFlybyInput);
    }

    let vinf_in = v_spacecraft_in_km_s - v_body_km_s;
    let vinf = vinf_in.norm();
    if vinf <= 0.0 {
        return Err(ThreeBodyError::InvalidFlybyInput);
    }

    let ecc = 1.0 + periapsis_radius_km * vinf * vinf / mu_body_km3_s2;
    let delta = 2.0 * (1.0 / ecc).asin();

    let s_hat = vinf_in / vinf;
    let z_hat = Vector3::new(0.0, 0.0, 1.0);
    let x_hat = Vector3::new(1.0, 0.0, 0.0);
    let ref_hat = if s_hat.cross(&z_hat).norm() > 1e-14 {
        z_hat
    } else {
        x_hat
    };
    let t_hat = s_hat.cross(&ref_hat).normalize();
    let r_hat = t_hat.cross(&s_hat).normalize();
    // Curtis 3rd Ed, Eq 8.35: B = R cos(theta) + T sin(theta)
    let b_hat = r_hat * b_plane_angle_rad.cos() + t_hat * b_plane_angle_rad.sin();
    let n_hat = b_hat.cross(&s_hat).normalize();
    let vinf_out_hat = s_hat * delta.cos() + n_hat * delta.sin();

    Ok((v_body_km_s + vinf * vinf_out_hat, delta))
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    use crate::bodies::VENUS;
    use super::compute_flyby;

    #[test]
    fn flyby_venus_curtis_example() {
        // Data from Curtis, Example 8.6
        let v_sc = Vector3::new(37.51, 2.782, 0.0);
        let v_body = Vector3::new(35.02, 0.0, 0.0);
        let mu = VENUS.mu_km3_s2;
        let r_venus = VENUS.mean_radius_km;
        let h_p = 300.0;
        let rp = r_venus + h_p;
        let theta = 0.0; // b_plane_angle

        let (v_out, delta) = compute_flyby(v_sc, v_body, mu, rp, theta).unwrap();

        let expected_v_out = Vector3::new(31.73, 1.766, 0.0);
        let expected_delta_deg = 103.6;

        assert_relative_eq!(v_out, expected_v_out, epsilon = 1e-1); // 1e-1 km/s tolerance as in Python (rtol=1e-3 ~ 0.03 km/s)
        assert_relative_eq!(delta.to_degrees(), expected_delta_deg, epsilon = 1e-1);
    }

    #[test]
    fn flyby_preserves_hyperbolic_excess_speed_in_planet_frame() {
        let v_sc = Vector3::new(12.0, 0.5, 0.0);
        let v_body = Vector3::new(8.0, 0.0, 0.0);
        let mu = 398600.4418;
        let rp = 7000.0;

        let (v_out, _delta) = compute_flyby(v_sc, v_body, mu, rp, 0.0).unwrap();
        let vinf_in = (v_sc - v_body).norm();
        let vinf_out = (v_out - v_body).norm();
        assert_relative_eq!(vinf_in, vinf_out, max_relative = 1e-12);
    }

    #[test]
    fn flyby_turn_angle_matches_hyperbolic_relation() {
        let v_sc = Vector3::new(12.0, 0.5, 0.0);
        let v_body = Vector3::new(8.0, 0.0, 0.0);
        let mu = 398600.4418;
        let rp = 7000.0;

        let (_v_out, delta) = compute_flyby(v_sc, v_body, mu, rp, 0.0).unwrap();
        let vinf = (v_sc - v_body).norm();
        let ecc = 1.0 + rp * vinf * vinf / mu;
        let expected = 2.0 * (1.0 / ecc).asin();
        assert_relative_eq!(delta, expected, max_relative = 1e-12);
    }

    #[test]
    fn invalid_inputs_return_error() {
        let v = Vector3::new(1.0, 0.0, 0.0);
        assert!(compute_flyby(v, v, 398600.0, 7000.0, 0.0).is_err());
        assert!(compute_flyby(v, Vector3::zeros(), -1.0, 7000.0, 0.0).is_err());
    }
}
