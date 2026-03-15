use std::f64::consts::PI;

/// Calculates the minimum and maximum values of ground-range angles.
///
/// # Arguments
///
/// * `h_km` - Altitude over surface.
/// * `eta_fov_rad` - Angle of the total area that a sensor can observe.
/// * `eta_center_rad` - Center boresight angle.
/// * `r_km` - Attractor equatorial radius.
///
/// # Returns
///
/// * `(lambda_min_rad, lambda_max_rad)`
pub fn min_and_max_ground_range(h_km: f64, eta_fov_rad: f64, eta_center_rad: f64, r_km: f64) -> (f64, f64) {
    let r_sat = r_km + h_km;
    let eta_max = eta_center_rad + eta_fov_rad / 2.0;
    let eta_min = eta_center_rad - eta_fov_rad / 2.0;
    
    // Gamma max calculation
    let sin_val_max = (r_sat * eta_max.sin()) / r_km;
    // Clamping to avoid domain error if close to 1
    let clamped_sin_max = sin_val_max.clamp(-1.0, 1.0);
    let mut gamma_max = clamped_sin_max.asin();
    
    // Logic from python: if abs(gamma) <= pi/2, gamma = pi - gamma.
    // Since asin returns [-pi/2, pi/2], this is always true.
    gamma_max = PI - gamma_max;

    // Gamma min calculation
    let sin_val_min = (r_sat * eta_min.sin()) / r_km;
    let clamped_sin_min = sin_val_min.clamp(-1.0, 1.0);
    let mut gamma_min = clamped_sin_min.asin();
    
    gamma_min = PI - gamma_min;

    // Maximum and minimum slant ranges
    let rho_max = r_km * gamma_max.cos() + r_sat * eta_max.cos();
    let rho_min = r_km * gamma_min.cos() + r_sat * eta_min.cos();
    
    // Calculate lambdas
    let lambda_max = ((rho_max * eta_max.sin()) / r_km).clamp(-1.0, 1.0).asin();
    let lambda_min = ((rho_min * eta_min.sin()) / r_km).clamp(-1.0, 1.0).asin();
    
    (lambda_min, lambda_max)
}

/// Calculates the difference in ground-range angles.
///
/// # Arguments
///
/// * `h_km` - Altitude over surface.
/// * `eta_fov_rad` - Angle of the total area that a sensor can observe.
/// * `eta_center_rad` - Center boresight angle.
/// * `beta_rad` - Phase angle, azimuth.
/// * `phi_nadir_rad` - Latitude of nadir point.
/// * `lambda_nadir_rad` - Longitude of nadir point.
/// * `r_km` - Attractor equatorial radius.
///
/// # Returns
///
/// * `Result<(delta_lambda_rad, phi_tgt_rad, lambda_tgt_rad), &'static str>`
pub fn ground_range_diff_at_azimuth(
    h_km: f64,
    eta_fov_rad: f64,
    eta_center_rad: f64,
    beta_rad: f64,
    phi_nadir_rad: f64,
    lambda_nadir_rad: f64,
    r_km: f64
) -> Result<(f64, f64, f64), &'static str> {
    if !(0.0..PI).contains(&beta_rad) {
        return Err("beta must be between 0 and PI radians");
    }

    let r_sat = r_km + h_km;
    
    // Calculate Gamma for center (boresight)
    let sin_val = (r_sat * eta_center_rad.sin()) / r_km;
    let clamped_sin = sin_val.clamp(-1.0, 1.0);
    let mut gamma = clamped_sin.asin();
    
    // Always true for asin output
    gamma = PI - gamma;
    
    let rho = r_km * gamma.cos() + r_sat * eta_center_rad.cos();
    
    // Calculate Lambda (capital Lambda in Vallado)
    let capital_lambda = ((rho * eta_center_rad.sin()) / r_km).clamp(-1.0, 1.0).asin();
    
    // Calculate target latitude
    let sin_phi_tgt = beta_rad.cos() * phi_nadir_rad.cos() * capital_lambda.sin() + phi_nadir_rad.sin() * capital_lambda.cos();
    let phi_tgt = sin_phi_tgt.clamp(-1.0, 1.0).asin();
    
    // Calculate delta capital lambda (difference in longitude related angle?)
    // Vallado formula: sin(d_lambda) = sin(beta) * sin(Lambda) / cos(phi_tgt)
    let sin_delta_capital_lambda = (beta_rad.sin() * capital_lambda.sin()) / phi_tgt.cos();
    let delta_capital_lambda = sin_delta_capital_lambda.clamp(-1.0, 1.0).asin();
    
    let lambda_tgt = lambda_nadir_rad + delta_capital_lambda;
    
    // Calculate delta_lambda (difference in ground range from boresight?)
    // Note: Python code returns delta_lambda as (Lambda_max - Lambda_min) / 2
    // Which is the half-width of the FOV ground range?
    // Wait, check python code carefully.
    // delta_λ = (Λ_max - Λ_min) / 2
    // Yes.
    
    let (lambda_min, lambda_max) = min_and_max_ground_range(h_km, eta_fov_rad, eta_center_rad, r_km);
    let delta_lambda = (lambda_max - lambda_min) / 2.0;
    
    Ok((delta_lambda, phi_tgt, lambda_tgt))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use crate::bodies::EARTH;

    #[test]
    fn test_max_and_min_ground_range() {
        let altitude = 800.0; // km
        let fov = 25.0_f64.to_radians();
        let boresight = 40.0_f64.to_radians();
        let r_earth = EARTH.equatorial_radius_km;
        
        let (lat_lon_min, lat_lon_max) = min_and_max_ground_range(altitude, fov, boresight, r_earth);
        
        let expected_max = 10.73517_f64.to_radians();
        let expected_min = 3.80977_f64.to_radians();
        
        assert_relative_eq!(lat_lon_max, expected_max, epsilon = 1e-4);
        assert_relative_eq!(lat_lon_min, expected_min, epsilon = 1e-4);
    }
    
    #[test]
    fn test_max_and_min_ground_range_nadir() {
        let altitude = 800.0; // km
        let fov = 25.0_f64.to_radians();
        let boresight = 0.0_f64.to_radians();
        let r_earth = EARTH.equatorial_radius_km;
        
        let (lat_lon_min, lat_lon_max) = min_and_max_ground_range(altitude, fov, boresight, r_earth);
        
        let expected_val = 1.5984_f64.to_radians();
        
        assert_relative_eq!(lat_lon_max, expected_val, epsilon = 1e-4);
        assert_relative_eq!(lat_lon_min, -expected_val, epsilon = 1e-4);
    }

    #[test]
    fn test_ground_range_diff_at_azimuth() {
        let altitude = 800.0;
        let fov = 25.0_f64.to_radians();
        let boresight = 40.0_f64.to_radians();
        let azimuth = 140.0_f64.to_radians();
        let nadir_lat = 50.0_f64.to_radians();
        let nadir_lon = 40.0_f64.to_radians();
        let r_earth = EARTH.equatorial_radius_km;
        
        let result = ground_range_diff_at_azimuth(altitude, fov, boresight, azimuth, nadir_lat, nadir_lon, r_earth);
        assert!(result.is_ok());
        let (ground_range_diff, target_lat, target_lon) = result.unwrap();
        
        // expected_ground_range_diff = (6.9254 / 2) deg
        // Note: 6.9254 is Lambda_max - Lambda_min for the params in test_max_and_min_ground_range?
        // Let's check: 10.73517 - 3.80977 = 6.9254. Correct.
        let expected_diff = (6.9254 / 2.0_f64).to_radians();
        let expected_target_lat = 44.9926_f64.to_radians();
        let expected_target_lon = 45.7577_f64.to_radians();
        
        assert_relative_eq!(ground_range_diff, expected_diff, epsilon = 1e-5);
        assert_relative_eq!(target_lat, expected_target_lat, epsilon = 1e-5);
        assert_relative_eq!(target_lon, expected_target_lon, epsilon = 1e-5);
    }
    
    #[test]
    fn test_exception_ground_range_diff_at_azimuth() {
        let altitude = 800.0;
        let fov = 25.0_f64.to_radians();
        let boresight = 40.0_f64.to_radians();
        let azimuth = 190.0_f64.to_radians(); // > PI
        let nadir_lat = 50.0_f64.to_radians();
        let nadir_lon = 40.0_f64.to_radians();
        let r_earth = EARTH.equatorial_radius_km;
        
        let result = ground_range_diff_at_azimuth(altitude, fov, boresight, azimuth, nadir_lat, nadir_lon, r_earth);
        assert!(result.is_err());
    }
}
