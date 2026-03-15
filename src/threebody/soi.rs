use crate::bodies::Body;
use crate::threebody::ThreeBodyError;

/// Calculates the sphere of influence radius using the Laplace definition (m/M)^(2/5).
/// This is commonly used for patched conics.
pub fn sphere_of_influence_radius_km(
    semi_major_axis_km: f64,
    secondary_mass_param: f64,
    primary_mass_param: f64,
) -> Result<f64, ThreeBodyError> {
    if secondary_mass_param <= 0.0 || primary_mass_param <= 0.0 {
        return Err(ThreeBodyError::InvalidMassRatio);
    }
    Ok(semi_major_axis_km * (secondary_mass_param / primary_mass_param).powf(2.0 / 5.0))
}

/// Calculates the Laplace radius (Sphere of Influence) for a body relative to its primary.
pub fn laplace_radius_km(body: &Body, primary: &Body) -> Result<f64, ThreeBodyError> {
    let a = body.mean_semi_major_axis_km.ok_or(ThreeBodyError::MissingSemiMajorAxis)?;
    sphere_of_influence_radius_km(a, body.mu_km3_s2, primary.mu_km3_s2)
}

/// Calculates the Hill radius (m/3M)^(1/3).
/// This is the region where the secondary body dominates the gravitational field.
pub fn hill_radius_km(body: &Body, primary: &Body) -> Result<f64, ThreeBodyError> {
    let a = body.mean_semi_major_axis_km.ok_or(ThreeBodyError::MissingSemiMajorAxis)?;
    let mu_ratio = body.mu_km3_s2 / (3.0 * primary.mu_km3_s2);
    if mu_ratio <= 0.0 {
        return Err(ThreeBodyError::InvalidMassRatio);
    }
    Ok(a * mu_ratio.cbrt())
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::bodies::{EARTH, JUPITER, MARS, MERCURY, NEPTUNE, SATURN, SUN, URANUS, VENUS};
    use super::*;

    #[test]
    fn earth_soi_around_sun_is_close_to_reference() {
        let soi = laplace_radius_km(&EARTH, &SUN).unwrap();
        // 9.25e8 m = 9.25e5 km
        assert_relative_eq!(soi, 925_000.0, max_relative = 1e-2);
    }

    #[test]
    fn invalid_masses_return_error() {
        assert!(sphere_of_influence_radius_km(1.0, 0.0, 1.0).is_err());
        assert!(sphere_of_influence_radius_km(1.0, 1.0, -1.0).is_err());
    }

    #[test]
    fn test_laplace_radius_planets() {
        // Data from Table A.2., Curtis "Orbital Mechanics for Engineering Students"
        // Converted to km
        // Note: Python tests use 1.12e8 m = 1.12e5 km, etc.
        let cases = vec![
            (MERCURY, 1.12e5),
            (VENUS, 6.16e5),
            (EARTH, 9.25e5),
            (MARS, 5.77e5),
            (JUPITER, 4.82e7),
            (SATURN, 5.48e7),
            (URANUS, 5.18e7),
            (NEPTUNE, 8.66e7),
        ];

        for (body, expected_km) in cases {
            let r_soi = laplace_radius_km(&body, &SUN).unwrap();
            assert_relative_eq!(r_soi, expected_km, max_relative = 1e-1);
        }
    }

    #[test]
    fn test_hill_radius_planets() {
        // Data from Chebotarev "Gravitational Spheres of the Major Planets, Moon and Sun"
        // Converted to km
        // Note: Poliastro implementation uses perihelion distance for Hill radius calculation for Mercury?
        // We use mean semi-major axis, so we expect the standard Hill Radius value (Chebotarev).
        let cases = vec![
            (MERCURY, 2.21e5), // Matches Chebotarev and standard calculation
            (VENUS, 1.03e6),
            (EARTH, 1.49e6),
            (MARS, 1.07e6),
            (JUPITER, 5.28e7),
            (SATURN, 6.50e7),
            (URANUS, 7.01e7),
            (NEPTUNE, 1.16e8),
        ];

        for (body, expected_km) in cases {
            let r_hill = hill_radius_km(&body, &SUN).unwrap();
            assert_relative_eq!(r_hill, expected_km, max_relative = 1e-1);
        }
    }
}
