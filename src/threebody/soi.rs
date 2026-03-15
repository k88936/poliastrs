use crate::threebody::ThreeBodyError;

pub fn sphere_of_influence_radius_km(
    semi_major_axis_km: f64,
    secondary_mass: f64,
    primary_mass: f64,
) -> Result<f64, ThreeBodyError> {
    if secondary_mass <= 0.0 || primary_mass <= 0.0 {
        return Err(ThreeBodyError::InvalidMassRatio);
    }
    Ok(semi_major_axis_km * (secondary_mass / primary_mass).powf(2.0 / 5.0))
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::sphere_of_influence_radius_km;

    #[test]
    fn earth_soi_around_sun_is_close_to_reference() {
        let soi = sphere_of_influence_radius_km(149_597_870.7, 5.972e24, 1.9885e30).unwrap();
        assert_relative_eq!(soi, 924_650.0, max_relative = 5e-3);
    }

    #[test]
    fn invalid_masses_return_error() {
        assert!(sphere_of_influence_radius_km(1.0, 0.0, 1.0).is_err());
        assert!(sphere_of_influence_radius_km(1.0, 1.0, -1.0).is_err());
    }
}
