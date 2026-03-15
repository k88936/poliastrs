use crate::bodies::Body;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeanElements {
    pub a_km: f64,
    pub ecc: f64,
    pub inc_rad: f64,
}

pub fn get_mean_elements(body: Body) -> Result<MeanElements, String> {
    match body.name {
        "Earth" => Ok(MeanElements {
            a_km: 149_598_023.0,
            ecc: 0.0167086,
            inc_rad: 0.00005_f64.to_radians(),
        }),
        _ => Err(format!(
            "No available mean elements for {:?}. It must be an instance of `poliastro.bodies.SolarSystemPlanet`",
            body
        )),
    }
}

#[cfg(test)]
mod tests {
    use crate::bodies::SUN;

    use super::get_mean_elements;

    #[test]
    fn get_mean_elements_raises_for_invalid_body() {
        let err = get_mean_elements(SUN).unwrap_err();
        assert!(err.contains(
            "No available mean elements for Body { name: \"Sun\""
        ));
        assert!(err.contains("poliastro.bodies.SolarSystemPlanet"));
    }
}
