use crate::bodies::Body;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeanElements {
    pub a_km: f64,
    pub ecc: f64,
    pub inc_rad: f64,
}

pub fn get_mean_elements(body: Body) -> Result<MeanElements, String> {
    match body.name {
        "Mercury" => Ok(MeanElements {
            a_km: 57_909_050.0,
            ecc: 0.205630,
            inc_rad: 7.00487_f64.to_radians(),
        }),
        "Venus" => Ok(MeanElements {
            a_km: 108_208_000.0,
            ecc: 0.00677323,
            inc_rad: 3.39471_f64.to_radians(),
        }),
        "Earth" => Ok(MeanElements {
            a_km: 149_598_023.0,
            ecc: 0.0167086,
            inc_rad: 0.00005_f64.to_radians(),
        }),
        "Mars" => Ok(MeanElements {
            a_km: 227_939_200.0,
            ecc: 0.0934025,
            inc_rad: 1.85061_f64.to_radians(),
        }),
        "Jupiter" => Ok(MeanElements {
            a_km: 778_547_200.0,
            ecc: 0.048498,
            inc_rad: 1.30530_f64.to_radians(),
        }),
        "Saturn" => Ok(MeanElements {
            a_km: 1_433_449_370.0,
            ecc: 0.055546,
            inc_rad: 2.485240_f64.to_radians(),
        }),
        "Uranus" => Ok(MeanElements {
            a_km: 2_872_460_000.0,
            ecc: 0.047318,
            inc_rad: 0.772556_f64.to_radians(),
        }),
        "Neptune" => Ok(MeanElements {
            a_km: 4_495_060_000.0,
            ecc: 0.008606,
            inc_rad: 1.76917_f64.to_radians(),
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
