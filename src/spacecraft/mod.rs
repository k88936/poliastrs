#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spacecraft {
    pub area_km2: f64,
    pub drag_coeff: f64,
    pub mass_kg: f64,
}

impl Spacecraft {
    pub fn new(area_km2: f64, drag_coeff: f64, mass_kg: f64) -> Self {
        Self { area_km2, drag_coeff, mass_kg }
    }
    
    pub fn ballistic_coefficient(&self) -> f64 {
        // C_D * A / m (km^2 / kg)
        self.drag_coeff * self.area_km2 / self.mass_kg
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_ballistic_coefficient() {
        let cd = 2.2;
        let area_m2 = PI / 4.0;
        let area_km2 = area_m2 * 1e-6;
        let mass_kg = 100.0;
        
        let sc = Spacecraft::new(area_km2, cd, mass_kg);
        
        let expected = 1.7278759594743866e-08;
        assert_relative_eq!(sc.ballistic_coefficient(), expected, epsilon = 1e-12);
    }
}
