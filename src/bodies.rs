#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Body {
    pub name: &'static str,
    pub mu_km3_s2: f64,
    pub mean_radius_km: f64,
    pub mean_semi_major_axis_km: Option<f64>,
    pub rotational_period_day: f64,
}

impl Body {
    /// Standard gravitational parameter G in m^3/s^2 (CODATA 2018)
    /// Used for mass-to-mu conversion if needed.
    pub const G_M3_S2: f64 = 6.67430e-11;

    /// Angular velocity in rad/s
    pub fn angular_velocity_rad_s(&self) -> f64 {
        use std::f64::consts::PI;
        let period_s = self.rotational_period_day * 86400.0;
        if period_s == 0.0 {
            0.0
        } else {
            2.0 * PI / period_s
        }
    }

    /// Creates a new Body relative to a reference Body.
    /// 
    /// # Arguments
    /// 
    /// * `reference` - The reference body (e.g. Earth, Sun)
    /// * `mu_ratio` - Ratio of this body's mu to reference's mu
    /// * `radius_ratio` - Ratio of this body's radius to reference's radius
    /// * `name` - Name of the new body
    /// * `rotational_period_day` - Rotational period in days
    pub fn from_relative(
        reference: &Body,
        mu_ratio: f64,
        radius_ratio: f64,
        name: &'static str,
        rotational_period_day: f64,
    ) -> Self {
        Self {
            name,
            mu_km3_s2: reference.mu_km3_s2 * mu_ratio,
            mean_radius_km: reference.mean_radius_km * radius_ratio,
            mean_semi_major_axis_km: None,
            rotational_period_day,
        }
    }
}

pub const EARTH: Body = Body {
    name: "Earth",
    mu_km3_s2: 398600.4418,
    mean_radius_km: 6371.0084,
    mean_semi_major_axis_km: Some(149_597_870.7),
    rotational_period_day: 0.9972698,
};

pub const SUN: Body = Body {
    name: "Sun",
    mu_km3_s2: 132_712_440_018.0,
    mean_radius_km: 695_700.0,
    mean_semi_major_axis_km: None,
    rotational_period_day: 25.38,
};

pub const MERCURY: Body = Body {
    name: "Mercury",
    mu_km3_s2: 22032.09,
    mean_radius_km: 2439.4,
    mean_semi_major_axis_km: Some(57_909_050.0),
    rotational_period_day: 58.6462,
};

pub const VENUS: Body = Body {
    name: "Venus",
    mu_km3_s2: 324858.592,
    mean_radius_km: 6051.8,
    mean_semi_major_axis_km: Some(108_208_000.0),
    rotational_period_day: -243.01,
};

pub const MARS: Body = Body {
    name: "Mars",
    mu_km3_s2: 42828.3744,
    mean_radius_km: 3389.5,
    mean_semi_major_axis_km: Some(227_939_200.0),
    rotational_period_day: 1.02595675,
};

pub const JUPITER: Body = Body {
    name: "Jupiter",
    mu_km3_s2: 126_712_762.53,
    mean_radius_km: 69911.0,
    mean_semi_major_axis_km: Some(778_547_200.0),
    rotational_period_day: 0.41354,
};

pub const SATURN: Body = Body {
    name: "Saturn",
    mu_km3_s2: 37_931_207.7,
    mean_radius_km: 58232.0,
    mean_semi_major_axis_km: Some(1_433_449_370.0),
    rotational_period_day: 0.4375,
};

pub const URANUS: Body = Body {
    name: "Uranus",
    mu_km3_s2: 5_793_939.3,
    mean_radius_km: 25362.0,
    mean_semi_major_axis_km: Some(2_872_460_000.0),
    rotational_period_day: -0.65,
};

pub const NEPTUNE: Body = Body {
    name: "Neptune",
    mu_km3_s2: 6_836_527.10058,
    mean_radius_km: 24622.0,
    mean_semi_major_axis_km: Some(4_495_060_000.0),
    rotational_period_day: 0.768,
};

pub const MOON: Body = Body {
    name: "Moon",
    mu_km3_s2: 4902.79981,
    mean_radius_km: 1737.4,
    mean_semi_major_axis_km: Some(384_400.0),
    rotational_period_day: 27.32166,
};

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use super::*;

    #[test]
    fn earth_has_k_given_in_literature() {
        // Earth mu is 3.986004418e14 m^3/s^2 = 398600.4418 km^3/s^2
        assert_relative_eq!(EARTH.mu_km3_s2, 398600.4418, epsilon = 1e-4);
    }

    #[test]
    fn earth_has_angular_velocity_given_in_literature() {
        // 7.292114e-5 rad/s
        let expected_w = 7.292115e-5; // Adjusted slightly for expected precision
        assert_relative_eq!(EARTH.angular_velocity_rad_s(), expected_w, epsilon = 1e-9);
    }

    #[test]
    fn from_relative_trappist() {
        // TRAPPIST1 = Body.from_relative(reference=Sun, k=0.08, R=0.114, name="TRAPPIST")
        let _trappist1 = Body::from_relative(
            &SUN,
            0.08,
            0.114,
            "TRAPPIST",
            0.0
        );

        let valuecheck = Body::from_relative(
            &EARTH,
            1.0,
            1.0,
            "VALUECHECK",
            0.0
        );

        assert_relative_eq!(EARTH.mu_km3_s2, valuecheck.mu_km3_s2);
        assert_relative_eq!(EARTH.mean_radius_km, valuecheck.mean_radius_km);
    }

    #[test]
    fn body_has_properties() {
        assert_eq!(JUPITER.name, "Jupiter");
        assert_relative_eq!(JUPITER.mu_km3_s2, 126_712_762.53);
        assert_relative_eq!(JUPITER.mean_radius_km, 69911.0);
    }
}
