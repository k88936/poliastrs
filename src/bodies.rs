#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Body {
    pub name: &'static str,
    pub mu_km3_s2: f64,
    pub mean_radius_km: f64,
}

pub const EARTH: Body = Body {
    name: "Earth",
    mu_km3_s2: 398600.4418,
    mean_radius_km: 6378.1363,
};

pub const SUN: Body = Body {
    name: "Sun",
    mu_km3_s2: 132_712_440_018.0,
    mean_radius_km: 695_700.0,
};
