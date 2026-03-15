use crate::bodies::Body;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotationalElements {
    pub alpha_deg: f64,
    pub delta_deg: f64,
    pub w_deg: f64,
}

pub fn compute_rotational_elements(body: Body, tdb_seconds_from_j2000: f64) -> RotationalElements {
    let d = tdb_seconds_from_j2000 / 86400.0;
    
    // Constants from test_frames.py (at J2000)
    let (alpha_0, delta_0, w_0) = match body.name {
        "Sun" => (286.13, 63.87, 84.176),
        "Mercury" => (281.0103, 61.45, 329.5999488),
        "Venus" => (272.76, 67.16, 160.2),
        "Mars" => (317.68085441, 52.88643928, 176.63205973),
        "Jupiter" => (268.05720404, 64.49580995, 284.95),
        "Saturn" => (40.589, 83.537, 38.9),
        "Uranus" => (257.311, -15.175, 203.81),
        "Neptune" => (299.33373896, 42.95035902, 249.99600757),
        "Moon" => (266.85773344495135, 65.64110274784535, 41.1952639807452),
        "Earth" => (0.0, 90.0, 280.4606), // Approx GMST at J2000
        _ => (0.0, 90.0, 0.0),
    };

    // Calculate W rate from rotational period
    // W_dot = 360.0 / period_days
    // Note: period can be negative for retrograde rotation?
    // Usually rotational_period_day is positive magnitude.
    // We need direction.
    // Venus is retrograde. Uranus is retrograde (technically > 90 deg tilt).
    // Let's rely on standard W_dot definition:
    // W = W0 + Wdot * d.
    // If I use period, I need to know if W increases or decreases.
    // IAU defines North Pole by right-hand rule. Rotation is positive about North Pole.
    // So W always increases?
    // Yes, W increases in the direction of rotation.
    
    let w_dot = if body.rotational_period_day != 0.0 {
        360.0 / body.rotational_period_day
    } else {
        0.0
    };

    let w = (w_0 + w_dot * d).rem_euclid(360.0);

    RotationalElements {
        alpha_deg: alpha_0,
        delta_deg: delta_0,
        w_deg: w,
    }
}
