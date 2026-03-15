
/// Converts from given geometric altitude to geopotential one.
/// z: Geometric altitude [km]
/// r0: Planet radius [km]
pub fn geometric_to_geopotential(z: f64, r0: f64) -> f64 {
    r0 * z / (r0 + z)
}

/// Converts from given geopotential altitude to geometric one.
/// h: Geopotential altitude [km]
/// r0: Planet radius [km]
pub fn geopotential_to_geometric(h: f64, r0: f64) -> f64 {
    r0 * h / (r0 - h)
}

/// Relates Earth gravity field magnitude with the geometric height.
/// z: Geometric height [km]
/// g0: Gravity value at sea level [m/s^2]
/// r0: Planet radius [km]
pub fn gravity(z: f64, g0: f64, r0: f64) -> f64 {
    g0 * (r0 / (r0 + z)).powi(2)
}

/// Finds element in list and returns index.
/// x: Element to be searched.
/// x_levels: List for searching (must be sorted).
pub fn get_index(x: f64, x_levels: &[f64]) -> usize {
    for (i, &value) in x_levels.iter().enumerate() {
        if value > x {
            if i == 0 { return 0; }
            return i - 1;
        }
    }
    x_levels.len().saturating_sub(1)
}

/// Checks if altitude is inside valid range.
/// alt: Altitude to be checked [km]
/// r0: Attractor radius [km]
/// geometric: If `true`, assumes geometric altitude kind.
/// limits: (min_z, max_z) [km]
/// Returns: (z, h)
pub fn check_altitude(alt: f64, r0: f64, geometric: bool, min_z: f64, max_z: f64) -> Result<(f64, f64), String> {
    let (z, h) = if geometric {
        (alt, geometric_to_geopotential(alt, r0))
    } else {
        (geopotential_to_geometric(alt, r0), alt)
    };
    
    if z < min_z || z > max_z {
        return Err(format!("Geometric altitude must be in range [{:.1} km, {:.1} km]", min_z, max_z));
    }
    
    Ok((z, h))
}
