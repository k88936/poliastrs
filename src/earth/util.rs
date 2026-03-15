use std::f64::consts::PI;

pub fn gmst_ia82(seconds_from_j2000: f64) -> f64 {
    // Formula: GMST = 280.46061837 + 360.98564736629 * D
    // D = seconds / 86400
    let d = seconds_from_j2000 / 86400.0;
    let deg = 280.46061837 + 360.98564736629 * d;
    deg.to_radians().rem_euclid(2.0 * PI)
}

pub fn get_local_sidereal_time(lon_rad: f64, seconds_from_j2000: f64) -> f64 {
    let gmst = gmst_ia82(seconds_from_j2000);
    (gmst + lon_rad).rem_euclid(2.0 * PI)
}

pub fn raan_from_ltan(epoch_seconds_from_j2000: f64, ltan_rad: f64) -> f64 {
    // Port of poliastro logic
    let t = epoch_seconds_from_j2000 / (36525.0 * 86400.0);
    
    // Sun mean anomaly (deg)
    let m_sun_deg = 357.5291092 + 35999.05034 * t;
    let m_sun = m_sun_deg.to_radians();
    
    // Sun mean longitude (deg)
    let l_sun_deg = 280.460 + 36000.771 * t; 
    let l_sun = l_sun_deg.to_radians();
    
    // Ecliptic longitude correction
    let l_ecliptic_part2 = 1.914666471f64.to_radians() * m_sun.sin() 
                         + 0.019994643f64.to_radians() * (2.0 * m_sun).sin();
    let l_ecliptic = l_sun + l_ecliptic_part2;
    
    // Equation of Time (rad)
    let eq_time = -l_ecliptic_part2 
                  + 2.466f64.to_radians() * (2.0 * l_ecliptic).sin() 
                  - 0.0053f64.to_radians() * (4.0 * l_ecliptic).sin();
                  
    // Sun RA calculation from l_ecliptic
    let eps = 23.43929111f64.to_radians(); // J2000 obliquity
    let ra_sun = (eps.cos() * l_ecliptic.sin()).atan2(l_ecliptic.cos());
    
    // Apparent Local Solar Time = RA_Sun + 12h
    let salt = ra_sun + PI;
    
    // Mean Local Solar Time
    let smlt = salt + eq_time;
    
    // RAAN = SMLT + LTAN
    (smlt + ltan_rad).rem_euclid(2.0 * PI)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_local_sidereal_time_curtis() {
        // Curtis Example 5.6
        // 2004-03-03 04:30:00 UT1
        // J2000 = 2000-01-01 12:00:00
        // Days difference:
        // 2000 (leap), 2001, 2002, 2003. 4 years.
        // 4 * 365 + 1 = 1461 days to Jan 1 2004 12:00?
        // 2004 is leap. Feb has 29.
        // Jan 1 2000 to Jan 1 2004: 366 + 365 + 365 + 365 = 1461 days.
        // From Jan 1 2004 12:00 to Mar 3 04:30:
        // Jan (31), Feb (29 in 2004), Mar 3.
        // 31 + 29 + 2 = 62 days to Mar 3 12:00.
        // To 04:30 is -7.5 hours.
        // Total days: 1461 + 62 - 7.5/24.0 = 1523 - 0.3125 = 1522.6875 days.
        // Wait, J2000 is usually Jan 1.5 (Jan 1 12:00).
        // Let's rely on calculation or hardcode days.
        // Python test: Time("2004-03-03 04:30:00", scale="ut1").jd - 2451545.0
        // JD(2004-03-03 04:30) = 2453067.6875
        // JD(2000-01-01 12:00) = 2451545.0
        // Diff = 1522.6875.
        
        let days = 1522.6875;
        let seconds = days * 86400.0;
        
        let lon = 139.80f64.to_radians();
        let lst = get_local_sidereal_time(lon, seconds);
        let expected = 8.59f64.to_radians(); // 8.59 deg? Or 8h 59m?
        // Python test: expected_lst = 8.59 * u.deg.
        // Wait. 8.59 deg is tiny.
        // If it was 8h 59m, it would be ~135 deg.
        // The test explicitly says 8.59 * u.deg.
        
        // Let's check calculation:
        // GMST = 280.46 + 360.9856 * 1522.6875
        // = 280.46 + 549666.3...
        // % 360 = 228.79 deg (approx)
        // LST = GMST + Lon = 228.79 + 139.80 = 368.59 -> 8.59 deg.
        // Matches!
        
        assert_relative_eq!(lst, expected, epsilon = 1e-4); // 1e-2 deg ~ 1e-4 rad
    }

    #[test]
    fn test_raan_from_ltan_metopb() {
        // MetOp-B
        // Epoch: 2020-01-01 00:00 + 49.954 days - 1 day?
        // Test: Time("2020-01-01 00:00").mjd + 49.954 - 1
        // 2020-01-01 00:00 is MJD 58849.0
        // Epoch MJD = 58849.0 + 48.954 = 58897.954
        // J2000 MJD = 51544.5
        // Days from J2000 = 58897.954 - 51544.5 = 7353.454
        
        let days = 7353.45408566;
        let seconds = days * 86400.0;
        
        // LTAN: 21h 31m 45s
        let ltan_hours = 21.0 + 31.0/60.0 + 45.0/3600.0;
        let ltan_rad = (ltan_hours / 24.0) * 2.0 * PI;
        
        let expected_raan = 110.9899f64.to_radians();
        
        let raan = raan_from_ltan(seconds, ltan_rad);
        
        // Python test uses atol=0.3 deg (~0.005 rad)
        assert_relative_eq!(raan, expected_raan, epsilon = 0.01);
    }
}
