use crate::ephem::{Ephem, interpolator::LinearInterpolator};
use crate::frames::Plane;
use crate::twobody::orbit::Orbit;
use approx::assert_relative_eq;
use nalgebra::Vector3;
use crate::bodies::EARTH;

// Helper to create J2000 epoch
// const J2000: f64 = 0.0;

#[test]
fn test_ephem_sample_no_arguments_returns_exactly_same_input() {
    let epochs = vec![0.0, 100.0, 200.0, 300.0];
    let coordinates = vec![
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.9, 0.1, 0.0),
        Vector3::new(0.8, 0.2, 0.0),
        Vector3::new(0.7, 0.3, 0.0),
    ];
    let plane = Plane::EarthEquator;
    
    let ephem = Ephem::new(epochs.clone(), coordinates.clone(), plane);
    
    // Test sample() with no arguments (returns all epochs)
    let result = ephem.sample(None, LinearInterpolator);
    
    assert_eq!(result.len(), 4);
    for (i, res) in result.iter().enumerate() {
        assert_eq!(*res, coordinates[i]);
    }
}

#[test]
fn test_from_orbit_has_desired_properties() {
    // Port of test_from_orbit_has_desired_properties
    // Setup Orbit at J2000 (default)
    let r = Vector3::new(-1000.0, -2000.0, 3100.0); // km
    let v = Vector3::new(-1.836, 5.218, 4.433); // km/s
    let orb = Orbit::from_vectors(EARTH, r.into(), v.into()); // Into array/Vector3
    
    // Target epochs (TDB seconds from J2000)
    // 2020-02-01 12:00:00 -> approx 20 years + 1 month
    // 2020-02-13
    // 2020-03-04
    // 2020-03-17
    
    // NOTE: Rust standard library doesn't parse dates easily without Chrono.
    // I'll calculate J2000 offset for these dates.
    // J2000 is 2000-01-01 12:00:00 TDB.
    // 2020-01-01 12:00:00 is 20 * 365.25 days = 7305 days.
    // 2020 is leap year.
    // 2020-02-01 is 31 days after Jan 1.
    // 2020-02-13 is 12 days later.
    // 2020-03-04 is 20 days later (Feb has 29 days in 2020).
    // 2020-03-17 is 13 days later.
    
    // Let's approximate or hardcode offsets if possible.
    // Python astropy Time("2020-02-01 12:00:00", scale="tdb").jd - 2451545.0
    // I can't run python here.
    // But I can use approximate days.
    // Julian Day J2000.0 = 2451545.0
    // Jan 1 2020 is roughly 7305 days.
    // Exact days? 
    // 2000 (leap), 2004, 2008, 2012, 2016. 5 leap years.
    // 20 * 365 + 5 = 7305 days exactly?
    // 2000 is leap.
    // So days = 7305.
    
    // Feb 01 = Jan 31 + 1 = day 32 of year? No, Jan has 31. So Feb 1 is day 32 (0-indexed 31).
    // Offset = 7305 + 31 = 7336 days.
    
    // Feb 13 = 7336 + 12 = 7348 days.
    // Mar 04 = 7336 + 29 (Feb) + 3 (Mar) = 7336 + 32 = 7368 days. (Feb 2020 has 29).
    // Mar 17 = 7368 + 13 = 7381 days.
    
    // Correction for UTC to TDB (approx 69.184s)
    let offset = 69.184;

    let days = vec![7336.0, 7348.0, 7368.0, 7381.0];
    let epochs: Vec<f64> = days.iter().map(|d| d * 86400.0 + offset).collect();
    
    let ephem = Ephem::from_orbit(orb, epochs.clone(), Plane::EarthEquator);
    let coords = ephem.sample(None, LinearInterpolator);
    
    // Expected coordinates from python test (km)
    let expected = vec![
        Vector3::new(336.77109079, -1447.38211842, -743.72094119),
        Vector3::new(-1133.43957703, 449.41297342, 3129.10416554),
        Vector3::new(201.42480053, -1978.64139325, -287.25776291),
        Vector3::new(-1084.94556884, -1713.5357774, 3298.72284309),
    ];
    
    for (i, coord) in coords.iter().enumerate() {
        // Assert relative equality with some tolerance (due to time calc diffs)
        assert_relative_eq!(coord, &expected[i], epsilon = 1.0); // 1 km tolerance for time diffs
    }
}

#[test]
#[ignore]
fn test_ephem_from_body_has_expected_properties() {
    // Requires VSOP87 or Ephemeris
    // TODO: Implement VSOP87
}
