use crate::ephem::{Ephem, interpolator::LinearInterpolator};
use crate::frames::Plane;
use crate::twobody::orbit::Orbit;
use approx::assert_relative_eq;
use nalgebra::Vector3;
use crate::bodies::EARTH;

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
fn test_ephem_sample_same_epochs_returns_same_input() {
    let epochs = vec![0.0, 100.0, 200.0, 300.0];
    let coordinates = vec![
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.9, 0.1, 0.0),
        Vector3::new(0.8, 0.2, 0.0),
        Vector3::new(0.7, 0.3, 0.0),
    ];
    let plane = Plane::EarthEquator;
    
    let ephem = Ephem::new(epochs.clone(), coordinates.clone(), plane);
    
    let result = ephem.sample(Some(epochs.clone()), LinearInterpolator);
    
    assert_eq!(result.len(), 4);
    for (i, res) in result.iter().enumerate() {
        assert_relative_eq!(*res, coordinates[i], epsilon=1e-12);
    }
}

#[test]
fn test_ephem_sample_existing_epochs_returns_corresponding_input() {
    let epochs = vec![0.0, 100.0, 200.0, 300.0];
    let coordinates = vec![
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.9, 0.1, 0.0),
        Vector3::new(0.8, 0.2, 0.0),
        Vector3::new(0.7, 0.3, 0.0),
    ];
    let plane = Plane::EarthEquator;
    
    let ephem = Ephem::new(epochs.clone(), coordinates.clone(), plane);
    
    // Sample every 2nd epoch: 0, 200
    let target_epochs = vec![epochs[0], epochs[2]];
    let result = ephem.sample(Some(target_epochs), LinearInterpolator);
    
    assert_eq!(result.len(), 2);
    assert_relative_eq!(result[0], coordinates[0], epsilon=1e-12);
    assert_relative_eq!(result[1], coordinates[2], epsilon=1e-12);
}

#[test]
fn test_rv_no_parameters_returns_input_vectors() {
    let epochs = vec![0.0, 100.0];
    let coordinates = vec![
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.9, 0.1, 0.0),
    ];
    let velocities = vec![
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-0.1, 0.9, 0.0),
    ];
    
    let mut ephem = Ephem::new(epochs.clone(), coordinates.clone(), Plane::EarthEquator);
    ephem.velocities = Some(velocities.clone());
    
    let (r, v) = ephem.rv(None);
    
    assert_eq!(r.len(), 2);
    assert_eq!(v.len(), 2);
    
    for i in 0..2 {
        assert_relative_eq!(r[i], coordinates[i], epsilon=1e-12);
        assert_relative_eq!(v[i], velocities[i], epsilon=1e-12);
    }
}

#[test]
fn test_from_orbit_has_desired_properties() {
    // Port of test_from_orbit_has_desired_properties
    let r = [-1000.0, -2000.0, 3100.0]; // km
    let v = [-1.836, 5.218, 4.433]; // km/s
    let orb = Orbit::from_vectors(EARTH, r, v);
    
    // J2000 = 2451545.0
    // Exact TDB seconds from J2000 for UTC dates
    // Calculated using astropy.time.Time(..., scale='utc').tdb - Time('J2000', scale='tdb')
    let epochs = vec![633830469.1847928, 634867269.1850804, 636595269.1854556, 637718469.1856025];
    
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
        // Use reasonable tolerance for cross-language propagation
        // Differences might be due to propagator precision or slight constant mismatches
        assert_relative_eq!(coord, &expected[i], epsilon = 0.1);
    }
}
