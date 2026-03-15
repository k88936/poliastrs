use crate::bodies::Body;
use nalgebra::Vector3;

/// Low level calculations for oblate spheroid locations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpheroidLocation {
    pub lon_rad: f64,
    pub lat_rad: f64,
    pub h_km: f64,
    pub body: Body,
}

impl SpheroidLocation {
    pub fn new(lon_rad: f64, lat_rad: f64, h_km: f64, body: Body) -> Self {
        Self { lon_rad, lat_rad, h_km, body }
    }

    /// Semi-major axis (a)
    pub fn a(&self) -> f64 {
        self.body.equatorial_radius_km
    }

    /// Semi-minor axis (c)
    pub fn c(&self) -> f64 {
        self.body.polar_radius_km()
    }

    /// First flattening (f)
    pub fn f(&self) -> f64 {
        self.body.flattening
    }

    /// Eccentricity squared (e^2)
    pub fn e2(&self) -> f64 {
        let f = self.f();
        f * (2.0 - f)
    }

    /// Normal vector of the ellipsoid at the given location (N)
    /// N = [2x/a^2, 2y/a^2, 2z/c^2] normalized?
    /// Wait, implementation in Python:
    /// N = np.array([2 * x / a**2, 2 * y / b**2, 2 * z / c**2]) normalized.
    /// Here b=a (equatorial radius). So 2x/a^2, 2y/a^2.
    pub fn n_vector(&self) -> Vector3<f64> {
        let (x, y, z) = self.cartesian_cords_tuple();
        let a = self.a();
        let c = self.c();
        let n = Vector3::new(
            2.0 * x / (a * a),
            2.0 * y / (a * a),
            2.0 * z / (c * c)
        );
        n.normalize()
    }

    /// Calculates cartesian coordinates.
    pub fn cartesian_cords(&self) -> Vector3<f64> {
        let (x, y, z) = self.cartesian_cords_tuple();
        Vector3::new(x, y, z)
    }

    fn cartesian_cords_tuple(&self) -> (f64, f64, f64) {
        let a = self.a();
        let e2 = self.e2();
        let lat = self.lat_rad;
        let lon = self.lon_rad;
        let h = self.h_km;
        
        let sin_lat = lat.sin();
        let cos_lat = lat.cos();
        let n_val = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
        
        let x = (n_val + h) * cos_lat * lon.cos();
        let y = (n_val + h) * cos_lat * lon.sin();
        let z = ((1.0 - e2) * n_val + h) * sin_lat;
        
        (x, y, z)
    }

    /// Returns orthonormal vectors tangential to the ellipsoid at the given location.
    pub fn tangential_vecs(&self) -> (Vector3<f64>, Vector3<f64>) {
        let n = self.n_vector();
        let mut u = Vector3::new(1.0, 0.0, 0.0);
        // Gram-Schmidt / projection
        u -= n * n.dot(&u);
        u = u.normalize();
        let v = n.cross(&u);
        (u, v)
    }

    /// Radius of curvature of the meridian at the latitude of the given location.
    pub fn radius_of_curvature(&self) -> f64 {
        let a = self.a();
        let e2 = self.e2();
        let lat = self.lat_rad;
        let sin_lat = lat.sin();
        
        a * (1.0 - e2) / (1.0 - e2 * sin_lat * sin_lat).powf(1.5)
    }

    /// Calculates the distance from an arbitrary point to the given location.
    pub fn distance(&self, point: Vector3<f64>) -> f64 {
        (self.cartesian_cords() - point).norm()
    }

    /// Determine whether an object located at a given point is visible from the given location.
    pub fn is_visible(&self, point: Vector3<f64>) -> bool {
        let c = self.cartesian_cords();
        let n = self.n_vector();
        
        // Plane equation: N . (P - C) >= 0 ?
        // Python: d = -(N @ c); p = (N @ u) + d; return p >= 0
        // p = N.u - N.c = N.(u - c)
        // So check N . (point - c) >= 0
        let vec = point - c;
        n.dot(&vec) >= 0.0
    }

    /// Converts cartesian coordinates to ellipsoidal coordinates for the given ellipsoid.
    /// Uses Bowring's approximation.
    pub fn cartesian_to_ellipsoidal(body: Body, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let a = body.equatorial_radius_km;
        let c = body.polar_radius_km();
        let e2 = body.flattening * (2.0 - body.flattening); // e^2
        let e2_prime = e2 / (1.0 - e2); // e'^2 = (a^2 - c^2) / c^2 ? Or e2 / (1-e2)
        
        let p = (x * x + y * y).sqrt();
        let th = (z * a / (p * c)).atan();
        let sin_th = th.sin();
        let cos_th = th.cos();
        
        let lon = y.atan2(x);
        let lat = (
            (z + e2_prime * c * sin_th * sin_th * sin_th) /
            (p - e2 * a * cos_th * cos_th * cos_th)
        ).atan();
        
        let sin_lat = lat.sin();
        let v = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
        
        let h = if lat.abs() < 1e-18 {
            p / lat.cos() - v
        } else {
            z / sin_lat - (1.0 - e2) * v
        };
        
        (lon, lat, h)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use crate::bodies::EARTH;

    #[test]
    fn test_cartesian_coordinates() {
        let expected_cords = Vector3::new(
            3764859.30127275 / 1000.0,
            2987201.67496698 / 1000.0,
            4179160.71540021 / 1000.0,
        );
        
        let el_cords = (38.43_f64.to_radians(), 41.2_f64.to_radians(), 0.0);
        
        let p = SpheroidLocation::new(el_cords.0, el_cords.1, el_cords.2, EARTH);
        let c_cords = p.cartesian_cords();
        
        // Use relative epsilon for large coordinates
        assert_relative_eq!(c_cords, expected_cords, max_relative = 1e-7);
    }

    #[test]
    fn test_tangential_vectors() {
        let el_cords = (38.43_f64.to_radians(), 41.2_f64.to_radians(), 0.0);
        let p = SpheroidLocation::new(el_cords.0, el_cords.1, el_cords.2, EARTH);
        
        let n = p.n_vector();
        let (v1, v2) = p.tangential_vecs();
        
        assert_relative_eq!(n.dot(&v1), 0.0, epsilon = 1e-7);
        assert_relative_eq!(n.dot(&v2), 0.0, epsilon = 1e-7);
    }

    #[test]
    fn test_visible() {
        let el_cords = (38.43_f64.to_radians(), 41.2_f64.to_radians(), 0.0);
        let p = SpheroidLocation::new(el_cords.0, el_cords.1, el_cords.2, EARTH);
        
        let cords = p.cartesian_cords();
        let n = p.n_vector();
        
        let p1 = cords + 10.0 / 1000.0 * n; // 10m up
        let p2 = cords - 10.0 / 1000.0 * n; // 10m down
        
        assert!(p.is_visible(p1));
        assert!(!p.is_visible(p2));
    }

    #[test]
    fn test_f() {
        let expected_f = 0.0033528131;
        let el_cords = (38.43_f64.to_radians(), 41.2_f64.to_radians(), 0.0);
        let p = SpheroidLocation::new(el_cords.0, el_cords.1, el_cords.2, EARTH);
        
        assert_relative_eq!(p.f(), expected_f, epsilon = 1e-9);
    }
    
    #[test]
    fn test_radius_of_curvature() {
        let expected_roc = 6363141.421601379 / 1000.0; // km
        let el_cords = (38.43_f64.to_radians(), 41.2_f64.to_radians(), 0.0);
        let p = SpheroidLocation::new(el_cords.0, el_cords.1, el_cords.2, EARTH);
        
        assert_relative_eq!(p.radius_of_curvature(), expected_roc, max_relative = 1e-7);
    }
    
    #[test]
    fn test_distance() {
        let expected_distance = 6368850.150294118 / 1000.0;
        let el_cords = (38.43_f64.to_radians(), 41.2_f64.to_radians(), 0.0);
        let point_cords = Vector3::new(10.5 / 1000.0, 35.5 / 1000.0, 45.5 / 1000.0);
        
        let p = SpheroidLocation::new(el_cords.0, el_cords.1, el_cords.2, EARTH);
        let dist = p.distance(point_cords);
        
        assert_relative_eq!(dist, expected_distance, max_relative = 1e-7);
    }

    #[test]
    fn test_cartesian_conversion_approximate() {
        let el_cords = (0.7190227, 0.670680, 0.0);
        let c_cords = (
            3764258.64785411 / 1000.0,
            3295359.33856106 / 1000.0,
            3942945.28570563 / 1000.0,
        );
        
        let (lon, lat, h) = SpheroidLocation::cartesian_to_ellipsoidal(EARTH, c_cords.0, c_cords.1, c_cords.2);
        
        assert_relative_eq!(lon, el_cords.0, epsilon = 1e-4);
        assert_relative_eq!(lat, el_cords.1, epsilon = 1e-4);
        assert_relative_eq!(h, el_cords.2, epsilon = 0.001);
    }

    #[test]
    fn test_h_calculation_near_lat_singularity() {
        // Test near equator (lat ~ 0)
        let lat = 1e-5_f64;
        let lon = 10.0_f64.to_radians();
        let h_expected = 5.0 / 1000.0; // 5m
        
        let p = SpheroidLocation::new(lon, lat, h_expected, EARTH);
        let coords = p.cartesian_cords();
        let (_lon_out, _lat_out, h_out) = SpheroidLocation::cartesian_to_ellipsoidal(EARTH, coords.x, coords.y, coords.z);
        
        assert_relative_eq!(h_out, h_expected, epsilon = 1e-6);
    }
}
