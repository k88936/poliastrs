use crate::bodies::{EARTH, SUN};
use crate::twobody::orbit::Orbit;
use crate::frames::Plane;

/// ISS orbit example
/// Taken from Plyades (c) 2012 Helge Eichhorn (MIT License)
pub fn iss() -> Orbit {
    let r_km = [8.59072560e2, -4.13720368e3, 5.29556871e3];
    let v_km_s = [7.37289205, 2.08223573, 4.39999794e-1];
    // "2013-03-18 12:00" UTC -> TDB seconds from J2000
    let epoch = 416880067.1855774;
    
    Orbit::from_vectors_at(EARTH, r_km, v_km_s, epoch, Plane::EarthEquator)
}

/// Molniya orbit example
pub fn molniya() -> Orbit {
    let a_km = 26600.0;
    let ecc = 0.75;
    let inc_deg: f64 = 63.4;
    let raan_deg: f64 = 0.0;
    let argp_deg: f64 = 270.0;
    let nu_deg: f64 = 80.0;
    
    Orbit::from_keplerian(
        EARTH,
        a_km,
        ecc,
        inc_deg.to_radians(),
        raan_deg.to_radians(),
        argp_deg.to_radians(),
        nu_deg.to_radians(),
        0.0, // J2000
        Plane::EarthEquator
    )
}

/// Soyuz geostationary transfer orbit (GTO) example
/// Taken from Soyuz User's Manual, issue 2 revision 0
pub fn soyuz_gto() -> Orbit {
    let r_earth = EARTH.equatorial_radius_km;
    
    let r_a = r_earth + 35950.0;
    let r_p = r_earth + 250.0;
    let a = (r_a + r_p) / 2.0;
    let ecc = r_a / a - 1.0;
    
    let inc_deg: f64 = 6.0;
    let raan_deg: f64 = 188.5;
    let argp_deg: f64 = 178.0;
    let nu_deg: f64 = 0.0;
    
    Orbit::from_keplerian(
        EARTH,
        a,
        ecc,
        inc_deg.to_radians(),
        raan_deg.to_radians(),
        argp_deg.to_radians(),
        nu_deg.to_radians(),
        0.0, // J2000
        Plane::EarthEquator
    )
}

/// Comet 67P/Churyumov–Gerasimenko orbit example
pub fn churi() -> Orbit {
    let au_km = 149_597_870.700; // Standard AU
    let a_km = 3.46250 * au_km;
    let ecc = 0.64;
    let inc_deg: f64 = 7.04;
    let raan_deg: f64 = 50.1350;
    let argp_deg: f64 = 12.8007;
    let nu_deg: f64 = 63.89;
    
    // "2015-11-05 12:00" UTC -> TDB seconds from J2000
    let epoch = 499996868.1825882;
    
    Orbit::from_keplerian(
        SUN,
        a_km,
        ecc,
        inc_deg.to_radians(),
        raan_deg.to_radians(),
        argp_deg.to_radians(),
        nu_deg.to_radians(),
        epoch,
        Plane::EarthEcliptic 
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_examples_creation() {
        let _iss = iss();
        let _molniya = molniya();
        let _soyuz = soyuz_gto();
        let _churi = churi();
    }
}
