use poliastrs::plotting::czml::CZMLExtractor;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::bodies::EARTH;
use poliastrs::core::elements::ClassicalElements;
use chrono::{Utc, TimeZone, Duration};

fn main() {
    println!("Visualizing orbital data with Cesium (CZML)");
    
    // Molniya orbit parameters
    // a = 26600 km, ecc = 0.74, inc = 63.4 deg, raan = 360 deg (0), argp = 270 deg, nu = 0
    let a = 26600.0;
    let ecc = 0.74;
    let inc_rad = 63.4_f64.to_radians();
    let raan_rad = 0.0_f64.to_radians(); // 360 is 0
    let argp_rad = 270.0_f64.to_radians();
    let nu_rad = 0.0_f64.to_radians();
    
    // p = a * (1 - e^2)
    let p_km = a * (1.0 - ecc * ecc);
    
    let coe = ClassicalElements {
        p_km,
        ecc,
        inc_rad,
        raan_rad,
        argp_rad,
        nu_rad,
    };
    
    // Create orbit. 
    // from_classical is available in Orbit
    let orbit = Orbit::from_classical(EARTH, coe);
    
    // Period
    let period_seconds = orbit.period_seconds().expect("Hyperbolic orbit or invalid parameters");
    println!("Orbit period: {:.2} minutes", period_seconds / 60.0);
    
    // Start epoch = J2000 (2000-01-01 12:00:00 UTC)
    let start_epoch = Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap();
    let end_epoch = start_epoch + Duration::milliseconds((period_seconds * 1000.0) as i64);
    
    println!("Simulating from {} to {}", start_epoch, end_epoch);
    
    let mut extractor = CZMLExtractor::new(start_epoch, end_epoch, 80);
    
    extractor.add_orbit(
        orbit, 
        "MolniyaOrbit", 
        2.0, 
        "Molniya", 
        [125, 80, 120, 255]
    );
    
    let packets = extractor.packets();
    
    let filename = "molniya.czml";
    let json = serde_json::to_string_pretty(&packets).unwrap();
    std::fs::write(filename, json).unwrap();
    
    println!("Generated CZML file: {}", filename);
    println!("You can load this file into Cesium to visualize the orbit.");
}
