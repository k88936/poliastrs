use poliastrs::bodies::SUN;
use poliastrs::twobody::orbit::Orbit;
use poliastrs::twobody::maneuver::Maneuver;
use poliastrs::core::elements::ClassicalElements;
use poliastrs::core::angles::{m_to_e, e_to_nu};
use chrono::{TimeZone, Utc, DateTime};

fn mean_to_true(mean: f64, ecc: f64) -> f64 {
    let e_anom = m_to_e(mean, ecc);
    e_to_nu(e_anom, ecc)
}

fn main() {
    println!("Going to Jupiter with Rust and poliastrs");
    
    // Dates
    let date_launch = Utc.with_ymd_and_hms(2011, 8, 5, 16, 25, 0).unwrap();
    let date_flyby = Utc.with_ymd_and_hms(2013, 10, 9, 19, 21, 0).unwrap();
    let date_arrival = Utc.with_ymd_and_hms(2016, 7, 5, 3, 18, 0).unwrap();
    
    let j2000 = Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap();
    
    // Helper to get seconds from J2000
    let get_seconds = |dt: DateTime<Utc>| (dt - j2000).num_milliseconds() as f64 / 1000.0;
    
    let t_launch = get_seconds(date_launch);
    let t_flyby = get_seconds(date_flyby);
    let t_arrival = get_seconds(date_arrival);
    
    // Approximate Ephemerides (J2000 Elements propagated)
    // Earth J2000: a=1.00000011 AU, e=0.01671022, i=0.00005 deg, L=100.46435 deg, w_bar=102.94719 deg, raan=-11.26064 deg
    // M = L - w_bar = 100.46435 - 102.94719 = -2.48284 deg
    // w = w_bar - raan = 102.94719 - (-11.26064) = 114.20783 deg
    
    let au_km = 149_597_870.7;
    let earth_a = 1.00000011 * au_km;
    let earth_e = 0.01671022;
    let earth_i = 0.00005_f64.to_radians();
    let earth_raan = (-11.26064_f64).to_radians();
    let earth_w = 114.20783_f64.to_radians();
    let earth_m = (-2.48284_f64).to_radians();
    
    let earth_nu = mean_to_true(earth_m, earth_e);
    
    let earth_j2000 = Orbit::from_classical(
        SUN,
        ClassicalElements {
            p_km: earth_a * (1.0 - earth_e * earth_e),
            ecc: earth_e,
            inc_rad: earth_i,
            raan_rad: earth_raan,
            argp_rad: earth_w,
            nu_rad: earth_nu,
        }
    );
    
    // Jupiter J2000: a=5.20336301 AU, e=0.04839266, i=1.30530 deg, L=34.40438 deg, w_bar=14.75385 deg, raan=100.55615 deg
    // M = L - w_bar = 34.40438 - 14.75385 = 19.65053 deg
    // w = w_bar - raan = 14.75385 - 100.55615 = -85.8023 deg
    
    let jup_a = 5.20336301 * au_km;
    let jup_e = 0.04839266;
    let jup_i = 1.30530_f64.to_radians();
    let jup_raan = 100.55615_f64.to_radians();
    let jup_w = (-85.8023_f64).to_radians();
    let jup_m = 19.65053_f64.to_radians();
    
    let jup_nu = mean_to_true(jup_m, jup_e);
    
    let jupiter_j2000 = Orbit::from_classical(
        SUN,
        ClassicalElements {
            p_km: jup_a * (1.0 - jup_e * jup_e),
            ecc: jup_e,
            inc_rad: jup_i,
            raan_rad: jup_raan,
            argp_rad: jup_w,
            nu_rad: jup_nu,
        }
    );
    
    // Propagate to dates
    let earth_launch = earth_j2000.propagate_seconds(t_launch).unwrap();
    let earth_flyby = earth_j2000.propagate_seconds(t_flyby).unwrap();
    let jupiter_arrival = jupiter_j2000.propagate_seconds(t_arrival).unwrap();
    
    println!("Earth Launch (J2000+{}s): {:?}", t_launch, earth_launch.state.r_km);
    
    // 1. Launch Maneuver
    let c3: f64 = 31.1; // km^2/s^2
    let v_inf = c3.sqrt();
    let v_e0 = earth_launch.state.v_km_s;
    let dv_launch = v_inf * v_e0.normalize();
    
    let man_launch = Maneuver::impulse(dv_launch);
    let ic1 = earth_launch.apply_maneuver(&man_launch);
    
    println!("IC1 Period: {:.2} days", ic1.period_seconds().unwrap() / 86400.0);
    
    // Propagate to aphelion (approximate transfer end)
    let ic1_end = ic1.propagate_to_anomaly(180.0_f64.to_radians()).unwrap();
    
    // 2. Deep Space Maneuver (Lambert to Earth Flyby)
    let man_flyby = Maneuver::lambert(&ic1_end, &earth_flyby).unwrap();
    let (dt1, dv1) = man_flyby.impulses[0];
    let (_dt2, dv2) = man_flyby.impulses[1];
    
    println!("DSM (at Aphelion) Delta-V: {:.4} km/s", dv1.norm());
    println!("Flyby Arrival Delta-V: {:.4} km/s (relative to Earth velocity)", dv2.norm());
    
    // Apply DSM
    let ic2 = ic1_end.apply_maneuver(&Maneuver::impulse(dv1));
    let ic2_end = ic2.propagate_seconds(man_flyby.impulses[1].0).unwrap();
    
    // 3. Jupiter Transfer (Lambert from Earth Flyby to Jupiter)
    let man_jupiter = Maneuver::lambert(&ic2_end, &jupiter_arrival).unwrap();
    let (dt_j1, dv_j1) = man_jupiter.impulses[0];
    let (_dt_j2, dv_j2) = man_jupiter.impulses[1];
    
    println!("Earth Flyby Departure Delta-V: {:.4} km/s", dv_j1.norm());
    println!("Jupiter Arrival Delta-V: {:.4} km/s (Insertion)", dv_j2.norm());
    
    // Total Delta-V
    let total_dv = v_inf + dv1.norm() + dv_j1.norm() + dv_j2.norm(); // Approximate
    
    println!("Total Delta-V (approx): {:.4} km/s", total_dv);
}
