use poliastrs::bodies::{EARTH, SUN, VENUS, Body};
use poliastrs::twobody::orbit::Orbit;
use poliastrs::core::elements::ClassicalElements;
use poliastrs::threebody::flyby::compute_flyby;
use poliastrs::twobody::maneuver::Maneuver;
use poliastrs::frames::Plane;
use nalgebra::Vector3;
use std::f64::consts::PI;

fn main() {
    println!("Analyzing the Parker Solar Probe flybys");
    println!("=======================================");
    
    // 1. Modulus of the exit velocity, some features of Orbit #2
    println!("\n## 1. Orbit #2 Features");
    
    let t_ref = 150.0 * 86400.0; // 150 days in seconds
    let k = SUN.mu_km3_s2;
    
    // T = 2 * PI * sqrt(a^3 / mu) => a = cbrt(mu * T^2 / (4 * PI^2))
    let a_ref = (k * t_ref * t_ref / (4.0 * PI * PI)).cbrt();
    let au_km = 149_597_870.7;
    println!("Reference Semimajor Axis (a_ref): {:.3} km ({:.3} AU)", a_ref, a_ref / au_km);
    
    let energy_ref = -k / (2.0 * a_ref);
    println!("Reference Energy: {:.3} km^2/s^2", energy_ref);
    
    // 2. Lambert arc between #0 and #1
    println!("\n## 2. Lambert arc between Launch and Flyby #1");
    
    // Dates
    // Launch: 2018-08-11
    // Flyby 1: 2018-09-28
    let jd_launch = julian_date(2018, 8, 11, 0.0);
    let jd_flyby1 = julian_date(2018, 9, 28, 0.0);
    
    // Epoch TDB seconds from J2000 (2451545.0)
    let epoch_launch = (jd_launch - 2451545.0) * 86400.0;
    let epoch_flyby1 = (jd_flyby1 - 2451545.0) * 86400.0;
    
    // Get approximate planetary positions
    // Earth at Launch
    let earth_launch = get_planet_orbit(EARTH, epoch_launch);
    let r0 = earth_launch.rv().0; // vector
    let v0_earth = earth_launch.rv().1;
    
    // Venus at Flyby 1
    let venus_flyby1 = get_planet_orbit(VENUS, epoch_flyby1);
    let r1 = venus_flyby1.rv().0; // vector
    let v_venus = venus_flyby1.rv().1; // Velocity of Venus
    
    let tof = epoch_flyby1 - epoch_launch;
    println!("Time of Flight: {:.2} days", tof / 86400.0);
    
    // Create Orbits for Lambert solver
    let orb_launch = Orbit::from_vectors_at(SUN, r0.into(), [0.0; 3], epoch_launch, Plane::EarthEquator);
    let orb_flyby1 = Orbit::from_vectors_at(SUN, r1.into(), [0.0; 3], epoch_flyby1, Plane::EarthEquator);
    
    let lambert = Maneuver::lambert(&orb_launch, &orb_flyby1).unwrap();
    // Lambert returns delta_v relative to initial orbit velocity.
    // orb_launch velocity was [0,0,0], so delta_v = v_departure.
    let v_dep = lambert.impulses[0].1; 
    
    // Second impulse is delta_v relative to final orbit velocity.
    // orb_flyby1 velocity was [0,0,0], so delta_v = v_final - v_arrival
    // => v_arrival = v_final - delta_v = 0 - delta_v = -delta_v
    let v_arr_vec = -lambert.impulses[1].1; // Arrival velocity on transfer orbit (relative to 0 target)
    
    let v1_pre = Vector3::new(v_arr_vec.x, v_arr_vec.y, v_arr_vec.z);
    
    println!("Velocity at arrival (pre-flyby): {:.3} km/s", v1_pre.norm());
    println!("Velocity of Venus: {:.3} km/s", Vector3::new(v_venus[0], v_venus[1], v_venus[2]).norm());
    
    // 3. Flyby #1 around Venus
    println!("\n## 3. Flyby #1 around Venus");
    
    let h = 2548.0; // Altitude km
    let d_flyby = VENUS.mean_radius_km + h;
    
    let v_venus_vec = Vector3::new(v_venus[0], v_venus[1], v_venus[2]);
    
    // Compute flyby with 0 theta (default B-plane angle)
    match compute_flyby(v1_pre, v_venus_vec, VENUS.mu_km3_s2, d_flyby, 0.0) {
        Ok((v_out, delta)) => {
            println!("Flyby (theta=0):");
            println!("  Exit Velocity: {:.3} km/s", v_out.norm());
            println!("  Turn Angle: {:.2} deg", delta.to_degrees());
        },
        Err(e) => println!("Flyby calculation failed: {:?}", e),
    }

    // 4. Optimization
    println!("\n## 4. Optimization");
    
    // We want period of post-flyby orbit to be T_ref (150 days)
    // Function to minimize: Period(theta) - T_ref
    
    // Simple grid search to find bracket
    let mut best_theta = 0.0;
    let mut min_diff = 1e9;
    
    // Search 0 to 2PI
    let steps = 100;
    for i in 0..steps {
        let theta = 2.0 * PI * (i as f64) / (steps as f64);
        
        // Calculate period diff
        let diff = match compute_flyby(v1_pre, v_venus_vec, VENUS.mu_km3_s2, d_flyby, theta) {
            Ok((v_out, _)) => {
                let r1_vec = Vector3::new(r1[0], r1[1], r1[2]);
                let orb = Orbit::from_vectors(SUN, r1_vec.into(), [v_out.x, v_out.y, v_out.z]);
                (orb.period_seconds().unwrap_or(0.0) - t_ref).abs()
            },
            Err(_) => 1e9,
        };
        
        if diff < min_diff {
            min_diff = diff;
            best_theta = theta;
        }
    }
    
    println!("Best theta from grid search: {:.2} rad ({:.2} deg), Diff: {:.2} s", 
             best_theta, best_theta.to_degrees(), min_diff);
             
    // Refine with simple bisection if we can find a bracket
    // Let's just use the grid search result as "close enough" for this demo
    
    // Recalculate with best theta
    if let Ok((v_out_opt, _)) = compute_flyby(v1_pre, v_venus_vec, VENUS.mu_km3_s2, d_flyby, best_theta) {
        let r1_vec = Vector3::new(r1[0], r1[1], r1[2]);
        let orb_opt = Orbit::from_vectors(SUN, r1_vec.into(), [v_out_opt.x, v_out_opt.y, v_out_opt.z]);
        
        println!("Optimized Orbit #2:");
        if let Some(period) = orb_opt.period_seconds() {
            println!("  Period: {:.2} days (Target: 150.0)", period / 86400.0);
        } else {
            println!("  Period: Infinite");
        }
        println!("  Semimajor Axis: {:.3} AU", orb_opt.a_km() / au_km);
        println!("  Eccentricity: {:.4}", orb_opt.ecc());
        println!("  Inclination: {:.2} deg", orb_opt.inc_rad().to_degrees());
    } else {
        println!("Optimization failed to produce valid flyby.");
    }

}

fn julian_date(year: i32, month: i32, day: i32, hour: f64) -> f64 {
    let (y, m) = if month > 2 {
        (year as f64, month as f64)
    } else {
        ((year - 1) as f64, (month + 12) as f64)
    };
    let d = day as f64 + hour / 24.0;
    let a = (y / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();
    
    (365.25 * (y + 4716.0)).floor() + (30.6001 * (m + 1.0)).floor() + d + b - 1524.5
}

// Simple analytical ephemeris for Earth and Venus (Mean Elements J2000)
// Very approximate!
fn get_planet_orbit(body: Body, epoch_tdb_seconds: f64) -> Orbit {
    let au = 149_597_870.7;
    // Mean elements (approximate, from Wikipedia/JPL)
    // a (km), e, i (rad), Omega (rad), w (rad), M0 (rad), n (rad/day)
    let (a, e, i, raan, argp, m0, n) = match body.name {
        "Earth" => (
            1.00000 * au, 
            0.01671, 
            0.00005_f64.to_radians(), 
            -11.26_f64.to_radians(), 
            102.94_f64.to_radians(), 
            357.517_f64.to_radians(), 
            0.9856_f64.to_radians()
        ), // deg/day
        "Venus" => (
            0.72333 * au, 
            0.00677, 
            3.39471_f64.to_radians(), 
            76.68_f64.to_radians(), 
            54.89_f64.to_radians(), 
            48.005_f64.to_radians(), 
            1.6021_f64.to_radians()
        ),
        _ => (1.0 * au, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    };

    let dt_days = epoch_tdb_seconds / 86400.0;
    // Mean anomaly at epoch
    let mean_anomaly = (m0 + n * dt_days).rem_euclid(2.0 * PI);
    
    // Solve Kepler's Equation for Eccentric Anomaly (E)
    // M = E - e sin E
    let mut eccentric_anomaly = mean_anomaly;
    for _ in 0..10 {
        eccentric_anomaly = mean_anomaly + e * eccentric_anomaly.sin();
    }
    
    // True Anomaly (nu) from Eccentric Anomaly (E)
    // tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2)
    let nu_rad = 2.0 * ((1.0 + e).sqrt() / (1.0 - e).sqrt() * (eccentric_anomaly / 2.0).tan()).atan();
    
    // Semilatus rectum p = a * (1 - e^2)
    let p_km = a * (1.0 - e * e);
    
    Orbit::from_classical_at(SUN, ClassicalElements {
        p_km,
        ecc: e,
        inc_rad: i,
        raan_rad: raan,
        argp_rad: argp,
        nu_rad,
    }, epoch_tdb_seconds, Plane::EarthEquator).unwrap()
}