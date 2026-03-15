use poliastrs::bodies::{SUN, Body};
use poliastrs::twobody::orbit::Orbit;
use poliastrs::core::elements::ClassicalElements;
use poliastrs::core::angles::{m_to_e, e_to_nu};
use poliastrs::frames::Plane;
use poliastrs::plotting::porkchop::PorkchopPlotter;
use std::path::Path;
use chrono::{TimeZone, Utc};

// Constants
const AU_KM: f64 = 149_597_870.7;

fn mean_to_true(mean: f64, ecc: f64) -> f64 {
    let e_anom = m_to_e(mean, ecc);
    e_to_nu(e_anom, ecc)
}

fn get_j2000_orbit(
    attractor: Body,
    a_au: f64,
    ecc: f64,
    inc_deg: f64,
    raan_deg: f64,
    argp_deg: f64,
    mean_anom_deg: f64,
) -> Orbit {
    let a_km = a_au * AU_KM;
    let mean_rad = mean_anom_deg.to_radians();
    let nu_rad = mean_to_true(mean_rad, ecc);
    
    Orbit::from_classical_at(
        attractor,
        ClassicalElements {
            p_km: a_km * (1.0 - ecc * ecc),
            ecc,
            inc_rad: inc_deg.to_radians(),
            raan_rad: raan_deg.to_radians(),
            argp_rad: argp_deg.to_radians(),
            nu_rad,
        },
        0.0, // J2000 epoch (approx)
        Plane::EarthEcliptic,
    ).unwrap()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Generating Porkchop Plot for Earth-Mars 2005");

    // Earth J2000
    let earth_j2000 = get_j2000_orbit(
        SUN, 
        1.00000011, 
        0.01671022, 
        0.00005, 
        -11.26064, 
        114.20783, 
        -2.48284
    );
    
    // Mars J2000
    let mars_j2000 = get_j2000_orbit(
        SUN,
        1.523679,
        0.093405,
        1.8497,
        49.558,
        286.502,
        19.412
    );

    // Time spans
    let j2000 = Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap();
    
    let launch_start = Utc.with_ymd_and_hms(2005, 4, 30, 0, 0, 0).unwrap();
    let launch_end = Utc.with_ymd_and_hms(2005, 10, 7, 0, 0, 0).unwrap();
    
    let arr_start = Utc.with_ymd_and_hms(2005, 11, 16, 0, 0, 0).unwrap();
    let arr_end = Utc.with_ymd_and_hms(2006, 12, 21, 0, 0, 0).unwrap();
    
    let t_launch_start = (launch_start - j2000).num_milliseconds() as f64 / 1000.0;
    let t_launch_end = (launch_end - j2000).num_milliseconds() as f64 / 1000.0;
    
    let t_arr_start = (arr_start - j2000).num_milliseconds() as f64 / 1000.0;
    let t_arr_end = (arr_end - j2000).num_milliseconds() as f64 / 1000.0;

    // Create Plotter
    let plotter = PorkchopPlotter::new(
        earth_j2000,
        mars_j2000,
        (t_launch_start, t_launch_end),
        (t_arr_start, t_arr_end),
    );
    
    // Plot
    // The implementation in src/plotting/porkchop.rs seems to take a Path and do everything inside.
    // Need to verify if `plot` method is public.
    // Wait, the struct fields are private? No, defined in `src/plotting/porkchop.rs` as:
    // pub struct PorkchopPlotter { ... }
    // but fields are not pub. So I must use `new`.
    // `new` is public.
    // `plot` is public.
    
    // Wait, check src/plotting/porkchop.rs again for `pub fn plot`.
    // Yes, `pub fn plot(&self, path: &Path)`.
    
    // Wait, one issue: `PorkchopPlotter` inside `src/plotting/porkchop.rs` is:
    // `use crate::iod::lambert::izzo;`
    // And inside `plot`: `izzo(...)`.
    // It seems it does not expose the raw data calculation, just plots it.
    
    plotter.plot(Path::new("porkchop_mars.png"))?; // Changed filename to match example request

    println!("Generated porkchop_mars.png");

    Ok(())
}
