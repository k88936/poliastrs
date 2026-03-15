use crate::plotting::orbit_plotter::OrbitPlotter;
use crate::bodies::{MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, SUN};
use crate::twobody::mean_elements::get_mean_elements;
use crate::twobody::orbit::Orbit;
use crate::core::elements::ClassicalElements;
use crate::frames::Plane;

pub fn plot_solar_system(outer: bool, epoch: f64) -> OrbitPlotter {
    let mut plotter = OrbitPlotter::new();
    plotter.set_attractor(SUN);

    let planets = if outer {
        vec![MARS, JUPITER, SATURN, URANUS, NEPTUNE]
    } else {
        vec![MERCURY, VENUS, EARTH, MARS]
    };

    for body in planets {
        if let Ok(els) = get_mean_elements(body) {
            let p = els.a_km * (1.0 - els.ecc * els.ecc);
            let orbit = Orbit::from_classical_at(
                SUN,
                ClassicalElements {
                    p_km: p,
                    ecc: els.ecc,
                    inc_rad: els.inc_rad,
                    raan_rad: 0.0,
                    argp_rad: 0.0,
                    nu_rad: 0.0,
                },
                epoch,
                Plane::EarthEcliptic
            ).unwrap();
            
            plotter.plot(&orbit, Some(body.name));
        }
    }
    
    plotter
}
