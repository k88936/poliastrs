use crate::twobody::orbit::Orbit;
use plotters::prelude::*;
use std::path::Path;

pub struct GabbardPlotter {
    orbits: Vec<(Orbit, String)>,
    pub title: String,
    pub dark_mode: bool,
}

impl GabbardPlotter {
    pub fn new(dark_mode: bool) -> Self {
        Self {
            orbits: Vec::new(),
            title: "Gabbard Plot".to_string(),
            dark_mode,
        }
    }

    pub fn plot_orbits(&mut self, orbits: &[Orbit], label: Option<&str>) {
        let l = label.unwrap_or("Orbits");
        for orbit in orbits {
            self.orbits.push((*orbit, l.to_string()));
        }
    }

    pub fn save(&self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let root = BitMapBackend::new(path, (800, 600)).into_drawing_area();
        
        let bg_color = if self.dark_mode { BLACK } else { WHITE };
        let fg_color = if self.dark_mode { WHITE } else { BLACK };
        
        root.fill(&bg_color)?;

        if self.orbits.is_empty() {
            return Ok(());
        }

        // Find ranges
        let mut min_period = f64::MAX;
        let mut max_period = f64::MIN;
        let mut min_alt = f64::MAX;
        let mut max_alt = f64::MIN;

        for (orbit, _) in &self.orbits {
             if let Some(period) = orbit.period_seconds() {
                 let p = period / 60.0; // minutes
                 let r_p = orbit.r_p_km();
                 let r_a = orbit.r_a_km().unwrap_or(r_p);
                 let mean_r = orbit.attractor.mean_radius_km;
                 
                 let alt_p = r_p - mean_r;
                 let alt_a = r_a - mean_r;

                 min_period = min_period.min(p);
                 max_period = max_period.max(p);
                 min_alt = min_alt.min(alt_p).min(alt_a);
                 max_alt = max_alt.max(alt_p).max(alt_a);
             }
        }
        
        if min_period == f64::MAX {
             // No valid orbits
             return Ok(());
        }

        // Add some margin
        let p_margin = (max_period - min_period).max(1.0) * 0.1;
        min_period -= p_margin;
        max_period += p_margin;
        
        let a_margin = (max_alt - min_alt).max(1.0) * 0.1;
        min_alt -= a_margin;
        max_alt += a_margin;

        let mut chart = ChartBuilder::on(&root)
            .caption(&self.title, ("sans-serif", 30).into_font().color(&fg_color))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(min_period..max_period, min_alt..max_alt)?;

        chart.configure_mesh()
            .x_desc("Period (min)")
            .y_desc("Altitude (km)")
            .axis_style(fg_color)
            .label_style(("sans-serif", 15).into_font().color(&fg_color))
            .draw()?;

        let mut unique_labels: Vec<String> = self.orbits.iter().map(|(_, l)| l.clone()).collect();
        unique_labels.sort();
        unique_labels.dedup();
        
        for (idx, label) in unique_labels.iter().enumerate() {
            let color = Palette99::pick(idx);
            let style = ShapeStyle {
                color: color.to_rgba(),
                filled: true,
                stroke_width: 1,
            };
            
            let points: Vec<(f64, f64)> = self.orbits.iter()
                .filter(|(_, l)| l == label)
                .filter_map(|(o, _)| {
                    if let Some(period) = o.period_seconds() {
                        let p = period / 60.0;
                        let r_p = o.r_p_km();
                        let r_a = o.r_a_km().unwrap_or(r_p);
                        let mean_r = o.attractor.mean_radius_km;
                        Some(vec![(p, r_p - mean_r), (p, r_a - mean_r)])
                    } else {
                        None
                    }
                })
                .flatten()
                .collect();
                
             chart.draw_series(
                points.iter().map(|(x, y)| Circle::new((*x, *y), 3, style.clone()))
             )?.label(label).legend(move |(x, y)| Circle::new((x, y), 3, style.clone()));
        }
        
        chart.configure_series_labels()
            .background_style(&bg_color.mix(0.8))
            .border_style(&fg_color)
            .label_font(("sans-serif", 15).into_font().color(&fg_color))
            .draw()?;

        Ok(())
    }
}
