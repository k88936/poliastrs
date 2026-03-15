use crate::bodies::Body;
use crate::twobody::mean_elements::get_mean_elements;
use plotters::prelude::*;
use std::f64::consts::PI;
use std::path::Path;

#[derive(Clone, Copy)]
pub enum TisserandKind {
    Apsis,
    Energy,
    Period,
}

pub struct TisserandPlotter {
    kind: TisserandKind,
    lines: Vec<(Vec<(f64, f64)>, String, RGBAColor)>,
}

impl TisserandPlotter {
    pub fn new(kind: TisserandKind) -> Self {
        Self { kind, lines: Vec::new() }
    }

    pub fn plot(&mut self, body: Body, vinf_span: (f64, f64), num_contours: usize, label: Option<&str>) {
        if let Ok(mean_el) = get_mean_elements(body) {
            let mu_sun = 132_712_440_018.0; // SUN.mu_km3_s2
            let r_body_val = mean_el.a_km;
            let v_body_val = (mu_sun / r_body_val).sqrt();
            
            let (v_min, v_max) = vinf_span;
            for i in 0..num_contours {
                let v_inf = if num_contours == 1 {
                    v_min
                } else {
                    v_min + (v_max - v_min) * (i as f64 / (num_contours as f64 - 1.0))
                };
                
                let v_inf_norm = v_inf / v_body_val;
                
                let mut points = Vec::new();
                let n_alpha = 100;
                for j in 0..n_alpha {
                    let alpha = PI * (j as f64 / (n_alpha as f64 - 1.0));
                    
                    let denom = (1.0 - v_inf_norm * v_inf_norm - 2.0 * v_inf_norm * alpha.cos()).abs();
                    if denom < 1e-6 { continue; }
                    let a_sc = 1.0 / denom; // dimensionless
                    
                    let term = (3.0 - 1.0 / a_sc - v_inf_norm * v_inf_norm) / 2.0;
                    let ecc_sc = (1.0 - 1.0 / a_sc * term * term).max(0.0).sqrt();
                    
                    let r_p = a_sc * r_body_val * (1.0 - ecc_sc);
                    let r_a = a_sc * r_body_val * (1.0 + ecc_sc);
                    let energy = -mu_sun / (2.0 * a_sc * r_body_val);
                    let period = 2.0 * PI * ((a_sc * r_body_val).powi(3) / mu_sun).sqrt();

                    let x = match self.kind {
                        TisserandKind::Apsis => r_a,
                        _ => r_p,
                    };
                    
                    let y = match self.kind {
                        TisserandKind::Apsis => r_p,
                        TisserandKind::Energy => energy,
                        TisserandKind::Period => period,
                    };
                    
                    points.push((x, y));
                }
                
                let l_str = label.unwrap_or(body.name).to_string();
                let color = Palette99::pick(self.lines.len() % 99).to_rgba();
                self.lines.push((points, l_str, color));
            }
        }
    }

    pub fn save(&self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let root = BitMapBackend::new(path, (800, 600)).into_drawing_area();
        root.fill(&WHITE)?;

        if self.lines.is_empty() {
            return Ok(());
        }

        // Find ranges
        let mut min_x = f64::MAX;
        let mut max_x = f64::MIN;
        let mut min_y = f64::MAX;
        let mut max_y = f64::MIN;

        for (points, _, _) in &self.lines {
            for (x, y) in points {
                if x.is_finite() && y.is_finite() {
                    min_x = min_x.min(*x);
                    max_x = max_x.max(*x);
                    min_y = min_y.min(*y);
                    max_y = max_y.max(*y);
                }
            }
        }
        
        // Log scale needs positive values
        if min_x <= 0.0 { min_x = 1e-3; } // fallback

        let mut chart = ChartBuilder::on(&root)
            .caption("Tisserand Plot", ("sans-serif", 30))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(
                (min_x..max_x).log_scale(),
                min_y..max_y // Energy is negative, so log scale won't work for Y if Energy
            )?;

        chart.configure_mesh()
            .x_desc("Rp (km)")
            .y_desc(match self.kind {
                TisserandKind::Energy => "Energy (km^2/s^2)",
                TisserandKind::Period => "Period (s)",
                TisserandKind::Apsis => "Ra (km)",
            })
            .draw()?;

        for (points, label, color) in &self.lines {
            chart.draw_series(LineSeries::new(points.clone(), color.clone()))?
                .label(label)
                .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.clone()));
        }

        chart.configure_series_labels()
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw()?;

        Ok(())
    }
}
