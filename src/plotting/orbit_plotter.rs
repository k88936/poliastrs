use crate::bodies::Body;
use crate::twobody::orbit::Orbit;
// use crate::ephem::Ephem; // Not used yet
use crate::twobody::sampling::TrueAnomalyBounds;
use plotters::prelude::*;
use std::path::Path;

pub struct OrbitPlotter {
    pub attractor: Option<Body>,
    pub trajectories: Vec<(Vec<(f64, f64)>, String, RGBAColor)>,
    pub trajectories_3d: Vec<(Vec<(f64, f64, f64)>, String, RGBAColor)>,
    pub dark_mode: bool,
}

impl OrbitPlotter {
    pub fn new() -> Self {
        Self {
            attractor: None,
            trajectories: Vec::new(),
            trajectories_3d: Vec::new(),
            dark_mode: false,
        }
    }

    pub fn set_attractor(&mut self, body: Body) {
        if let Some(attractor) = self.attractor {
            if attractor != body {
                panic!("Attractor already set to different body");
            }
        } else {
            self.attractor = Some(body);
        }
    }

    pub fn plot(&mut self, orbit: &Orbit, label: Option<&str>) {
        if self.attractor.is_none() {
            self.attractor = Some(orbit.attractor);
        } else if self.attractor != Some(orbit.attractor) {
            panic!("Orbit attractor does not match plotter attractor");
        }

        let sampler = TrueAnomalyBounds::default();
        let (coords, _) = sampler.sample(*orbit);
        
        let label = label.unwrap_or("Orbit").to_string();
        let color = Palette99::pick(self.trajectories.len()).to_rgba();

        let points_2d: Vec<(f64, f64)> = coords.iter().map(|p| (p[0], p[1])).collect();
        self.trajectories.push((points_2d, label.clone(), color));
        
        let points_3d: Vec<(f64, f64, f64)> = coords.iter().map(|p| (p[0], p[1], p[2])).collect();
        self.trajectories_3d.push((points_3d, label, color));
    }

    pub fn save_2d(&self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let root = BitMapBackend::new(path, (800, 600)).into_drawing_area();
        
        let bg_color = if self.dark_mode { BLACK } else { WHITE };
        let fg_color = if self.dark_mode { WHITE } else { BLACK };
        
        root.fill(&bg_color)?;

        if self.trajectories.is_empty() {
            return Ok(());
        }

        // Find ranges
        let mut min_x = f64::MAX;
        let mut max_x = f64::MIN;
        let mut min_y = f64::MAX;
        let mut max_y = f64::MIN;

        for (points, _, _) in &self.trajectories {
            for (x, y) in points {
                min_x = min_x.min(*x);
                max_x = max_x.max(*x);
                min_y = min_y.min(*y);
                max_y = max_y.max(*y);
            }
        }
        
        // Add margin
        let margin_x = (max_x - min_x).max(1.0) * 0.1;
        let margin_y = (max_y - min_y).max(1.0) * 0.1;
        min_x -= margin_x;
        max_x += margin_x;
        min_y -= margin_y;
        max_y += margin_y;
        
        // Make aspect ratio 1:1 approximately
        let range_x = max_x - min_x;
        let range_y = max_y - min_y;
        if range_x > range_y {
            let diff = range_x - range_y;
            min_y -= diff / 2.0;
            max_y += diff / 2.0;
        } else {
            let diff = range_y - range_x;
            min_x -= diff / 2.0;
            max_x += diff / 2.0;
        }

        let mut chart = ChartBuilder::on(&root)
            .caption("Orbit Plot", ("sans-serif", 30).into_font().color(&fg_color))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(min_x..max_x, min_y..max_y)?;

        chart.configure_mesh()
            .x_desc("x (km)")
            .y_desc("y (km)")
            .axis_style(fg_color)
            .label_style(("sans-serif", 15).into_font().color(&fg_color))
            .draw()?;

        // Draw attractor if set
        if let Some(attractor) = self.attractor {
            let color = &fg_color; // or specific body color
            chart.draw_series(std::iter::once(
                Circle::new((0.0, 0.0), 5, ShapeStyle::from(color).filled())
            ))?.label(attractor.name).legend(|(x, y)| Circle::new((x, y), 5, ShapeStyle::from(&BLACK).filled()));
        }

        for (points, label, color) in &self.trajectories {
            chart.draw_series(
                LineSeries::new(points.clone(), color)
            )?.label(label).legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.clone()));
        }
        
        chart.configure_series_labels()
            .background_style(&bg_color.mix(0.8))
            .border_style(&fg_color)
            .label_font(("sans-serif", 15).into_font().color(&fg_color))
            .draw()?;

        Ok(())
    }
}
