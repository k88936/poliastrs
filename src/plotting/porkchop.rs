use crate::twobody::orbit::Orbit;
use crate::iod::lambert::izzo;
use plotters::prelude::*;
use std::path::Path;
use nalgebra::Vector3;

pub struct PorkchopPlotter {
    dep_orbit: Orbit,
    arr_orbit: Orbit,
    launch_span: (f64, f64), // Start, End seconds
    arrival_span: (f64, f64), // Start, End seconds
}

impl PorkchopPlotter {
    pub fn new(dep: Orbit, arr: Orbit, launch_span: (f64, f64), arrival_span: (f64, f64)) -> Self {
        Self {
            dep_orbit: dep,
            arr_orbit: arr,
            launch_span,
            arrival_span,
        }
    }

    pub fn plot(&self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let root = BitMapBackend::new(path, (800, 600)).into_drawing_area();
        root.fill(&WHITE)?;

        let n_launch = 50;
        let n_arrival = 50;

        let t_launch_start = self.launch_span.0;
        let t_launch_end = self.launch_span.1;
        let t_arr_start = self.arrival_span.0;
        let t_arr_end = self.arrival_span.1;

        let mut data_c3 = Vec::new();
        // let mut data_tof = Vec::new();

        let mu = self.dep_orbit.attractor.mu_km3_s2; // Assumes same attractor (Sun)

        for i in 0..n_launch {
            let t1 = t_launch_start + (t_launch_end - t_launch_start) * (i as f64 / (n_launch - 1) as f64);
            // Propagate departure body to t1
            // Use dep_orbit epoch as reference.
            let dt1 = t1 - self.dep_orbit.epoch_tdb_seconds;
            let dep_state = self.dep_orbit.propagate_seconds(dt1).unwrap();
            let r1 = Vector3::new(dep_state.state.r_km.x, dep_state.state.r_km.y, dep_state.state.r_km.z);
            let v1_body = Vector3::new(dep_state.state.v_km_s.x, dep_state.state.v_km_s.y, dep_state.state.v_km_s.z);

            for j in 0..n_arrival {
                let t2 = t_arr_start + (t_arr_end - t_arr_start) * (j as f64 / (n_arrival - 1) as f64);
                
                if t2 <= t1 {
                    continue;
                }

                let dt2 = t2 - self.arr_orbit.epoch_tdb_seconds;
                let arr_state = self.arr_orbit.propagate_seconds(dt2).unwrap();
                let r2 = Vector3::new(arr_state.state.r_km.x, arr_state.state.r_km.y, arr_state.state.r_km.z);
                let v2_body = Vector3::new(arr_state.state.v_km_s.x, arr_state.state.v_km_s.y, arr_state.state.v_km_s.z);

                let tof = t2 - t1;
                
                // Solve Lambert
                // izzo(mu, r1, r2, tof, 0, false, false, ...) ?
                // Prograde usually.
                // Assuming prograde = true (direct transfer)
                match izzo(mu, r1, r2, tof, 0, true, false, 35, 1e-8) {
                    Ok((v1_trans, v2_trans)) => {
                         let v_inf_dep = (v1_trans - v1_body).norm();
                         let c3 = v_inf_dep * v_inf_dep;
                         
                         // Store for contour
                         // We map (t1, t2) -> c3
                         data_c3.push((t1, t2, c3));
                    },
                    Err(_) => {
                        // Failed to solve
                    }
                }
            }
        }
        
        // Simple scatter plot for now as plotters contour support is manual or complex
        // We will plot points colored by C3
        
        let mut chart = ChartBuilder::on(&root)
            .caption("Porkchop Plot (C3 km^2/s^2)", ("sans-serif", 30))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(t_launch_start..t_launch_end, t_arr_start..t_arr_end)?;

        chart.configure_mesh()
            .x_desc("Launch Time (s)")
            .y_desc("Arrival Time (s)")
            .draw()?;

        // Color scale
        let max_c3 = data_c3.iter().map(|(_,_,c)| *c).fold(0.0/0.0, f64::max);
        let min_c3 = data_c3.iter().map(|(_,_,c)| *c).fold(0.0/0.0, f64::min);
        
        chart.draw_series(
            data_c3.iter().map(|(x, y, c)| {
                let v = (c - min_c3) / (max_c3 - min_c3);
                let color = HSLColor(0.7 - v * 0.7, 1.0, 0.5); // Blue to Red
                Circle::new((*x, *y), 3, ShapeStyle::from(color.filled()))
            })
        )?;

        Ok(())
    }
}
