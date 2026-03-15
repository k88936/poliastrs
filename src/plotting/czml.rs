use serde::Serialize;
use crate::twobody::orbit::Orbit;
use chrono::{DateTime, Utc, TimeZone};

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct DocumentPacket {
    pub id: String,
    pub version: String,
    pub name: String,
    pub clock: Clock,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Clock {
    pub interval: String,
    pub current_time: String,
    pub multiplier: f64,
    pub range: String,
    pub step: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Packet {
    pub id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub availability: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub position: Option<Position>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub path: Option<Path>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<Label>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Position {
    pub epoch: String,
    pub interpolation_algorithm: String,
    pub interpolation_degree: i32,
    pub reference_frame: String,
    pub cartesian: Vec<f64>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Path {
    pub width: f64,
    pub show: bool,
    pub material: Material,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Label {
    pub text: String,
    pub show: bool,
    pub fill_color: Color,
    pub outline_color: Color,
    pub font: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Material {
    pub solid_color: SolidColor,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SolidColor {
    pub color: Rgba,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Rgba {
    pub rgba: [u8; 4],
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Color {
    pub rgba: [u8; 4],
}

pub struct CZMLExtractor {
    start_epoch: DateTime<Utc>,
    end_epoch: DateTime<Utc>,
    n_samples: usize,
    packets: Vec<Packet>,
    document_packet: DocumentPacket,
}

impl CZMLExtractor {
    pub fn new(start_epoch: DateTime<Utc>, end_epoch: DateTime<Utc>, n_samples: usize) -> Self {
        let interval = format!("{}/{}", start_epoch.to_rfc3339(), end_epoch.to_rfc3339());
        
        let document = DocumentPacket {
            id: "document".to_string(),
            version: "1.0".to_string(),
            name: "document_packet".to_string(),
            clock: Clock {
                interval: interval.clone(),
                current_time: start_epoch.to_rfc3339(),
                multiplier: 60.0,
                range: "LOOP_STOP".to_string(),
                step: "SYSTEM_CLOCK_MULTIPLIER".to_string(),
            },
        };
        
        Self {
            start_epoch,
            end_epoch,
            n_samples,
            packets: vec![],
            document_packet: document,
        }
    }
    
    pub fn add_orbit(&mut self, mut orbit: Orbit, id_name: &str, path_width: f64, label_text: &str, color: [u8; 4]) {
        let duration = self.end_epoch - self.start_epoch;
        let total_seconds = duration.num_milliseconds() as f64 / 1000.0;
        let dt = total_seconds / (self.n_samples as f64 - 1.0);
        
        let mut cartesian = Vec::new();
        
        // Assume J2000 reference epoch for Orbit.epoch_tdb_seconds
        // J2000: 2000-01-01T12:00:00Z
        let j2000 = Utc.with_ymd_and_hms(2000, 1, 1, 12, 0, 0).unwrap();
        let start_seconds_from_j2000 = (self.start_epoch - j2000).num_milliseconds() as f64 / 1000.0;
        
        // Propagate orbit to start time
        let delta_to_start = start_seconds_from_j2000 - orbit.epoch_tdb_seconds;
        orbit = orbit.propagate_seconds(delta_to_start).unwrap();
        
        for i in 0..self.n_samples {
            let current_dt = i as f64 * dt;
            
            // We propagate from the initial state (at start_epoch) by current_dt
            let state_at_t = orbit.propagate_seconds(current_dt).unwrap();
            
            // Convert to meters for CZML
            let x = state_at_t.state.r_km[0] * 1000.0;
            let y = state_at_t.state.r_km[1] * 1000.0;
            let z = state_at_t.state.r_km[2] * 1000.0;
            
            cartesian.push(current_dt);
            cartesian.push(x);
            cartesian.push(y);
            cartesian.push(z);
        }
        
        let position = Position {
            epoch: self.start_epoch.to_rfc3339(),
            interpolation_algorithm: "LAGRANGE".to_string(),
            interpolation_degree: 5,
            reference_frame: "INERTIAL".to_string(),
            cartesian,
        };
        
        let path = Path {
            width: path_width,
            show: true,
            material: Material {
                solid_color: SolidColor {
                    color: Rgba { rgba: color },
                },
            },
        };
        
        let label = Label {
            text: label_text.to_string(),
            show: true,
            fill_color: Color { rgba: color },
            outline_color: Color { rgba: [0, 0, 0, 255] },
            font: "12pt Lucida Console".to_string(),
        };
        
        let packet = Packet {
            id: id_name.to_string(),
            name: Some(id_name.to_string()),
            availability: Some(format!("{}/{}", self.start_epoch.to_rfc3339(), self.end_epoch.to_rfc3339())),
            position: Some(position),
            path: Some(path),
            label: Some(label),
        };
        
        self.packets.push(packet);
    }

    pub fn packets(&self) -> Vec<serde_json::Value> {
        let mut result = Vec::new();
        result.push(serde_json::to_value(&self.document_packet).unwrap());
        for p in &self.packets {
            result.push(serde_json::to_value(p).unwrap());
        }
        result
    }
}
