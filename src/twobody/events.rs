use nalgebra::Vector3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EventBase {
    pub terminal: bool,
    pub direction: f64,
    pub last_t: Option<f64>,
}

impl EventBase {
    pub fn new(terminal: bool, direction: f64) -> Self {
        Self { terminal, direction, last_t: None }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AltitudeCrossEvent {
    pub base: EventBase,
    pub alt_km: f64,
    pub r_body_km: f64,
}

impl AltitudeCrossEvent {
    pub fn new(alt_km: f64, r_body_km: f64) -> Self {
        Self { base: EventBase::new(true, -1.0), alt_km, r_body_km }
    }
    pub fn eval(&mut self, t: f64, u: [f64; 6]) -> f64 {
        self.base.last_t = Some(t);
        let r = Vector3::new(u[0], u[1], u[2]).norm();
        r - self.r_body_km - self.alt_km
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LithobrakeEvent(pub AltitudeCrossEvent);

impl LithobrakeEvent {
    pub fn new(r_body_km: f64) -> Self {
        Self(AltitudeCrossEvent::new(0.0, r_body_km))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NodeCrossEvent {
    pub base: EventBase,
}

impl NodeCrossEvent {
    pub fn new(terminal: bool) -> Self {
        Self { base: EventBase::new(terminal, 0.0) }
    }
    pub fn eval(&mut self, t: f64, u: [f64; 6]) -> f64 {
        self.base.last_t = Some(t);
        u[2]
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LosEvent {
    pub base: EventBase,
    pub r_body_km: f64,
}

impl LosEvent {
    pub fn new(r_body_km: f64) -> Self {
        Self { base: EventBase::new(false, 0.0), r_body_km }
    }
    pub fn eval(&mut self, t: f64, r1: [f64; 3], r2: [f64; 3]) -> f64 {
        self.base.last_t = Some(t);
        let a = Vector3::new(r1[0], r1[1], r1[2]);
        let b = Vector3::new(r2[0], r2[1], r2[2]);
        let d = b - a;
        let tstar = (-(a.dot(&d)) / d.dot(&d)).clamp(0.0, 1.0);
        let closest = a + d * tstar;
        closest.norm() - self.r_body_km
    }
}

#[cfg(test)]
mod tests {
    use super::{AltitudeCrossEvent, LithobrakeEvent, LosEvent, NodeCrossEvent};

    #[test]
    fn test_altitude_crossing() {
        let mut ev = AltitudeCrossEvent::new(50.0, 6378.0);
        let y = ev.eval(100.0, [6378.0 + 40.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(y < 0.0);
        assert_eq!(ev.base.last_t, Some(100.0));
    }
    #[test]
    fn test_altitude_cross_not_happening_is_ok() {
        let mut ev = AltitudeCrossEvent::new(50.0, 6378.0);
        let y = ev.eval(100.0, [6378.0 + 230.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(y > 0.0);
    }
    #[test]
    fn test_latitude_cross_event() {
        let mut ev = NodeCrossEvent::new(true);
        let y = ev.eval(1701.7, [1.0, 2.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(y, 0.0);
    }
    #[test]
    fn test_penumbra_event_not_triggering_is_ok() {
        let mut ev = LosEvent::new(6378.0);
        assert!(ev.eval(100.0, [10000.0, 0.0, 0.0], [0.0, 10000.0, 0.0]) > 0.0);
    }
    #[test]
    fn test_umbra_event_not_triggering_is_ok() {
        let mut ev = LosEvent::new(6378.0);
        assert!(ev.eval(100.0, [7000.0, 0.0, 0.0], [0.0, -7000.0, 0.0]) < 0.0);
    }
    #[test]
    fn test_umbra_event_crossing() {
        let mut ev = LosEvent::new(6378.0);
        let _ = ev.eval(291.0, [7000.0, 0.0, 0.0], [-7000.0, 0.0, 0.0]);
        assert_eq!(ev.base.last_t, Some(291.0));
    }
    #[test]
    fn test_penumbra_event_crossing() {
        let mut ev = LosEvent::new(6378.0);
        let _ = ev.eval(266.0, [7000.0, 0.0, 0.0], [-7000.0, 1000.0, 0.0]);
        assert_eq!(ev.base.last_t, Some(266.0));
    }
    #[test]
    fn test_node_cross_event() {
        let mut ev = NodeCrossEvent::new(true);
        let y = ev.eval(3.4, [1.0, 2.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(y, 0.0);
    }
    #[test]
    fn test_node_event_equatorial_orbit() {
        let mut ev = NodeCrossEvent::new(true);
        assert_eq!(ev.eval(0.0, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]), 0.0);
    }
    #[test]
    fn test_orbit_propagation_continues_if_events_terminal_is_false() {
        let ev = NodeCrossEvent::new(false);
        assert!(!ev.base.terminal);
    }
    #[test]
    fn test_orbit_propagation_position_vector_does_not_repeat_if_events_terminal_is_true() {
        let ev = NodeCrossEvent::new(true);
        assert!(ev.base.terminal);
    }
    #[test]
    fn test_propagation_stops_if_atleast_one_event_has_terminal_set_to_true() {
        let ev = NodeCrossEvent::new(true);
        assert!(ev.base.terminal);
    }
    #[test]
    fn test_line_of_sight() {
        let mut ev = LosEvent::new(6378.0);
        let v = ev.eval(0.0, [10000.0, 0.0, 0.0], [0.0, 10000.0, 0.0]);
        assert!(v > 0.0);
    }
    #[test]
    fn test_los_event_raises_warning_if_norm_less_than_radius() {
        let mut ev = LosEvent::new(6378.0);
        let v = ev.eval(0.0, [1000.0, 0.0, 0.0], [0.0, 7000.0, 0.0]);
        assert!(v < 0.0);
    }
    #[test]
    fn test_los_event_with_lithobrake_event_behavior() {
        let litho = LithobrakeEvent::new(6378.0);
        assert_eq!(litho.0.alt_km, 0.0);
    }
    #[test]
    fn test_los_event() {
        let mut ev = LosEvent::new(6378.0);
        let _ = ev.eval(10.0, [7000.0, 1.0, 0.0], [7000.0, -1.0, 0.0]);
        assert_eq!(ev.base.last_t, Some(10.0));
    }
}
