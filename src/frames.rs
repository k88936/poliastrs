#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Plane {
    EarthEquator,
    EarthEcliptic,
}

#[cfg(test)]
mod tests {
    use super::Plane;

    #[test]
    fn planes_are_distinct() {
        assert_ne!(Plane::EarthEquator, Plane::EarthEcliptic);
    }
}
