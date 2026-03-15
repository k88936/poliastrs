/// Generates a range of times (f64 seconds).
/// 
/// Arguments:
/// * `start` - Start time in seconds
/// * `end` - End time in seconds (optional)
/// * `spacing` - Spacing in seconds (optional)
/// * `num_values` - Number of values (optional)
/// 
/// Returns a vector of times.
pub fn time_range(
    start: f64,
    end: Option<f64>,
    spacing: Option<f64>,
    num_values: Option<usize>,
) -> Result<Vec<f64>, String> {
    match (end, spacing, num_values) {
        (Some(e), None, Some(n)) => {
            if n < 2 {
                return Ok(vec![start]); // Or error?
            }
            let step = (e - start) / (n as f64 - 1.0);
            let mut res = Vec::with_capacity(n);
            for i in 0..n {
                res.push(start + i as f64 * step);
            }
            Ok(res)
        }
        (None, Some(s), Some(n)) => {
            let mut res = Vec::with_capacity(n);
            for i in 0..n {
                res.push(start + i as f64 * s);
            }
            Ok(res)
        }
        (Some(e), Some(s), None) => {
            let n = ((e - start) / s).floor() as usize + 1;
            let mut res = Vec::with_capacity(n);
            for i in 0..n {
                let t = start + i as f64 * s;
                if t > e + 1e-9 { break; } // Floating point safety
                res.push(t);
            }
            Ok(res)
        }
        _ => Err("Invalid arguments: need (end, num_values), (spacing, num_values), or (end, spacing)".to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_time_range_spacing_num_values() {
        let start = 0.0;
        let spacing = 60.0; // 1 minute
        let num_values = 5;
        let expected_duration = 240.0; // 4 minutes

        // Case 1: spacing + num_values
        let res1 = time_range(start, None, Some(spacing), Some(num_values)).unwrap();
        assert_eq!(res1.len(), num_values);
        assert_relative_eq!(res1.last().unwrap() - res1.first().unwrap(), expected_duration);

        // Case 2: end + num_values
        let end = 240.0;
        let res2 = time_range(start, Some(end), None, Some(num_values)).unwrap();
        assert_eq!(res2.len(), num_values);
        assert_relative_eq!(res2.last().unwrap() - res2.first().unwrap(), expected_duration);
    }
    
    #[test]
    fn test_time_range_errors() {
        assert!(time_range(0.0, None, None, None).is_err());
        assert!(time_range(0.0, Some(10.0), Some(1.0), Some(10)).is_err()); // Ambiguous? 
        // Python implementation raises ValueError if too many args provided?
        // My implementation returns Err for the catch-all case.
    }
}
