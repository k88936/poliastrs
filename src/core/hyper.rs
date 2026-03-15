pub fn hyp2f1b(x: f64) -> f64 {
    if x >= 1.0 {
        f64::INFINITY
    } else {
        let mut res = 1.0;
        let mut term = 1.0;
        let mut ii = 0;
        loop {
            term = term * (3.0 + ii as f64) * (1.0 + ii as f64)
                / (2.5 + ii as f64)
                * x
                / (ii as f64 + 1.0);
            let res_old = res;
            res += term;
            if (res_old - res).abs() < 1e-15 {
                return res;
            }
            ii += 1;
            if ii > 10000 {
                 // Safety break to prevent infinite loops in edge cases
                 break;
            }
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_hyp2f1b_values() {
        let test_cases = vec![
            (0.0, 1.0),
            (0.1, 1.13542436662003),
            (0.2, 1.310435208989098),
            (0.3, 1.5444243078411781),
            (0.4, 1.8713706568239126),
            (0.5, 2.356194490192345),
            (0.6, 3.138589617245872),
            (0.7, 4.576573199723625),
            (0.8, 7.893449518324751),
            (0.9, 20.68119128330911),
        ];

        for (x, expected) in test_cases {
            let result = hyp2f1b(x);
            assert_relative_eq!(result, expected, epsilon = 1e-12);
        }
    }
}
