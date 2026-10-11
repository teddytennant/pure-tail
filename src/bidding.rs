//! Online bidding, Section 3 and the proof of Lemma 7.
//!
//! `BID_r` draws ξ uniformly from [0, 1) and submits r^ξ, r^{ξ+1}, r^{ξ+2}, ...
//! until the first bid that meets the unknown threshold T ≥ 1. The cost is the
//! sum of the submitted bids. OPT is T.

/// γ(r) = r^2 / (r - 1), for r in [2, e].
pub fn gamma(r: f64) -> f64 {
    r * r / (r - 1.0)
}

/// Theorem 1: expected competitive ratio at most r / ln r.
pub fn expected_ratio_bound(r: f64) -> f64 {
    r / r.ln()
}

/// Run `BID_r`. Returns (cost, last bid, number of bids).
pub fn bid(r: f64, xi: f64, threshold: f64) -> (f64, f64, usize) {
    assert!(r > 1.0, "r must be > 1");
    assert!((0.0..1.0).contains(&xi), "xi must lie in [0, 1)");
    assert!(threshold >= 1.0, "threshold must be >= 1");
    let mut exponent = xi;
    let mut cost = 0.0;
    let mut last;
    let mut n = 0usize;
    loop {
        last = r.powf(exponent);
        cost += last;
        n += 1;
        if last >= threshold {
            return (cost, last, n);
        }
        exponent += 1.0;
        if n > 10_000 {
            panic!("bid sequence did not reach the threshold");
        }
    }
}

pub fn competitive_ratio(r: f64, xi: f64, threshold: f64) -> f64 {
    let (cost, _, _) = bid(r, xi, threshold);
    cost / threshold
}

/// ζ = log_r(last bid) - log_r(T), which lies in [0, 1).
///
/// The proof of Lemma 7 uses this, not the fractional part of log_r(last bid)
/// alone. The main text's wording and the proof differ; the proof is the one
/// that produces r^{ζ+1} / (r - 1).
pub fn zeta(r: f64, last_bid: f64, threshold: f64) -> f64 {
    last_bid.log(r) - threshold.log(r)
}

/// Equation (9): quality is at most r^{ζ+1} / (r - 1).
pub fn quality_bound(r: f64, zeta: f64) -> f64 {
    r.powf(zeta + 1.0) / (r - 1.0)
}
