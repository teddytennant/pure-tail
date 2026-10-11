//! Line search (the cow-path problem), Section 4 and Section E.
//!
//! `SMARTCOW_r` (Kao, Reif, Tate 1996) draws ξ uniformly from [0, 1) and a
//! fair starting direction, then walks to r^ξ, r^{ξ+1}, r^{ξ+2}, ... on
//! alternating rays until it reaches the target.

/// γ(r) = 1 + 2 r^2 / (r - 1), for r in [2, r0].
pub fn gamma(r: f64) -> f64 {
    1.0 + 2.0 * r * r / (r - 1.0)
}

/// Lemma 16: expected competitive ratio 1 + (1 + r) / ln r.
pub fn expected_ratio_bound(r: f64) -> f64 {
    1.0 + (1.0 + r) / r.ln()
}

/// r0 ≈ 3.59112, the unique solution of r = (1 + r) / ln r above 1.
pub fn r0() -> f64 {
    let mut r: f64 = 3.591121476668622;
    for _ in 0..8 {
        let ln = r.ln();
        r -= (r * ln - r - 1.0) / ln;
    }
    r
}

/// Distance traveled before the target at `t` on `target_positive` is found.
///
/// A trip that does not contain the target costs twice the turning distance
/// (out and back). The successful trip costs `t`.
pub fn cost(r: f64, xi: f64, start_positive: bool, target_positive: bool, t: f64) -> f64 {
    assert!(r > 1.0);
    assert!((0.0..1.0).contains(&xi));
    assert!(t >= 1.0);
    let mut side = start_positive;
    let mut traveled = 0.0;
    let mut i = 0u32;
    loop {
        let dist = r.powf(xi + f64::from(i));
        if side == target_positive && dist >= t {
            return traveled + t;
        }
        traveled += 2.0 * dist;
        side = !side;
        i += 1;
        if i > 10_000 {
            panic!("line search did not reach the target");
        }
    }
}

pub fn competitive_ratio(
    r: f64,
    xi: f64,
    start_positive: bool,
    target_positive: bool,
    t: f64,
) -> f64 {
    cost(r, xi, start_positive, target_positive, t) / t
}
