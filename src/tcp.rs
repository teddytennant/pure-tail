//! TCP acknowledgment, Sections 5 and 6.
//!
//! `DET_z` acknowledges at the first time T for which some earlier time τ in
//! the current interval satisfies P(last, τ) * (T - τ) = z. On a finite arrival
//! sequence that minimum is the min over pending packets of
//! arrival_j + z / j, with j counted from 1 in arrival order.
//!
//! `RAND_r` draws z from the distribution in Section 5.2 and runs `DET_z`.

/// γ(r) = 1 + 1/r, for r in (0, 1].
pub fn gamma(r: f64) -> f64 {
    1.0 + 1.0 / r
}

fn denom(r: f64) -> f64 {
    (1.0 + r) * std::f64::consts::E - r.exp()
}

/// Point mass q_r at z = r.
pub fn mass_at_r(r: f64) -> f64 {
    r * r.exp() / denom(r)
}

/// Density p(z) on (r, 1].
pub fn density(r: f64, z: f64) -> f64 {
    (r + 1.0) * z.exp() / denom(r)
}

/// Lemma 12: expected competitive ratio 1 + e^r / ((1 + r) e - e^r).
pub fn expected_ratio_bound(r: f64) -> f64 {
    1.0 + r.exp() / denom(r)
}

/// Theorem 4 lower bound. r is in [(√5 - 1) / 2, 1].
pub fn lower_bound(r: f64) -> f64 {
    let pow = r.powf(2.0 + r);
    1.0 + (2.0 + r) * pow / ((1.0 + r) * (1.0 + r) - pow)
}

/// Inverse CDF of the Section 5.2 distribution. `u` is in [0, 1).
pub fn sample_z(r: f64, u: f64) -> f64 {
    let q = mass_at_r(r);
    if u <= q || (1.0 - q) <= f64::EPSILON {
        return r;
    }
    let v = (u - q) / (1.0 - q);
    let er = r.exp();
    (er + v * (std::f64::consts::E - er)).ln()
}

fn ack_candidate(pending: &[f64], z: f64) -> f64 {
    let mut best = f64::INFINITY;
    for (j, &arrival) in pending.iter().enumerate() {
        let cand = arrival + z / (j as f64 + 1.0);
        if cand < best {
            best = cand;
        }
    }
    best
}

/// Cost of `DET_z`: number of acknowledgments plus total latency.
///
/// `arrivals` must be sorted ascending. Packets that arrive at the
/// acknowledgment time are included in that acknowledgment.
pub fn det_cost(arrivals: &[f64], z: f64) -> f64 {
    assert!(z > 0.0);
    let n = arrivals.len();
    if n == 0 {
        return 0.0;
    }
    let mut start = 0usize;
    let mut acks = 0.0;
    let mut latency = 0.0;
    while start < n {
        let mut end = start + 1;
        loop {
            let cand = ack_candidate(&arrivals[start..end], z);
            // Include a packet that arrives as the acknowledgment fires.
            if end < n && arrivals[end] <= cand {
                end += 1;
                continue;
            }
            acks += 1.0;
            for &arrival in &arrivals[start..end] {
                latency += cand - arrival;
            }
            start = end;
            break;
        }
    }
    acks + latency
}

/// Offline optimum. Acknowledgments only at arrival times; the last packet
/// is always acknowledged. Waiting between arrivals only adds latency.
pub fn opt_cost(arrivals: &[f64]) -> f64 {
    let n = arrivals.len();
    if n == 0 {
        return 0.0;
    }
    let mut prefix = vec![0.0; n + 1];
    for i in 0..n {
        prefix[i + 1] = prefix[i] + arrivals[i];
    }
    let mut best = vec![f64::INFINITY; n + 1];
    best[0] = 0.0;
    for i in 1..=n {
        let t = arrivals[i - 1];
        for j in 0..i {
            let count = (i - j) as f64;
            let lat = count * t - (prefix[i] - prefix[j]);
            let cost = best[j] + 1.0 + lat;
            if cost < best[i] {
                best[i] = cost;
            }
        }
    }
    best[n]
}

/// Expected `DET_z` cost under the `RAND_r` distribution, by a trapezoid
/// over the continuous part plus the point mass at r.
pub fn expected_det_cost(arrivals: &[f64], r: f64, steps: usize) -> f64 {
    let q = mass_at_r(r);
    let mut total = q * det_cost(arrivals, r);
    if r >= 1.0 || steps == 0 {
        return total;
    }
    let width = (1.0 - r) / steps as f64;
    let mut prev_z = r;
    let mut prev_w = density(r, r) * det_cost(arrivals, r);
    for s in 1..=steps {
        let z = r + width * s as f64;
        let w = density(r, z) * det_cost(arrivals, z);
        total += 0.5 * (prev_w + w) * (z - prev_z);
        prev_z = z;
        prev_w = w;
    }
    total
}
