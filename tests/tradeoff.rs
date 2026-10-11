use pure_tail::{bidding, euler_ratio, line_search, tcp};

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-9 * (1.0 + b.abs()), "{a} != {b}");
}

/// Midpoint integral of the ratio, split at the jump ξ = {log_r T}.
/// A uniform left grid overshoots r/ln r by about 1/N even when the
/// integral itself is strictly below the bound.
fn mean_bid_ratio(r: f64, t: f64, n: usize) -> f64 {
    let log_t = t.log(r);
    let frac = log_t - log_t.floor();
    let mut acc = 0.0;
    if frac > 0.0 {
        for i in 0..n {
            let xi = frac * (i as f64 + 0.5) / n as f64;
            acc += bidding::competitive_ratio(r, xi, t) * frac;
        }
    }
    let right = 1.0 - frac;
    if right > 0.0 {
        for i in 0..n {
            let xi = frac + right * (i as f64 + 0.5) / n as f64;
            if xi >= 1.0 {
                continue;
            }
            acc += bidding::competitive_ratio(r, xi, t) * right;
        }
    }
    acc / n as f64
}

#[test]
fn bidding_matches_lemma_7() {
    for r in [2.0, 2.4, std::f64::consts::E] {
        let gamma = bidding::gamma(r);
        let bound = bidding::expected_ratio_bound(r);
        let mut worst = 0.0_f64;
        for s in 0..48 {
            let xi = s as f64 / 48.0;
            let mut thresholds = vec![1.0, 1.0 + 1e-9, 1e6];
            let mut prev = 1.0_f64;
            for k in 0..8 {
                let x = r.powf(xi + k as f64);
                thresholds.push(x);
                thresholds.push((prev + x) * 0.5);
                // Just above the previous bid is the dependability worst case.
                if k > 0 {
                    thresholds.push(prev * (1.0 + 1e-12));
                }
                prev = x;
            }
            for t in thresholds {
                if t < 1.0 {
                    continue;
                }
                let ratio = bidding::competitive_ratio(r, xi, t);
                let (cost, last, _) = bidding::bid(r, xi, t);
                assert!((cost / t - ratio).abs() <= 1e-12);
                assert!(
                    ratio <= gamma * (1.0 + 1e-9),
                    "r={r} xi={xi} T={t} ratio={ratio} gamma={gamma}"
                );
                let z = bidding::zeta(r, last, t);
                assert!(z >= -1e-9 && z < 1.0 + 1e-9, "zeta {z}");
                assert!(
                    ratio <= bidding::quality_bound(r, z.max(0.0)) * (1.0 + 1e-9),
                    "r={r} xi={xi} T={t} ratio={ratio}"
                );
                worst = worst.max(ratio);
            }
        }
        assert!(worst > gamma * 0.5);

        for t in [1.0, 1.0 + 1e-6, r.powf(3.2), 1e6] {
            let mean = mean_bid_ratio(r, t, 2000);
            assert!(
                mean <= bound + 1e-6,
                "mean {mean} > {bound} at r={r} T={t}"
            );
        }
        let large = mean_bid_ratio(r, 1e6, 2000);
        assert!(
            (large - bound).abs() < 1e-3,
            "large-T mean {large} bound {bound}"
        );
    }
}

#[test]
fn geometric_bids_obey_the_step_constraint() {
    let r: f64 = 2.5;
    let xi = 0.3;
    let mut sum = 0.0;
    let mut prev = 1.0; // x_0 in the proof
    for k in 0..12 {
        let x = r.powf(xi + k as f64);
        sum += x;
        assert!(sum <= bidding::gamma(r) * prev + 1e-8);
        prev = x;
    }
}

#[test]
fn smartcow_is_dependable_and_near_its_expectation() {
    let r0 = line_search::r0();
    assert!((r0 * r0.ln() - r0 - 1.0).abs() < 1e-9);
    assert!((3.59112..=3.59113).contains(&r0));
    assert!((line_search::gamma(2.0) - 9.0).abs() < 1e-12);

    for r in [2.0, 2.8, r0] {
        let gamma = line_search::gamma(r);
        let bound = line_search::expected_ratio_bound(r);
        for s in 0..24 {
            let xi = s as f64 / 24.0;
            let mut thresholds = vec![1.0, 1.0 + 1e-9];
            for k in 0..6 {
                let dist = r.powf(xi + k as f64);
                if dist >= 1.0 {
                    thresholds.push(dist);
                    thresholds.push(dist * (1.0 + 1e-12));
                }
            }
            thresholds.push(r.powf(18.0));
            for t in thresholds {
                for target_positive in [true, false] {
                    for start_positive in [true, false] {
                        let ratio = line_search::competitive_ratio(
                            r,
                            xi,
                            start_positive,
                            target_positive,
                            t,
                        );
                        assert!(
                            ratio <= gamma * (1.0 + 1e-9),
                            "r={r} xi={xi} T={t} ratio={ratio} gamma={gamma}"
                        );
                    }
                }
            }
        }

        // Large T: midpoint grid over ξ and both starts. The mean approaches
        // 1 + (1+r)/ln r from below; a short T never gets within 15%.
        let t = r.powf(24.0);
        let n = 4000usize;
        let mut sum = 0.0;
        for i in 0..n {
            let xi = (i as f64 + 0.5) / n as f64;
            for start_positive in [true, false] {
                sum += line_search::competitive_ratio(r, xi, start_positive, true, t);
            }
        }
        let mean = sum / (2.0 * n as f64);
        assert!(
            mean <= bound + 1e-6,
            "mean {mean} exceeds {bound} at r={r}"
        );
        assert!(
            (mean - bound).abs() < 1e-3,
            "large-T mean {mean} bound {bound}"
        );
    }
}

#[test]
fn tcp_distribution_and_single_packet_meet_lemma_12() {
    for r in [0.2, 0.5, 0.9, 1.0] {
        let d = (1.0 + r) * std::f64::consts::E - r.exp();
        let integral = (r + 1.0) / d * (std::f64::consts::E - r.exp());
        close(tcp::mass_at_r(r) + integral, 1.0);
        close(tcp::det_cost(&[0.0], 0.4_f64.max(r)), 1.0 + 0.4_f64.max(r));
        // E[z] = e^r / ((1+r) e - e^r), so a packet at time 0 meets the bound.
        let e_z = tcp::mass_at_r(r) * r
            + (r + 1.0) / d * (1.0 - r) * r.exp();
        close(1.0 + e_z, tcp::expected_ratio_bound(r));
    }
    close(tcp::expected_ratio_bound(1.0), 2.0);
    close(tcp::lower_bound(1.0), 2.0);
    let tiny = 1e-6;
    assert!(tcp::expected_ratio_bound(tiny) < euler_ratio() + 1e-4);
    assert!(tcp::expected_ratio_bound(tiny) > euler_ratio() - 1e-3);
}

#[test]
fn det_z_is_dependable_on_finite_instances() {
    let instances: &[&[f64]] = &[
        &[0.0],
        &[0.0, 0.1],
        &[0.0, 1.0, 1.0, 3.5],
        &[0.2, 0.2, 0.2, 5.0, 5.1, 9.0],
        &[1.0, 1.01, 1.02, 1.5, 8.0],
    ];
    for arrivals in instances {
        let opt = tcp::opt_cost(arrivals);
        assert!(opt >= 1.0);
        for z in [0.05, 0.2, 0.5, 1.0] {
            let cost = tcp::det_cost(arrivals, z);
            let ratio = cost / opt;
            assert!(
                ratio <= tcp::gamma(z) + 1e-8,
                "z={z} ratio={ratio} gamma={} arrivals={arrivals:?}",
                tcp::gamma(z)
            );
        }
        for r in [0.25, 0.6, 1.0] {
            let expected = tcp::expected_det_cost(arrivals, r, 80);
            let ratio = expected / opt;
            assert!(
                ratio <= tcp::expected_ratio_bound(r) + 1e-3,
                "r={r} ratio={ratio} bound={}",
                tcp::expected_ratio_bound(r)
            );
        }
    }
}

#[test]
fn sample_z_respects_the_mass() {
    let r = 0.3;
    let q = tcp::mass_at_r(r);
    assert_eq!(tcp::sample_z(r, 0.0), r);
    assert_eq!(tcp::sample_z(r, q), r);
    let z = tcp::sample_z(r, 0.5 * (q + 1.0));
    assert!(z > r && z <= 1.0 + 1e-12);
    close(tcp::sample_z(r, 1.0 - 1e-12), 1.0);
}

#[test]
fn section6_identities_and_a_burst() {
    let phi = 0.5 * (5.0_f64.sqrt() - 1.0);
    for r in [phi, 0.75, 1.0] {
        let g = tcp::gamma(r);
        close(g, 1.0 + 1.0 / r);
        let fixed = 1.0 + r;
        close(g * (fixed - 1.0), fixed);
        close(fixed, g / (g - 1.0));
        let c = r * r + r - 1.0;
        assert!((0.0..=1.0).contains(&c) || c.abs() < 1e-12);
    }
    for n in [1usize, 4, 17] {
        let arrivals = vec![0.0; n];
        for z in [0.2, 0.5, 1.0] {
            close(tcp::det_cost(&arrivals, z), 1.0 + z);
            close(tcp::opt_cost(&arrivals), 1.0);
        }
    }
}

#[test]
fn lower_bound_stays_under_the_upper_bound() {
    let phi = 0.5 * (5.0_f64.sqrt() - 1.0);
    let mut r = phi;
    while r <= 1.0 {
        let lower = tcp::lower_bound(r);
        let upper = tcp::expected_ratio_bound(r);
        assert!(lower <= upper + 1e-9, "r={r} lower={lower} upper={upper}");
        assert!(lower >= 1.0);
        r += 0.05;
    }
}
