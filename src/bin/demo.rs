//! Print a few points on the three tradeoff curves.

use pure_tail::{bidding, euler_ratio, line_search, tcp};

fn main() {
    let r0 = line_search::r0();
    println!("online bidding  (gamma = r^2/(r-1),  rho = r/ln r)");
    for r in [2.0, std::f64::consts::E] {
        println!(
            "  r = {r:.4}  gamma = {:.4}  expected ratio <= {:.4}",
            bidding::gamma(r),
            bidding::expected_ratio_bound(r)
        );
    }

    println!("line search  (r0 = {r0:.5})");
    for r in [2.0, r0] {
        println!(
            "  r = {r:.4}  gamma = {:.4}  expected ratio <= {:.4}",
            line_search::gamma(r),
            line_search::expected_ratio_bound(r)
        );
    }

    println!(
        "TCP acknowledgment  (unconstrained randomized ratio e/(e-1) = {:.4})",
        euler_ratio()
    );
    let phi = 0.5 * (5.0_f64.sqrt() - 1.0);
    for r in [phi, 0.8, 1.0] {
        println!(
            "  r = {r:.4}  gamma = {:.4}  upper = {:.4}  lower = {:.4}",
            tcp::gamma(r),
            tcp::expected_ratio_bound(r),
            tcp::lower_bound(r)
        );
    }
}
