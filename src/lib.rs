//! Dependability-competitiveness tradeoffs from arXiv:2609.35337.
//!
//! The simulators are the upper-bound algorithms: `BID_r` for online bidding,
//! `SMARTCOW_r` for line search, and `RAND_r` / `DET_z` for TCP acknowledgment.
//! Lower bounds are closed forms. They are not proof checkers.

pub mod bidding;
pub mod line_search;
pub mod tcp;

/// e / (e - 1), the unconstrained randomized TCP / ski-rental ratio.
pub fn euler_ratio() -> f64 {
    std::f64::consts::E / (std::f64::consts::E - 1.0)
}
