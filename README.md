# Pure tail constraints

Dependable randomized algorithms for online bidding, line search, and TCP acknowledgment (arXiv:2609.35337).

```
cargo test
cargo run --bin demo
```

`BID_r` is the bidding algorithm from Section 3. `SMARTCOW_r` is the line-search algorithm from Section E (Kao, Reif, Tate 1996). `RAND_r` draws the Section 5.2 distribution and runs `DET_z`. The TCP optimum is the dynamic program that acknowledges only at arrival times.

This is not a TCP stack and not a proof of the lower bounds. Arrivals are a finite list, not a network. The lower bounds are the closed forms in Theorems 1, 2, and 4, not proof checkers. ζ in the bidding quality bound follows the proof of Lemma 7, log_r(last bid) - log_r(T), which is what makes equation (9) true. The main text calls this the fractional part of log_r of the last bid.
