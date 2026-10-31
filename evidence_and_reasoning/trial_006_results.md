# Trial 006 Results — E8 automorphism basis changes for triple Okubo

Created: 2026-04-09.  Updated: 2026-10-28.

Date: 2026-04-09

## Algebra tested

Same as trial 005 (para-octonion + Okubo_τ + Okubo_τ²), but with E8 lattice
automorphisms applied as basis changes within each block.

## Results

| Search | Tested | Best closure |
|---|---|---|
| Baseline (identity) | 1 | 0.0% (0/80) |
| Single-block random | 300 | 2.5% |
| Three-block random | 300 | 3.8% |

## Interpretation

The baseline drops to 0.0% (from 6.9% in trial 005's larger sample) because
this trial's test set includes neither type1×type1 nor type1×type3 pairs, the
only combinations with products in Λ in trial 005 (53.5% and 4.0%).  Even with automorphisms, the best rate is 3.8%,
essentially zero.

The irrational √3 factors in the Petersson construction's structure constants
cannot be absorbed by E8 automorphisms, which are rational matrices (the ones
tested are signed permutations).

## Verdict

**FAIL.** Confirms trial 005's conclusion: the fundamental problem is
irrational structure constants, not the choice of basis.
