# Trial 006 Results — E8 automorphism basis changes for triple Okubo

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
this trial's test set does not include type1×type1 pairs (the only partially
successful combination).  Even with automorphisms, the best rate is 3.8% —
essentially zero.

The irrational √3 factors in the Petersson construction's structure constants
cannot be absorbed by E8 automorphisms (which are signed permutations with
integer/half-integer entries).

## Verdict

**FAIL.** Confirms trial 005's conclusion: the fundamental problem is
irrational structure constants, not the choice of basis.
