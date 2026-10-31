# Trial 001 Results — Triple-octonion algebra (O₁ ⊕ O₂ ⊕ O₃)

Created: 2026-04-09.  Updated: 2026-10-28.

Date: 2026-04-09

## Algebra tested

Three copies of the standard octonion algebra on R²⁴:
- O₁ at indices 0–7, O₂ at indices 8–15, O₃ at indices 16–23.
- Same-block products: Oα × Oα → Oα (standard octonion multiplication).
- Cross-block products: Oα × Oβ → Oγ where {α,β,γ} = {1,2,3}.
- The octonion structure constants are identical for all block pairs (no sign
  changes or conjugations on cross terms).

## Closure test results

Tested against the 196,560 minimal-shell vectors of Λ (squared norm 8).

| Factor types | Pairs tested | Failures | Failure rate |
|-------------|-------------|----------|-------------|
| type1 × type1 | 518,400 (exhaustive) | 0 | 0% |
| type1 × type2 | 5,000 (random) | 0 | 0% |
| type1 × type3 | 5,000 (random) | 0 | 0% |
| type2 × type2 | 5,000 (random) | 0 | 0% |
| type2 × type3 | 5,000 (random) | 0 | 0% |
| type3 × type3 | 5,000 (random) | 3,740 | 74.8% |

No program in the repository reproduces this run or the norm² table below.
`python_project/src/trial_001_triple_octonion.py` tests 1,300 pairs (200
type1 × type1, 200 type1 × type2, 500 type2 × type2, 100 type1 × type3,
100 type2 × type3 and 200 type3 × type3, the type-3 vectors from a fixed
sample of 2,000) and finds 160 failures, all type3 × type3 (80.0% of 200)
and all at condition 3, with norm² from 32 to 88.  On 5,000 type3 × type3
pairs, `python_project/src/consistency_checks.py`, check 8, finds 3,669
failures (73.4%), all at condition 3.

## Key finding: failures are exclusively type3 × type3

All five other factor-type combinations show **zero failures** (including the
exhaustive 518,400-pair test of type1 × type1). Only type3 × type3 products
fail, and they fail at a very high rate (≈75%).

## Failure mode analysis

**Every single failure violates exactly one Wilson condition: `x+y+z ∈ Ls`
(condition 3).** The other conditions, components in L (condition 1) and
pairwise sums in Ls̄ (condition 2), are satisfied in all failures.

This is a remarkably clean failure signature:
- Wilson condition 1 (x, y, z ∈ L): always satisfied ✓
- Wilson condition 2 (x+y, x+z, y+z ∈ Ls̄): always satisfied ✓
- Wilson condition 3 (x+y+z ∈ Ls): violated in 100% of failures ✗

## Norm² distribution of failing products

| ||prod||² | Count |
|-----------|-------|
| 24 | 6 |
| 32 | 73 |
| 40 | 129 |
| 48 | 582 |
| 56 | 503 |
| 64 | 1,135 |
| 72 | 486 |
| 80 | 591 |
| 88 | 127 |
| 96 | 83 |
| 104 | 8 |
| 112 | 14 |
| 120 | 3 |

Products are NOT confined to the minimal shell (norm² = 8). They range from
norm² = 24 to 120, with a peak at 64.

## Sanity check: product norms

- Type1 × type1: all products have norm² = 64 (= 8 × 8, consistent with
  composition algebra norm-multiplicativity acting on norm-8 vectors).
- Type2 × type2: norms² ∈ {32, 48, 64, 80, 96}.
- Products are much larger than the minimal shell norm. No rescaling was applied
  for this trial because the test is whether products stay in Λ (not just on the
  minimal shell).

## Interpretation

The triple-octonion construction preserves the E8 lattice structure (condition 1)
and the Ls̄ sublattice structure (condition 2) perfectly. It fails only on the
finest condition, membership in Ls, and only when both factors are type-3
vectors.

Type-3 vectors are the family ((λs)j, λk, (λj)k). The failure has an exact
reason (paper Proposition 5.3 with the identity twist; key claim 007): the
block sum of a product is the product of the block sums, so condition 3 holds
for all products of vectors of Λ exactly when Ls is closed under the standard
octonion product, and it is not (paper Remark 4.5).

## Verdict

**FAIL.** The triple-octonion algebra is not closed on Λ. The failure is
confined to type3 × type3 products and is caused exclusively by Wilson
condition 3 (x+y+z ∈ Ls). Conjugation, sign and routing changes in the
cross-block terms do not remove it (trial 003). The triple product with the
σ-twisted product in all nine block products closes on Λ (trial 007; paper
Theorem 1.1): σ(Ls̄) and σ(Ls) are closed under the standard octonion product
(paper Lemma 4.4 and Proposition 5.3).
