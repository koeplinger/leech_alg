# Trial 002 Results — Per-block scaled triple-octonion algebra

Created: 2026-04-09.  Updated: 2026-10-28.

Date: 2026-04-09

## Algebra tested

Same triple-octonion multiplication as trial 001 (O₁ ⊕ O₂ ⊕ O₃ with Z₃
cross-block mixing), but with per-block scaling (s₁, s₂, s₃).  The product
of two R²⁴ vectors has block k equal to:

    Block k = (1/sₖ) · same_k + sₖ/(sᵢsⱼ) · cross_k

where same_k is the within-block octonion product and cross_k is the sum of
cross-block octonion products targeting block k.

Since overall scaling factors out (all sₖ × c gives product / c), the
search space reduces to two ratios: u = s₂/s₁, v = s₃/s₁, with the overall
scale determined by the norm constraint ||product||² = 8.

## Search methodology

For each pair (a, b) from Min(Λ) and each (u, v):
1. Compute the unnormalized product shape q(u, v).
2. Normalize: p = (√8/||q||) · q, so that ||p|| = √8.
3. Check whether p ∈ Λ (Wilson's conditions).

Three search strategies:
- **Z₃-symmetric** (u = v = 1): equivalent to normalizing the trial-001 product.
- **Grid search**: 60 × 60 log-spaced grid over (u, v) ∈ [0.1, 10]².
- **Random search**: 50,000 random (u, v) samples, log-uniform on [0.1, 10]².

## Results

### Search outcomes

| Search | (u, v) tested | Best Min(Λ) hits | Out of 300 pairs |
|--------|--------------|-------------------|-----------------|
| Z₃-symmetric | 1 | 1 (0.3%) | 300 |
| Grid (3,600 pts) | 3,600 | 0 | 300 |
| Random (50,000 pts) | 50,000 | 0 | 300 |

**No tested per-block scaling closes Min(Λ).** One type2 × type2 product
lands on Min(Λ), at u = v = 1; none does at the grid and random points.

### Root cause analysis

The two computations in this section, a numerical minimization of ||q(u,v)||²
and a targeted search for products of norm 8, are not part of
`python_project/src/trial_002_scaled_triple_octonion.py`, which runs the three
searches above; no program in the repository reproduces them.

**Norm bound at s₁ = 1 (type3 × type3)**

For type3 × type3 products, the minimum achievable ||q(u,v)||² over ALL (u, v)
was computed via numerical optimization for 100 random pairs:

    Minimum ||q||²: 14.6
    First 20 sorted values: [14.6, 14.9, 15.3, 17.9, 19.8, 20.0, ...]

This minimum is taken at the overall scale s₁ = 1.  The overall scale factors
out (all sₖ × c gives product / c), so a suitable s₁ gives any one product
norm 8; the minimum is no obstruction to reaching the minimal shell.

**Obstruction: lattice membership (all types)**

For type2 × type2 products, the norm CAN reach 8 at specific (u, v) values.
A targeted search found 1,188 products with ||q||² = 8 (to within 0.001).
None of them satisfy Wilson's lattice conditions:

    Products at ||q||² = 8 checked: 1,188
    Of those in Λ: 0

The normalization factor √(8/||q||²) is generically irrational (since ||q||²
is a sum of products of E8 root coordinates and half-integers). Multiplying a
lattice vector by an irrational factor produces a nonlattice vector. Even when
the norm happens to be exactly 8, the individual coordinates are irrational
combinations of the building blocks and do not satisfy the D8+ integrality
conditions or Wilson's sublattice conditions.

## Interpretation

Two observations bear on Min(Λ) closure under per-block scaling:

1. **Norm**: at the overall scale s₁ = 1, type3 × type3 products have
   ||q||² ≥ 14.6 on the pairs sampled; since the overall scale is free, this
   does not prevent a product from reaching norm² = 8.

2. **Lattice obstruction**: the per-block scaling introduces irrational mixing
   of lattice coordinates, destroying the D8+ integrality structure.

These observations are specific to per-block (diagonal) scaling.  A more general
basis change (nondiagonal, mixing coordinates within blocks) might overcome
the lattice obstruction.

## Verdict

**FAIL.** Under per-block scaling the normalized products miss Min(Λ): one
type2 × type2 product lands on it at u = v = 1, and none at the grid and random
points tested.  The obstruction is lattice membership: the normalization makes
the coordinates irrational in general.

## What this rules out

Any algebra of the form "three copies of octonions with per-block scaling
and Z₃-symmetric cross-block mixing (product always goes to the third block)"
is ruled out as a candidate for a product that maps Min(Λ) × Min(Λ) into
Min(Λ).  That is shell-closure, which an order does not need (paper Remark
1.3).  The family does contain products closed on Λ: at s₁ = s₂ = s₃ = ½ the
product is twice the trial-001 product, which maps Λ × Λ into 2L³ ⊆ Λ, since
L · L ⊆ L and 2L ⊆ Ls̄ ∩ Ls.

## Other modifications

1. **Different cross-block structure constants**: conjugation, sign flips
   and routing changes in the cross-block terms, tested in trial 003 (all
   1,536 variants fail).
2. **Nondiagonal basis changes**: a general GL(8) transformation within each
   block, not just scaling.  Trial 004 samples signed permutations, and none
   of those sampled closes.  The same transposition σ in all three blocks,
   which trial 004 did not sample, closes on Λ (trial 007; paper
   Theorem 1.1).
3. **Different block decomposition**: not O₁ ⊕ O₂ ⊕ O₃ but some other
   partition of 24 dimensions; not tested.
