# Trial 002 Results — Per-block scaled triple-octonion algebra

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

**No per-block scaling produces any closure on Min(Λ).** The single hit at
u = v = 1 was a coincidental type2 × type2 product.

### Root cause analysis: two independent obstructions

**Obstruction 1 — Norm bound (type3 × type3)**

For type3 × type3 products, the minimum achievable ||q(u,v)||² over ALL (u, v)
was computed via numerical optimization for 100 random pairs:

    Minimum ||q||²: 14.6
    First 20 sorted values: [14.6, 14.9, 15.3, 17.9, 19.8, 20.0, ...]

The minimum is nearly twice the target (8). **No per-block scaling can bring
type3 × type3 products to the minimal shell.** This is a hard geometric
obstruction: the same-block and cross-block contributions cannot cancel
sufficiently because all three blocks carry nonzero components in type-3
vectors.

**Obstruction 2 — Lattice membership (all types)**

For type2 × type2 products, the norm CAN reach 8 at specific (u, v) values.
A targeted search found 1,188 products with ||q||² = 8 (to within 0.001).
None of them satisfy Wilson's lattice conditions:

    Products at ||q||² = 8 checked: 1,188
    Of those in Λ: 0

The normalization factor √(8/||q||²) is generically irrational (since ||q||²
is a sum of products of E8 root coordinates and half-integers). Multiplying a
lattice vector by an irrational factor produces a non-lattice vector. Even when
the norm happens to be exactly 8, the individual coordinates are irrational
combinations of the building blocks and do not satisfy the D8+ integrality
conditions or Wilson's sublattice conditions.

## Interpretation

The triple-octonion algebra with per-block scaling faces **two independent
and fundamental obstructions** to Min(Λ) closure:

1. **Norm obstruction**: For the most complex vector pairs (type3 × type3),
   the product norm is bounded below by ≈ 14.6, preventing it from ever
   reaching the minimal shell at norm² = 8.

2. **Lattice obstruction**: Even for simpler pairs where norm 8 is achievable,
   the per-block scaling introduces irrational mixing of lattice coordinates,
   destroying the D8+ integrality structure.

These obstructions are specific to per-block (diagonal) scaling.  A more general
basis change (non-diagonal, mixing coordinates within blocks) might overcome
obstruction 2. A fundamentally different multiplication rule (not the plain
octonion product on each block pair) would be needed to overcome obstruction 1.

## Verdict

**FAIL — structural obstruction.** Per-block scaling of the triple-octonion
algebra cannot achieve Min(Λ) closure. The obstruction is not a search failure
(insufficiently fine grid) but a mathematical impossibility for the type3 × type3
sector. The minimum achievable product norm (≈ 14.6) is bounded well above the
target (8) for all (u, v).

## What this rules out

Any algebra of the form "three copies of octonions with per-block scaling
and Z₃-symmetric cross-block mixing (product always goes to the third block)"
is ruled out as a candidate for an order on Λ closed on Min(Λ).

## What remains open

1. **Different cross-block structure constants** — using conjugation, sign flips,
   or other modifications in the cross-block terms (not just the plain octonion
   product).
2. **Non-diagonal basis changes** — a general GL(8) transformation within each
   block, not just scaling.
3. **Different block decomposition** — not O₁ ⊕ O₂ ⊕ O₃ but some other
   partition of 24 dimensions.
