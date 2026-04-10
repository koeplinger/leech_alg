# Trial 005 Results — Triple symmetric-composition algebra (para-octonion + Okubo_τ + Okubo_τ²)

Date: 2026-04-09

## Algebra tested

Three different 8-dimensional symmetric composition algebras on R²⁴, one per
block, forming a Z₃ orbit under the powers of an order-3 octonion automorphism τ:

- Block 0 (indices 0–7): **para-octonion** x *₀ y = x̄ · ȳ (Petersson with τ⁰)
- Block 1 (indices 8–15): **Okubo from τ** x *₁ y = τ(x̄) · τ²(ȳ)
- Block 2 (indices 16–23): **Okubo from τ²** x *₂ y = τ²(x̄) · τ(ȳ)

τ fixes the quaternion subalgebra {e₀, e₁, e₃, e₇} and rotates (e₂, e₅)
and (e₄, e₆) by 2π/3.  Reference: [MarraniCorradettiZucconi2025] eq. (1.5).

Cross-block products use the target block's algebra: Bα × Bβ → Bγ using *_γ.

### Why not three shifted-mediator Okubo algebras?

The user's original intent was to cyclically rotate which element of a Fano
line serves as mediator, giving three different Okubo algebras from the same
line.  Computational verification showed that for EVERY Fano line, only 1 or
2 of the 3 cyclic mediators yield valid automorphisms.  No line admits all
three.  The Z₃ orbit {τ⁰, τ¹, τ²} is the closest valid construction.

## Results

### Base algebra: catastrophic failure

| Factor types | Tested | In Λ | Closure |
|---|---|---|---|
| t1×t1 | 200 | 107 | 53.5% |
| t1×t2 | 200 | 0 | 0.0% |
| t2×t2 | 500 | 0 | 0.0% |
| t1×t3 | 100 | 4 | 4.0% |
| t2×t3 | 100 | 0 | 0.0% |
| t3×t3 | 500 | 0 | 0.0% |
| **Total** | **1600** | **111** | **6.9%** |

**Contrast with triple-octonion (trial 001):** the triple-octonion had ~68.8%
closure with failures confined to type3×type3.  This algebra fails on ALL type
combinations except a fraction of type1×type1.

### Failure mode: products leave the E8 lattice

| Wilson condition | Violations | Of 1,489 failures |
|---|---|---|
| x+y+z ∈ Ls | 1,467 | 98.5% |
| y+z ∈ Ls̄ | 1,457 | 97.9% |
| z ∈ L | 1,362 | 91.5% |
| x+z ∈ Ls̄ | 1,362 | 91.5% |
| y ∈ L | 1,362 | 91.5% |
| x+y ∈ Ls̄ | 1,362 | 91.5% |
| x ∈ L | 104 | 7.0% |

**91.5% of failures violate condition 1 (y ∈ L, z ∈ L).** The products don't
even stay in the E8 lattice.  This is a much more severe failure than the
triple-octonion, which always satisfied conditions 1 and 2.

### Root cause: irrational structure constants

The Petersson construction with τ introduces structure constants involving
cos(2π/3) = −1/2 and sin(2π/3) = √3/2.  When applied to E8 lattice vectors
(which have integer or half-integer coordinates), the products acquire
irrational √3 components that cannot lie in the E8 lattice.

Block 0 (para-octonion) is less affected because x̄ · ȳ has the same
structure constants as the octonion product (just with sign changes from
conjugation).  This explains why type1×type1 partially succeeds (53.5%):
type-1 vectors have only one nonzero block, so the same-block product in
block 0 uses the para-octonion (integer structure constants).

### Discrete variants: all 1,536 fail

| Metric | Value |
|---|---|
| Variants tested | 1,536 (exhaustive) |
| Best closure rate | 9.4% |
| Best variant | route=left, conj=[ā*b̄,a*b̄,a*b̄], sign=−−− |

No discrete modification (conjugation, sign, routing) can fix the
irrational structure constant problem.

## Verdict

**FAIL — fundamental obstruction.** The Okubo/para-octonion structure
constants involve √3, which makes products of E8 lattice vectors land
outside the E8 lattice.  This is more severe than the triple-octonion
failure (which preserved E8 membership but failed on condition 3).  No
conjugation, sign, or routing variant can fix irrational structure constants.

## What this rules out

Any triple algebra built from Petersson isotopes of the octonion algebra
(para-octonion and/or Okubo algebras) using the standard 2π/3-rotation
automorphism τ, with the same block structure and cross-block routing as
the triple-octonion trials.
