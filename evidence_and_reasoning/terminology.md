# Terminology Reference

This file defines terminology used in this project. Where established
mathematical terminology exists, it is cited. Where constructions are
specific to this project, terms are marked as **project-specific**.

---

## Established terms

### Symmetric composition algebra

An 8-dimensional composition algebra (V, *, n) satisfying
(x * y) * x = x * (y * x) = n(x) y for all x, y.  Over a field of
characteristic ≠ 2, 3, these are exactly the **para-Hurwitz algebras** and the
**Okubo algebras**.  Reference: [Elduque2000_Triality], [KMRT, Chapter VIII].

### Para-Hurwitz algebra (para-octonion algebra)

The symmetric composition algebra obtained from an octonion algebra (O, ·)
by the **Petersson construction** with the identity automorphism (τ = id):
x * y = x̄ · ȳ, where x̄ denotes octonion conjugation.  This is the simplest
symmetric composition algebra.  Reference: [Elduque2000_Triality] §3.

### Okubo algebra (pseudo-octonion algebra)

The symmetric composition algebra obtained from an octonion algebra (O, ·)
via the **Petersson construction** with a non-trivial order-3 automorphism τ:
x * y = τ(x̄) · τ²(ȳ).  Named after Susumu Okubo, who first constructed
it on traceless 3×3 matrices over C.  Over R, there is (up to isomorphism)
one Okubo algebra; over other fields there may be several.
References: [MarraniCorradettiZucconi2025], [Elduque2000_Triality].

### Petersson construction

Given an octonion algebra (O, ·) and an order-3 automorphism τ of O, the
**Petersson isotope** is the algebra (O, *_τ) with product
x *_τ y = τ(x̄) · τ²(ȳ).  The resulting algebra is always a symmetric
composition algebra.  When τ = id, this gives the para-Hurwitz algebra;
when τ ≠ id, it gives an Okubo algebra.
Reference: [MarraniCorradettiZucconi2025] eq. (1.5); originally due to
H. Petersson (1969).

### Maximal order (of integral octonions)

A lattice Γ ⊂ O (the real octonion algebra) that is closed under
multiplication, contains the identity, and is maximal with respect to
inclusion among such lattices.  There are exactly 7 maximal orders in the
octonions, each isometric to the E8 root lattice.  They are obtained by
choosing one of the 7 imaginary basis elements as a distinguished unit and
applying the **index-doubling permutation** (correcting Kirmse's original
non-closed construction).  Reference: Conway and Smith, "On Quaternions and
Octonions" (2003), Chapter 9; Coxeter (1946).

### Kirmse integers / index-doubling permutation

Kirmse (1924) proposed a set of 240 "integral octonions" forming the E8
lattice, but this set is NOT closed under multiplication.  Coxeter (1946)
showed that applying an **index-doubling permutation** (a specific
permutation of the 7 imaginary basis elements) to Kirmse's set produces a
genuine maximal order.  In this project, we use **Kirmse twist** as shorthand
for this correction.  Reference: Conway and Smith (2003), §9.3.

### Triality (D₄ triality)

The outer automorphism of the Lie group Spin(8) (or equivalently, the
S₃-symmetry of the D₄ Dynkin diagram) that permutes the three 8-dimensional
representations (vector, left-spinor, right-spinor).  For composition
algebras, the **Principle of Local Triality** (Elduque) relates the three
symmetric composition algebras {para-O, Okubo_τ, Okubo_τ²} obtained from
a single octonion algebra O and an order-3 automorphism τ.
Reference: [Elduque2000_Triality] §4.

---

## Project-specific terms

### Triple octonion product

(Trials 001, 002.)  A 24-dimensional algebra R²⁴ = O₁ ⊕ O₂ ⊕ O₃ built
from three copies of the standard octonion algebra with **Z₃ cross-block
routing**: products of vectors from blocks α and β land in block γ, where
{α, β, γ} = {1, 2, 3}, using the standard octonion multiplication.
Same-block products Oα × Oα → Oα also use the standard octonion product.

Trial 001 tested this product without rescaling.  Trial 002 tested it with
per-block scaling.  Both failed, but trial 001 exhibited a distinctive
**2+1 closure pattern** (see below).

### Petersson triality triple

(Trial 005, prompt 023.)  A 24-dimensional algebra R²⁴ = B₀ ⊕ B₁ ⊕ B₂
where the three blocks carry different symmetric composition algebras
forming the **Z₃ orbit** under powers of an order-3 octonion automorphism τ:

- B₀: para-octonion (Petersson with τ⁰ = id)
- B₁: Okubo_τ   (Petersson with τ)
- B₂: Okubo_τ²  (Petersson with τ²)

Cross-block products use the target block's algebra.  This construction
arose when the user's original request (cyclically rotating the mediator
of a Fano line) was shown to be mathematically impossible for all 7 Fano
lines; the Petersson triality triple is the closest valid alternative.

The construction failed due to a **√3 obstruction**: the Petersson
construction with non-trivial τ introduces irrational structure constants
(cos 2π/3, sin 2π/3), so products of E8 lattice vectors leave the E8
lattice entirely.

### 2+1 closure pattern

A recurring phenomenon in this project's trials where closure tests exhibit
an asymmetry between three nominally symmetric components, with two
succeeding (or partially succeeding) and one failing:

- **Trial 001** (triple octonion product): type1×type1, type1×type2,
  type1×type3, type2×type2, type2×type3 all closed at 100%.  Only
  type3×type3 failed (≈75% failure rate).  Since type-1 and type-2 vectors
  each populate only two of the three blocks while type-3 vectors populate
  all three, this is effectively a 2+1 split based on block occupancy.

- **Trial 005** (Petersson triality triple): block 0 (para-octonion)
  partially succeeded (53.5% on type1×type1) while blocks 1 and 2
  (Okubo algebras) failed catastrophically.  Here the 2+1 split is between
  the two Okubo blocks (which share irrational structure constants) and the
  para-octonion block (which has integer structure constants).

The 2+1 pattern may be a structural signal reflecting the underlying
asymmetry between the para-Hurwitz and Okubo components of the triality
orbit, or between the block-occupancy patterns of the three vector types
in Wilson's Leech lattice construction.

### √3 obstruction

The failure mode specific to the Petersson triality triple (trial 005):
the Petersson construction with a non-trivial order-3 automorphism τ
introduces structure constants involving cos(2π/3) = −1/2 and
sin(2π/3) = √3/2.  Products of E8 lattice vectors (integer or
half-integer coordinates) then acquire irrational √3 components and
cannot lie in the E8 lattice.  This is a more fundamental obstruction
than the Wilson-condition-3 failure of the triple octonion product.

---

Last updated: 2026-04-11
