# Terminology Reference

This file defines terminology used in this project. Where established
mathematical terminology exists, it is cited. Where constructions are
specific to this project, terms are marked as **(project-specific)**.

---

## Index

**Established terms:**
Isomorphism (of algebras) — Isotopy (of algebras) — Kirmse integers —
Kirmse twist — Leech lattice Λ — Maximal order (of integral octonions) —
Okubo algebra — Para-Hurwitz algebra — Petersson construction —
Symmetric composition algebra — Triality (D₄ triality)

**Project-specific terms:**
2+1 closure pattern — Petersson triality triple —
Transposition-twisted triple octonion product — Triple octonion product —
√3 obstruction

---

## Established terms

### Isomorphism (of algebras)

An **isomorphism** between two algebras (A, ·) and (B, *) is a bijective
linear map φ: A → B satisfying φ(x · y) = φ(x) * φ(y) for all x, y ∈ A.
A single map governs both input and output.

Isomorphism is a special case of **isotopy** (see below): it corresponds to
the isotopy triple (φ, φ, φ) where all three maps are the same.

In this project: the transposition-twisted algebra swap(s,t) is isomorphic
to the standard octonion algebra.  The map σ: e_s ↔ e_t (extended linearly)
satisfies σ(x ·_std y) = σ(x) ·_swap σ(y) for all basis pairs.  Verified
computationally: 0/64 mismatches.  However, σ does not preserve the E8
lattice L (it moves half-integer basis vectors), so L is closed under the
standard product but the swap product acts differently on L.

### Isotopy (of algebras)

An **isotopy** between two algebras (A, ·) and (B, *) is a triple (f, g, h)
of bijective linear maps A → B satisfying h(x · y) = f(x) * g(y) for all
x, y ∈ A.  Unlike isomorphism, the maps applied to the two inputs and the
output may all differ.  Two algebras related by an isotopy are called
**isotopic**; the new algebra is an **isotope** of the original.

When f = g = h = φ, this reduces to an **isomorphism** (see above).

In this project: the **Petersson construction** x *_τ y = τ(x̄) · τ²(ȳ)
is an isotopy of the octonion algebra (O, ·), with f: x ↦ τ(x̄),
g: y ↦ τ²(ȳ), and h = id.  Since f ≠ g (when τ ≠ id), it is a genuine
isotopy that is not an isomorphism.  The resulting algebra (a symmetric
composition algebra) has fundamentally different properties — for example,
it has no identity element.

Reference: McCrimmon, K., "A Taste of Jordan Algebras" (2004), §5.3;
Albert, A.A., "Non-associative algebras I", Ann. of Math. 43 (1942),
685–707 (introduced the concept of isotopy for non-associative algebras).

### Leech lattice Λ — norm conventions

The **Leech lattice** Λ is the unique even unimodular lattice of rank 24
with no vectors of squared norm 2.  Its automorphism group is 2·Co₁, where
Co₁ is Conway's first sporadic group.

**Norm convention used in this project.**  Throughout all code and
documentation, we use the **standard squared Euclidean norm** in R²⁴:

  N(v) = v · v = Σᵢ vᵢ²

This is the norm returned by `np.dot(v, v)` in the codebase.  Key values:

| Object                  | N(v) = v·v |
|-------------------------|------------|
| E8 roots                | 2          |
| Min(Λ) (minimal shell)  | 8          |
| Second shell of Λ       | 12         |

Wilson (2009) uses a **halved norm** N_W = ½ v·v, under which the minimal
shell has N_W = 4 and E8 roots have N_W = 1.  When citing Wilson's results,
we translate to the standard convention above.

The Leech lattice is constructed in this project via Wilson's
characterisation:  Λ = {(x, y, z) ∈ L³ : conditions 1–3}, where L is the
E8 lattice (= D₈⁺, the Coxeter–Dickson maximal order of integral
octonions).  See `leech_wilson.py` for the implementation.

References:
- Conway, J.H. and Sloane, N.J.A., "Sphere Packings, Lattices and Groups"
  (3rd ed., 1999), Chapter 4.
- Wilson, R.A., "Octonions and the Leech lattice", J. Algebra 322 (2009),
  2186–2190.

### Kirmse integers — Kirmse twist

Johannes Kirmse (born 1894, Schmölln, Thuringia) received his PhD from the
University of Leipzig in 1923 under Gustav Herglotz and published his sole
known paper on integral octonions in 1924, while working as a high school
teacher in Apolda:

> J. Kirmse, "Über die Darstellbarkeit natürlicher ganzer Zahlen als Summen
> von acht Quadraten und über ein mit diesem Problem zusammenhängendes
> nichtkommutatives und nichtassoziatives Zahlensystem,"
> Ber. Verh. Sächs. Akad. Wiss. Leipzig, Math.-Phys. Kl. **76** (1924),
> 63–82.

In this paper Kirmse accomplished the following:

1. **Kirmse's identities**: x(x̄y) = n(x)y = (yx̄)x — appearing in print
   presumably for the first time (Petersson 2018).
2. He exhibited the E8 lattice as a positive definite unimodular integral
   quadratic lattice of rank 8 inside the real octonion algebra, by a
   judicious choice of basis vectors.
3. He may have inadvertently initiated the study of alternative algebras;
   Petersson (2018) speculates that Emil Artin's lifelong interest in
   alternative algebras was inspired by Kirmse's paper, since both Kirmse
   and Artin were students of Herglotz at Leipzig (Artin received his PhD
   in 1921, Kirmse in 1923).

Kirmse then claimed — without proof — that his lattice is closed under
octonion multiplication.  Coxeter (1946), upon studying Kirmse's paper,
was unable to verify this and eventually showed the claim is false: certain
products of several factors, where the multiplication is non-associative,
land outside the lattice.  Coxeter, together with Bruck, then remedied the
defect by a modification of Kirmse's construction (an **index-doubling
permutation** of the imaginary basis elements), producing a genuine maximal
order.

In this project, **Kirmse twist** is used as shorthand for the
index-doubling permutation that corrects Kirmse's construction to a closed
one.

**A note on attribution and fairness.**  The term "Kirmse integers" in
modern usage refers exclusively to the non-closed system — Kirmse's error.
This practice reduces the man to his one mistake and erases his actual
contributions (the identities, the lattice construction, and the
pioneering investigation of octonion arithmetic).  Petersson (2018), who
gives the most detailed modern account, treats Kirmse with respect, calling
his work "the bold attempt of a young mathematician" tackling "a weird and
bizarre topic even by today's standards, let alone the ones of a hundred
years ago."  This project follows Petersson's example and records Kirmse's
positive contributions alongside the error.

References:
- Kirmse (1924), as cited above.
- Coxeter, H.S.M., "Integral Cayley numbers," Duke Math. J. **13** (1946),
  561–578.
- Petersson, H.P., "Integral octonions," lecture at the Málaga Workshop on
  Non-Associative Algebras, September 6, 2018.  Available at:
  https://www.fernuni-hagen.de/mi/fakultaet/emeriti/docs/petersson/ass.-rem.-int.-oct.pdf
- Conway, J.H. and Smith, D.A., "On Quaternions and Octonions" (2003),
  Chapter 9.

### Maximal order (of integral octonions)

A lattice Γ ⊂ O (the real octonion algebra) that is closed under
multiplication, contains the identity, and is maximal with respect to
inclusion among such lattices.  There are exactly 7 maximal orders in the
octonions, each isometric to the E8 root lattice.  They are obtained by
choosing one of the 7 imaginary basis elements as a distinguished unit and
applying the Kirmse twist (index-doubling permutation).

Note: Petersson (2018) documents that Dickson (1923) had constructed
integral octonions isomorphic to Coxeter's more than twenty years before
Coxeter's 1946 paper; Coxeter himself acknowledged this in a postscript.

References: Conway and Smith (2003), Chapter 9; Coxeter (1946);
Dickson, L.E., "A new simple theory of hypercomplex integers,"
J. Math. Pures Appl. (1923).

### Okubo algebra (pseudo-octonion algebra)

The symmetric composition algebra obtained from an octonion algebra (O, ·)
via the **Petersson construction** with a non-trivial order-3 automorphism τ:
x * y = τ(x̄) · τ²(ȳ).  Named after Susumu Okubo, who first constructed
it on traceless 3×3 matrices over C.  Over R, there is (up to isomorphism)
one Okubo algebra; over other fields there may be several.
References: [MarraniCorradettiZucconi2025], [Elduque2000_Triality].

### Para-Hurwitz algebra (para-octonion algebra)

The symmetric composition algebra obtained from an octonion algebra (O, ·)
by the **Petersson construction** with the identity automorphism (τ = id):
x * y = x̄ · ȳ, where x̄ denotes octonion conjugation.  This is the simplest
symmetric composition algebra.  Reference: [Elduque2000_Triality] §3.

### Petersson construction

Given an octonion algebra (O, ·) and an order-3 automorphism τ of O, the
**Petersson isotope** is the algebra (O, *_τ) with product
x *_τ y = τ(x̄) · τ²(ȳ).  The resulting algebra is always a symmetric
composition algebra.  When τ = id, this gives the para-Hurwitz algebra;
when τ ≠ id, it gives an Okubo algebra.
Reference: [MarraniCorradettiZucconi2025] eq. (1.5); originally due to
H. Petersson (1969).

### Symmetric composition algebra

An 8-dimensional composition algebra (V, *, n) satisfying
(x * y) * x = x * (y * x) = n(x) y for all x, y.  Over a field of
characteristic ≠ 2, 3, these are exactly the **para-Hurwitz algebras** and
the **Okubo algebras**.
Reference: [Elduque2000_Triality], [KMRT, Chapter VIII].

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

### Transposition-twisted triple octonion product

(Trial 007, prompt 025.)  The triple octonion product (see below) where the
standard Fano-plane multiplication is replaced by a **transposition-twisted**
multiplication: apply any transposition (s ↔ t) of two imaginary basis
elements {1,...,7} to the standard Fano triples before building the
OctonionAlgebra.  This changes the signs of certain structure constants,
producing a different but isomorphic octonion algebra on the same R⁸.

All 21 transpositions achieve **100% Leech lattice closure** on 593,400
tested pairs (including 518,400 exhaustive type1×type1 and 50,000
type3×type3).  This fixes the Wilson condition 3 failure of the
(un-twisted) triple octonion product.

Since all transpositions are in the same orbit under the Fano-plane
automorphism group GL(3, F₂), they all produce the same multiplication
table up to basis relabeling.  The transposition-twisted product is
therefore essentially unique (up to Fano-plane automorphism).

### Triple octonion product

(Trials 001, 002.)  A 24-dimensional algebra R²⁴ = O₁ ⊕ O₂ ⊕ O₃ built
from three copies of the standard octonion algebra with **Z₃ cross-block
routing**: products of vectors from blocks α and β land in block γ, where
{α, β, γ} = {1, 2, 3}, using the standard octonion multiplication.
Same-block products Oα × Oα → Oα also use the standard octonion product.

Trial 001 tested this product without rescaling.  Trial 002 tested it with
per-block scaling.  Both failed, but trial 001 exhibited a distinctive
**2+1 closure pattern** (see above).

### √3 obstruction

The failure mode specific to the Petersson triality triple (trial 005):
the Petersson construction with a non-trivial order-3 automorphism τ
introduces structure constants involving cos(2π/3) = −1/2 and
sin(2π/3) = √3/2.  Products of E8 lattice vectors (integer or
half-integer coordinates) then acquire irrational √3 components and
cannot lie in the E8 lattice.  This is a more fundamental obstruction
than the Wilson-condition-3 failure of the triple octonion product.

---

Last updated: 2026-04-11 (isotopy/isomorphism, Leech lattice norm conventions added)
