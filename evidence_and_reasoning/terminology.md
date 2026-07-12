# Terminology Reference

This file defines terminology used in this project. Where established
mathematical terminology exists, it is cited. Where constructions are
specific to this project, terms are marked as **(project-specific)**.

---

## Index

**Established terms:**
Idempotent — Isomorphism (of algebras) — Isotopy (of algebras) —
Kirmse integers — Kirmse twist — Leech lattice Λ —
Maximal order (of integral octonions) — Okubo algebra —
Para-Hurwitz algebra — Petersson construction —
Square-zero element (square nilpotent) —
Symmetric composition algebra — Triality (D₄ triality)

**Project-specific terms:**
2+1 closure pattern — Automorphism group Aut(Λ, +, ⋆) —
Ideal decomposition R²⁴ = D ⊕ T — Petersson triality triple —
Triple octonion product — Z₃-symmetric triple-octonion product
(historically: transposition-twisted triple octonion product) —
√3 obstruction

---

## Established terms

### Idempotent

An element u of an algebra (A, ·) with u · u = u.  (0 is always one; a
*nonzero* idempotent is the interesting case.)

In this project (main paper §5.2, Remark 5.3;
`verify_idempotent_classification.py`): the ambient algebra
(R²⁴, +, ⋆) has **exactly eight** idempotents.  Throughout this file, ⋆
is the Z₃-symmetric triple-octonion product (see the project-specific
section below), in the paper's notation.  With ε₁ = (e₀, 0, 0),
ε₂ = (0, e₀, 0), ε₃ = (0, 0, e₀) and ω = ε₁ + ε₂ + ε₃, they are

    0;   ε₁, ε₂, ε₃;   ω/3;   εᵢ − ω/3  (i = 1, 2, 3).

Every block of every one is a real multiple of e₀: there is no
idempotent with a nonzero imaginary part, and no positive-dimensional
family.  **None lies in Λ**, because e₀ ∉ L; the order (Λ, +, ⋆)
therefore has *no idempotent at all*.  Each nonzero idempotent does
span a rational ray meeting Λ, the least positive lattice multiples
being 4εᵢ (norm 16), 6(εᵢ − ω/3) (norm 24) and 4ω (norm 48), each
satisfying u ⋆ u = n u for its own multiplier n
(`verify_idempotent_lattice_rescaling.py`).

Not to be confused with an **identity element**: (R²⁴, +, ⋆) has none,
so the order is non-unital (main paper §5.1).

### Isomorphism (of algebras)

An **isomorphism** between two algebras (A, ·) and (B, *) is a bijective
linear map φ: A → B satisfying φ(x · y) = φ(x) * φ(y) for all x, y ∈ A.
A single map governs both input and output.

Isomorphism is a special case of **isotopy** (see below): it corresponds to
the isotopy triple (φ, φ, φ) where all three maps are the same.

In this project: the octonion product ·_σ (main paper Definition 3.2)
is isomorphic to the standard octonion product.  The map σ: e_s ↔ e_t
(extended linearly)
satisfies σ(x ·_std y) = σ(x) ·_σ σ(y) for all basis pairs.  Verified
computationally: 0/64 mismatches.  For the coordinate-symmetric placement
L = D₈⁺ used throughout the paper, σ preserves L (it is a coordinate
permutation, which preserves both integer/half-integer parity and the
coordinate sum); σ moves the sublattices Ls and Ls̄ of L, and that is
where the twist does its work.  See "Kirmse integers — Kirmse twist"
below for the contrast with the axis-asymmetric placement of E₈ used
by Kirmse (1924), where the twist moves the lattice itself.

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
characterization:  Λ = {(x, y, z) ∈ L³ : conditions 1–3}, where L is the
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
octonion multiplication.  In fact Kirmse stated that there are *eight*
such maximal domains (his p. 70); the true number is **seven**, and the
one he exhibited explicitly, $J_1$, is the spurious entry.  Coxeter
(1946, §4) re-derived the failure with an explicit one-line witness,
$\tfrac12(1{+}i_1{+}i_2{+}i_3)\cdot\tfrac12(1{+}i_1{+}i_4{+}i_5)
=\tfrac12(i_1{+}i_2{+}i_4{+}i_7)\notin J_1$, and recorded the remedy
supplied by **R. H. Bruck**: "Bruck's domain $J$ can be derived from
$J_1$ by transposing two of the $i$'s" — a **transposition of two
imaginary basis units**, a linear involution of $\mathbb{R}^8$.  The
$\binom{7}{2}=21$ such transpositions yield the seven maximal orders
in seven sets of three.  Verified from the originals in this project
(\texttt{python\_project/src/verify\_kirmse\_1924.py},
\texttt{verify\_coxeter\_1946.py}; notes in \texttt{paper/reviews/}).

**Kirmse twist = the transposition $\sigma$ of the main paper.** The
Kirmse twist is shorthand for the transposition of two imaginary basis
units described above; it is the **same map** as the $\sigma$ used
throughout the main paper — a linear involution of $\mathbb{R}^8$.
There is no difference of operation, only of name.

What differs is the lattice on which the twist is applied. Kirmse's
$J_1$ (and the seven Kirmse–Coxeter–Bruck orders that the twist relates
to it) are **axis-asymmetric** placements of $E_8$ in $\mathbb{R}^8$:
their coordinate stabilizers in $S_8$ are of order $1344 \cong
\mathrm{AGL}(3,2)$ and contain no transposition, so every transposition
moves them. It is exactly this axis-asymmetry that lets the Kirmse
twist turn a lattice not closed under the octonion product
(Kirmse's $J_1$) into one that does (Bruck's $J$).

The main paper, and the Leech-lattice construction of Wilson on which
it builds, take the opposite representation: $L = D_8^+$ is the
**axis-symmetric** placement of $E_8$, invariant under every
coordinate permutation. The Kirmse twist therefore *fixes* $L$
($\sigma(L) = L$), and $L$ closes under the octonion product directly.
The twist's effect in the main paper lies entirely on sublattices
of $L$ that are *not* $\sigma$-invariant — chiefly $Ls$ — where it
plays the same kind of role (turning a non-closed sublattice into a
closed one) that the Kirmse twist plays on Kirmse's $J_1$ one
stratum up.

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
inclusion among such lattices.  There are exactly **seven** maximal
orders in the octonions, each isometric to the $E_8$ root lattice.  The
$\binom{7}{2}=21$ transpositions of pairs of imaginary basis units act
transitively on the seven (Coxeter 1946, §4); they fall into seven sets
of three, with each set fixing one common imaginary unit (Coxeter's
"special unit" — one of $i_1,\ldots,i_7$).

Note: Dickson (1923, Theorem XV) gave a correct construction of one
maximal order more than twenty years before Coxeter; his theorem
nonetheless states three maximal orders rather than seven, an undercount
arising from an extra assumption he states openly in the derivation
(p. 321: that the orders contain the Hurwitz quaternions
$\tfrac12(\pm1\pm i\pm j\pm k)$).  Coxeter (1946, §14, postscript)
identifies this cause independently.  Verified from the originals in
this project (\texttt{verify\_dickson\_1923.py},
\texttt{verify\_coxeter\_1946.py}).

References: Dickson (1923); Kirmse (1924); Coxeter (1946); Conway and
Smith (2003), Chapter 9.

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

### Square-zero element (square nilpotent)

An element u of an algebra (A, ·) with u · u = 0 and u ≠ 0.  In a
non-associative algebra the general word "nilpotent" is ambiguous
(powers must be bracketed), so this project says **square-zero** and
means exactly u ⋆ u = 0.

In this project (main paper §5.2, Remark 5.4;
`verify_square_zero_classification.py`):

- In the ambient algebra (R²⁴, +, ⋆), u = (x, y, z) satisfies u ⋆ u = 0
  **iff** x, y, z are purely imaginary, of equal norm, and sum to zero:
  a hexagonal triple in Im(O) = R⁷.  The solution set is a cone of
  dimension 12; it lies in T and meets D only in 0.
- **The order (Λ, +, ⋆) does contain square-zero elements: 4,032 of
  them, all of norm 12.**  They arise from the hexagonal triples of the
  E₇ root system formed by the 126 norm-4 vectors of Ls̄ ∩ Im(O).
- **Min(Λ) contains none**, for an arithmetic reason: N(u) = 3 N(x)
  with N(x) ∈ 2Z, and every norm of Λ lies in 4Z, so every square-zero
  vector of Λ has norm in 12Z, while the minimal norm is 8.

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

### Automorphism group Aut(Λ, +, ⋆)

The group of additive bijections of Λ that preserve the product ⋆
(`verify_aut_lambda_star.py`, `verify_aut_octonion_crosscheck.py`,
`gap_project/aut_lambda_star.g`; main paper §5.3, Remark 5.5;
standalone note `paper/automorphism_group_2026-07-12.pdf`):

- Aut(Λ, +, ⋆) is **finite of order 36**, of structure C₆ × S₃.
  Center C₆, derived subgroup C₃; −I₂₄ is *not* in it.  Element orders:
  1 (×1), 2 (×7), 3 (×8), 6 (×20).
- Every element is a blockwise octonion automorphism composed with a
  permutation of the three blocks.  The S₃ permutes the three blocks;
  the C₆ is generated by the signed permutation of imaginary units
  fixing e₀:
  e₁ ↦ −e₅, e₂ ↦ −e₇, e₃ ↦ e₂, e₄ ↦ −e₄, e₅ ↦ e₆, e₆ ↦ e₁, e₇ ↦ −e₃.
- **Ambient group**: Aut(R²⁴, ⋆) = G₂ × G₂ × S₃, compact, of
  dimension 28 (one G₂ on each of the ideals D and T; see below).
- **Aut(Λ, +, ⋆) ⊊ Co₀** is *proved*, from the ambient algebra alone:
  an automorphism of the *real* octonions preserves the
  positive-definite octonion norm, and D ⊥ T, so every ⋆-automorphism
  of R²⁴ is orthogonal.  The inclusion is strict because −I₂₄ ∈ Co₀
  does not preserve ⋆.  (The trace-form route fails: tr(L_u L_v) has
  signature (3, 21), not the Leech form.)
- **Completeness caveat**, as stated in the paper: the order-36 claim
  and the Co₀ containment rest on one non-enumerative step, the
  classification of the automorphisms of the complexification of
  (T, ⋆).  The machine-verified lower bound C₆ × S₃ ⊆ Aut(Λ, +, ⋆) and
  the strictness witness do not depend on it.

Related, and not to be confused with the above: the octonion-automorphism
stabilizers of the Wilson sublattices, |Stab(L)| = 1344 = 2³·L₃(2),
|Stab(Ls)| = 168 = 2³:(7:3), |Stab(Ls̄)| = 12096 = U₃(3).2 = G₂(2)
(`gap_project/octonion_stabilisers.g`).

Do not confuse either with Aut(Λ) = 2·Co₁, the isometry group of the
*lattice alone* (no product); see "Leech lattice Λ" above.

### Ideal decomposition R²⁴ = D ⊕ T

(Main paper §5.3; `verify_ideal_decomposition.py`,
`verify_star_algebra_structure.py`.)  As a real algebra,
R²⁴ = O³ under the Z₃-symmetric triple-octonion product ⋆ splits as a
direct sum of two **mutually annihilating two-sided ideals**:

- **D** = {(a, a, a) : a ∈ O} — the **diagonal copy of the octonions**,
  of dimension 8.  Here a is an *octonion*, not a real number.
- **T** = {(p, q, r) ∈ O³ : p + q + r = 0} — of dimension 16.

D ⋆ T = T ⋆ D = 0, and D ⊥ T for the Euclidean form.  This
decomposition is what settles the automorphism group and the Co₀
containment (see above).  It also organizes the idempotents: ω/3 lies
in D, the three εᵢ − ω/3 lie in T, and the eight idempotents are
exactly the sums d + t of an idempotent d ∈ {0, ω/3} of D and an
idempotent t ∈ {0, ε₁ − ω/3, ε₂ − ω/3, ε₃ − ω/3} of T.  The
square-zero cone lies entirely in T.

Intersections with the lattice: Λ ∩ D = {(a, a, a) : a ∈ Ls};
Λ ∩ T = {(p, q, r) ∈ (Ls̄)³ : p + q + r = 0}; the glue group between
them is (Z/3)⁸, of order 6,561.

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

### Z₃-symmetric triple-octonion product
(historically: transposition-twisted triple octonion product)

(Trial 007, prompt 025.)  The triple
octonion product (see below) built from the octonion product ·_σ in
all three blocks, where ·_σ is defined by x ·_σ y := σ(σ(x) · σ(y))
for a transposition σ = (s t) of two imaginary basis elements
{1,...,7}.  The product ·_σ is simply another octonion product on R⁸,
isomorphic to the standard one via σ itself.  The transposition σ is
the device that exhibits the fit of this octonion product with
Wilson's representation of the Leech lattice: σ leaves L = D₈⁺
invariant but moves Wilson's sublattices Ls and Ls̄.

All 21 transpositions close on Λ.  The evidence is the symbolic proof of
paper §4 for σ = (1 2), together with the exact Z-basis test of all 21
(`verify_consecutive_twists_exact.py`); random-pair sampling
(`trial_007_kirmse_twist.py`, 900 pairs per transposition, and the 12M+ pairs
of the scaled harnesses) agrees.  This fixes the Wilson condition 3 failure of
the triple octonion product taken over Wilson's own multiplication
convention.

Since all transpositions are in the same orbit under the Fano-plane
automorphism group GL(3, F₂), they all produce the same multiplication
table up to basis relabeling.  The construction is therefore
essentially unique (up to Fano-plane automorphism).

Transpositions are not the only permutations of the imaginary axes
whose twist closes.  Writing x ·_π y := π(π(x) · π(y)) for any
permutation π of e₁, …, e₇ and testing exactly on a Z-basis of Λ, the
census over all cycle types of S₇ is: transposition 21/21 (100%),
3-cycle 35/70 (50%), 4-cycle 0/210 (0%), 5-cycle 252/504 (50%),
6-cycle 210/840 (25%), 7-cycle 336/720 (46.7%); total 854/2,365
(36.1%).  Exhaustive, not sampled (`verify_all_cycles_exact.py`,
`consistency_checks.py`; main paper §5.5, Remark 5.7).  Of the 105
(2,2)-double transpositions, 42 close (main paper Remark 4.9,
`verify_consecutive_twists_exact.py`).

The historical name "transposition-twisted triple octonion product"
(used in trial 007 and key claim 008, both preserved per Manifesto
§12) over-privileges the twist: the product is not an exotic twisted
object but an ordinary octonion product; the twist concerns its fit
with a chosen representation.

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
