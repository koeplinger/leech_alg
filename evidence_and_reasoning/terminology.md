# Terminology Reference

Created: 2026-04-11.  Updated: 2026-10-28.

This file defines terminology used in this project. Where established
mathematical terminology exists, it is cited. Where constructions are
specific to this project, the terms are listed under **Project-specific
terms**.

---

## Index

**Established terms:**
Idempotent — Isomorphism (of algebras) — Isotopy (of algebras) —
Kirmse integers — Kirmse twist (project shorthand; the paper's Coxeter–Bruck transposition) — Leech lattice Λ —
Maximal order (of integral octonions; the paper's 𝒪) — Moufang identities;
alternative and flexible laws — Okubo algebra — Order —
Para-Hurwitz algebra — Petersson construction —
Square-zero element (square nilpotent) — Standard octonion product —
Symmetric composition algebra — Triality (D₄ triality)

**Project-specific terms:**
2+1 closure pattern — Automorphism group Aut(Λ, +, ⋆) —
Ideal decomposition R²⁴ = D ⊕ T — Petersson triality triple —
Polarization of L/2L —
Reflection ρ — Span index — Triple octonion product — Twist —
Z₃-symmetric triple-octonion product
(triple product for short; in the trial files: transposition-twisted
triple octonion product) —
√3 obstruction

---

## Established terms

### Idempotent

An element u of an algebra (A, ·) with u · u = u.  (0 is always one; a
*nonzero* idempotent is the interesting case.)

In this project (main paper §5.2;
`verify_idempotent_classification.py`): the ambient algebra
(R²⁴, +, ⋆) has **exactly eight** idempotents.  Throughout this file, ⋆
is the Z₃-symmetric triple-octonion product (see the project-specific
section below), in the paper's notation.  With ε₁ = (e₀, 0, 0),
ε₂ = (0, e₀, 0), ε₃ = (0, 0, e₀) and ω = ε₁ + ε₂ + ε₃, they are

    0;   ε₁, ε₂, ε₃;   ω/3;   εᵢ − ω/3  (i = 1, 2, 3).

Every block of every one is a real multiple of e₀: there is no
idempotent with a nonzero imaginary part, and no positive-dimensional
family.  **None of the seven nonzero idempotents lies in Λ**, because
e₀ ∉ L; the order (Λ, +, ⋆)
therefore has *no nonzero idempotent*.  Each nonzero idempotent does
span a rational ray meeting Λ, the least positive lattice multiples
being 4εᵢ (norm 16), 6(εᵢ − ω/3) (norm 24) and 4ω (norm 48), each
satisfying u ⋆ u = n u for its own multiplier n
(`verify_idempotent_lattice_rescaling.py`).

Not to be confused with an **identity element**: (R²⁴, +, ⋆) has none,
so the order is nonunital (main paper §5.1).

### Isomorphism (of algebras)

An **isomorphism** between two algebras (A, ·) and (B, *) is a bijective
linear map φ: A → B satisfying φ(x · y) = φ(x) * φ(y) for all x, y ∈ A.
A single map governs both input and output.

Isomorphism is a special case of **isotopy** (see below): it corresponds to
the isotopy triple (φ, φ, φ) where all three maps are the same.

In this project: the σ-twisted octonion product ·_σ (main paper
Definition 3.2) is isomorphic to the standard octonion product · of eq.
(2.1).  The map σ: e_i ↔ e_j (extended linearly, fixing e₀)
satisfies σ(x · y) = σ(x) ·_σ σ(y) for all x, y, by definition
(σ² = id); `consistency_checks.py` (check 2) confirms it for σ = (1 2)
on the 64 basis pairs.  For the coordinate-symmetric placement
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
composition algebra) has fundamentally different properties: for example,
it has no identity element.

Reference: McCrimmon, K., "A Taste of Jordan Algebras" (2004), §5.3;
Albert, A.A., "Non-associative algebras. I. Fundamental concepts and
isotopy", Ann. of Math. (2) 43 (1942),
685–707 (introduced the concept of isotopy for nonassociative algebras).

### Leech lattice Λ — norm conventions

The **Leech lattice** is the unique even unimodular lattice of rank 24
with no vectors of squared length 2.  Its automorphism group is the Conway
group Co₀ = 2·Co₁, a double cover of the sporadic simple group Co₁.  As in
the paper (§1), Λ denotes the √2-scaled copy of this lattice, at the scale
of Wilson's construction; its minimal shell Min(Λ) consists of the 196,560
vectors of squared length 8.

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
characterization (paper §2.3), which starts from the E8 lattice L = D₈⁺ of
paper §2.2, a √2-scaled copy of the maximal order 𝒪 of the same section
(L = (e₁ + e₂)𝒪; see "Maximal order" below).  With
s = ½(−e₀ + e₁ + ⋯ + e₇) ∈ L and s̄ = ½(−e₀ − e₁ − ⋯ − e₇) ∉ L,
**Wilson's sublattices** are Ls̄ = {ℓ·s̄ : ℓ ∈ L} and Ls = {ℓ·s : ℓ ∈ L},
each of index 16 in L, with Ls̄ + Ls = L and Ls̄ ∩ Ls = 2L.  Wilson's three
conditions define the lattice (paper Definition 2.1):
Λ = {(x, y, z) ∈ L³ : (1) x, y, z ∈ L; (2) x + y, x + z, y + z ∈ Ls̄;
(3) x + y + z ∈ Ls}.  Its 196,560 minimal vectors fall into three types by
the number of nonzero blocks: Type 1 (720), Type 2 (11,520) and Type 3
(184,320); paper §2.3 lists their forms.  See `leech_wilson.py` for the
implementation.

References:
- Conway, J.H. and Sloane, N.J.A., "Sphere Packings, Lattices and Groups"
  (3rd ed., 1999), Chapter 4.
- Wilson, R.A., "Octonions and the Leech lattice", J. Algebra 322 (2009),
  2186–2190.

### Kirmse integers — Kirmse twist

Johannes Kirmse (born 1894, Schmölln, Thuringia) received his PhD from the
University of Leipzig in 1923 under Gustav Herglotz and published his sole
known paper on integral octonions in 1924, under the by-line "Johannes
Kirmse in Hannover" (Kirmse 1924, p. 63):

> J. Kirmse, "Über die Darstellbarkeit natürlicher ganzer Zahlen als Summen
> von acht Quadraten und über ein mit diesem Problem zusammenhängendes
> nichtkommutatives und nichtassoziatives Zahlensystem,"
> Ber. Verh. Sächs. Akad. Wiss. Leipzig, Math.-Phys. Kl. **76** (1924),
> 63–82.

In this paper Kirmse accomplished the following:

1. **Kirmse's identities**: x(x̄y) = n(x)y = (yx̄)x, appearing in print
   presumably for the first time (Petersson 2018).
2. He exhibited the E8 lattice as a positive definite unimodular integral
   quadratic lattice of rank 8 inside the real octonion algebra, by a
   judicious choice of basis vectors.
3. He may have inadvertently initiated the study of alternative algebras;
   Petersson (2018) speculates that Emil Artin's lifelong interest in
   alternative algebras was inspired by Kirmse's paper, since both Kirmse
   and Artin were students of Herglotz at Leipzig (Artin received his PhD
   in 1921, Kirmse in 1923).

Kirmse then claimed, without proof, that his lattice is closed under
octonion multiplication.  In fact Kirmse stated that there are *eight*
such maximal domains (his p. 70); the true number is **seven**, and the
one he exhibited explicitly, $J_1$, is the spurious entry.  Coxeter
(1946, §4) identified the failure with an explicit one-line witness,
$\tfrac12(1{+}i_1{+}i_2{+}i_3)\cdot\tfrac12(1{+}i_1{+}i_4{+}i_5)
=\tfrac12(i_1{+}i_2{+}i_4{+}i_7)\notin J_1$, and recorded the remedy
supplied by **R. H. Bruck**: "Bruck's domain $J$ can be derived from
$J_1$ by transposing two of the $i$'s": a **transposition of two
imaginary basis units**, a linear involution of $\mathbb{R}^8$.  The
$\binom{7}{2}=21$ such transpositions yield the seven maximal orders
in seven sets of three.  Verified from the originals in this project
(`python_project/src/verify_kirmse_1924.py`,
`verify_coxeter_1946.py`; notes in `paper/reviews/`).

**Kirmse twist (this project's shorthand, used in the trial record) = the Coxeter–Bruck transposition of the main paper's Appendix B = the transposition $\sigma$ of Definition 3.1.** The
Kirmse twist is shorthand for the transposition of two imaginary basis
units described above; it is the **same map** as the $\sigma$ used
throughout the main paper, a linear involution of $\mathbb{R}^8$.
There is no difference of operation, only of name.

What differs is the lattice on which the twist is applied. Kirmse's
$J_1$ (and the seven Kirmse–Coxeter–Bruck orders that the twist relates
to it) are **axis-asymmetric** placements of $E_8$ in $\mathbb{R}^8$:
their coordinate stabilizers in $S_8$ are isomorphic to
$\mathrm{AGL}(3,2)$, of order $1344$, and contain no transposition, so every transposition
moves them. It is exactly this axis-asymmetry that lets the Kirmse
twist turn a lattice that does not close under the octonion product
(Kirmse's $J_1$) into one that does (Bruck's $J$).

The main paper, and the Leech-lattice construction of Wilson on which
it builds, take the opposite representation: $L = D_8^+$ is the
**axis-symmetric** placement of $E_8$, invariant under every
coordinate permutation. The Kirmse twist therefore *fixes* $L$
($\sigma(L) = L$), and $L$ closes under the octonion product directly.
The twist's effect in the main paper lies entirely on sublattices
of $L$ that are *not* $\sigma$-invariant (chiefly $Ls$), where it
plays the same kind of role (turning a nonclosed sublattice into a
closed one) that the Kirmse twist plays on Kirmse's $J_1$ one
stratum up.

On $L\bar s$ the twist acts on Kirmse's set itself, up to scale.  In the
main paper's coordinates $L\bar s = 2J_1$, with $J_1$ the analogue of
Kirmse's set for the multiplication of the paper's eq. (2.1) (Wilson's
$L\bar s = 2B$, with $B$ the set of "Kirmse's mistake"; Wilson 2009,
p. 2187), and $\sigma(J_1) = \mathcal O$, the maximal order of the main
paper's §2.2, so that $\sigma(L\bar s) = 2\mathcal O$ (claim 019).  Unlike
$J_1$, the lattice $L\bar s = 2J_1$ is already closed under that
multiplication, and $L \cdot L\bar s \subseteq L\bar s$ holds without the
twist (`symbolic_proof_checks.py`), so no closure on $L\bar s$ depends on
the twist; the twist is needed for $Ls$ (main paper Lemma 4.4).

**A note on attribution and fairness.**  The term "Kirmse integers" in
modern usage refers exclusively to the nonclosed system, Kirmse's error.
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

For a fixed octonion multiplication, the orders of the octonion algebra
over Q (the Q-span of e₀, …, e₇) are the lattices Γ that are closed under
multiplication and contain the identity; their elements are integral
octonions.  A maximal order
is one maximal with respect to inclusion among them.  There are infinitely
many maximal orders (G₂(Q) is infinite, and only finitely many of its
elements preserve a given maximal order; the external review of v8,
§6.1(a)); exactly **seven** of them contain all eight basis units (the
lattice $\mathbb{Z}\langle e_0, \ldots, e_7\rangle$ they span), each an $E_8$ lattice (at
minimal norm 1, so a 1/√2-scaled copy of the norm-2 root lattice) (Coxeter
1946; Conway–Smith 2003, Ch. 9; for the multiplication of the paper's eq.
(2.1), exhaustively by `verify_sigma_sublattice_maximal_order.py`).  The
$\binom{7}{2}=21$ transpositions of pairs of imaginary basis units, applied
to Kirmse's $J_1$, yield the seven (Coxeter 1946, §4); they fall into seven sets
of three, with each set fixing one common imaginary basis unit (Coxeter's
"special unit", one of $i_1,\ldots,i_7$).

Note: Dickson (1923, Theorem XV) gave a correct construction of one
maximal order more than twenty years before Coxeter; his theorem
nonetheless states three maximal orders rather than seven, an undercount
arising from an extra assumption he states openly in the derivation
(p. 321: that the orders contain the Hurwitz quaternions
$\tfrac12(\pm1\pm i\pm j\pm k)$).  Coxeter (1946, §14, postscript)
identifies Dickson's undercount with this cause.  Verified from the originals in
this project (`verify_dickson_1923.py`,
`verify_coxeter_1946.py`).

In this project (main paper §2.2, eq. (2.2); `verify_structural_proof.py`):
the main paper uses one maximal order throughout, 𝒪 := ½(e₁ + e₂)L.  It
contains the eight basis units, is closed under the standard octonion
product, and is one of the seven maximal orders above, the ring A₃ of
Wilson (2009, §2); the ledger (claims 019, 021) writes it O⁽⁴⁾.  Since
(e₁ + e₂)² = −2 and L = −L, the left alternative law x(xy) = x²y gives L = (e₁ + e₂)𝒪, so
L ⊆ 𝒪 and L is a √2-scaled copy of 𝒪.  Lemmas 4.2 to 4.4 of the main
paper (§4.1) rest on three facts: σ(Ls̄) = 2𝒪 and σ(Ls) = p𝒪 with
p = e₂ − e₄ + e₅ − e₆ (eq. (4.1)); (xa)(xb) = x((ax)b) for every purely
imaginary octonion x and all octonions a, b (eq. (4.2)); hence
(xM)(xM) ⊆ xM, where xM := {x·m : m ∈ M}, for every lattice M closed under
the standard octonion product and every purely imaginary x ∈ M (eq. (4.3)).
With M = 𝒪, x = e₁ + e₂ gives Lemma 4.2 and x = p gives Lemma 4.4; Lemma
4.3 follows from L ⊆ 𝒪 and the closure of 𝒪.  The lattice σ(Ls̄) = 2𝒪 is
a two-sided ideal of L; σ(Ls) = p𝒪 is not an ideal of L on either side
(Remark 4.5).

References: Dickson (1923); Kirmse (1924); Coxeter (1946); Conway and
Smith (2003), Chapter 9.

### Moufang identities; alternative and flexible laws

An algebra is **alternative** when it satisfies the left alternative law
x(xy) = x²y and the right alternative law (yx)x = yx² for all x, y.  An
alternative algebra is **flexible**, x(yx) = (xy)x, and satisfies the
**Moufang identities**, among them the left Moufang identity
(xyx)z = x(y(xz)).  The octonions are alternative.

In this project (main paper §4.1; `verify_structural_proof.py`, item (D)):
from these laws and x² = −N(x), the paper derives (4.2),
(xa)(xb) = x((ax)b), for purely imaginary x and all octonions a, b.  For
general x the two sides differ by 2 Re(x)[x, a, b], where
[x, a, b] = (xa)b − x(ab).  Ledger claim 021 records the literature of (4.2).

References: Schafer (1966), Chapter III; Conway and Smith (2003).

### Okubo algebra (pseudo-octonion algebra)

The symmetric composition algebra obtained from an octonion algebra (O, ·)
via the **Petersson construction** with a nontrivial order-3 automorphism τ:
x * y = τ(x̄) · τ²(ȳ).  Named after Susumu Okubo, who first constructed
it on traceless 3×3 matrices over C.  Over R, there is (up to isomorphism)
one Okubo algebra with positive definite norm, the pseudo-octonion
algebra (a division algebra); over other fields there may be several.
References: [MarraniCorradettiZucconi2025], [Elduque2000_Triality].

### Order

(Main paper Corollary 1.2 with its footnote, and Remark 1.3.)  An
**order** in a finite-dimensional real algebra is a subring that is a free
Z-module of rank equal to the algebra's dimension: a full-rank lattice
closed under the product.  The classical definition also asks for the
identity element (see "Maximal order" above); the ambient algebra
(R²⁴, +, ⋆) has none, and the paper uses the term "in the nonunital
sense" (Corollary 1.2, footnote).  (Λ, +, ⋆) is such an order, and so is L
under the standard octonion product (L · L ⊆ L, Lemma 4.2; e₀ ∉ L).  An
order is **maximal** when no larger order of the same algebra contains it;
(Λ, +, ⋆) is not maximal, since L³ is a ⋆-closed overlattice of index 2¹²
(claim 016; main paper §5.3).

### Para-Hurwitz algebra (para-octonion algebra)

The symmetric composition algebra obtained from an octonion algebra (O, ·)
by the **Petersson construction** with the identity automorphism (τ = id):
x * y = x̄ · ȳ, where x̄ denotes octonion conjugation.  This is the simplest
symmetric composition algebra.  Reference: [Elduque2000_Triality] §3.

### Petersson construction

Given an octonion algebra (O, ·) and an automorphism τ of O with τ³ = id, the
**Petersson isotope** is the algebra (O, *_τ) with product
x *_τ y = τ(x̄) · τ²(ȳ).  The resulting algebra is always a symmetric
composition algebra.  When τ = id, this gives the para-Hurwitz algebra;
when τ ≠ id, it gives an Okubo algebra.
Reference: [MarraniCorradettiZucconi2025] eqs. (1.5)–(1.6); originally due to
H. Petersson [Petersson1969].

### Square-zero element (square nilpotent)

An element u of an algebra (A, ·) with u · u = 0 and u ≠ 0.  In a
nonassociative algebra the general word "nilpotent" is ambiguous
(powers must be bracketed), so this project says **square-zero** and
means exactly u ⋆ u = 0.

In this project (main paper §5.5;
`verify_square_zero_classification.py`):

- In the ambient algebra (R²⁴, +, ⋆), u = (x, y, z) satisfies u ⋆ u = 0
  **iff** x, y, z are purely imaginary, of equal norm, and sum to zero:
  a hexagonal triple in Im(O) = R⁷.  The solution set is a cone of
  dimension 12; it lies in T and meets D only in 0.
- **The order (Λ, +, ⋆) does contain square-zero elements, with norms
  in 12Z** (N(u) = 3 N(x) with N(x) ∈ 2Z, and every norm of Λ lies in
  4Z).  The norm-12 stratum is exactly the 4,032 ordered hexagonal triples of
  the E₇ root system formed by the 126 norm-4 vectors of Ls̄ ∩ Im(O);
  norm 24 occurs as well (witness:
  `verify_square_zero_norm24_witness.py`).
- **Min(Λ) contains none**: the minimal norm is 8, not a multiple
  of 12.

### Standard octonion product

(Main paper §2.1, eq. (2.1).)  The multiplication of the real octonion
algebra O with basis e₀, …, e₇, e₀ the identity, determined by the seven
**Fano triples** (1,2,4), (2,3,5), (3,4,6), (4,5,7), (5,6,1), (6,7,2),
(7,1,3): for each triple (a, b, c), e_a · e_b = e_c, cyclic permutations
carrying a positive sign and anti-cyclic ones a negative sign [Baez2002];
equivalently, e_a · e_{a+1} = e_{a+3} for a = 1, …, 7 (indices mod 7 in
{1, …, 7}).  Each **imaginary basis unit** e₁, …, e₇ squares to −e₀.  An
octonion x = Σ xᵢeᵢ is **purely imaginary** when x₀ = 0, that is, when
x̄ = −x.  The σ-twisted octonion product ·_σ (main paper Definition 3.2)
and the twisted octonion products ·_g (see "Twist" below) are defined
from it.

### Symmetric composition algebra

An 8-dimensional composition algebra (V, *, n) satisfying
(x * y) * x = x * (y * x) = n(x) y for all x, y.  Over a field of
characteristic ≠ 2, 3, these are exactly the **para-Hurwitz algebras** and
the **Okubo algebras**.
Reference: [Elduque2000_Triality]; [KnusMerkurjevRostTignol1998], Chapter VIII.

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
  populate one and two of the three blocks, respectively, while type-3
  vectors populate all three, this is effectively a 2+1 split based on
  block occupancy.

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
`gap_project/aut_lambda_star.g`; main paper §5.4;
standalone note `paper/automorphism_group_2026-07-12.pdf`):

- Aut(Λ, +, ⋆) is **finite of order 36**, of structure C₆ × S₃.
  Center C₆, derived subgroup C₃; −I₂₄ is *not* in it.  Element orders:
  1 (×1), 2 (×7), 3 (×8), 6 (×20).
- Every element is a blockwise octonion automorphism composed with a
  permutation of the three blocks.  The S₃ permutes the three blocks;
  the C₆ is generated by the signed permutation of imaginary basis units
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
- **Completeness**: the order-36 count and the Co₀ containment are
  established by two independent routes in the standalone note
  (`paper/automorphism_group_2026-07-12.pdf`): Route 1 uses the
  classification of the automorphisms of the complexification of
  (T, ⋆) (nonenumerative); Route 2 is entirely finite and
  machine-checkable, so nothing load-bearing is left open.  The
  machine-verified lower bound C₆ × S₃ ⊆ Aut(Λ, +, ⋆) and
  the strictness witness are independent of both.

Related, and not to be confused with the above: the stabilizers of the
Wilson sublattices in Aut(O, ·_σ), the automorphism group of the σ-twisted
product (the setting of the standalone note
`paper/automorphism_group_2026-07-12.tex`, §5.2): |Stab(L)| = 1344 =
2³·L₃(2), |Stab(Ls)| = 168 = 2³:(7:3), |Stab(Ls̄)| = 12096 = U₃(3).2 =
G₂(2) (`gap_project/octonion_stabilisers.g`).  Under the standard octonion
product the orders are 1344, 21 and 1344.

Do not confuse either with Aut(Λ) = 2·Co₁, the isometry group of the
*lattice alone* (no product); see "Leech lattice Λ" above.

### Ideal decomposition R²⁴ = D ⊕ T

(Main paper §5.4; `verify_ideal_decomposition.py`,
`verify_star_algebra_structure.py`.)  As a real algebra,
R²⁴ = O³ under the Z₃-symmetric triple-octonion product ⋆ splits as a
direct sum of two **mutually annihilating two-sided ideals**:

- **D** = {(a, a, a) : a ∈ O}, the **diagonal copy of the octonions**,
  of dimension 8.  Here a is an *octonion*, not a real number.
- **T** = {(p, q, r) ∈ O³ : p + q + r = 0}, of dimension 16.

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

Cross-block products use the target block's algebra.  Rotating the
mediator of a Fano line cyclically gives three order-3 automorphisms: on
each of the 7 Fano lines, each of the 3 elements mediates one, for one of
the two relative senses of the two rotations (worked examples in
`python_project/src/okubo_samples.py`).  Trial 005 tests the Petersson
triality triple, not the algebras these three give.

The construction failed due to a **√3 obstruction**: the Petersson
construction with trial 005's τ (a rotation by 2π/3 in two planes)
introduces irrational structure constants
(cos 2π/3, sin 2π/3), so products of E8 lattice vectors leave the E8
lattice entirely.

### Polarization of L/2L

(Main paper Appendix A, "The closure criterion modulo 2", label
`app:mod2quotient`; `verify_mod2_quotient.py`,
`verify_polarization_lagrangians.py`.)  Since L · L ⊆ L, the standard
octonion product descends to an F₂-bilinear product μ̄ on L̄ := L/2L ≅ F₂⁸.
A lattice between 2L and L is closed exactly when its image in L̄ is a
subalgebra, so Proposition 5.3 (`prop:criterion`) reads modulo 2: ⋆_g
closes on Λ if and only if g(Ls̄)/2L and g(Ls)/2L are subalgebras of
(L̄, μ̄).  With q(v + 2L) := N(v)/2 mod 2, a subspace is **totally
isotropic** when q vanishes on it; the 4-dimensional ones are the
**Lagrangians**, and an ordered pair of complementary Lagrangians is a
**polarization**.  Wilson's pair (Ls̄/2L, Ls/2L) is a polarization, and σ
maps it to V := σ(Ls̄)/2L = 2𝒪/2L and W := σ(Ls)/2L = p𝒪/2L (display
(A.3)).  Lemmas 4.3 and 4.4, read modulo 2L, say that V is a left ideal and
W a subalgebra; V is also a two-sided ideal, while Ls/2L is not a
subalgebra, so the identity twist does not close.  Since W(E₈) acts on L̄ as
the full orthogonal group of q (Wilson's book, §3.12.4), every polarization
is (g(Ls̄)/2L, g(Ls)/2L) for an isometry g of L, and ⋆_g closes exactly when
both halves are subalgebras.  Credits (the appendix's "Relation to prior
work"): the orthogonal space (L̄, q) is classical [WilsonBook2009];
Lepowsky and Meurman (1982, §2) identify the lattices between 2L and L
isometric to √2 L with the preimages of the Lagrangians and construct a
lattice isometric to Λ from any polarization, noting Tits; Silverman found
the identification independently (part 12 of [Baez2014]), read Egan's
construction modulo 2 and found that Wilson's pair is a polarization
[Silverman2016]; Calderbank related Egan's construction to Lepowsky and
Meurman's [Calderbank2023].

### Reflection ρ

(Main paper §5.6, label `sec:census`; `verify_twist_span_index_census.py`;
claims 017, 018.  Written s_r in the external review of v8.)  The
reflection of R⁸ in the root r = ½(e₀ + e₁ + ⋯ + e₆ − e₇) of L, a twist (see
"Twist" below).  The triple product ⋆_ρ closes on Λ with span index 1.  The
transposition σ is itself the reflection in the root e₁ − e₂ of L.  Of the
120 reflections in roots of L, 56 close on Λ, 28 with span index 2¹⁶ and
28 with span index 1 (claim 017).  Span index 1 does not make (Λ, +, ⋆_ρ)
a maximal order: L³ is a ⋆_ρ-closed overlattice of Λ of index 2¹² (claim
018).

### Span index

(Main paper §5.6, label `sec:census`, and Proposition 5.4, `prop:span`;
`verify_product_span_index.py`, `verify_product_span_structure.py`,
`verify_product_span_proof.py`, `verify_twist_span_index_census.py`; claims
015, 017, 018.)  When a triple product ⋆_g (see "Twist" below)
closes on Λ, the integer combinations of its products u ⋆_g v (u, v ∈ Λ)
form a sublattice of Λ, the **span** of ⋆_g.  The **span index** of ⋆_g is
the index of the span in Λ: the number of disjoint translates of the span
that make up Λ.  Span index 1 means that every vector of Λ is an integer
combination of products.  An order with a multiplicative identity has span
index 1; the ambient algebra (R²⁴, +, ⋆) has none, and the index can be
larger.  By bilinearity, closure and span index are decided by the 576
products of pairs from a Z-basis of Λ; the index is computed exactly by
Hermite normal form.  Values: the span of ⋆ = ⋆_σ consists of the vectors
(w, w, w) + 2u with w ∈ Ls and u ∈ Λ, so it contains 2Λ and has span index
2¹⁶ = 65,536 (Proposition 5.4, proved by hand; claim 015); every odd
closing permutation twist has span index 2¹⁶ and every even one 2¹⁰; the
reflection ρ has span index 1 (claim 017); the induced product φ of the
Baez–Egan closure, at the paper's scale (main paper §6, `sec:related`), has
span index 2²⁵ (`verify_phi_span_index.py`).  The span index is invariant under
ring isomorphism, so orders with different span indices are not isomorphic:
an order from an odd closing permutation (2¹⁶) is not isomorphic to one
from an even closing permutation (2¹⁰), and neither is isomorphic to
(Λ, +, ⋆_ρ) (span index 1).  It is independent of maximality
(claim 018).

### Z₃-symmetric triple-octonion product
(triple product for short; in the trial files: transposition-twisted
triple octonion product)

(Trial 007, prompt 025; main paper Definition 3.3.)  The triple
octonion product (see below) built from the σ-twisted octonion product
·_σ in all three blocks, where ·_σ is defined by
x ·_σ y := σ(σ(x) · σ(y)) for a transposition σ = (i j) on {1, …, 7},
extended to the linear involution of R⁸ that fixes e₀ and swaps the
imaginary basis units e_i and e_j (main paper Definitions 3.1 and 3.2).
The product ·_σ is simply another octonion product on R⁸,
isomorphic to the standard one via σ itself.  The transposition σ is
the device that exhibits the fit of this octonion product with
Wilson's representation of the Leech lattice: σ leaves L = D₈⁺
invariant but moves Wilson's sublattices Ls and Ls̄.

All 21 transpositions close on Λ.  The evidence is the symbolic proof of
paper §4 for σ = (1 2), together with the exact Z-basis test of all 21
(`verify_consecutive_twists_exact.py`); random-pair sampling
(`trial_007_kirmse_twist.py`, 900 pairs per transposition, and the 12M+ pairs
of the scaled harnesses) agrees.  This fixes the failure of Wilson's condition
(3) for the triple octonion product taken over Wilson's own multiplication
convention.

The 21 transpositions form a single orbit under conjugation by the
order-21 group of coordinate permutations preserving the oriented Fano
triples.  Each such permutation is an octonion automorphism that fixes s
and s̄ and preserves L, so conjugation carries the paper's §4 proof for
σ = (1 2) to every transposition (paper §4, first footnote).

Transpositions are not the only permutations of the imaginary basis units
for which the assembled triple product closes.  Writing
x ·_π y := π⁻¹(π(x) · π(y)) for any
permutation π of e₁, …, e₇ (the paper's twisted product ·_g with
g = π; see "Twist" below) and testing exactly on a Z-basis of Λ, the
census over all 5,040 permutations of S₇ is: transposition 21/21
(100%), 3-cycle 35/70 (50%), (2,2) 42/105 (40%), 4-cycle 0/210 (0%),
(3,2) 126/420 (30%), 5-cycle 252/504 (50%), (2,2,2) 21/105 (20%),
(3,3) 112/280 (40%), (4,2) 294/630 (46.7%), 6-cycle 210/840 (25%),
(3,2,2) 105/210 (50%), (4,3) 168/420 (40%), (5,2) 42/504 (8.3%),
7-cycle 336/720 (46.7%); total 1,764/5,039 (35.0%); the identity does
not close.  Exhaustive, not sampled
(`verify_all_permutations_exact.py`; partial checks in
`verify_all_cycles_exact.py` (the 3- to 7-cycles, exact) and
`consistency_checks.py` (three transpositions, sampled); main paper
§5.6).

The name "transposition-twisted triple octonion product", used in the
trial 007 files, over-privileges the twist: the
block product ·_σ is not an exotic twisted object but an ordinary octonion
product; the twist concerns its fit with a chosen representation.

### Triple octonion product

(Trials 001, 002.)  A 24-dimensional algebra R²⁴ = O₁ ⊕ O₂ ⊕ O₃ built
from three copies of the standard octonion algebra with **Z₃ cross-block
routing**: products of vectors from blocks α and β land in block γ, where
{α, β, γ} = {1, 2, 3}, using the standard octonion multiplication.
Same-block products Oα × Oα → Oα also use the standard octonion product.

Trial 001 tested this product without rescaling.  Trial 002 tested it with
per-block scaling, each product normalized to norm 8 (closure of the minimal
shell).  Both failed, but trial 001 exhibited a distinctive
**2+1 closure pattern** (see above).

In the main paper's notation this product is ⋆_g with g the identity
(see "Twist" below), which the paper calls the untwisted product
(Appendix C).  By Proposition 5.3 (`prop:criterion`) it closes on Egan's
lattice Λ(L₁, L₂) (main paper §6, display (6.2), `eq:egan`) exactly when L₁
and L₂ are closed under the standard product: on σ(Λ) = Λ(2𝒪, p𝒪) it
closes, on Wilson's Λ = Λ(Ls̄, Ls) it does not, since Ls is not closed
(Remark 4.5).  The paper's *triple product* (Definition 3.3) is the
Z₃-symmetric triple-octonion product ⋆ = ⋆_σ (see above).

### Twist

(Main paper §5.6, label `sec:census`, and Proposition 5.3,
`prop:criterion`.)  A **twist** is an invertible linear map g of R⁸ with
g(L) = L.  It defines the **twisted octonion product** x ·_g y :=
g⁻¹(g(x) · g(y)) on O, another octonion product, isomorphic to the standard
one via g; the triple product ⋆_g is that of main paper Definition 3.3 with
·_g in all nine block products.  The twist g = σ gives the σ-twisted product
of Definition 3.2 (σ⁻¹ = σ), and ⋆_σ = ⋆.  The **closure criterion**
(Proposition 5.3, proved in the paper): for every twist g, ⋆_g closes on Λ
if and only if Ls̄ and Ls are closed under ·_g, that is, if and only if
g(Ls̄) and g(Ls) are closed under the standard octonion product.  For g = σ
it gives Theorem 1.1 again; the identity does not close, since Ls is not
closed (Remark 4.5).  The abstract and Section 7 (`sec:conclusion`) name
it, Appendix A (`app:mod2quotient`) reads it modulo 2, and Section 6
(display (6.2), `eq:egan`) applies it to Egan's lattices.  The same proof
holds for every Z-bilinear product μ with μ(L, L) ⊆ L (records only:
`evidence_and_reasoning/2026-10-20_closed_form_plan.md`, F1).  It is
Theorem A of the external review of 13 October 2026, stated there for
twists in W(E₈) (claim 017); `verify_twist_span_index_census.py` checks it on
all 5,040 permutations.  The twists tested are the permutations π of
e₁, …, e₇ (the S₇ census, under "Z₃-symmetric triple-octonion product"
above) and the reflections in roots of L (see "Reflection ρ" above).

### √3 obstruction

The failure mode specific to the Petersson triality triple (trial 005):
the Petersson construction with the order-3 automorphism τ of trial 005
(a rotation by 2π/3 in two planes) introduces structure constants involving cos(2π/3) = −1/2 and
sin(2π/3) = √3/2.  Products of E8 lattice vectors (integer or
half-integer coordinates) then acquire irrational √3 components and
cannot lie in the E8 lattice.  This is a more fundamental obstruction
than the failure of Wilson's condition (3) for the triple octonion product.
