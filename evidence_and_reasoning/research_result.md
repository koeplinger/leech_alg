# Research Result

Published: J. Köplinger, *An order on the Leech lattice from a Z₃-symmetric
triple-octonion product*, [doi:10.13140/RG.2.2.22093.19686](https://doi.org/10.13140/RG.2.2.22093.19686)
(ResearchGate); this is v7 of the paper.  v8 is in progress in
`paper/main.tex`.

## Summary

The Leech lattice Λ admits an order under a bilinear product derived from
the octonion algebra.

### The product

Let O denote the real octonion algebra with standard Fano-plane
multiplication, and let σ: O → O be the linear map that swaps two
imaginary basis units (any transposition of {e₁, …, e₇}).  Define a
second multiplication ·_σ on O by

    x ·_σ y  =  σ( σ(x) · σ(y) )        (†)

where · is the standard octonion product.  The map σ is an isomorphism
from (O, ·) to (O, ·_σ): σ(x · y) = σ(x) ·_σ σ(y).  So (O, ·_σ) is
again an octonion algebra, not a new kind of algebra; the two
multiplication tables differ in 30 of 64 structure constants, and the
twist is simply how this octonion product fits Wilson's representation
of Λ.  In particular σ is *not* an octonion automorphism; it fixes
L = D₈⁺ setwise but moves the sublattices Ls and Ls̄.

Now set R²⁴ = O₁ ⊕ O₂ ⊕ O₃, three copies of (O, ·_σ), all using the
**same** swapped multiplication, with Z₃ cross-block routing:

    Same block:      O_α × O_α → O_α
    Cross block:     O_α × O_β → O_γ      where {α,β,γ} = {1,2,3}

All nine block-pair products use the single multiplication ·_σ.  The
resulting product ⋆ on R²⁴ is the paper's Z₃-symmetric triple-octonion
product (Definition 3.3).

### The claim

For every u, v ∈ Min(Λ) (the 196,560 vectors of norm 8), the
product u ⋆ v lies in Λ.

Since Min(Λ) generates Λ over Z, bilinearity extends closure to all of Λ:
for every u, v ∈ Λ, the product u ⋆ v lies in Λ.  Therefore (Λ, +, ⋆) is
an order in the 24-dimensional R-algebra (R²⁴, +, ⋆).

### What makes this simple

The entire construction has three ingredients:

1. **One octonion product**, related to Wilson's standard convention
   by a single transposition σ of two imaginary basis units.
2. **Three identical copies** of that product on R²⁴.
3. **Z₃ cross-block routing**: products of vectors from different blocks
   land in the third block.

No rescaling, no conjugation, no distinct algebras on different blocks.
The same product on all three copies.

### Verification status

| Test | Pairs / samples | Failures |
|---|---|---|
| Trial 007 base (per transposition) | 900 | 0 |
| Scaled tests (three harnesses: `trial_007_{scaled_test,fast,exhaust}.py`) | 12,000,000+ | 0 |
| All 21 transpositions (small sample) | ~15,000 | 0 |
| Symbolic proof (Lemmas 4.1–4.4, paper §4): its finite content, the closure of one maximal order (paper §2.2) and two lattice identities (paper §4.1) | 64 + 8 + 8 products, plus covolumes | 0 |
| Lemma inclusions on Z-bases, in Python and GAP (validation) | 8 basis vectors (Lemma 4.1) + 192 basis products (Lemmas 4.2–4.4) | 0 |
| §5 algebraic-properties test (N = 10⁶) | 1,000,000 per property | n/a (rates) |

Exhaustive verification of all 196,560² ≈ 3.86 × 10¹⁰ pairs is
computationally feasible (~2 hours with 16 CPU cores) but is **not
required**: the symbolic proof in paper §4 establishes closure
independently of sampling.

All 21 transpositions of imaginary basis units produce the same result.
The 21 transpositions form a single orbit under conjugation by the order-21
group of coordinate permutations preserving the oriented Fano triples; each
such permutation is an octonion automorphism fixing s and s̄ and preserving
L, so conjugation carries the proof for σ = (1 2) to every transposition
(paper §4, footnote).

Transpositions are not the only permutations of the imaginary basis units for
which the assembled triple product closes.  The census over all 5,040
permutations of S₇ is exhaustive
(`verify_all_permutations_exact.py`; earlier partial runs
`verify_all_cycles_exact.py`, `consistency_checks.py`; paper §5.6):

| Cycle type of π | Permutations | Close | Rate |
|---|---|---|---|
| transposition | 21 | 21 | 100% |
| 3-cycle | 70 | 35 | 50% |
| (2,2) | 105 | 42 | 40% |
| 4-cycle | 210 | 0 | 0% |
| (3,2) | 420 | 126 | 30% |
| 5-cycle | 504 | 252 | 50% |
| (2,2,2) | 105 | 21 | 20% |
| (3,3) | 280 | 112 | 40% |
| (4,2) | 630 | 294 | 46.7% |
| 6-cycle | 840 | 210 | 25% |
| (3,2,2) | 210 | 105 | 50% |
| (4,3) | 420 | 168 | 40% |
| (5,2) | 504 | 42 | 8.3% |
| 7-cycle | 720 | 336 | 46.7% |
| **total** | **5,039** | **1,764** | **35.0%** |

The identity (the untwisted product) does not close.

### Symbolic proof

Closure is proved symbolically via four lemmas (paper §4):

- **Lemma 4.1**: σ(L) = L (a coordinate permutation preserves D₈⁺).
- **Lemma 4.2**: L · L ⊆ L.
- **Lemma 4.3**: L · σ(Ls̄) ⊆ σ(Ls̄).
- **Lemma 4.4**: σ(Ls) · σ(Ls) ⊆ σ(Ls).

Lemmas 4.2 to 4.4 follow from the closure of one maximal order of integral
octonions, 𝒪 := ½(e₁ + e₂)L (paper §2.2; Wilson's ring A₃).  Three lattice
identities hold: L = (e₁ + e₂)𝒪, σ(Ls̄) = 2𝒪, and σ(Ls) = p𝒪 with
p = e₂ − e₄ + e₅ − e₆ (claims 019, 021).  For purely imaginary x and all
octonions a, b, the alternative, flexible and left Moufang laws give
(xa)(xb) = x((ax)b), so x𝒪 is closed for every purely imaginary x ∈ 𝒪.
Both e₁ + e₂ and p are such elements, which gives Lemmas 4.2 and 4.4;
Lemma 4.3 reads L · 2𝒪 ⊆ 2𝒪, which holds because L ⊆ 𝒪.

The finite content (the closure of 𝒪 on a Z-basis; σ(Ls̄) = 2𝒪 and
σ(Ls) = p𝒪 on eight products each, plus the covolume count) is checked in
exact arithmetic by `verify_structural_proof.py`.  The lemma inclusions
themselves are checked on Z-bases by `symbolic_proof_checks.py`
(Lemmas 4.1, 4.3, 4.4), `verify_structural_proof.py` (Lemmas 4.2 to 4.4)
and `gap_project/tests/test_lemmas.g` (all four).

The standard product fails condition 3 because Ls·Ls ⊄ Ls.  The twist
maps the condition-3 sublattice from Ls to σ(Ls), where closure holds.
The accompanying observation σ(Ls) ≠ Ls (verified by an explicit
witness in `symbolic_proof_checks.py`) is recorded as a remark in
Section 4 of the paper rather than a lemma: it confirms the construction
is nontrivial, but is not used in the closure argument.

One candidate explanation for Lemma 4.4 is **excluded**
(`verify_sigma_Ls_ideal_exclusion.py`; paper Remark 4.5; claim 013):
σ(Ls) = p𝒪 is *not* an ideal of L on either side.  On the script's
Z-bases, only 32 of 64 basis-pair products satisfy L · σ(Ls) ⊆ σ(Ls), and
only 21 of 64 satisfy σ(Ls) · L ⊆ σ(Ls).  (Control: L · Ls ⊆ Ls holds on
24 of 64.)  By contrast σ(Ls̄) = 2𝒪 *is* a two-sided ideal of L, 64/64 on
both sides; Lemma 4.3 is its left-handed half.  So the closure of σ(Ls)
under its own product is not inherited from an ideal property.

### Algebraic properties

The order (Λ, +, ⋆) is nonunital, noncommutative, nonassociative,
not alternative, not flexible, and not power-associative.  The product
is not norm-multiplicative, and it is not a symmetric composition
algebra.  Sampled rates, N = 10⁶ per property, exact integer
arithmetic (`verify_section5_properties.py`; paper §5.1).

### Structure of the ambient algebra and of the order

Below, ⋆ denotes the triple product of the Summary above (the paper's
notation); scripts are in `python_project/src/` unless stated.

**Ideal decomposition.**  As a real algebra, R²⁴ = O³ splits as
D ⊕ T, where D = {(a, a, a) : a ∈ O} is the 8-dimensional diagonal
copy of the octonions (a is an *octonion*, not a real number) and
T = {(p, q, r) ∈ O³ : p + q + r = 0} is 16-dimensional.  D and T are
mutually annihilating two-sided ideals: D ⋆ T = T ⋆ D = 0.
(`verify_ideal_decomposition.py`, `verify_star_algebra_structure.py`.)

**Automorphism group** (`verify_aut_lambda_star.py`, `verify_aut_octonion_crosscheck.py`,
`gap_project/aut_lambda_star.g`; standalone note
`paper/automorphism_group_2026-07-12.pdf`; paper §5.4).

- Aut(Λ, +, ⋆) is **finite of order 36**, with structure C₆ × S₃.
  Center C₆, derived subgroup C₃; −I₂₄ is *not* in it.  Element
  orders: 1 (×1), 2 (×7), 3 (×8), 6 (×20).
- Every automorphism is a blockwise octonion automorphism followed by
  a permutation of the three blocks.  The C₆ is generated by the
  signed permutation fixing e₀: e₁ ↦ −e₅, e₂ ↦ −e₇, e₃ ↦ e₂,
  e₄ ↦ −e₄, e₅ ↦ e₆, e₆ ↦ e₁, e₇ ↦ −e₃.
- Ambient: Aut(R²⁴, ⋆) = G₂ × G₂ × S₃, compact, of dimension 28.
- **Aut(Λ, +, ⋆) ⊊ Co₀**, and this follows from the ambient algebra alone:
  an automorphism of the *real* octonions preserves the
  positive-definite octonion norm, and D ⊥ T, so every
  ⋆-automorphism of R²⁴ is orthogonal.  The inclusion is strict
  because −I₂₄ ∈ Co₀ does not preserve ⋆.  (The trace-form route
  fails: tr(L_u L_v) has signature (3, 21), not the Leech form.)
- Completeness: the order-36 count and the Co₀ containment are
  established by two independent routes in the companion note
  (`paper/automorphism_group_2026-07-12.pdf`): Route 1 uses the
  classification of the automorphisms of the complexification of
  (T, ⋆) (nonenumerative); Route 2 is entirely finite and
  machine-checkable, so nothing load-bearing is left open.  The
  machine-verified lower bound C₆ × S₃ ⊆ Aut and the strictness
  witness are independent of both.
- Octonion-automorphism stabilizers of the Wilson sublattices:
  |Stab(L)| = 1344 = 2³·L₃(2); |Stab(Ls)| = 168 = 2³:(7:3);
  |Stab(Ls̄)| = 12096 = U₃(3).2 = G₂(2)
  (`gap_project/octonion_stabilisers.g`).

**Idempotents: complete classification** (paper §5.2;
`verify_idempotent_classification.py`,
`verify_idempotent_lattice_rescaling.py`).  The ambient algebra
(R²⁴, +, ⋆) has **exactly eight** idempotents.  With ε₁ = (e₀, 0, 0),
ε₂ = (0, e₀, 0), ε₃ = (0, 0, e₀) and ω = ε₁ + ε₂ + ε₃ they are

    0;   ε₁, ε₂, ε₃;   ω/3;   εᵢ − ω/3  (i = 1, 2, 3).

Every block of every one is a multiple of e₀: no imaginary part, no
positive-dimensional family, and the algebra has no identity element.
**None of the seven nonzero idempotents lies in Λ** (because e₀ ∉ L), so
the order has no nonzero idempotent.
Each nonzero idempotent spans a rational ray meeting Λ; the least
positive lattice multiples are 4εᵢ (norm 16), 6(εᵢ − ω/3) (norm 24),
and 4ω (norm 48), satisfying u ⋆ u = 4u, 6u and 12u respectively.  The
rescaling is already visible on the minimal shell: for each of the 84 purely
imaginary roots λ of L, (2λ, 0, 0) ⋆ (2λ, 0, 0) = −8ε₁ ∈ Λ.  The order is
not maximal: L³ is ⋆-closed (the proof of Proposition 4.6 uses only Wilson's
condition (1)), an even integral overlattice of index 2¹²; and ε₁ acts on
either side as the swap of blocks 2 and 3, which preserves Λ, with
ε₁ ⋆ ε₁ = ε₁, so Λ + Zε₁ is a second ⋆-closed overlattice, of index 4, not
contained in L³ and not integral (paper §5.3; key claim 016).

**Square-zero elements** (paper §5.5;
`verify_square_zero_classification.py`).

- In the ambient algebra, u = (x, y, z) satisfies u ⋆ u = 0 iff x, y, z
  are purely imaginary, of equal norm, and sum to zero: a hexagonal
  triple in Im(O) = R⁷.  The solution set is a 12-dimensional cone
  inside T, meeting D only in 0.
- **Λ does contain nonzero square-zero elements, with norms in 12Z**
  (N(u) = 3 N(x) with N(x) ∈ 2Z, and every Λ-norm lies in 4Z).  The
  norm-12 stratum is exactly the 4,032 hexagonal triples of the E₇
  root system formed by the 126 norm-4 vectors of Ls̄ ∩ Im(O); higher
  strata are realized (norm-24 witness:
  `verify_square_zero_norm24_witness.py`).
- The minimal shell sees none: the minimal norm is 8, not a multiple
  of 12.

Min(Λ) itself contains no idempotent and no square-zero vector.

**Span of the image** (paper §5.6;
`verify_product_span_index.py`, `verify_product_span_structure.py`).
S := Z-span{u ⋆ v : u, v ∈ Λ} has index [Λ : S] = 65,536 = 2¹⁶, with
2Λ ⊆ S and Λ/S ≅ (Z/2)¹⁶.  Closure is not surjectivity.

The index depends on the twist (paper §5.6, v8; claims 017, 018): among the closing
permutation twists, every even permutation gives 2¹⁰ and every odd one 2¹⁶;
the reflection in the root (e₀+e₁+⋯+e₆−e₇)/2 of L closes on Λ with index 1,
its products Z-spanning Λ.  L³ is a closed overlattice for every twist by an
isometry of L, so no such twist makes the order maximal.
