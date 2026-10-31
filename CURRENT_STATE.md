# Current State of the Research

Created: 2026-04-09.  Updated: 2026-10-31.

The result is published: J. Köplinger, *An order on the Leech lattice from a
Z₃-symmetric triple-octonion product*,
[doi:10.13140/RG.2.2.22093.19686](https://doi.org/10.13140/RG.2.2.22093.19686)
(ResearchGate), the released v7 of the paper.  A further revision, v8, is
in progress in `paper/main.tex`; its plans are the dated `2026-10-*_plan.md`
files in `evidence_and_reasoning/`, listed with their status in the repository
map below.  References to the paper (sections, equations, and LaTeX labels in
code font) are to `paper/main.tex` unless stated otherwise; `paper/main.aux`
resolves the labels to numbers.

This document is the technical summary of the research, for a human or an AI,
with or without prior context.  It records what was established, what was ruled
out, and what remains open.  All claims are backed by evidence in the files
referenced below.

---

## Research goal

Find a bilinear multiplication · : R²⁴ × R²⁴ → R²⁴ such that the Leech
lattice Λ is closed under ·, making (Λ, +, ·) an **order** in the R-algebra
(R²⁴, +, ·).

Since Min(Λ) (the 196,560 minimal vectors, of norm 8) generates Λ as a
lattice and · is bilinear, it suffices to verify closure on Min(Λ).

The order should be highly symmetric.  The Leech lattice's automorphism group
(the Conway group Co₀) is the largest context for symmetry; any order structure
would ideally respect or interact with this symmetry in a meaningful way.

**Status: FOUND.**  Trial 007 identifies a bilinear product that closes on Λ:
the *Z₃-symmetric triple-octonion product*.  It is simply an octonion product
in all three blocks; the transposition σ is the device that exhibits its fit
with Wilson's representation.  (The trial files call it the "transposition-twisted triple octonion product".)
See "The finding" below.
Two other candidate families are comprehensively ruled out: the (untwisted)
triple-octonion algebra and the triple Okubo/para-octonion algebra
(Petersson isotopes).  See below.

---

## The finding: the Z₃-symmetric triple-octonion product closes on Λ (trial 007)

### The construction

Given a transposition σ = (i j) on {1,…,7}, the product
·_σ defined by x ·_σ y := σ(σ(x) · σ(y)) is another octonion product
on R⁸, isomorphic to the standard one via σ itself:
σ(x · y) = σ(x) ·_σ σ(y).  In Wilson's representation of the Leech
lattice (which fixes the standard product), the transposition σ is
what exhibits the fit: σ leaves the E₈ lattice L = D₈⁺ invariant but
moves Wilson's sublattices Ls and Ls̄.

The **Z₃-symmetric triple-octonion product** ⋆ on R²⁴ = O₁ ⊕ O₂ ⊕ O₃
uses ·_σ in all three blocks with Z₃-symmetric cross-block routing:
- Same-block: Oα × Oα → Oα using ·_σ
- Cross-block: Oα × Oβ → Oγ using ·_σ, where {α,β,γ} = {1,2,3}

All 21 transpositions close on Λ, verified exactly on a Z-basis
(`verify_consecutive_twists_exact.py`); they are equivalent up to relabeling of
the Fano triples: the 21 transpositions form a single orbit under conjugation
by the order-21 group of coordinate permutations preserving the oriented Fano
triples, and conjugation carries the proof for σ = (1 2) to every
transposition (paper §4, first footnote).  They are not the only twists that close: see
the cycle census below and the closure criterion (paper Proposition 5.3,
`prop:criterion`).

### Verification status

| Test | Pairs / samples | Failures | File |
|------|-----------------|----------|------|
| Initial (trial 007 base) | 900 per transposition | 0 | `trial_007_kirmse_twist.py` |
| Scaled (4M random) | 4,000,000 | 0 | `trial_007_scaled_test.py` |
| Fast (4M random) | 4,000,000 | 0 | `trial_007_fast.py` |
| Multiprocessor (4M random) | 4,000,000 | 0 | `trial_007_exhaust.py` |
| Symbolic proof (Lemmas 4.1–4.4, paper §4) | exact: the proof's finite content (the closure of the maximal order 𝒪 = ½(e₁+e₂)L of paper §2.2 on its 64 basis products; σ(Ls̄) = 2𝒪 and σ(Ls) = p𝒪, with p = e₂ − e₄ + e₅ − e₆, on 8 products each plus covolumes) and the identity (xa)(xb) = x((ax)b) for purely imaginary x, symbolically, in `verify_structural_proof.py`; the four lemmas on Z-bases (Lemma 4.1 on the 8 basis vectors of L, Lemmas 4.2–4.4 on 64 basis products each): all four in GAP, Lemmas 4.1, 4.3, 4.4 in `symbolic_proof_checks.py`, Lemmas 4.2–4.4 in `verify_structural_proof.py` | 0 | `verify_structural_proof.py`, `symbolic_proof_checks.py` (Python) + `gap_project/tests/test_lemmas.g` (GAP) |
| Identity witnesses (paper §5.1, `sec:identities`) | one exact witness per identity | n/a (each identity fails) | `verify_identity_witnesses.py` |
| Identity-test suite (n = 10⁶, paper Appendix D, `app:samples`) | 1,000,000 samples per property | n/a (rates reported) | `verify_section5_properties.py` |
| Exhaustive (all 38.6B pairs) | **not required** | — | `trial_007_exhaust.py --exhaustive` |

**Zero failures across 12+ million tested pairs**, plus a symbolic proof
of closure (paper §4, Theorem 1.1) that does not depend on sampling.

### How it differs from the ruled-out triple-octonion (trial 001)

The standard (untwisted) triple-octonion product fails on roughly three
quarters of type3×type3 products, all due to Wilson condition 3 (x+y+z ∉ Ls).
The rate is a sample estimate: 73.4% on the 5,000-pair sample of
`consistency_checks.py` check 8, 74.8% on the 5,000-pair sample recorded in
`evidence_and_reasoning/trial_001_results.md`.
The transposition twist changes 30 of 64 multiplication table entries.
This fixes condition 3 precisely: all other conditions were already
satisfied by the untwisted product (`consistency_checks.py`, check 8).  By the
closure criterion (paper Proposition 5.3, `prop:criterion`), the untwisted
product fails because Ls is not closed under it (paper Remark 4.5).

### Algebraic properties (paper §5.1, `sec:identities`, and Appendix D, `app:samples`)

The order (Λ, +, ⋆) is:
- **Nonunital**: the ambient algebra (R²⁴, +, ⋆) has *no* multiplicative
  identity at all (paper Remark 5.2, `rem:no-identity`).  In particular
  (e₀, e₀, e₀) is not one: it sends (x', y', z') to
  (x'+y'+z', x'+y'+z', x'+y'+z').
- **Nonassociative, noncommutative, not alternative, not flexible, not
  power-associative (cube and quartic), not norm-multiplicative, and not a
  symmetric composition algebra**: each failure has a witness in the
  ambient algebra that is checked by hand, and holds on the order by
  scaling (paper §5.1; `verify_identity_witnesses.py`).  ⋆ does not preserve
  the minimal shell (paper Remark 5.1, `rem:product-norms`).
- **Sampled on minimal vectors** (paper Appendix D; n = 10⁶ per identity;
  `verify_section5_properties.py`, and for the product norms
  `verify_pair_norm_histogram.py`): <0.1% of Min(Λ) pairs commute; product
  norms distribute across {16, 32, 48, 64, 80, 96, 112, 128}, with
  64 = N(u)·N(v) as the mode (~47%).

The key property is closure, not any classical algebraic identity.

### Structure of the algebra (paper §2.2, §4.1, §5.2–§5.6, §6 and Appendix A)

These results are exact and exhaustive.

**Ideal decomposition** (`verify_ideal_decomposition.py`,
`verify_star_algebra_structure.py`).  As a real algebra, R²⁴ = O³ splits as
D ⊕ T, where D = {(a, a, a) : a ∈ O} is the 8-dimensional diagonal copy of
the octonions (a is an *octonion*, not a real number) and
T = {(p, q, r) ∈ O³ : p + q + r = 0} is 16-dimensional.  Both are two-sided
ideals and they annihilate each other: D ⋆ T = T ⋆ D = 0.

**The automorphism group** (`verify_aut_lambda_star.py`,
`verify_aut_octonion_crosscheck.py`, `gap_project/aut_lambda_star.g`;
full derivation in `paper/automorphism_group_2026-07-12.pdf`).

- Aut(Λ, +, ⋆) is **finite of order 36**, with structure **C₆ × S₃**
  (center C₆, derived subgroup C₃).  Element orders: 1 (×1), 2 (×7),
  3 (×8), 6 (×20).  −I₂₄ is not in it.
- Every automorphism is a blockwise octonion automorphism followed by a
  permutation of the three blocks.  The S₃ is the block permutation; the
  C₆ is generated by the signed permutation fixing e₀ and sending
  e₁ ↦ −e₅, e₂ ↦ −e₇, e₃ ↦ e₂, e₄ ↦ −e₄, e₅ ↦ e₆, e₆ ↦ e₁, e₇ ↦ −e₃.
- Ambient: Aut(R²⁴, ⋆) = G₂ × G₂ × S₃, compact, of dimension 28.
- **Aut(Λ, +, ⋆) ⊊ Co₀**, proved from the ambient algebra and before
  the lattice is imposed: an automorphism of the *real* octonions preserves
  the positive-definite norm, and D ⊥ T, so every ⋆-automorphism of R²⁴ is
  orthogonal.  The inclusion is strict because −I₂₄ ∈ Co₀ does not
  preserve ⋆.  (The trace-form route fails: tr(L_u L_v) has signature
  (3, 21), not the Leech form.)
- Stabilizers of the Wilson sublattices in Aut(O, ·_σ), the automorphism
  group of the σ-twisted product (the setting of the automorphism-group note,
  `paper/automorphism_group_2026-07-12.pdf`, its §5.2):
  |Stab(L)| = 1344 = 2³·L₃(2); |Stab(Ls)| = 168 = 2³:(7:3);
  |Stab(Ls̄)| = 12096 = U₃(3).2 = G₂(2) (`gap_project/octonion_stabilisers.g`).
  Under the standard octonion product the orders are 1344, 21 and 1344.
- Completeness: the order-36 count and the Co₀ containment are established
  by two independent routes in that note: Route 1 uses the
  classification of the automorphisms of the complexification of (T, ⋆)
  (nonenumerative); Route 2 is entirely finite and machine-checkable, so
  nothing load-bearing is left open.  The machine-verified lower bound
  C₆ × S₃ ⊆ Aut(Λ, +, ⋆) and the strictness witness are independent of
  both.

**Idempotents: complete classification** (`verify_idempotent_classification.py`,
`verify_idempotent_lattice_rescaling.py`).  (R²⁴, +, ⋆) has exactly eight
idempotents.  With ε₁ = (e₀,0,0), ε₂ = (0,e₀,0), ε₃ = (0,0,e₀) and
ω = ε₁+ε₂+ε₃, they are 0; ε₁, ε₂, ε₃; ω/3; and ε_i − ω/3 (i = 1,2,3).
Every block of every one of them is a multiple of e₀: no imaginary part, no
positive-dimensional family.  **None of the seven nonzero idempotents lies
in Λ** (because e₀ ∉ L), so the order has no nonzero idempotent.  Each
nonzero idempotent spans a rational ray that
meets Λ; the least positive multiples that do are 4ε_i (N = 16),
6(ε_i − ω/3) (N = 24) and 4ω (N = 48), satisfying u ⋆ u = 4u, 6u and 12u
respectively.  The
rescaling is already visible on the minimal shell: for each of the 84 purely
imaginary roots λ of L, (2λ, 0, 0) ⋆ (2λ, 0, 0) = −8ε₁ ∈ Λ.

**The order is not maximal** (paper §5.3; claim 016).  L³ itself is
⋆-closed: the proof of Proposition 4.6 uses only Wilson's condition (1)
(L closed under ·_σ, by Lemmas 4.1 and 4.2), so L³ is an even integral
overlattice of Λ of index [L³ : Λ] = 2¹².  A second ⋆-closed overlattice is
not contained in L³ and not integral: ε₁ = (e₀, 0, 0) acts on either side as
the swap of blocks 2 and 3, which preserves Λ, and ε₁ ⋆ ε₁ = ε₁; since
Λ ∩ Rε₁ = 4Zε₁, the lattice Λ + Zε₁ is ⋆-closed, of index 4 over Λ
(verification record `evidence_and_reasoning/2026-10-12_plan.md`, §4).  The
same holds for ε₂ and ε₃, and Λ + Zε₁ + Zε₂ + Zε₃ is ⋆-closed of index 16
(single-pass machine check, recorded in the same place).

**Square-zero elements** (`verify_square_zero_classification.py`).  In the ambient algebra,
u = (x, y, z) satisfies u ⋆ u = 0 exactly when x, y, z are purely imaginary,
of equal norm, and sum to zero: a hexagonal triple in Im(O) = R⁷, giving a
12-dimensional cone inside T.  **Λ does contain nonzero square-zero
elements**, with norms in 12Z (N(u) = 3·N(x) with N(x) ∈ 2Z, and every
Λ-norm lies in 4Z).  The norm-12 stratum is exactly the 4,032 ordered
hexagonal triples of the E₇ root system formed by the 126 norm-4 vectors of
Ls̄ ∩ Im(O); norm 24 occurs as well (witness:
`verify_square_zero_norm24_witness.py`).  Min(Λ) itself contains no
idempotent and no square-zero vector (minimal norm 8, not a multiple
of 12).

**σ(Ls) is not an ideal of L** (paper Remark 4.5, without the counts, which
depend on the bases; claim 013; `verify_sigma_Ls_ideal_exclusion.py`).  A
candidate structural explanation for Lemma 4.4 is excluded: of the 64
basis-pair products, only 32 of L · σ(Ls) and only 21 of σ(Ls) · L land
inside σ(Ls).  By contrast σ(Ls̄) *is* a two-sided ideal of L (64/64 on both
sides), of which Lemma 4.3 is the left-handed half.  Control: L · Ls ⊆ Ls
holds for only 24 of 64.

**Closure of permutation cycles: the full S₇ census**
(`verify_all_permutations_exact.py`, all 5,040 permutations, reusing the
per-permutation test of `verify_consecutive_twists_exact.py`; the
3- to 7-cycles also in `verify_all_cycles_exact.py`; paper §5.6,
`sec:census`).  Writing
x ·_π y := π⁻¹(π(x) · π(y)) for a permutation π of e₁,…,e₇ (for the
involution σ this is the paper's twisted product) and testing closure of
the assembled triple product exactly on the 576 basis-pair products of a
Z-basis of Λ:

| Cycle type of π | Permutations | Close | Rate |
|---|---|---|---|
| Transposition | 21 | 21 | 100% |
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
| **Total** | **5,039** | **1,764** | **35.0%** |

Exhaustive, not sampled.  The identity (the untwisted product) does not
close, since Ls is not closed under the standard product (paper §5.6,
`sec:census`, by Proposition 5.3).  For the 3-cycles the closing half admits a description: one of the
two cyclic orientations per three-element index subset.

**Which twists close; the span index depends on the twist** (paper §5.6,
`sec:census`, Propositions 5.3 and 5.4, `prop:criterion` and `prop:span`;
claims 015, 017, 018; the external review of v8 in
`paper/reviews/2026-10-13_claude_opus_5_5_review_v8/`, second pass
`verify_twist_span_index_census.py` and the note
`paper/reviews/2026-10-14_claude_fable_5_1_index_one_verification.md`).  For
a twist g, an invertible linear map of R⁸ with g(L) = L, write
x ·_g y := g⁻¹(g(x)·g(y)), ⋆_g for the triple product of Definition 3.3
built from ·_g, and S_g for the Z-span of the products u ⋆_g v.  The closure
criterion is the paper's Proposition 5.3, proved there for every twist: ⋆_g
closes on Λ iff Ls̄ and Ls are both closed under ·_g.  The same proof holds
for every Z-bilinear product μ with μ(L, L) ⊆ L (records only: the
closed-form plan's F1).  It is Theorem A of the review (lines 125–169,
stated there for g in W(E₈)), checked against direct closure on all 5,040
permutations (`verify_twist_span_index_census.py`), on all 120 root
reflections and all 17,280 polarizations (the review's package), and on
twists in GL(L) that are not isometries (prompts 315 and 317).  Among the 1,764 closing permutation twists above, every
even permutation gives [Λ : S_π] = 2¹⁰ and every odd one (σ among them) 2¹⁶.
The reflection ρ in the root r = (e₀+e₁+⋯+e₆−e₇)/2 of L (the paper's symbol;
s_r in the review) gives a product that closes on Λ with [Λ : S_ρ] = 1: its
products Z-span Λ.  Of the 120 root
reflections, 56 close, 28 with index 2¹⁶ and 28 with index 1.  The span of
⋆ = ⋆_σ consists of the vectors (w, w, w) + 2u with w ∈ Ls and u ∈ Λ,
proved by hand in the paper (Proposition 5.4; `verify_product_span_proof.py`).
The review's index formula [Λ : S_g] = [Ls : ⟨Ls ·_g Ls⟩]·[Ls̄ : ⟨Ls̄ ·_g Ls̄⟩]²
(its Theorem B, with a proof sketch there; not in the paper) holds on every
closing twist checked.  Index 1 concerns the
span of the products, not maximality: L ·_g L ⊆ L for every such g (g(L) = L
and L · L ⊆ L), so L³ is a ⋆_g-closed overlattice of Λ of index 2¹² for every
twist in this family, and (Λ, ⋆_g) is never a maximal order; for the index-1
reflection, Λ + Z(1_ρ, 0, 0) is ⋆_ρ-closed of index 4 over Λ and, like
Λ + Zε₁ for σ, not contained in L³ (1_ρ = ρ(e₀), the identity of ·_ρ, is not in L).  Above L³ lies an explicit maximal order, of
index 2²⁴·3⁸ over Λ, decomposable and nonintegral (claim 020; two machine
passes).  Independently of the twist, the closed
lattice proj_D(Λ) ⊕ proj_T(Λ) contains Λ with index 3⁸, so Wilson's 3-adic
gluing alone rules out maximality (claim 018, single-pass strengthening).  Open: a structural proof of the index values other
than σ's (the two sublattice factors take only the values 1, 2⁵ and 2⁸ in
the review's census over all of W(E₈)), and the ring structure of the
index-1 orders.

**The closure criterion modulo 2** (paper Appendix A, `app:mod2quotient`;
claim 017; `verify_mod2_quotient.py`, `verify_polarization_lagrangians.py`).
The standard product descends to μ̄ on L̄ = L/2L, and a lattice between 2L
and L is closed exactly when its image is a subalgebra, so ⋆_g closes on Λ
iff g(Ls̄)/2L and g(Ls)/2L are subalgebras of (L̄, μ̄).  Under
q(v + 2L) = N(v)/2 mod 2, Wilson's pair (Ls̄/2L, Ls/2L) is a polarization, an
ordered pair of complementary Lagrangians; σ maps it to (2𝒪/2L, p𝒪/2L), a
two-sided ideal and a subalgebra.  Every polarization is
(g(Ls̄)/2L, g(Ls)/2L) for an isometry g of L (W(E₈) acts on L̄ as the full
orthogonal group of q), and ⋆_g closes exactly when both halves are
subalgebras.  Not in the paper: 4,704 of the 17,280 polarizations close,
and of the 270 Lagrangians 142 are subalgebras, 16 right ideals and 8
two-sided ideals (the review's §3.3 and its package,
`paper/reviews/2026-10-13_claude_opus_5_5_review_v8/`; recomputed at
prompts 317 and 339).  With σ applied in each block, σ(Λ) = Λ(2𝒪, p𝒪) is
Egan's lattice for the pair (2𝒪, p𝒪) (paper §6, display (6.2), `eq:egan`;
claim 023).

**Wilson's sublattices are Kirmse's set and Bruck's order, scaled** (claims
019 and 021; `verify_sigma_sublattice_maximal_order.py`,
`verify_structural_proof.py`; the verification note above, §4).  With
J₁ = Z⁸ + C_K/2, the analogue of Kirmse's set for (2.1), and
O⁽⁴⁾ = Z⁸ + (0 4)C_K/2 the maximal order obtained from it by transposing the
coordinates 0 and 4: Ls̄ = 2J₁, L = J₁·s, Ls = J₁·s², σ(Ls̄) = 2O⁽⁴⁾ and
σ(Ls) = p·O⁽⁴⁾ with p = e₂ − e₄ + e₅ − e₆.  O⁽⁴⁾ is the paper's maximal
order 𝒪 = ½(e₁ + e₂)L (paper §2.2; Wilson's ring A₃), and
L = (e₁ + e₂)𝒪 ⊆ 𝒪.  The paper's proof of Lemmas 4.2 to 4.4 (§4.1) uses
σ(Ls̄) = 2𝒪 and σ(Ls) = p𝒪 (its eq. (4.1)): Lemma 4.3 (L · 2𝒪 ⊆ 2𝒪) and its
right-handed counterpart (2𝒪 · L ⊆ 2𝒪) follow from L ⊆ 𝒪 and the closure
of 𝒪, and Lemmas 4.2 and 4.4 from (xa)(xb) = x((ax)b) for purely imaginary
x, which makes x𝒪 closed for x = e₁ + e₂ and x = p.  The factor 2⁸ in σ's
span index is [2𝒪 : 4𝒪].  Ls̄ = 2J₁ and L = J₁·s are Wilson's (2009,
pp. 2187, 2189: Ls̄ = 2B with B the set of "Kirmse's mistake", and Bs = L);
σ(Ls̄) = 2𝒪 and σ(Ls) = p𝒪 go beyond Wilson.

### Evidence files

- `evidence_and_reasoning/research_result.md`: condensed summary
- `evidence_and_reasoning/trial_007_results.md`: detailed trial results
- `evidence_and_reasoning/references/prior_art_orders_on_leech.txt`:
  literature search confirming novelty
- `python_project/src/consistency_checks.py`: consistency checks of the construction

---

## What "comprehensively ruled out" means

An algebra is **comprehensively ruled out** when:

1. The base algebra has been tested for Min(Λ) closure and found to fail.
2. The failure mode has been identified precisely (which vector types fail,
   which lattice membership conditions are violated, and what the norm
   distribution of failing products looks like).
3. All natural modifications of the algebra (discrete variants: conjugation,
   sign flips, routing; continuous variants: basis scaling, basis changes via
   lattice automorphisms) have been tested and shown to fail.
4. The failure has been shown to be **structural**: it persists across all
   tested modifications and is traceable to a specific architectural property
   of the algebra that no parameter change can fix.

This is not a proof of impossibility in the mathematical sense.  It is
systematic computational evidence that rules out a well-defined family of
algebras and identifies the structural reason for failure.

---

## Established foundations

These are the project's tools and verified constructions.  The trials build
on the octonion algebra and Wilson's E8 and Leech constructions; trials 005
and 006 also on the Okubo algebra.

### Verified lattice constructions

| Construction | File | Vectors | Status |
|---|---|---|---|
| Wilson's E8 (L) | `python_project/src/e8_wilson.py` | 240 roots | ESTABLISHED (key claim 002) |
| Dixon's E8 | `python_project/src/e8_dixon.py` | 240 roots | ESTABLISHED (key claim 003) |
| Wilson's Leech (Λ) | `python_project/src/leech_wilson.py` | 196,560 minimal vectors | ESTABLISHED (key claim 004) |
| Dixon's Leech | `python_project/src/leech_dixon.py` | 196,560 minimal vectors | ESTABLISHED (key claim 005) |

**Wilson's construction is the primary tool for testing.**  Membership in Λ is
checked via `is_in_leech(x, y, z)` in `leech_wilson.py`, which implements
Wilson's three conditions:
1. x, y, z ∈ L (E8 lattice)
2. x+y, x+z, y+z ∈ Ls̄
3. x+y+z ∈ Ls

The 196,560 minimal vectors decompose as:
- Type 1 (720): one nonzero block, two zero blocks
- Type 2 (11,520): two nonzero blocks, one zero block
- Type 3 (184,320): all three blocks nonzero

Key claim 006 records a systematic double-check of all Wilson/Dixon
conclusions.  Neither Wilson [Wilson2009] nor Dixon [Dixon2010] answers the
central research question (whether Λ admits an order); Dixon poses it as open
[Dixon2010, §7].

### Octonion algebra

`python_project/src/octonions.py` implements the standard octonion algebra
(Fano-plane rule: e_a · e_{a+1} = e_{a+3}) and Dixon's X-product variant.
Key claim 001 establishes that Wilson and Dixon use the same multiplication
table (up to index relabeling).

### Okubo algebra

`python_project/src/okubo.py` implements the Okubo algebra by the Petersson
construction from octonions via an order-3 automorphism τ, and on traceless
3×3 Hermitian matrices; `okubo_samples.py` gives worked examples of the
Petersson construction.  Trials 005–006 test the Petersson construction
against the Leech lattice (see below).

### Test suite

198 tests verify the foundations.  Run with:
```
cd python_project && python3 -m pytest tests/ -v
```

A parallel **GAP / LOOPS** verification suite re-derives the paper's and
companion's key claims independently: same mathematical content, different
language and runtime.  110 checks, all arithmetic exact (rational); the
LOOPS package is used to cross-verify the octonion multiplication as a
nonassociative Moufang loop of order 16, isomorphic to the standard
octonion loop MoufangLoop(16, 3).  Run with:
```
gap -q gap_project/run_all.g
```
See [gap_project/README.md](gap_project/README.md) for the test-by-test
mapping to paper sections and to the Python suite.

---

## Ruled-out algebra: triple-octonion (O₁ ⊕ O₂ ⊕ O₃)

**Key claim 007 records the full ruling-out argument.**

### The algebra

Three copies of the standard octonion algebra on R²⁴:
- O₁ at indices 0–7, O₂ at indices 8–15, O₃ at indices 16–23
- Same-block: Oα × Oα → Oα (standard octonion product)
- Cross-block: Oα × Oβ → Oγ where {α,β,γ} = {1,2,3} (Z₃-symmetric routing)

### What was tested (4 trials)

| Trial | What it tests | File | Result |
|---|---|---|---|
| 001 | Base algebra, no modifications | `trial_001_triple_octonion.py` | FAIL: ~3/4 of t3×t3 products leave Λ |
| 002 | Per-block scaling (s₁, s₂, s₃), products normalized to norm 8 | `trial_002_scaled_triple_octonion.py` | FAIL: no closure on Min(Λ) (lattice obstruction) |
| 003 | Conjugation × sign × routing (1,536 variants, exhaustive) | `trial_003_discrete_variants.py` | FAIL: all 1,536 variants fail; best = baseline |
| 004 | E8 automorphism basis changes (1,128 tested) | `trial_004_basis_automorphisms.py` | FAIL: no sampled change beats the identity; most make it worse |

### The failure mode

In the base algebra (trial 001) and in the best discrete variant (trial 003),
the failure is:
- **Exclusively type3 × type3 products** (all other type combinations pass)
- **Exclusively Wilson condition 3** (x+y+z ∈ Ls; conditions 1 and 2 always pass)
- **Persistent**: type3 × type3 products fail under every tested modification

### Why it fails (structural diagnosis)

Type-3 vectors are the most complex family: all three blocks carry nonzero
components.  The triple-octonion product distributes each block pair's
octonion product into a target block.  For type3 × type3, this creates 9
nonzero contributions across 3 target blocks.  The x+y+z sum (used in
Wilson's condition 3) aggregates information from all three blocks.  The
triple-octonion architecture does not preserve the delicate Ls sublattice
structure that condition 3 requires.

None of the tested modifications fixes it:
- Conjugation doesn't fix it (all 4³ = 64 patterns fail; uniform patterns do best)
- Sign flips leave the failure count unchanged (all 2³ = 8 patterns give the
  same number of failures; the rates differ where products vanish)
- Routing changes don't help (third-block is already optimal)
- The sampled basis changes do not help (none beats the identity; most make it
  worse, and some add new failures); the basis change σ in all three blocks,
  not among them, gives the product of the finding, which closes
- Per-block scaling with each product rescaled to norm 8 fails (trial 002: no
  tested t3×t3 product lands in Min(Λ))

### What this rules out

The base algebra and the variants tested for closure on Λ in trials 001, 003
and 004 (key claim 007), and the per-block scalings of trial 002 for closure of
the minimal shell only:
- Three copies of the standard octonion algebra
- Block structure O₁ ⊕ O₂ ⊕ O₃ mapped to indices {0–7, 8–15, 16–23}
- Any conjugation pattern on cross-block terms
- Any sign pattern on cross-block terms
- Any of three routing rules (third/left/right block)
- The E8 automorphism basis changes sampled in trial 004
- The per-block scalings of trial 002 (products normalized to the minimal shell)

### What this does NOT rule out

- Other octonion products in the blocks: one of them, the σ-twisted product
  (trial 004's form with the basis change σ in all three blocks), closes on Λ
  (the finding above)
- Algebras with a different block decomposition (not 8+8+8)
- Algebras where the cross-block product is not an octonion product at all
  (e.g., Okubo product, or a product defined by different structure constants)
- Algebras where the 24-dimensional product is not built from 8-dimensional
  block products (e.g., a product defined directly on R²⁴)
- Nonbilinear multiplication rules

---

## Ruled-out algebra: triple Okubo/para-octonion (Petersson isotopes)

### The algebra

A Z₃ orbit of Petersson isotopes on R²⁴, using the order-3 octonion
automorphism τ (fixes {e₀, e₁, e₃, e₇}, rotates (e₂, e₅) and (e₄, e₆)
by 2π/3):

- Block 0 (indices 0–7): **para-octonion** x * y = x̄ · ȳ (Petersson with τ⁰)
- Block 1 (indices 8–15): **Okubo_τ** x * y = τ(x̄) · τ²(ȳ)
- Block 2 (indices 16–23): **Okubo_τ²** x * y = τ²(x̄) · τ(ȳ)

Cross-block products use the target block's algebra: Bα × Bβ → Bγ using *_γ.

**Three Okubo algebras from one Fano line:** rotating the mediator of a Fano
line cyclically needs each of its 3 elements to mediate an order-3
automorphism, and on every Fano line each does, for one of the two relative
senses of the two rotations (worked examples in
`python_project/src/okubo_samples.py`).  Trial 005 tests the Z₃ orbit above,
not these three algebras.

### What was tested (2 trials)

| Trial | What it tests | File | Result |
|---|---|---|---|
| 005 | Base algebra + 1,536 discrete variants (exhaustive) | `trial_005_triple_okubo.py` | FAIL: 91.5% of failures leave the E8 lattice |
| 006 | E8 automorphism basis changes (600 tested) | `trial_006_triple_okubo_automorphisms.py` | FAIL: 3.8% best; √3 not absorbable |

### The failure mode

This failure is **more severe** than the triple-octonion:

- The triple-octonion preserved E8 membership (conditions 1 & 2) and failed
  only on condition 3 for type3×type3.
- The triple Okubo/para-octonion fails on **ALL type combinations** (not just
  type3×type3) and **91.5% of the failing products leave the E8 lattice**
  (condition 1 fails).

### Why it fails (structural diagnosis)

The Petersson construction with the order-3 automorphism τ above
introduces structure constants involving cos(2π/3) = −1/2 and sin(2π/3) = √3/2.
When applied to E8 lattice vectors (integer or half-integer coordinates), the
products acquire irrational √3 components that cannot lie in the E8 lattice.

Block 0 (para-octonion) is less affected because x̄ · ȳ has the same structure
constants as the octonion product.  This explains why type1×type1 partially
succeeds (53.5%): type-1 vectors have only one nonzero block, so the
same-block product in block 0 uses rational structure constants.

No conjugation, sign, routing, or E8 automorphism change can eliminate
irrational structure constants.  This is a fundamental algebraic obstruction.

### What this rules out

The triple algebra tested (the para-octonion algebra and the Okubo algebras
from τ and τ² in blocks 0, 1 and 2, with the order-3 automorphism τ above, the
2π/3-rotation), its 1,536 discrete variants (trial 005) and the 600 sampled E8
automorphism basis changes (trial 006), with the same 8+8+8 block structure and
cross-block routing.  The trials do not test an algebra with the para-octonion
product in every block.

### What this does NOT rule out

- Okubo-type algebras constructed without the Petersson 2π/3-rotation
  (i.e., with a different order-3 automorphism that has rational matrix entries)
- Products where the cross-block rule is not "use the target block's algebra"
- Non-block-decomposed products on R²⁴

---

## How to read the evidence

### Trial results

Each trial's results are in `evidence_and_reasoning/trial_NNN_results.md`.
These record: what was tested, how it was tested, what passed, what failed,
the failure mode analysis, and the verdict.

### Key claims

`evidence_and_reasoning/key_claims/NNN_*.txt` and the ledger
`key_claims/README.md` (rows 009–024): the logical backbone.  Each claim has a
status and its evidence.

### Prompt logs

`prompt_logs/NNN_*.txt`: chronological record of every research step.
An entry may be appended to while its work is open; once closed it is
immutable, per the Manifesto (§12) and `DOCUMENT_GENRES.md`.
The README in that folder has a summary table.

### Trial code

`python_project/src/trial_NNN_*.py`: self-contained, readable top-to-bottom.
Each file has a docstring explaining the algebra, a main function that runs
the trial, and prints detailed results.

### Methodology

`TRIAL_METHODOLOGY.md`: defines the structure and philosophy for all trials.

---

## The paper and its supporting material

1. **The symbolic proof.**  A symbolic proof of closure on Wilson's
   three sublattice conditions is established in `paper/main.tex`
   (Section 4, "Proof of closure").  The proof rests on four lemmas
   about the interaction between the transposition σ, the order L
   (a √2-scaled copy of a maximal order), and the sublattices Ls and Ls̄:
   Lemma 4.1 (σ(L) = L);
   Lemma 4.2 (L · L ⊆ L); Lemma 4.3 (L · σ(Ls̄) ⊆ σ(Ls̄));
   Lemma 4.4 (σ(Ls) · σ(Ls) ⊆ σ(Ls)).  Lemma 4.1 is a parity argument;
   Lemmas 4.2 to 4.4 follow from the closure of the maximal order
   𝒪 = ½(e₁ + e₂)L of Section 2.2 (eq. (2.2)), with L = (e₁ + e₂)𝒪, the
   lattice identities σ(Ls̄) = 2𝒪 and σ(Ls) = p𝒪, p = e₂ − e₄ + e₅ − e₆
   (eq. (4.1)), and the identity (xa)(xb) = x((ax)b) for purely imaginary x
   (eq. (4.2)) with its consequence for lattices (eq. (4.3)).  The finite
   content of the proof (the closure of 𝒪 on its 64 basis products; the
   identities (4.1) on eight products each, with covolumes) and the
   identity (4.2), symbolically, are checked in
   `python_project/src/verify_structural_proof.py`.  The lemmas themselves
   are checked on Z-bases in exact rational arithmetic: Lemmas 4.1, 4.3
   and 4.4 in `python_project/src/symbolic_proof_checks.py`, Lemmas 4.2
   to 4.4 in `verify_structural_proof.py`, and all four in
   `gap_project/tests/test_lemmas.g`.  The three coefficient tables of v7's
   Appendix A (192 rows) are recomputed and compared entry by entry with
   the printed ones by `python_project/src/verify_appendix_tables_v7.py`.
   The observation σ(Ls) ≠ Ls is recorded as a remark in the paper
   (Remark 4.5): nontriviality of the construction, not a step of the
   closure proof.

2. **The paper.**  v7, dated 19 July 2026, is the last released version,
   `paper/main_2026-07-19.tex` (24 pages), published on ResearchGate:
   doi:10.13140/RG.2.2.22093.19686.  The working source `paper/main.tex`
   is the v8 draft (22 pages, dated 31 October 2026, marked "v8 (WIP)").
   The draft contains abstract (standalone, with
   keywords and MSC 2020 codes; it names the closure criterion),
   preliminaries, construction, symbolic proof, algebraic properties
   (Section 5, in six subsections: §5.1 identities, each failure
   with a witness checked by hand (`sec:identities`); §5.2 idempotents;
   §5.3 maximality (L³ and Λ + Zε₁ are ⋆-closed overlattices); §5.4 the
   automorphism group; §5.5 square-zero elements; §5.6 other twists,
   closure and span of the image (`sec:census`): the closure criterion,
   Proposition 5.3 (`prop:criterion`), for every twist g with g(L) = L;
   the full S₇ census read through it, with the span index of each
   closing permutation, 2¹⁶ for σ and the odd ones and 2¹⁰ for the even
   ones; the span of ⋆ by hand, Proposition 5.4 (`prop:span`); and a
   reflection twist with span index 1), related work (`sec:related`:
   Egan's construction, display (6.2) (`eq:egan`), with Λ and σ(Λ) among
   its lattices, and the span index 2²⁵ of the Baez–Egan induced product
   φ), conclusion (`sec:conclusion`: a summary that names the criterion,
   and three open questions), outlook, declarations, and four appendices:
   the closure criterion modulo 2, with the answer to which polarizations
   of L/2L close (Appendix A, `app:mod2quotient`); a historical note
   synthesizing the 1923–1946 integral-octonion literature from the
   primary sources, Dickson (1923), Kirmse (1924), Mahler (1942), and
   Coxeter (1946) (Appendix B, `app:history`); the research methodology
   (Appendix C, `sec:methodology`); and the identity frequencies sampled
   on minimal vectors and the norms of products of minimal vectors
   (Appendix D, `app:samples`).  The released versions are
   `paper/main_2026-04-29.tex` (v3), `paper/main_2026-05-25.tex` (v4),
   `paper/main_2026-06-07.tex` (v5), `paper/main_2026-07-18.tex` (v6)
   and `paper/main_2026-07-19.tex` (v7); the changes from one version to
   the next are in `paper/v3_to_v4_summary.md`,
   `paper/v4_to_v5_changelog.md`, `paper/v5_to_v6_changelog.md`,
   `paper/v6_to_v7_changelog.md` and `paper/v7_to_v8_changelog.md` (in
   progress).

3. **Computational verification.**  Over 12,000,000 random
   minimal-vector pairs tested across multiple test harnesses with
   zero failures.  In addition, eight classical algebraic identities
   (commutativity, norm multiplicativity, left and right alternativity, flexibility,
   cube and quartic power-associativity, symmetric composition) have
   been tested on n = 10⁶ samples; rates are tabulated in the paper's
   Appendix D (`app:samples`), and §5.1 (`sec:identities`) gives a witness
   checked by hand for each failure.

4. **σ and the Kirmse twist.**  The transposition σ used in
   the paper is the same kind of linear involution that Bruck applied
   (recorded in Coxeter 1946) to repair Kirmse's nonclosed candidate
   J₁.  In the paper's coordinate-symmetric placement L = D₈⁺, σ
   leaves L invariant (σ(L) = L) and is an algebra isomorphism
   (O, ·) → (O, ·_σ) between two octonion products, not an octonion
   automorphism; it does its work by moving Wilson's sublattices Ls
   and Ls̄.

5. **Structure of the algebra.**  The ideal
   decomposition, the automorphism group, the idempotents, the
   square-zero elements, the exclusion of the σ(Ls)-as-ideal
   explanation, the exhaustive cycle census, the closure criterion (paper
   Proposition 5.3, `prop:criterion`) and the span of ⋆ by hand
   (Proposition 5.4, `prop:span`) are recorded under "Structure of the
   algebra" above.  Scripts:
   `verify_ideal_decomposition.py`, `verify_star_algebra_structure.py`,
   `verify_aut_lambda_star.py`, `verify_aut_octonion_crosscheck.py`,
   `verify_idempotent_classification.py`,
   `verify_idempotent_lattice_rescaling.py`,
   `verify_square_zero_classification.py`,
   `verify_sigma_Ls_ideal_exclusion.py`, `verify_all_permutations_exact.py`,
   `verify_twist_span_index_census.py`, `verify_product_span_proof.py`
   (all in `python_project/src/`); `aut_lambda_star.g`,
   `aut_lambda_star_gens.g`, `octonion_stabilisers.g`,
   `octonion_stabilisers_gens.g` (in `gap_project/`).

6. **Standalone note on the automorphism group.**
   `paper/automorphism_group_2026-07-12.tex` / `.pdf` ("The automorphism
   group of the Z₃-symmetric triple-octonion product", 11 pages, 12 July
   2026) gives the full derivation, the complete enumerations, and an independent
   second computational pass.  It is self-contained and records what
   failed verification (the trace-form route to Co₀) as well as what
   passed.  It is companion material, not part of the paper.

## What remains open

These questions are the ones the paper leaves open (Sections 7 and 8), together
with further detail recorded here.

1. **Where Aut(Λ, +, ⋆) ≅ C₆ × S₃ sits inside Co₀ up to conjugacy**, and
   whether the order-6 octonion automorphism generating the C₆ admits an
   intrinsic description.  The group itself and the strict containment
   Aut(Λ, +, ⋆) ⊊ Co₀ are known (see "Structure of the algebra" above).

2. **A structural proof of the span-index values** other than σ's: 2¹⁶
   for the other odd closing permutations, 2¹⁰ for the even ones, 1 for the
   reflection ρ in the root (e₀+e₁+⋯+e₆−e₇)/2 of L (claim 017; "Structure
   of the algebra" above), and the ring structure of the index-1 orders.
   σ's value is proved in the paper (Proposition 5.4, `prop:span`).  Not
   posed in the paper, which states the other values in §5.6.

3. **Classification** of the algebra (R²⁴, +, ⋆): it is
   nonassociative, nonalternative, nonflexible, non-power-associative,
   noncommutative and nonunital.  Its
   idempotents, square-zero elements, ideal decomposition and
   automorphism group are known (above), but no classification follows
   from them.  A more natural classification may emerge from a
   ternary rather than binary viewpoint (Elduque); §8 of the paper
   sketches the direction.

4. **A reason for the closure pattern of the census.**  All 5,040
   permutations of S₇ have been enumerated exactly (the census above), and
   the closure criterion (paper Proposition 5.3, `prop:criterion`) decides
   closure for every twist: π closes exactly when π(Ls̄) and π(Ls) are
   closed under the standard octonion product.  The 2-cycles close by
   Theorem 1.1, and for the 3-cycles the closing half is described (for
   each three-element subset of imaginary indices, one of the two cyclic
   orientations closes).  Open (paper §7, `sec:conclusion`, its third
   question): a reason for the failure of every 4-cycle and for the exact
   one-half rates.

5. **Exhaustive verification** of all 38.6 billion minimal-vector
   pairs.  Not required given the symbolic proof; available as an
   independent computational check.

---

## Repository map

```
leech_alg/
├── MANIFESTO.md                   # Operating rules for AI-assisted research
├── DOCUMENT_GENRES.md             # Immutable / public record / current state: how each artifact may change
├── TRIAL_METHODOLOGY.md           # Structure and philosophy for trial files
├── tools/                         # lint_docs.py: enforces the document genres mechanically
├── CURRENT_STATE.md               # THIS FILE: technical summary of the research
├── README.md                      # General repository description
├── LICENSE.md, LICENSE-CODE       # CC BY 4.0 for text, evidence and paper; MIT for code
│
├── evidence_and_reasoning/
│   ├── key_claims/                # Logical backbone: one file per major claim
│   │   ├── 001–006                # Foundations (Wilson, Dixon, octonions)
│   │   ├── 007                    # Untwisted triple-octonion ruled out
│   │   ├── 008                    # Z₃-symmetric triple-octonion product: order on Λ
│   │   └── README.md              # The ledger; claims 009–024 are ledger-only rows
│   ├── references/                # Central reference registry
│   ├── 2026-05-22_plan.md         # Plan of the v3 → v4 revision (executed)
│   ├── 2026-10-12_plan.md         # Plan of the v7 → v8 revision (executed)
│   ├── 2026-10-14_plan.md         # Plan: the span-index findings in v8 (executed)
│   ├── 2026-10-20_plan.md         # Plan: the symbolic proof replaces the tables (executed)
│   ├── 2026-10-20_closed_form_plan.md  # Plan: closure criterion, Egan's construction, §5.1 by witnesses (executed)
│   ├── 2026-10-21_plan.md         # Item 1: §5.1 by witnesses, the sampled data to Appendix D (executed)
│   ├── 2026-10-21_criterion_plan.md  # Item 2: the closure criterion in §5.6 (executed)
│   ├── 2026-10-22_appendix_a_plan.md  # Item 3: Appendix A around the criterion mod 2 (executed)
│   ├── 2026-10-22_section6_plan.md  # Item 4: Section 6, the Baez–Egan block (executed)
│   ├── 2026-10-23_span_plan.md  # Item 5: the span of the triple product by hand, §5.6 (executed)
│   ├── 2026-10-24_consistency_plan.md  # Consistency before the fine-comb read (executed)
│   ├── trial_001_results.md       # Base triple-octonion: t3×t3 fails cond. 3
│   ├── trial_002_results.md       # Per-block scaling, normalized products: lattice obstruction
│   ├── trial_003_results.md       # Discrete variants: all 1,536 fail
│   ├── trial_004_results.md       # Sampled automorphism basis changes: none beats the identity
│   ├── trial_005_results.md       # Triple Okubo/para-octonion: √3 leaves E8
│   ├── trial_006_results.md       # Okubo + automorphisms: √3 not absorbable
│   ├── trial_007_results.md       # The σ twist: closure on Λ (THE FINDING)
│   ├── editorial_standards.md     # The five prose standards (US English)
│   ├── terminology.md             # Project-specific and established terms
│   ├── research_statement.md      # The research goal as posed, with its resolution
│   └── research_result.md         # Condensed summary of the finding
│
├── python_project/
│   ├── src/
│   │   ├── octonions.py           # Octonion algebra (shared tool)
│   │   ├── e8_wilson.py           # Wilson's E8 lattice (shared tool)
│   │   ├── e8_dixon.py            # Dixon's E8 lattice (shared tool)
│   │   ├── leech_wilson.py        # Wilson's Leech lattice (shared tool)
│   │   ├── leech_dixon.py         # Dixon's Leech lattice (shared tool)
│   │   ├── okubo.py               # Okubo algebra (shared tool)
│   │   ├── okubo_samples.py       # Okubo examples
│   │   ├── trial_001_*.py         # Trial 001: base triple-octonion
│   │   ├── trial_002_*.py         # Trial 002: per-block scaling
│   │   ├── trial_003_*.py         # Trial 003: discrete variants
│   │   ├── trial_004_*.py         # Trial 004: automorphism basis changes
│   │   ├── trial_005_*.py         # Trial 005: triple Okubo/para-octonion
│   │   ├── trial_006_*.py         # Trial 006: Okubo + E8 automorphisms
│   │   ├── trial_007_*.py         # Trial 007: the σ twist fit (THE FINDING)
│   │   ├── symbolic_proof_checks.py  # Lemmas 4.1, 4.3, 4.4, exact arithmetic
│   │   ├── verify_structural_proof.py  # Lemmas 4.2 to 4.4: the maximal order 𝒪, the lattice identities, the Moufang identity (paper §2.2, §4.1)
│   │   ├── verify_appendix_tables_v7.py  # v7's Appendix A tables recomputed (runs verify_appendix_tables.py on main_2026-07-19.tex)
│   │   ├── verify_aut_*.py        # Aut(Λ, +, ⋆) = C₆ × S₃, order 36
│   │   ├── verify_idempotent_*.py # The eight idempotents; lattice rescaling
│   │   ├── verify_square_zero_classification.py  # Square-zero classification (norm-12 stratum)
│   │   ├── verify_ideal_decomposition.py         # R²⁴ = D ⊕ T
│   │   ├── verify_star_algebra_structure.py      # Ambient algebra structure
│   │   ├── verify_sigma_Ls_ideal_exclusion.py    # σ(Ls) is not an ideal of L
│   │   ├── verify_all_permutations_exact.py      # Closure census over all of S₇ (3- to 7-cycles: verify_all_cycles_exact.py)
│   │   ├── verify_*.py            # Further verifications (span and the twist span-index census,
│   │   │                          #   the seven maximal orders and Wilson's sublattices, mod 2,
│   │   │                          #   historical sources, consecutive twists, the §5.1 witnesses,
│   │   │                          #   Egan's and Dixon's Leech constructions; sampled: the Appendix D identities and norms)
│   │   ├── egan_baez_count.py     # Egan's count 17,280 = 270 × 64 (paper §6)
│   │   ├── probe_aut_lambda_star.py  # Candidate Co₀ elements tested against ⋆ (a probe)
│   │   └── consistency_checks.py  # Consistency checks of the construction (checks 1-3, 5-10)
│   └── tests/                     # 198 tests verifying foundations
│
├── gap_project/                   # Independent GAP / LOOPS re-derivation
│   ├── aut_lambda_star*.g         # Automorphism group, second pass
│   ├── octonion_stabilisers*.g    # Stabilizers of L, Ls, Ls̄ in Aut(O, ·_σ)
│   ├── src/                       # Shared modules (octonions, E8, sublattices, twist, Leech, triple product)
│   ├── tests/                     # 8 test files, one CHECK per claim
│   └── run_all.g                  # 110 checks, exact rational arithmetic
│
├── prompt_logs/                   # Chronological AI interaction log
├── paper/                         # main.tex (v8, in progress); v7 of 2026-07-19, released (24 pp.);
│                                  #   bibliography inline via \bibitem;
│                                  #   released versions main_2026-{04-29,05-25,06-07,07-18,07-19}.tex
│                                  #   and main_2026-07-19.pdf; changelogs v3→v4 … v7→v8 (in progress);
│                                  #   companion and standalone notes, reviews/
└── source_documents/              # Primary source PDFs (freely redistributable)
```
