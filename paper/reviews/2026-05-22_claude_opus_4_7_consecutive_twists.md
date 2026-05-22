# Item (1c) — consecutive σ-twists and closure on the Leech lattice

**Date:** 2026-05-22
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 113 (Phase A of the v4 update plan).
**Scripts:** `python_project/src/verify_consecutive_twists_exact.py`
(closure test) and `verify_twist_characterization.py` (the Fano analysis).

---

## The question

The paper's product is built from the octonion product twisted by a single
transposition σ of two imaginary basis indices. What happens under **two
consecutive σ-twists**?

Twisting the product by σ₁ and then by σ₂ is a single twist by the product
permutation π = σ₂σ₁ (the conjugation twist x ·_π y = π⁻¹(π(x)·π(y))). A
product of two transpositions is exactly one of three things:

- the **identity** (σ₁ = σ₂ — "two identical twists");
- a **3-cycle** (σ₁, σ₂ share one index);
- a **(2,2)** double transposition (σ₁, σ₂ disjoint).

## Method, and a corrected misstep

A first attempt used a lemma criterion — ⋆_π closes iff
π(Ls)·π(Ls) ⊆ π(Ls) **and** L·π(Ls̄) ⊆ π(Ls̄). The second lemma is
**sufficient but not necessary**: the three terms of a pairwise block sum are
coupled by membership in Λ, so closure can hold even when L·π(Ls̄) ⊄ π(Ls̄).
A direct cross-check on Λ caught a 3-cycle that closes although the criterion
said it should not. The criterion was discarded.

The reliable test, used here: build a genuine **Z-basis of Wilson's Λ** (24
vectors; verified — covolume 2¹² = [L³:Λ] in L³-coordinates, every basis
vector in Λ). By bilinearity, **⋆_π closes on Λ iff all 24×24 = 576
basis-pair products lie in Λ** — exact, no sampling. Sanity: the identity
gives no closure, every single transposition closes (the paper's theorem) —
both reproduced.

## Result

| π (two consecutive σ-twists) | closes on Λ |
|---|---|
| identity (σ₁ = σ₂) | **no** |
| *single* transposition (one twist — baseline) | all 21 |
| 3-cycle (σ₁, σ₂ share an index) | **35 of 70** |
| (2,2) double transposition (σ₁, σ₂ disjoint) | **42 of 105** |

**Consecutive σ-twists do close on Λ — for a substantial fraction.** Not "all,
of course" (half the 3-cycles), and not "none."

## Characterisation — the Fano plane governs it

- **3-cycles.** For *every one* of the 35 three-element subsets {a,b,c},
  **exactly one** of the two 3-cycles on it closes and the other does not
  (35 subsets, uniformly: 0 with both orientations closing, 35 with exactly
  one, 0 with none). This holds identically for the 7 Fano-line subsets and
  the 28 non-line subsets. So each 3-subset contributes exactly one closing
  3-cycle — 35 in all.
- **(2,2) double transpositions.** A double transposition (a b)(c d) fixes
  three of the seven imaginary indices. If those three **form a Fano line**,
  it **never** closes (0 of 21). If they do **not** form a Fano line,
  **exactly half** close (42 of 84). Total 42 of 105.

So closure of a consecutive-twist product is not arbitrary: it is controlled
by the Fano-plane structure — orientation for 3-cycles, the Fano-line status
of the fixed triple for (2,2)'s.

## Significance

The paper's Remark "rem:unique" establishes that all 21 transposition-twists
are equivalent (one GL(3,𝔽₂) orbit) and notes that "longer-cycle coordinate
permutations are not covered by this argument and may yield inequivalent
products." This item fills that in for the next case: 3-cycle and (2,2)
twists **do** yield products that close on Λ, for definite, Fano-structured
fractions — and the identity (the untwisted product) does not, consistent
with the paper.

Whether to fold this into the paper (a remark extending rem:unique, or an
outlook line) is decision **DP7** — to be taken when §6/Conclusion are
revised.

## Status

Verified facts: the closure counts (identity 0; transpositions 21/21;
3-cycles 35/70; (2,2) 42/105) and the Fano characterisation are exact
(24-vector Z-basis of Λ, 576 basis-pair products, exact arithmetic). The
finer split among the 84 non-line (2,2)'s — why 42 of them and not the other
42 — is not resolved here and is left as an open question.
