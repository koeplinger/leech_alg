# Key Claims

This folder contains one file per major claim, lemma, or reasoning step in the research. Together these files form the logical backbone of the investigation.

## Purpose

Breaking the argument into discrete, documented claims allows:
- Independent verification of each step
- Clear identification of which steps are established facts vs. conjectural
- A structured foundation for writing the formal paper

## Naming Convention

Files are named descriptively:

```
NNN_short_claim_name.txt
```

The three-digit prefix preserves logical ordering where a claim depends on earlier ones.

## Status Labels

Each claim file should carry one of the following status labels:

- `ESTABLISHED` — Follows directly from cited references with no novel reasoning required
- `DERIVED` — Follows from established claims by reasoning documented in this project; second-pass verified
- `CONJECTURE` — Stated as a hypothesis; not yet proven
- `UNCERTAIN` — Flagged for review; do not build further claims on this until resolved
- `REFUTED` — Was asserted in this project and has since been shown false; the
  entry is kept, with the date and the counter-evidence, so the record of what
  was believed and when is not lost (Manifesto §12, forward-evolving corrections)

## Contents

| File | Claim | Status |
|---|---|---|
| [001_wilson_dixon_multiplication_tables_compared.txt](001_wilson_dixon_multiplication_tables_compared.txt) | Wilson and Dixon use identical octonion multiplication tables (up to index relabeling) | ESTABLISHED |
| [002_wilson_e8_verified.txt](002_wilson_e8_verified.txt) | Wilson's 240-element octonionic construction is a valid E8 lattice | ESTABLISHED |
| [003_dixon_e8_verified.txt](003_dixon_e8_verified.txt) | Dixon's two 240-element constructions (Ξ^even, A^odd) are valid E8 minimal-vector shells | ESTABLISHED |
| [004_wilson_leech_verified.txt](004_wilson_leech_verified.txt) | Wilson's three families (720+11520+184320=196560) are the Leech lattice minimal vectors | ESTABLISHED |
| [005_dixon_leech_verified.txt](005_dixon_leech_verified.txt) | Dixon's three families (720+11520+184320=196560) form a different embedding of the Leech minimal shell; 17,232 vectors shared with Wilson's | ESTABLISHED |
| [006_double_check_wilson_dixon.txt](006_double_check_wilson_dixon.txt) | Double-check: Wilson fully established; Dixon structural verification only (even/self-dual gaps); order question explicitly OPEN *(that open question was answered by claim 008: Λ does admit an order)* | ESTABLISHED |
| [007_triple_octonion_ruled_out.txt](007_triple_octonion_ruled_out.txt) | Triple-octonion algebra (O₁ ⊕ O₂ ⊕ O₃) comprehensively ruled out as order on Λ — 4 trials, all discrete/continuous variants fail on type3×type3 condition 3 | ESTABLISHED |
| [008_transposition_twist_order.txt](008_transposition_twist_order.txt) | Z₃-symmetric triple-octonion product (·_σ in all three blocks with Z₃ routing; historical file title "transposition-twisted") is an order on Λ — symbolic proof via 4 lemmas, plus 12M+ random-pair verifications | ESTABLISHED |

## Structure of the order

Claims 009–015 have no dedicated `NNN_*.txt` file; the rows below are their
ledger entries.  They are recorded in `paper/main.tex` and in the standalone
note `paper/automorphism_group_2026-07-12.tex` / `.pdf`.  Scripts are in
`python_project/src/` and `gap_project/` unless stated.

| # | Claim | Status | Verified by | Recorded in |
|---|---|---|---|---|
| 009 | **Ideal decomposition.** R²⁴ = O³ splits under ⋆ as D ⊕ T with D = {(a,a,a) : a ∈ O} (dim 8, the diagonal copy of the *octonions*) and T = {(p,q,r) ∈ O³ : p+q+r = 0} (dim 16); D and T are mutually annihilating two-sided ideals (D ⋆ T = T ⋆ D = 0) and are orthogonal. Λ ∩ D = {(a,a,a) : a ∈ Ls}; Λ ∩ T = {(p,q,r) ∈ (Ls̄)³ : p+q+r = 0}; glue (Z/3)⁸ of order 6,561. | ESTABLISHED | `verify_ideal_decomposition.py`, `verify_star_algebra_structure.py`, `verify_aut_lambda_star.py` (the Λ-slices and the (Z/3)⁸ glue) | paper §5.3; note §2 (ideal decomposition), §5 (lattice slices) |
| 010 | **Automorphism group.** Aut(Λ, +, ⋆) is finite of order 36, structure C₆ × S₃ (center C₆, derived subgroup C₃, −I₂₄ ∉ Aut; element orders 1×1, 2×7, 3×8, 6×20). Every element is a blockwise octonion automorphism composed with a block permutation. Ambient: Aut(R²⁴, ⋆) = G₂ × G₂ × S₃, compact, dim 28. **Aut(Λ, +, ⋆) ⊊ Co₀**, proved (real octonion automorphisms preserve the positive-definite norm and D ⊥ T, so every ⋆-automorphism of R²⁴ is orthogonal; strict because −I₂₄ does not preserve ⋆). The trace-form route fails: tr(L_u L_v) has signature (3,21). Completeness is established by two independent routes in the companion note: Route 1 via the classification of Aut of the complexification of (T, ⋆) (non-enumerative), Route 2 entirely finite and machine-checkable; the lower bound C₆ × S₃ ⊆ Aut and the strictness witness are independent of both. | ESTABLISHED | `verify_aut_lambda_star.py`, `verify_aut_octonion_crosscheck.py`, `gap_project/aut_lambda_star.g` (+ `_gens.g`) | paper §5.3; note §4 (ambient group), §5 (lattice stabilizer), §6 (containment in Co₀), §7 (what is proved vs. computed) |
| 010a | **Octonion-automorphism stabilizers of the Wilson sublattices.** \|Stab(L)\| = 1344 = 2³·L₃(2); \|Stab(Ls)\| = 168 = 2³:(7:3); \|Stab(Ls̄)\| = 12096 = U₃(3).2 = G₂(2). | ESTABLISHED | `gap_project/octonion_stabilisers.g` (+ `_gens.g`) | note §5.2 (the three octonion-automorphism stabilizers) |
| 011 | **Idempotents — complete classification.** (R²⁴, +, ⋆) has exactly eight idempotents: 0; ε₁, ε₂, ε₃; ω/3; εᵢ − ω/3 (ε₁ = (e₀,0,0) etc., ω = ε₁+ε₂+ε₃). All real (every block a multiple of e₀); no imaginary part, no positive-dimensional family; the algebra has no identity element. **None of the seven nonzero idempotents lies in Λ** (e₀ ∉ L), so the order has no nonzero idempotent. Least positive lattice multiples: 4εᵢ (norm 16), 6(εᵢ − ω/3) (norm 24), 4ω (norm 48), each with u ⋆ u = n u; minimal-shell witness (2λ,0,0) ⋆ (2λ,0,0) = −8ε₁ for each of the 84 purely imaginary roots λ of L. | ESTABLISHED | `verify_idempotent_classification.py`, `verify_idempotent_lattice_rescaling.py` | paper §5.2 (classification and least multiples; the −8ε₁ shell witness is in the script only) |
| 012 | **Square-zero elements.** u = (x,y,z) satisfies u ⋆ u = 0 iff x, y, z are purely imaginary, of equal norm, and sum to zero (a hexagonal triple in Im(O) = R⁷): a 12-dimensional cone inside T. Every square-zero norm of Λ lies in 12Z (N(u) = 3N(x) with N(x) ∈ 2Z, and every Λ-norm in 4Z); the norm-12 stratum is exactly the **4,032** hexagonal triples of the E₇ root system formed by the 126 norm-4 vectors of Ls̄ ∩ Im(O), and higher strata are realized (a norm-24 witness exists). Min(Λ) contains none (minimal norm 8). | ESTABLISHED | `verify_square_zero_classification.py`, `verify_square_zero_norm24_witness.py` | paper §5.4 |
| 012a | **Min(Λ) contains no idempotent and no square-zero vector.** This is the exhaustively checked statement, and it is scoped to the minimal shell; the unscoped "the order has no nilpotents" is false (claim 012). | ESTABLISHED (scoped to Min(Λ)) | `verify_idempotents_min_shell.py`; the scope is forced by `verify_square_zero_classification.py` | paper §5.1 Remark 5.1, §5.4 |
| 013 | **σ(Ls) is NOT an ideal of L** (a candidate explanation of Lemma D, excluded). L · σ(Ls) ⊆ σ(Ls) holds on only 32 of 64 basis-pair products; σ(Ls) · L ⊆ σ(Ls) on only 21 of 64. (Control: L · Ls ⊆ Ls on 24 of 64.) By contrast σ(Ls̄) **is** a two-sided ideal of L (64/64 on both sides); Lemma 4.3 of the paper is its left-handed half. So the closure σ(Ls) · σ(Ls) ⊆ σ(Ls) is not inherited from an ideal property, and the question of an integer-level structural reason for Lemma 4.4 stays open. | ESTABLISHED (exclusion) | `verify_sigma_Ls_ideal_exclusion.py` | paper §7 (Conclusion), the open-question paragraph |
| 014 | **Closure census — which permutation cycles close.** Exhaustive over **all 5,040 permutations** of S₇, testing the triple product assembled from x ·_π y := π⁻¹(π(x)·π(y)) on a Z-basis of Λ (for involutions this is the paper's twisted product; the two conventions coincide exactly there): transposition 21/21 (100%), 3-cycle 35/70 (50%), (2,2) 42/105 (40%), 4-cycle 0/210 (0%), (3,2) 126/420 (30%), 5-cycle 252/504 (50%), (2,2,2) 21/105 (20%), (3,3) 112/280 (40%), (4,2) 294/630 (46.7%), 6-cycle 210/840 (25%), (3,2,2) 105/210 (50%), (4,3) 168/420 (40%), (5,2) 42/504 (8.3%), 7-cycle 336/720 (46.7%); total 1,764/5,039 (35.0%). The identity does not close. Not sampled. | ESTABLISHED | `verify_all_permutations_exact.py` (all of S₇; reuses the per-permutation test of `verify_consecutive_twists_exact.py`); earlier partial runs `verify_all_cycles_exact.py` | paper §5.6 |
| 015 | **Span of the image.** S := Z-span{u ⋆ v : u, v ∈ Λ} has index [Λ : S] = 65,536 = 2¹⁶ in Λ, with 2Λ ⊆ S and Λ/S ≅ (Z/2)¹⁶. Closure is not surjectivity. Identifying S/2Λ inside Λ/2Λ structurally is **open**. | ESTABLISHED / OPEN (structural identification) | `verify_product_span_index.py`, `verify_product_span_structure.py` | paper §5.5 |
