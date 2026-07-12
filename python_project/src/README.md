# Source Code

This folder contains the Python source code for the project: lattice constructions, octonion / Okubo algebra implementations, trial experiments, and the exact-arithmetic checks behind the paper's symbolic proof.

## Shared tools (foundations)

These modules are verified by the 197-test suite in `../tests/` and used by all trials.

| Module | Purpose |
|---|---|
| `octonions.py` | Octonion algebra: standard product and Dixon X-product |
| `e8_wilson.py` | Wilson's E8 lattice L (240 roots, membership test, structural element s) |
| `e8_dixon.py` | Dixon's E8 constructions (Ξ^even, A^odd) |
| `leech_wilson.py` | Wilson's Leech lattice Λ (196,560 minimal vectors, membership test) |
| `leech_dixon.py` | Dixon's Leech lattice construction |
| `okubo.py` | Okubo algebra (Petersson construction via order-3 automorphism τ) |
| `okubo_samples.py` | Instructive examples of Petersson τ constructions |

## Trial files

Each trial is a self-contained Python script testing a specific algebra against the Leech lattice.  Run any trial directly:

```
cd python_project/src && python3 trial_NNN_*.py
```

| Trial | File | What it tests | Verdict |
|---|---|---|---|
| 001 | `trial_001_triple_octonion.py` | Base triple-octonion algebra | FAIL (cond. 3) |
| 002 | `trial_002_scaled_triple_octonion.py` | Per-block scaling search | FAIL (norm + lattice obstruction) |
| 003 | `trial_003_discrete_variants.py` | Conjugation × sign × routing (1,536 variants) | FAIL (all variants) |
| 004 | `trial_004_basis_automorphisms.py` | E8 automorphism basis changes | FAIL (identity already optimal) |
| 005 | `trial_005_triple_okubo.py` | Triple Okubo/para-octonion (base + discrete variants) | FAIL (√3 leaves E8) |
| 006 | `trial_006_triple_okubo_automorphisms.py` | Triple Okubo/para-octonion + E8 automorphisms | FAIL (√3 not absorbable) |
| 007 | `trial_007_kirmse_twist.py` | Z₃-symmetric triple-octonion product (the order on Λ) | **PASS** |
| 007 | `trial_007_explanation.py` | Worked-example exposition of the twist | — |
| 007 | `trial_007_fast.py` | 4M-pair random closure check (vectorised) | PASS |
| 007 | `trial_007_scaled_test.py` | 4M-pair scaled test harness | PASS |
| 007 | `trial_007_exhaust.py` | Multiprocessor harness, supports `--exhaustive` (38.6B pairs) | PASS |

All trials use fixed random seeds for reproducibility.  Results are recorded in `../../evidence_and_reasoning/trial_NNN_results.md`.

## Verification scripts

| File | Purpose |
|---|---|
| `symbolic_proof_checks.py` | Exact-rational-arithmetic verification of the four lemmas in Section 4 of the paper (σ(L) = L; L · L ⊆ L; L · σ(Lš) ⊆ σ(Lš); σ(Ls) · σ(Ls) ⊆ σ(Ls)), plus the accompanying non-triviality check σ(Ls) ≠ Ls now recorded as a remark in Section 4.  No floating-point. |
| `consistency_checks.py` | Ten pre-paper consistency checks: construction well-definedness, isomorphism claim, table differences, exhaustive verification harness, generation arguments, claimed algebraic properties, transposition independence, untwisted-vs-twisted comparison, cross-reference with Wilson's paper, code correctness. |
| `verify_ideal_decomposition.py` | The ideal decomposition R²⁴ = D ⊕ T of (R²⁴, ⋆): both are two-sided ideals, D ⋆ T = T ⋆ D = 0, and (D, ⋆) ≅ (O, ⋆ₛ) rescaled by 3. |
| `verify_star_algebra_structure.py` | The ambient algebraic group Aut(R²⁴, ⋆): annihilator characterisation of D and T; the C₃-Fourier normal form over Q(ζ₃); h_{A,B} = P⊗A + Q⊗B and the block permutations; dim Der = 28; the ⋆-intrinsic trace form (signature (3,21), **not** the Leech form). |
| `verify_aut_lambda_star.py` | **Aut(Λ, +, ⋆): exact order 36, structure C₆ × S₃.** Λ∩D and Λ∩T; complete enumeration of the octonion-automorphism stabilisers of L, Ls, Ls̄ (orders 1344, 168, 12096); the surviving pairs (A,B) by two independent methods; all 36 elements verified against Λ, ⋆ and the Leech Gram matrix (so the group lies inside Co₀). Writes the GAP generators. |
| `verify_aut_octonion_crosscheck.py` | Independent cross-checks of those three stabiliser enumerations: re-run from a different generating triple, and an assumption-free brute force over all 645 120 signed permutations for Stab(L). Writes 8×8 GAP generators. |

## Guidelines

- Each module has a single, clear responsibility.
- All non-trivial algorithms reference the mathematical source they implement (paper or textbook section).
- Trial files are readable top-to-bottom by a human without requiring context from prior trials.
- No hardcoded personal information or environment-specific paths.
