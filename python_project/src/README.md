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
| 007 | `trial_007_fast.py` | 4M-pair random closure check (vectorized) | PASS |
| 007 | `trial_007_scaled_test.py` | 4M-pair scaled test harness | PASS |
| 007 | `trial_007_exhaust.py` | Multiprocessor harness, supports `--exhaustive` (38.6B pairs) | PASS |

All trials use fixed random seeds for reproducibility.  Results are recorded in `../../evidence_and_reasoning/trial_NNN_results.md`.

## Verification scripts

| File | Purpose |
|---|---|
| `symbolic_proof_checks.py` | Exact-rational-arithmetic verification of the four lemmas in Section 4 of the paper (σ(L) = L; L · L ⊆ L; L · σ(Ls̄) ⊆ σ(Ls̄); σ(Ls) · σ(Ls) ⊆ σ(Ls)), plus the accompanying non-triviality check σ(Ls) ≠ Ls recorded as a remark in Section 4.  No floating-point. |
| `consistency_checks.py` | Pre-paper consistency checks 1-3 and 5-10 (check 4, the exhaustive harness, is `trial_007_exhaust.py`): construction well-definedness, isomorphism claim, table differences, generation arguments, transposition independence, untwisted-vs-twisted comparison, cross-reference with Wilson's paper, code correctness. |
| `verify_all_cycles_exact.py` | The cycle census of paper Section 5.5: which twists close.  Exhaustive (not sampled) over the 3- to 7-cycles; the transpositions and the (2,2)-doubles are in `verify_consecutive_twists_exact.py`.  Of the 2,365 permutations of the seven imaginary units that are a **single cycle**, 854 close (36.1%): transpositions 21/21, 3-cycles 35/70, 4-cycles 0/210, 5-cycles 252/504, 6-cycles 210/840, 7-cycles 336/720.  The (2,2)-double transpositions are counted separately (42/105 close, paper Remark 4.9). |
| `verify_all_permutations_exact.py` | The full-S₇ census behind the paper's §5.6 table: all 5,040 permutations of the imaginary units, grouped by cycle type, each tested exactly on the 576 basis-pair products (multiprocessing; reproduces the per-type counts of the single-cycle scripts as a built-in cross-check). 1,764 of the 5,039 non-identity permutations close (35.0%). |
| `verify_sigma_Ls_ideal_exclusion.py` | Excludes the "ideal" explanation of Lemma 4.4: σ(Ls) is **not** an ideal of L (32/64 basis-pair products on the left, 21/64 on the right), while σ(Ls̄) **is** a two-sided ideal (64/64 both sides). Control: L · Ls ⊆ Ls is 24/64. |
| `verify_section5_properties.py` | The Section-5 identity table: commutativity, norm multiplicativity, alternativity, flexibility, cube and quartic power-associativity, symmetric composition, on N = 10⁶ samples each, exact integer arithmetic. The paper's table is the sum of two runs: `python3 verify_section5_properties.py 20260601 600 500000` and `python3 verify_section5_properties.py 20260602 600 500000` (the no-argument default reproduces the earlier 802k single run). |
| `verify_pair_norm_histogram.py` | Sampled histogram of the pair-product norms N(u ⋆ v) on Min(Λ) × Min(Λ) (two-seed 10⁶ convention as the identity table): observed values {16,...,128}, all multiples of 16, mode 64; supports the sampled half of the §5.1 product-norm remark. |
| `verify_consecutive_twists_exact.py` | Which consecutive σ-twists close on Λ (paper Remark 4.9): exact test on a Z-basis. Transpositions 21/21, 3-cycles 35/70, (2,2)-doubles 42/105. |
| `verify_twist_characterization.py` | Characterizes *which* of them close, against the Fano-line structure. |
| `verify_product_span_index.py` | The span index [Λ : S] = 65,536 = 2¹⁶ for S = Z-span(Λ ⋆ Λ), by Hermite normal form (paper Remark 5.6). |
| `verify_product_span_structure.py` | The structure of that span: 2Λ ⊆ S and Λ/S ≅ (Z/2)¹⁶. |
| `verify_phi_span_index.py` | The same index for the Baez–Egan projected Jordan product φ: [Λ : S_φ] = 2²⁵, which does not contain 2Λ (paper Section 6). |
| `verify_mod2_quotient.py` | Appendix Tables A.2/A.3 as statements about the F₂-algebra L/2L: V a two-sided ideal, W a subalgebra. |
| `verify_discovery1_W_subalgebra.py` | The multiplication table of W = σ(Ls)/2L inside L/2L. |
| `verify_discovery2_V_isotropic.py` | V = σ(Ls̄)/2L totally isotropic; V and W as a complementary pair of Lagrangians. |

### Historical verification: paper Appendix B

Each of the four primary papers re-derived from the original, in exact arithmetic.

| File | Purpose |
|---|---|
| `verify_dickson_1923.py` | Dickson 1923, §19: his 8-element basis and the closure of all 64 basis products; his Theorem XV undercount (three maximal orders; the true count is seven). |
| `verify_kirmse_1924.py` | Kirmse 1924: his multiplication table (1), p. 64, and his order J₁, p. 70, which is not closed. |
| `verify_kirmse_1924_forensic.py` | Reconstruction of Kirmse's count of eight, using only his own data: the necessary form (thirty candidates) and the unit-stability narrowing. |
| `verify_mahler_1942.py` | Mahler 1942, on ideals in the Cayley–Dickson algebra. |
| `verify_coxeter_1946.py` | Coxeter 1946, §4, "Kirmse's mistake": Bruck's transposition, checked against Kirmse's own data. |

### Prior art

| File | Purpose |
|---|---|
| `egan_baez_count.py` | Independent check of Egan's count of 17,280 Jordan rings (17,280 = 270 × 64) in the Baez–Egan construction. |
| `verify_egan_lower_bound.py` | Reproduces Egan's lower bound 244,035,421 on Leech copies in (E₈)³ as the orbit-stabilizer ceiling ⌈3!·\|W(E₈)\|³ / \|Co₀\|⌉ (paper §6). Egan's derivation, not his enumeration. |

### The structure of (R²⁴, +, ⋆): paper Section 5

| File | Purpose |
|---|---|
| `verify_ideal_decomposition.py` | The ideal decomposition R²⁴ = D ⊕ T of (R²⁴, ⋆): both are two-sided ideals, D ⋆ T = T ⋆ D = 0, and (D, ⋆) ≅ (O, ⋆ₛ) rescaled by 3. Here D = {(a,a,a) : a ∈ O} is the 8-dimensional diagonal copy of the **octonions** (a is an octonion, not a real number) and T = {(p,q,r) : p+q+r = 0} is 16-dimensional. |
| `verify_idempotent_classification.py` | **Complete** classification of the idempotents of (R²⁴, +, ⋆): there are exactly **eight**: 0; ε₁, ε₂, ε₃; ω/3; and εᵢ − ω/3 (i = 1,2,3), where εᵢ is e₀ in block i and ω = ε₁+ε₂+ε₃. All are real (every block a multiple of e₀); there is no positive-dimensional family, and no identity element. |
| `verify_idempotent_lattice_rescaling.py` | How those eight sit relative to Λ: **none of them lies in Λ** (e₀ ∉ L), so the order has no idempotent. Least positive lattice multiples: 4εᵢ (norm 16), 6(εᵢ − ω/3) (norm 24), 4ω (norm 48), each satisfying u ⋆ u = n u. Minimal-shell witness: for each of the 84 purely imaginary roots λ of L, (2λ,0,0) ⋆ (2λ,0,0) = −8ε₁ ∈ Λ. |
| `verify_square_zero_classification.py` | **Complete** classification of the square-zero elements. u = (x,y,z) has u ⋆ u = 0 iff x, y, z are purely imaginary, of equal norm, and sum to zero (a hexagonal triple in Im(O) = R⁷): a 12-dimensional cone inside T. **Λ does contain nonzero square-zero elements**, with norms in 12Z: the norm-12 stratum is exactly the 4,032 hexagonal triples of the E₇ root system formed by the 126 norm-4 vectors of (Ls̄) ∩ Im(O), and higher strata are realized (norm-24 witness: `verify_square_zero_norm24_witness.py`). The minimal shell sees none: the minimal norm is 8, not a multiple of 12. |
| `verify_square_zero_norm24_witness.py` | Constructs a square-zero vector of Λ of norm 24 (blockwise sum of two blockwise orthogonal hexagonal triples): the norm-12 stratum is not the whole square-zero locus. |
| `verify_idempotents_min_shell.py` | The exhaustive minimal-shell search: Min(Λ) contains no idempotent and no square-zero vector. The shell search supports that scoped statement only; Λ itself does contain square-zero vectors (see `verify_square_zero_classification.py` above). |

### The automorphism group: paper Section 5.3

| File | Purpose |
|---|---|
| `verify_star_algebra_structure.py` | The ambient group Aut(R²⁴, ⋆) = G₂ × G₂ × S₃, compact, of dimension 28. D and T are two-sided ideals that annihilate each other and are **simple**, which is what forces g(D) = D for every automorphism (a dimension count on ideals of T; an annihilator argument alone would be circular, as the module docstring explains); the C₃-Fourier normal form over Q(ζ₃); h_{A,B} = P⊗A + Q⊗B and the block permutations; dim Der = 28; the ⋆-intrinsic trace form (signature (3,21), **not** the Leech form, so the trace-form route to Co₀ fails). |
| `verify_aut_lambda_star.py` | **Aut(Λ, +, ⋆): exact order 36, structure C₆ × S₃.** Center C₆, derived subgroup C₃, −I₂₄ ∉ it; element orders 1(×1), 2(×7), 3(×8), 6(×20). Λ∩D = {(a,a,a) : a ∈ Ls} and Λ∩T = {(p,q,r) ∈ (Ls̄)³, sum zero}, with glue (Z/3)⁸ of order 6,561; complete enumeration of the octonion-automorphism stabilizers of L, Ls, Ls̄ (orders 1344, 168, 12096); the surviving pairs (A,B) by two independent methods; all 36 elements re-verified against Λ and ⋆. Containment in Co₀ is **proved**, not assumed: every octonion automorphism preserves the positive-definite norm and D ⊥ T, so every ⋆-automorphism of R²⁴ is orthogonal (confirmed elementwise on the Leech Gram matrix). Writes the GAP generators. |
| `verify_aut_octonion_crosscheck.py` | Independent cross-checks of those three stabilizer enumerations: re-run from a different generating triple, and an assumption-free brute force over all 645,120 signed permutations for Stab(L). Writes 8×8 GAP generators. |
| `probe_aut_lambda_star.py` | Tests candidate Co₀ elements against ⋆ on all 576 basis-pair products (a sampling probe; the complete answer is `verify_aut_lambda_star.py`). |

The group-theoretic half of the automorphism calculation (`Size`, `StructureDescription`) is done independently in GAP: see [`../../gap_project/README.md`](../../gap_project/README.md). The derivation is written up in `paper/automorphism_group_2026-07-12.tex`.

**Completeness caveat** (also stated in the paper): the order-36 claim and the Co₀ containment rest on one non-enumerative step, the classification of Aut of the complexification of (T, ⋆). The machine-verified lower bound C₆ × S₃ ≤ Aut and the strictness witness do not depend on it.

## Guidelines

- Each module has a single, clear responsibility.
- All non-trivial algorithms reference the mathematical source they implement (paper or textbook section).
- Trial files are readable top-to-bottom by a human without requiring context from prior trials.
- No hardcoded personal information or environment-specific paths.
