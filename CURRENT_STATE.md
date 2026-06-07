# Current State of the Research

Last updated: 2026-05-23

This document is the entry point for anyone continuing this research — human or
AI, with or without prior context.  It summarises what has been established,
what has been ruled out, and what remains open.  All claims are backed by
evidence in the files referenced below.

---

## Research goal

Find a bilinear multiplication · : R²⁴ × R²⁴ → R²⁴ such that the Leech
lattice Λ is closed under · — making (Λ, +, ·) an **order** in the R-algebra
(R²⁴, +, ·).

Since Min(Λ) (the 196,560 minimal vectors at squared norm 8) generates Λ as a
lattice and · is bilinear, it suffices to verify closure on Min(Λ).

The order should be highly symmetric.  The Leech lattice's automorphism group
(the Conway group Co₀) is the largest context for symmetry; any order structure
would ideally respect or interact with this symmetry in a meaningful way.

**Status: FOUND.**  Trial 007 identifies a bilinear product that closes on Λ:
the *transposition-twisted triple octonion product*.  See "The finding" below.
Two earlier candidate families were comprehensively ruled out before this
discovery: the (untwisted) triple-octonion algebra and the triple
Okubo/para-octonion algebra (Petersson isotopes).  See below.

---

## The finding: transposition-twisted triple octonion product (trial 007)

### The construction

Given a transposition σ = (s t) on {1,…,7}, define a new multiplication
·_σ on R⁸ by applying σ to the standard Fano triples.  This gives an
octonion algebra isomorphic to the standard one via σ itself:
σ(x · y) = σ(x) ·_σ σ(y).

The **transposition-twisted triple octonion product** on R²⁴ = O₁ ⊕ O₂ ⊕ O₃
uses ·_σ in all three blocks with Z₃-symmetric cross-block routing:
- Same-block: Oα × Oα → Oα using ·_σ
- Cross-block: Oα × Oβ → Oγ using ·_σ, where {α,β,γ} = {1,2,3}

All 21 transpositions give the same result up to GL(3,F₂) relabelling
(consistency check 1c, check 7).

### Verification status

| Test | Pairs / samples | Failures | File |
|------|-----------------|----------|------|
| Initial (trial 007 base) | 593,412 | 0 | `trial_007_triple_octonion_swap.py` |
| Scaled (4M random) | 4,000,000 | 0 | `trial_007_scaled_test.py` |
| Fast (4M random) | 4,000,000 | 0 | `trial_007_fast.py` |
| Multiprocessor (4M random) | 4,000,000 | 0 | `trial_007_exhaust.py` |
| Symbolic proof (Lemmas 4.1–4.4, paper §4) | 192 basis products (exact) | 0 | `symbolic_proof_checks.py` (Python) + `gap_project/tests/test_lemmas.g` (GAP) |
| §5 algebraic-properties test (N = 10⁶, paper §5) | 1,000,000 samples per property | n/a (rates reported) | `verify_section5_properties.py` |
| Exhaustive (all 38.6B pairs) | **not required** | — | `trial_007_exhaust.py --exhaustive` |

**Zero failures across 12+ million tested pairs**, plus a symbolic proof
of closure (paper §4, Theorem 1.1) that does not depend on sampling.

### How it differs from the ruled-out triple-octonion (trial 001)

The standard (untwisted) triple-octonion product fails on 73.4% of
type3×type3 products, all due to Wilson condition 3 (x+y+z ∉ Ls).
The transposition twist changes 30 of 64 multiplication table entries.
This fixes condition 3 precisely — all other conditions were already
satisfied by the untwisted product (consistency check 8).

### Algebraic properties (consistency check 6)

The order (Λ, +, ·_σ) is:
- **Non-unital**: the ambient identity (1,0..0, 1,0..0, 1,0..0) has N=3 and
  is not in Λ.
- **Non-commutative**: <0.1% of Min(Λ) pairs commute.
- **Not norm-multiplicative**: product norms distribute across
  {16, 32, 48, 64, 80, 96, 112, 128}, with 64 = N(u)·N(v) as the mode
  (~47%).
- **Not alternative, not flexible, not power-associative** (numerically).

The key property is closure, not any classical algebraic identity.

### Evidence files

- `evidence_and_reasoning/research_result.md` — condensed summary
- `evidence_and_reasoning/trial_007_results.md` — detailed trial results
- `evidence_and_reasoning/references/prior_art_orders_on_leech.txt` —
  literature search confirming novelty
- `python_project/src/consistency_checks.py` — all pre-paper checks

---

## What "comprehensively ruled out" means

An algebra is **comprehensively ruled out** when:

1. The base algebra has been tested for Min(Λ) closure and found to fail.
2. The failure mode has been identified precisely (which vector types fail,
   which lattice membership conditions are violated, and what the norm
   distribution of failing products looks like).
3. All natural modifications of the algebra — discrete variants (conjugation,
   sign flips, routing), continuous variants (basis scaling, basis changes via
   lattice automorphisms) — have been tested and shown to fail.
4. The failure has been shown to be **structural**: it persists across all
   tested modifications and is traceable to a specific architectural property
   of the algebra that no parameter change can fix.

This is not a proof of impossibility in the mathematical sense.  It is
systematic computational evidence that rules out a well-defined family of
algebras and identifies the structural reason for failure.

---

## Established foundations (prompts 001–014)

These are the tools and verified constructions that all trials build on.

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
conclusions.  The central research question (whether Λ admits an order) is
confirmed as OPEN by both Wilson [Wilson2009] and Dixon [Dixon2010].

### Octonion algebra

`python_project/src/octonions.py` implements the standard octonion algebra
(Fano-plane rule: e_a · e_{a+1} = e_{a+3}) and Dixon's X-product variant.
Key claim 001 establishes that Wilson and Dixon use the same multiplication
table (up to index relabelling).

### Okubo algebra

`python_project/src/okubo.py` and `okubo_samples.py` implement the Okubo
algebra (Petersson construction from octonions via order-3 automorphism τ).
These have been tested against the Leech lattice in trials 005–006 (see below).

### Test suite

197 tests verify the foundations.  Run with:
```
cd python_project && python3 -m pytest tests/ -v
```

A parallel **GAP / LOOPS** verification suite re-derives the paper's and
companion's key claims independently — same mathematical content, different
language and runtime.  110 checks, all arithmetic exact (rational); the
LOOPS package is used to cross-verify the octonion multiplication as the
unique non-associative Moufang loop of order 16.  Run with:
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
- Cross-block: Oα × Oβ → Oγ where {α,β,γ} = {0,1,2} (Z₃-symmetric routing)

### What was tested (4 trials)

| Trial | What it tests | File | Result |
|---|---|---|---|
| 001 | Base algebra, no modifications | `trial_001_triple_octonion.py` | FAIL: 74.8% of t3×t3 products leave Λ |
| 002 | Per-block scaling (s₁, s₂, s₃) | `trial_002_scaled_triple_octonion.py` | FAIL: norm obstruction (min ≈ 14.6 > 8) + lattice obstruction |
| 003 | Conjugation × sign × routing (1,536 variants, exhaustive) | `trial_003_discrete_variants.py` | FAIL: all 1,536 variants fail; best = baseline |
| 004 | E8 automorphism basis changes (1,128 tested) | `trial_004_basis_automorphisms.py` | FAIL: identity is already optimal; changes make it worse |

### The failure mode

Across all four trials, the failure is:
- **Exclusively type3 × type3 products** (all other type combinations pass)
- **Exclusively Wilson condition 3** (x+y+z ∈ Ls) — conditions 1 and 2 always pass
- **Structurally invariant** under all tested modifications

### Why it fails (structural diagnosis)

Type-3 vectors are the most complex family: all three blocks carry nonzero
components.  The triple-octonion product distributes each block pair's
octonion product into a target block.  For type3 × type3, this creates 9
nonzero contributions across 3 target blocks.  The x+y+z sum (used in
Wilson's condition 3) aggregates information from all three blocks.  The
triple-octonion architecture does not preserve the delicate Ls sublattice
structure that condition 3 requires.

This is not a parameter-tuning problem:
- Conjugation doesn't fix it (all 4³ = 64 patterns fail identically)
- Sign flips have literally zero effect (all 2³ = 8 patterns tie)
- Routing changes don't help (third-block is already optimal)
- Basis changes make it worse (identity is optimal; automorphisms add new failures)
- Scaling can't fix it (norm is bounded below ≈ 14.6 for t3×t3)

### What this rules out

The entire family of algebras with:
- Three copies of any octonion algebra (standard, X-product, or related)
- Block structure O₁ ⊕ O₂ ⊕ O₃ mapped to indices {0–7, 8–15, 16–23}
- Any conjugation pattern on cross-block terms
- Any sign pattern on cross-block terms
- Any of three routing rules (third/left/right block)
- Any E8 automorphism basis change within blocks
- Any per-block scaling

### What this does NOT rule out

- Algebras with a different block decomposition (not 8+8+8)
- Algebras where the cross-block product is not an octonion product at all
  (e.g., Okubo product, or a product defined by different structure constants)
- Algebras where the 24-dimensional product is not built from 8-dimensional
  block products (e.g., a product defined directly on R²⁴)
- Non-bilinear multiplication rules

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

**Note on the original intent:** The user's idea was to take three Okubo
algebras from the same Fano line by cyclically rotating the mediator.
Computational verification showed that no Fano line admits all 3 cyclic
mediators as valid automorphisms — only 1 or 2 work.  The Z₃ orbit above
is the closest valid construction.

### What was tested (2 trials)

| Trial | What it tests | File | Result |
|---|---|---|---|
| 005 | Base algebra + 1,536 discrete variants (exhaustive) | `trial_005_triple_okubo.py` | FAIL: 91.5% leave E8 lattice |
| 006 | E8 automorphism basis changes (600 tested) | `trial_006_triple_okubo_automorphisms.py` | FAIL: 3.8% best; √3 not absorbable |

### The failure mode

This failure is **more severe** than the triple-octonion:

- The triple-octonion preserved E8 membership (conditions 1 & 2) and failed
  only on condition 3 for type3×type3.
- The triple Okubo/para-octonion fails on **ALL type combinations** (not just
  type3×type3) and **91.5% of products leave the E8 lattice entirely**
  (condition 1 fails).

### Why it fails (structural diagnosis)

The Petersson construction with any non-trivial order-3 automorphism τ
introduces structure constants involving cos(2π/3) = −1/2 and sin(2π/3) = √3/2.
When applied to E8 lattice vectors (integer or half-integer coordinates), the
products acquire irrational √3 components that cannot lie in the E8 lattice.

Block 0 (para-octonion) is less affected because x̄ · ȳ has the same structure
constants as the octonion product.  This explains why type1×type1 partially
succeeds (53.5%): type-1 vectors have only one nonzero block, so the
same-block product in block 0 uses rational structure constants.

No conjugation, sign, routing, or E8 automorphism change can eliminate
irrational structure constants — this is a fundamental algebraic obstruction.

### What this rules out

Any 24-dimensional algebra built from Petersson isotopes of the octonion
algebra (para-octonion and/or Okubo algebras) using a non-trivial order-3
automorphism τ, with the same 8+8+8 block structure and cross-block routing.

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

`evidence_and_reasoning/key_claims/NNN_*.txt` — the logical backbone.
Each claim has a status (ESTABLISHED / OPEN / etc.) and detailed evidence.

### Prompt logs

`prompt_logs/NNN_*.txt` — chronological record of every research step.
These are forward-evolving (never altered after commit) per the Manifesto.
The README in that folder has a summary table.

### Trial code

`python_project/src/trial_NNN_*.py` — self-contained, readable top-to-bottom.
Each file has a docstring explaining the algebra, a main function that runs
the trial, and prints detailed results.

### Methodology

`TRIAL_METHODOLOGY.md` — defines the structure and philosophy for all trials.

---

## What has been done since the finding

1. **Symbolic proof complete.**  A symbolic proof of closure on Wilson's
   three sublattice conditions is established in `paper/main.tex`
   (Section 4, "Proof of closure").  The proof rests on four lemmas
   about the interaction between the transposition σ, the maximal
   order L, and the sublattices Ls and Ls̄: Lemma 4.1 (σ(L) = L);
   Lemma 4.2 (L · L ⊆ L); Lemma 4.3 (L · σ(Ls̄) ⊆ σ(Ls̄));
   Lemma 4.4 (σ(Ls) · σ(Ls) ⊆ σ(Ls)).  The finite-case verifications
   (Lemmas 4.2–4.4, 64 Z-basis products each) are executed with exact
   rational arithmetic in `python_project/src/symbolic_proof_checks.py`,
   and recomputed independently in `gap_project/tests/test_lemmas.g`.
   The observation σ(Ls) ≠ Ls is recorded as a remark in the paper
   (rem:Ls-not-closed): non-triviality of the construction, not a step
   of the closure proof.

2. **Formal paper at v5, dated 7 June 2026.**  `paper/main.tex`
   (19 pages) contains abstract, preliminaries, construction,
   symbolic proof, algebraic properties at N = 10⁶ samples, related
   work, conclusion, outlook, and three appendices: explicit basis
   tables with a mod-2-quotient reading and a Lagrangian-polarization
   framing (Appendix A); a historical note synthesising the 1923–1946
   integral-octonion literature from the primary sources, Dickson
   (1923), Kirmse (1924), Mahler (1942), and Coxeter (1946)
   (Appendix B); and the research methodology (Appendix C).  v4 was
   assembled per `evidence_and_reasoning/2026-05-22_plan.md` and
   frozen as `paper/main_2026-05-25.tex` after two full referee-style
   review rounds (paper/reviews/); v5 is a minimal revision
   incorporating `paper/update_brief_2026-05-25_v4.pdf` and is
   frozen as `paper/main_2026-06-07.tex`.  Top-level changelogs:
   v3 → v4 in `paper/v3_to_v4_summary.md`, v4 → v5 in
   `paper/v4_to_v5_changelog.md`.

3. **Computational verification extended.**  Over 12,000,000 random
   minimal-vector pairs tested across multiple test harnesses with
   zero failures.  In addition, all classical algebraic identities
   (commutativity, norm multiplicativity, alternativity, flexibility,
   cube and quartic power-associativity, symmetric composition) have
   been tested on N = 10⁶ samples; rates are tabulated in §5.

4. **σ identified as the Kirmse twist.**  The transposition σ used in
   the paper is the same kind of linear involution that Bruck applied
   (recorded in Coxeter 1946) to repair Kirmse's non-closed candidate
   J₁.  In the paper's coordinate-symmetric placement L = D₈⁺, σ
   leaves L invariant and is an algebra isomorphism of the octonions;
   it does its work by moving Wilson's sublattices Ls and Ls̄.

## What remains open

1. **The structural reason for Lemma 4.4** (σ(Ls) closed under the
   standard octonion product).  The appendix tables certify the
   inclusion, but do not explain *why*.  Appendix Section A.1 lifts
   this to the mod-2 quotient: V := σ(Ls̄)/2L is a two-sided ideal
   and W := σ(Ls)/2L is a subalgebra of the descended product on
   L/2L, and — equipped with the plus-type quadratic form
   q(v + 2L) := N(v)/2 mod 2 — both V and W are Lagrangians (maximal
   totally isotropic subspaces). The decomposition L/2L = V ⊕ W is
   thus a Witt decomposition into a complementary pair of Lagrangians:
   σ produces a polarization of (L/2L, q). The remaining open part is
   the integer-level lift — an intrinsic identification of V and W
   that does not reference σ — and the natural research program:
   which polarizations of L/2L arise from this construction, and
   whether any lifts to a closed bilinear product on Λ itself.

2. **The automorphism group of (Λ, +, ⋆)** and its relationship to
   the Conway group Co₀ = Aut(Λ).  Currently open; a first-step probe
   found −I₂₄ ∈ Co₀ but −I₂₄ ∉ Aut(Λ, +, ⋆), so the inclusion
   Aut(Λ, +, ⋆) ⊆ Co₀ is strict, but the full group has not been
   characterised.

3. **Maximality** of the order in any appropriate sense.

4. **Classification** of the algebra (R²⁴, +, ⋆): it is
   non-associative, non-alternative, non-flexible, and non-unital
   (as an order).  A more natural classification may emerge from a
   ternary rather than binary viewpoint (Elduque); §8 of the paper
   sketches the direction.

5. **Other linear coordinate permutations of R⁸** — beyond simple
   transpositions of imaginary axes — that conjugate the octonion
   multiplication into a triple product closing on Λ.  The
   consecutive-twist footnote in §1 enumerates the result for
   compositions of two transpositions; the general picture is open.

6. **Exhaustive verification** of all 38.6 billion minimal-vector
   pairs.  Not required given the symbolic proof; available as an
   independent computational check.

---

## Repository map

```
leech_alg/
├── MANIFESTO.md                   # Operating rules for AI-assisted research
├── TRIAL_METHODOLOGY.md           # Structure and philosophy for trial files
├── CURRENT_STATE.md               # THIS FILE — the entry point
├── README.md                      # General repository description
│
├── evidence_and_reasoning/
│   ├── key_claims/                # Logical backbone: one file per major claim
│   │   ├── 001–006                # Foundations (Wilson, Dixon, octonions)
│   │   ├── 007                    # Untwisted triple-octonion ruled out
│   │   └── 008                    # Transposition-twisted triple octonion: order on Λ
│   ├── references/                # Central reference registry
│   ├── trial_001_results.md       # Base triple-octonion: t3×t3 fails cond. 3
│   ├── trial_002_results.md       # Per-block scaling: norm + lattice obstruction
│   ├── trial_003_results.md       # Discrete variants: all 1,536 fail
│   ├── trial_004_results.md       # Automorphism basis changes: identity optimal
│   ├── trial_005_results.md       # Triple Okubo/para-octonion: √3 leaves E8
│   └── trial_006_results.md       # Okubo + automorphisms: √3 not absorbable
│
├── python_project/
│   ├── src/
│   │   ├── octonions.py           # Octonion algebra (shared tool)
│   │   ├── e8_wilson.py           # Wilson's E8 lattice (shared tool)
│   │   ├── e8_dixon.py            # Dixon's E8 lattice (shared tool)
│   │   ├── leech_wilson.py        # Wilson's Leech lattice (shared tool)
│   │   ├── leech_dixon.py         # Dixon's Leech lattice (shared tool)
│   │   ├── okubo.py               # Okubo algebra (shared tool, not yet trialed)
│   │   ├── okubo_samples.py       # Okubo examples
│   │   ├── trial_001_*.py         # Trial 001: base triple-octonion
│   │   ├── trial_002_*.py         # Trial 002: per-block scaling
│   │   ├── trial_003_*.py         # Trial 003: discrete variants
│   │   ├── trial_004_*.py         # Trial 004: automorphism basis changes
│   │   ├── trial_005_*.py         # Trial 005: triple Okubo/para-octonion
│   │   ├── trial_006_*.py         # Trial 006: Okubo + E8 automorphisms
│   │   ├── trial_007_*.py         # Trial 007: transposition-twisted (THE FINDING)
│   │   └── consistency_checks.py  # Pre-paper verification (checks 1-10)
│   └── tests/                     # 157 tests verifying foundations
│
├── prompt_logs/                   # Chronological AI interaction log
├── paper/                         # Formal paper: main.tex, references.bib
└── source_documents/              # Primary source PDFs (freely redistributable)
```
