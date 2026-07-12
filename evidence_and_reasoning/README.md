# Evidence and Reasoning

This folder documents the structured reasoning behind the research: key claims,
logical steps, supporting arguments, and the evidence base for each.

## Purpose

This folder serves two functions:

1. **Transparency**: A reader can trace how each conclusion was reached, what
   assumptions were made, and what evidence supports them.
2. **Paper preparation**: When writing the formal paper, this folder provides
   the organized body of reasoning to draw from.

## Structure

```
evidence_and_reasoning/
├── key_claims/             # One file per major claim or lemma, with argument and references
│                           #   (README.md is the ledger, incl. the post-v5 claims 009–015)
├── references/             # Central registry of all cited works, organized by topic
├── trial_NNN_*.md          # Results of each computational trial (001–007)
├── research_statement.md   # The research goal as posed, with its resolution
├── research_result.md      # Condensed summary of the finding (Z₃-symmetric triple-octonion
│                           #   product) and of the structure of the resulting order
├── terminology.md          # Authoritative glossary of terms used across the corpus
├── editorial_standards.md  # The five prose standards; spelling and dash conventions
├── 2026-05-22_plan.md      # Planning record for the v3 → v4 paper revision (complete)
└── README.md               # This file
```

## Key claims (logical backbone)

See [key_claims/README.md](key_claims/README.md) for the full index.

| # | Claim | Status |
|---|---|---|
| 001 | Wilson and Dixon use identical octonion multiplication tables | ESTABLISHED |
| 002 | Wilson's 240-element octonionic construction is a valid E8 lattice | ESTABLISHED |
| 003 | Dixon's two 240-element constructions are valid E8 minimal-vector shells | ESTABLISHED |
| 004 | Wilson's three families give the 196,560 Leech lattice minimal vectors | ESTABLISHED |
| 005 | Dixon's three families give a different embedding of the Leech minimal shell | ESTABLISHED |
| 006 | Double-check of Wilson/Dixon; order question confirmed OPEN (answered by claim 008) | ESTABLISHED |
| 007 | Triple-octonion product over Wilson's own multiplication convention comprehensively ruled out | ESTABLISHED |
| 008 | **Z₃-symmetric triple-octonion product is an order on Λ** (historical file title: "transposition-twisted") | ESTABLISHED |

Claims 009–015 have no dedicated file yet; they are
recorded in the ledger at [key_claims/README.md](key_claims/README.md), each
with the script that verifies it:

| # | Claim | Status |
|---|---|---|
| 009 | Ideal decomposition R²⁴ = D ⊕ T (mutually annihilating two-sided ideals, dim 8 + 16) | ESTABLISHED |
| 010 | **Aut(Λ, +, ⋆) is finite of order 36, C₆ × S₃, and ⊊ Co₀** (was OPEN through v5) | ESTABLISHED (completeness caveat stated) |
| 011 | Idempotents: exactly eight in R²⁴, **none in Λ**; the order has no idempotent | ESTABLISHED |
| 012 | Square-zero elements: **Λ contains 4,032 of them, all of norm 12**; Min(Λ) contains none | ESTABLISHED |
| 012a | The unscoped claim "Λ has no nilpotents" | **REFUTED** (scoped to Min(Λ) it holds) |
| 013 | σ(Ls) is *not* an ideal of L (candidate explanation of Lemma D excluded) | ESTABLISHED |
| 014 | Cycle census: 854 of the 2,365 single-cycle permutations of the imaginary axes close (36.1%) | ESTABLISHED |
| 015 | Span of the image: [Λ : S] = 2¹⁶; structural identification of S/2Λ | ESTABLISHED / OPEN |

## Trial results

| Trial | Algebra | Verdict | Results file |
|---|---|---|---|
| 001 | Triple-octonion (base) | FAIL: t3×t3 on cond. 3 | [trial_001_results.md](trial_001_results.md) |
| 002 | Triple-octonion + per-block scaling | FAIL: norm + lattice obstruction | [trial_002_results.md](trial_002_results.md) |
| 003 | Triple-octonion + discrete variants (1,536) | FAIL: all variants fail | [trial_003_results.md](trial_003_results.md) |
| 004 | Triple-octonion + E8 automorphism basis changes | FAIL: identity is optimal | [trial_004_results.md](trial_004_results.md) |
| 005 | Triple Okubo/para-octonion (base + 1,536 discrete variants) | FAIL: √3 leaves E8 | [trial_005_results.md](trial_005_results.md) |
| 006 | Triple Okubo/para-octonion + E8 automorphism basis changes | FAIL: √3 not absorbable | [trial_006_results.md](trial_006_results.md) |
| 007 | Z₃-symmetric triple-octonion product (the order on Λ) | **PASS**: 12.5M+ random pairs, 0 failures; symbolic proof in [paper/main.tex](../paper/main.tex) | [trial_007_results.md](trial_007_results.md) |

## Guidelines

- Each key claim should be documented in its own file under `key_claims/`,
  linking to its references in `references/`.
- Claims must be stated precisely. Informal summaries are permitted as a
  companion to the formal statement, not as a substitute.
- Any step that is uncertain or unverified must be explicitly marked as such.
