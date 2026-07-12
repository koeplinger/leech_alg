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
│                           #   (README.md is the ledger; claims 009–015 are ledger-only rows, with no dedicated file)
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

The claims ledger is [key_claims/README.md](key_claims/README.md): one row per
claim, with its status, the script that verifies it, and where it is recorded.
Claims with a dedicated argument have their own `NNN_*.txt` file in that folder.

## Trial results

| Trial | Algebra | Verdict | Results file |
|---|---|---|---|
| 001 | Triple-octonion (base) | FAIL: t3×t3 on cond. 3 | [trial_001_results.md](trial_001_results.md) |
| 002 | Triple-octonion + per-block scaling | FAIL: norm + lattice obstruction | [trial_002_results.md](trial_002_results.md) |
| 003 | Triple-octonion + discrete variants (1,536) | FAIL: all variants fail | [trial_003_results.md](trial_003_results.md) |
| 004 | Triple-octonion + E8 automorphism basis changes | FAIL: identity is optimal | [trial_004_results.md](trial_004_results.md) |
| 005 | Triple Okubo/para-octonion (base + 1,536 discrete variants) | FAIL: √3 leaves E8 | [trial_005_results.md](trial_005_results.md) |
| 006 | Triple Okubo/para-octonion + E8 automorphism basis changes | FAIL: √3 not absorbable | [trial_006_results.md](trial_006_results.md) |
| 007 | Z₃-symmetric triple-octonion product (the order on Λ) | **PASS**: 12M+ random pairs, 0 failures; symbolic proof in [paper/main.tex](../paper/main.tex) | [trial_007_results.md](trial_007_results.md) |

## Guidelines

- Each key claim should be documented in its own file under `key_claims/`,
  linking to its references in `references/`.
- Claims must be stated precisely. Informal summaries are permitted as a
  companion to the formal statement, not as a substitute.
- Any step that is uncertain or unverified must be explicitly marked as such.
