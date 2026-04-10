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
├── key_claims/        # One file per major claim or lemma, with argument and references
├── references/        # Central registry of all cited works, organized by topic
├── trial_NNN_*.md     # Results of each computational trial
└── README.md          # This file
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
| 006 | Double-check of Wilson/Dixon; order question confirmed OPEN | ESTABLISHED / OPEN |
| 007 | **Triple-octonion algebra comprehensively ruled out** | ESTABLISHED |

## Trial results

| Trial | Algebra | Verdict | Results file |
|---|---|---|---|
| 001 | Triple-octonion (base) | FAIL: t3×t3 on cond. 3 | [trial_001_results.md](trial_001_results.md) |
| 002 | Triple-octonion + per-block scaling | FAIL: norm + lattice obstruction | [trial_002_results.md](trial_002_results.md) |
| 003 | Triple-octonion + discrete variants (1,536) | FAIL: all variants fail | [trial_003_results.md](trial_003_results.md) |
| 004 | Triple-octonion + E8 automorphism basis changes | FAIL: identity is optimal | [trial_004_results.md](trial_004_results.md) |

## Guidelines

- Each key claim should be documented in its own file under `key_claims/`,
  linking to its references in `references/`.
- Claims must be stated precisely. Informal summaries are permitted as a
  companion to the formal statement, not as a substitute.
- Any step that is uncertain or unverified must be explicitly marked as such.
