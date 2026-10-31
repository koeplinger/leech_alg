# Evidence and Reasoning

Created: 2026-04-06.  Updated: 2026-10-31.

This folder documents the structured reasoning behind the research: key claims,
logical steps, supporting arguments, and the evidence base for each.

## Purpose

This folder serves two functions:

1. **Transparency**: A reader can trace how each conclusion was reached, what
   assumptions were made, and what evidence supports them.
2. **The body of reasoning behind the paper**: the organized material the
   formal write-up draws on.

The result is published (doi:10.13140/RG.2.2.22093.19686), as v7 of the paper.
The plans of the v8 revision in progress are the dated `2026-10-*_plan.md`
files listed under Structure below, each with its status.

## Structure

```
evidence_and_reasoning/
├── key_claims/             # One file per major claim or lemma, with argument and references
│                           #   (README.md is the ledger; claims 009–024 are ledger-only rows, with no dedicated file)
├── references/             # Central registry of all cited works, organized by topic
├── trial_NNN_*.md          # Results of each computational trial (001–007)
├── research_statement.md   # The research goal as posed, with its resolution
├── research_result.md      # Condensed summary of the finding (Z₃-symmetric triple-octonion
│                           #   product) and of the structure of the resulting order
├── terminology.md          # Authoritative glossary of terms used across the corpus
├── editorial_standards.md  # The five prose standards; spelling and dash conventions
├── 2026-05-22_plan.md      # Planning record for the v3 → v4 paper revision (executed)
├── 2026-10-12_plan.md      # Planning record for the v7 → v8 paper revision (executed)
├── 2026-10-14_plan.md      # Plan: the span-index findings in v8, for the wrap-up (executed)
├── 2026-10-20_plan.md      # Plan: the symbolic proof replaces the table certificates (executed)
├── 2026-10-20_closed_form_plan.md  # Plan: the closure criterion, Egan's construction, Section 5.1 by witnesses (executed)
├── 2026-10-21_plan.md      # Item 1 of the closed-form plan: Section 5.1 by witnesses, the sampled data to Appendix D (executed)
├── 2026-10-21_criterion_plan.md  # Item 2 of the closed-form plan: the closure criterion in Section 5.6 (executed)
├── 2026-10-22_appendix_a_plan.md  # Item 3 of the closed-form plan: Appendix A around the criterion modulo 2 (executed)
├── 2026-10-22_section6_plan.md  # Item 4 of the closed-form plan: Section 6, the Baez–Egan block (executed)
├── 2026-10-23_span_plan.md  # Item 5 of the closed-form plan: the span of the triple product by hand, Appendix E (executed)
├── 2026-10-24_consistency_plan.md  # Plan: consistency of the paper and its supporting documents before the fine-comb read (executed)
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
| 002 | Triple-octonion + per-block scaling, products normalized to norm 8 | FAIL: lattice obstruction | [trial_002_results.md](trial_002_results.md) |
| 003 | Triple-octonion + discrete variants (1,536) | FAIL: all variants fail | [trial_003_results.md](trial_003_results.md) |
| 004 | Triple-octonion + E8 automorphism basis changes | FAIL: no sampled change beats the identity | [trial_004_results.md](trial_004_results.md) |
| 005 | Triple Okubo/para-octonion (base + 1,536 discrete variants) | FAIL: √3 leaves E8 | [trial_005_results.md](trial_005_results.md) |
| 006 | Triple Okubo/para-octonion + E8 automorphism basis changes | FAIL: √3 not absorbable | [trial_006_results.md](trial_006_results.md) |
| 007 | Z₃-symmetric triple-octonion product (the order on Λ) | **PASS**: 12M+ random pairs, 0 failures; symbolic proof in [paper/main.tex](../paper/main.tex) | [trial_007_results.md](trial_007_results.md) |

## Guidelines

- Each key claim is documented in its own file under `key_claims/` or, for
  the structure claims 009 onward, as a row of the ledger `key_claims/README.md`,
  linking to its references in `references/`.
- Claims must be stated precisely. Informal summaries are permitted as a
  companion to the formal statement, not as a substitute.
- Any step that is uncertain or unverified must be explicitly marked as such.
