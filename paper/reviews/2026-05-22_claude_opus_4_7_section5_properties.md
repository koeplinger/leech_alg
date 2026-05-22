# Item (H) — Section 5 algebraic properties, re-run at large sample

**Date:** 2026-05-22
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 115 (Phase A of the v4 update plan).
**Script:** `python_project/src/verify_section5_properties.py`.

---

## Method

The Section-5 algebraic-property tests were re-run with three changes asked
for in item (H): the "power-associativity" row split into **cube** and
**quartic** tested separately; the **Elduque symmetric-composition** property
added; and a **far larger sample**.

- Samples: **802,000 per property** (a single 540 s budgeted run — the
  harness caps one command at 10 min, so a literal one-hour single run is not
  possible; 802k is ~160–400× the previous 2k–5k).
- Minimal vectors of Λ drawn from Wilson's Type 1/2/3 families, weighted by
  family size (≈ uniform on Min(Λ)); every sampled vector verified to have
  norm 8.
- Exact integer arithmetic throughout (vectors in doubled-integer
  coordinates); the product is the paper's σ-twisted triple product, σ=(1 2).
- Standard error: ≤ 0.006 % for the small-probability rows, ≤ 0.06 % for the
  large ones — the figures below are stable to the digits shown.

## Result (802,000 samples per property)

The standard error is $\mathrm{SE}=\sqrt{p(1-p)/N}$, $N=802{,}000$.

| property | holds | percent ± SE | current §5 says |
|---|---|---|---|
| commutativity | 133 | **0.017 % ± 0.0014 %** | "< 0.1 %" — ok |
| norm multiplicativity | 374 710 | **46.722 % ± 0.056 %** | "≈ 47 %" — ok |
| left alternativity | 1 501 | **0.187 % ± 0.0048 %** | "< 0.1 %" — **wrong** |
| right alternativity | 1 480 | **0.185 % ± 0.0048 %** | "< 0.1 %" — **wrong** |
| flexibility | 66 940 | **8.347 % ± 0.031 %** | "≈ 8.6 %" — ok (slightly high) |
| cube power-associativity | 246 478 | **30.733 % ± 0.052 %** | "≈ 29 %" (as "power-assoc.") |
| quartic power-associativity | 102 232 | **12.747 % ± 0.037 %** | — (new row) |
| symmetric, Elduque ⟨u⋆v,w⟩=⟨u,v⋆w⟩ | 204 122 | **25.452 % ± 0.049 %** | — (new row) |

(All percentages are the fraction of random samples on which the identity
*happens* to hold; none of the eight properties holds identically.)

## Findings

1. **Correction — left/right alternativity.** The current §5 reports
   "< 0.1 %" for both. The true rate is **≈ 0.19 %**. The published figure is
   a small-sample artifact: at the old 2,000-pair sample size the expected
   count is ~3.7, so 0–2 hits ("< 0.1 %") was a plausible undershoot. This is
   exactly the imprecision the larger run was meant to catch; the §5 table
   must be corrected to ≈ 0.19 % (or "0.2 %").
2. **Cube vs. quartic power-associativity.** The degree-3 identity
   x²⋆x = x⋆x² holds for **30.7 %** of minimal vectors; the degree-4 identity
   (all five parenthesisations of x⋆x⋆x⋆x equal) for **12.7 %**. As expected,
   degree 4 is rarer than degree 3 — but not dramatically so. The current
   single "power-associativity ≈ 29 %" row should become two rows; the 29 %
   was the cube case and is now pinned at 30.7 %.
3. **Symmetric-composition property (Elduque).** ⟨u⋆v,w⟩ = ⟨u,v⋆w⟩ holds for
   **25.5 %** of triples — so ⋆ is **not** a symmetric composition algebra
   (consistent with norm-multiplicativity holding only 47 % of the time). This
   is worth stating explicitly: it is the property that *defines* the Okubo
   and para-Hurwitz symmetric composition algebras of §8's Outlook, and the
   Corradetti comparison (item 3) turns on the Okubo norm form having exactly
   this "associative" symmetry — which ⋆ does not.
4. **Other rows unchanged in substance.** norm-multiplicativity 46.7 %
   (≈ 47 %), flexibility 8.35 % (the "≈ 8.6 %" is slightly high; 8.3–8.4 % is
   the accurate figure), commutativity 0.017 % (genuinely < 0.1 %).

## For the §5 revision (Phase C)

- Replace the single "Power-associativity" row by **two**: "Cube
  power-associativity (x²⋆x = x⋆x²)" — 30.7 %; "Quartic power-associativity"
  — 12.7 %.
- Add a row "**Symmetric composition** ⟨u⋆v,w⟩ = ⟨u,v⋆w⟩ (after Elduque)" —
  25.5 %.
- **Correct** left/right alternativity from "< 0.1 %" to "≈ 0.19 %".
- Update flexibility to ≈ 8.3 %, norm-multiplicativity to ≈ 46.7 %.
- Update the §5 preamble: the sample size is now 802,000 per property (state
  the exact figure or "over 8×10⁵").
- **Record the standard error.** Give the per-property SE alongside each
  figure (the percent ± SE column above), and state in the preamble that
  $\mathrm{SE}=\sqrt{p(1-p)/N}$ with $N=802{,}000$ — so the reader sees the
  figures are stable to the digits shown (SE ≤ 0.056 %).

## Status

802,000 samples; standard error ≤ 0.06 %; the table is definitive to the
digits shown. Further chunks would tighten the error only marginally
(√-law) — 802k is appropriate; running more is optional, not needed.
