# v4 → v5 revision response

**v4** — frozen 25 May 2026, 19 pages. Source preserved at [paper/main_2026-05-25.tex](main_2026-05-25.tex).

**v5** — dated 7 June 2026, 20 pages. Source at [paper/main.tex](main.tex).

v5 is a minimal revision implementing the four issues recorded in [update_brief_2026-05-25_v4.pdf](update_brief_2026-05-25_v4.pdf). No mathematical content has been added or withdrawn; the four edits clarify three statements that were already in v4 and restore one citation that was dangling in the bibliography.

---

## Summary of changes

1. §2.2 — Petersson 2018 survey cited at the end of the integral-octonion-history sentence.
2. §6 — "cross-block $\mathbb{Z}_3$ symmetry" replaced by "cross-block $S_3$ extending the $\mathbb{Z}_3$ routing".
3. §7 — Automorphism-group bullet extended with the empirical bound $S_3 \subseteq \mathrm{Aut}(\Lambda,+,\star) \subsetneq \mathrm{Co}_0$ and a pointer to `probe_aut_lambda_star.py`.
4. Appendix B — Kirmse paragraph extended with one sentence on his ideal-theoretic treatment and the alternativity context picked up by Mahler.
5. Mechanical: subtitle date/version bump; `\DeclareMathOperator{\Aut}{Aut}` added to the preamble in support of change 3.

---

## Per-change record

### Change 1 — Petersson 2018 citation restored in §2.2

- **What changed.** Appended "; for a modern survey of the integral-octonion history, see Petersson~\cite{Petersson2018}." to the closing sentence of §2.2 that cites Coxeter and Conway–Smith.
- **Where.** §2.2, lines 227–228 of [paper/main.tex](main.tex), final clause of the paragraph.
- **Rationale.** Petersson 2018 was already in the bibliography (added in v4, item C.1) but had no body-text citation, leaving the reference dangling; update brief item 1 flagged this.

### Change 2 — $S_3$ vs $\mathbb{Z}_3$ in §6 Comparison paragraph

- **What changed.** "the cross-block $\mathbb{Z}_3$ symmetry is supplied by the routing" rewritten as "the cross-block $S_3$ extending the $\mathbb{Z}_3$ routing is supplied by the routing".
- **Where.** §6 "Comparison" paragraph, lines 698–701 of [paper/main.tex](main.tex).
- **Rationale.** The triple product is fixed by the full symmetric group $S_3$ on the three Λ-blocks, not just the cyclic $\mathbb{Z}_3$ rotation built into Definition 3.3; update brief item 2 corrects the understatement.

### Change 3 — §7 automorphism-group bullet extended

- **What changed.** The bullet listing the open question on $\mathrm{Aut}(\Lambda,+,\star)$ now records the empirical containment $S_3 \subseteq \mathrm{Aut}(\Lambda,+,\star) \subsetneq \mathrm{Co}_0$ with the reasoning ($S_3$ from Change 2 above; strict inclusion in $\mathrm{Co}_0$ because $\star$ is not $\mathrm{Co}_0$-invariant) and points to `probe_aut_lambda_star.py` for the exact-arithmetic check.
- **Where.** §7 conclusion, lines 806–812 of [paper/main.tex](main.tex).
- **Rationale.** Update brief item 3 asked for the bound to be stated explicitly so readers can see what is and is not known before the question is posed.

### Change 4 — Appendix B Kirmse paragraph

- **What changed.** One sentence appended to the Kirmse 1924 paragraph recording his ideal-theoretic framing, Mahler's continuation of that thread, and the alternativity context.
- **Where.** Appendix B (Historical note), lines 1357–1362 of [paper/main.tex](main.tex).
- **Rationale.** Update brief item 4, in line with the fair-attribution standard adopted during v4: record positive contributions alongside the $J_1$ closure error.

### Change 5 — Mechanical preamble and subtitle

- **What changed.** Subtitle "25 May 2026 (v4)" → "7 June 2026 (v5)"; `\DeclareMathOperator{\Aut}{Aut}` added to the preamble.
- **Where.** Preamble line 70, title block line 75 of [paper/main.tex](main.tex).
- **Rationale.** Version bookkeeping; the `\Aut` operator is needed to typeset Change 3.

---

## Page count and structural impact

Page count: v4 = 19, v5 = 20. The single added page is absorbed by the §7 bullet expansion (Change 3); §2.2, §6, and Appendix B each grow by one short sentence and produce no reflow at the page level. No sections were added, removed, renumbered, or reordered. The bibliography, figure/table inventory, theorem and lemma numbering, and cross-references are unchanged.

---

## Verification

The four content edits were checked against the wording of the update brief by an independent adversarial pass through the project workflow (`petersson-2018-citation-2.2`, `comparison-s3-z3-routing`, `conclusion-aut-bullet-s3`, `kirmse-paragraph-addition`); all four passed verbatim, with the one stylistic note that the brief's em-dashes appear as the LaTeX en-dash convention (`--`) in the source. An independent diff of [paper/main_2026-05-25.tex](main_2026-05-25.tex) against [paper/main.tex](main.tex) returned six hunks, all six classified as user-requested; no unintended differences were detected. The `probe_aut_lambda_star.py` reference added in Change 3 is a working script in the project's probe set.

---

*For background: the previous revision and the brief driving the present one are recorded in [paper/v3_to_v4_summary.md](v3_to_v4_summary.md) and [paper/update_brief_2026-05-25_v4.pdf](update_brief_2026-05-25_v4.pdf).*

— Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger, 2026-06-07
