# v4 → v5 revision response

**v4** — frozen 25 May 2026, 19 pages. Source preserved at [paper/main_2026-05-25.tex](main_2026-05-25.tex).

**v5** — dated 7 June 2026, 20 pages. Source at [paper/main.tex](main.tex).

v5 is a minimal revision implementing the four issues recorded in [update_brief_2026-05-25_v4.pdf](update_brief_2026-05-25_v4.pdf). No mathematical content has been added or withdrawn; the four edits clarify three statements that were already in v4 and restore one citation that was dangling in the bibliography.

---

## Summary of changes

1. §2.2 — Petersson 2018 survey cited at the end of the integral-octonion-history sentence.
2. §6 — "cross-block $\mathbb{Z}_3$ symmetry" replaced by "cross-block $S_3$ extending the $\mathbb{Z}_3$ routing".
3. §7 — Automorphism-group bullet extended with the empirical bound $S_3 \subseteq \mathrm{Aut}(\Lambda,+,\star) \subsetneq \mathrm{Co}_0$.
4. Appendix B — Kirmse paragraph extended with one sentence on his ideal-theoretic treatment and the alternativity context picked up by Mahler.
5. §8 — Outlook physics paragraph: "$S_3$ symmetry" added to the list of structural themes; Gresnigt–Gourlay–Varma 2023 (*Eur. Phys. J. C* 83 (2023), 747) added to the cite list with sorting by year of publication.
6. Acknowledgments — Holger P. Petersson added in alphabetical order.
7. Mechanical: subtitle date/version bump; `\DeclareMathOperator{\Aut}{Aut}` added to the preamble in support of change 3.

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

- **What changed.** The bullet listing the open question on $\mathrm{Aut}(\Lambda,+,\star)$ now records the empirical containment $S_3 \subseteq \mathrm{Aut}(\Lambda,+,\star) \subsetneq \mathrm{Co}_0$ with the reasoning ($S_3$ from Change 2 above; strict inclusion in $\mathrm{Co}_0$ because $-I_{24} \in \mathrm{Co}_0$ does not preserve $\star$).
- **Where.** §7 conclusion, the bullet on $\mathrm{Aut}(\Lambda,+,\star)$.
- **Rationale.** Update brief item 3 asked for the bound to be stated explicitly so readers can see what is and is not known before the question is posed.

### Change 4 — Appendix B Kirmse paragraph

- **What changed.** One sentence appended to the Kirmse 1924 paragraph recording his ideal-theoretic framing, Mahler's continuation of that thread, and the alternativity context.
- **Where.** Appendix B (Historical note), lines 1357–1362 of [paper/main.tex](main.tex).
- **Rationale.** Update brief item 4, in line with the fair-attribution standard adopted during v4: record positive contributions alongside the $J_1$ closure error.

### Change 5 — §8 outlook: $S_3$ symmetry token and Gresnigt 2023 citation

- **What changed.** The §8 outlook physics paragraph now reads "physics on normed-division-algebra, triality, $S_3$ symmetry, or symmetric-composition-algebra structures" (the inserted item being "$S_3$ symmetry"), and the accompanying citation list now sorts the four physics references by year of publication and adds the new one: Gresnigt–Gourlay–Varma 2023 — *Three generations of colored fermions with $S_3$ family symmetry from Cayley–Dickson sedenions*, Eur. Phys. J. C **83** (2023), 747 (arXiv:2306.13098).
- **Where.** §8 outlook physics paragraph (sole-sentence paragraph after §8.1); new bibitem inserted in alphabetical order between FureyHughes 2025 and Hall 2019 in the bibliography.
- **Rationale.** Change 3 establishes $S_3 \subseteq \mathrm{Aut}(\Lambda,+,\star)$; the Gresnigt–Gourlay–Varma 2023 paper is the closest physics-side reference using $S_3$ as a family-symmetry organising principle on a Cayley–Dickson algebra, and naturally belongs in the outlook list alongside the existing normed-division-algebra and triality references.

### Change 6 — Acknowledgments: Holger P. Petersson added

- **What changed.** The acknowledgments list now reads "Matthew Barley, Geoffrey M. Dixon, **Holger P. Petersson**, Petr Vojtěchovský, and Robert A. Wilson" — Petersson inserted in alphabetical order.
- **Where.** §Acknowledgments of [paper/main.tex](main.tex).
- **Rationale.** Petersson is the modern source whose 2018 lecture established the project's stance on Kirmse and the chain of 1923–1946 contributions; the citation added in Change 1 puts him in the bibliography, and the acknowledgment recognises his off-paper input alongside the other named collaborators.

### Change 7 — Mechanical preamble and subtitle

- **What changed.** Subtitle "25 May 2026 (v4)" → "7 June 2026 (v5)"; `\DeclareMathOperator{\Aut}{Aut}` added to the preamble.
- **Where.** Preamble line 70, title block line 75 of [paper/main.tex](main.tex).
- **Rationale.** Version bookkeeping; the `\Aut` operator is needed to typeset Change 3.

---

## Page count and structural impact

Page count: v4 = 19, v5 = 20. The single added page is absorbed by the §7 bullet expansion (Change 3); §2.2, §6, §8, and Appendix B each grow by one short sentence and produce no reflow at the page level. No sections were added, removed, renumbered, or reordered. The bibliography gains one entry (Gresnigt–Gourlay–Varma 2023 from Change 5, inserted in alphabetical order); the figure/table inventory, theorem and lemma numbering, and cross-references are unchanged.

---

## Verification

The four content edits driven by the update brief were checked against the wording of the brief by an independent adversarial pass through the project workflow (`petersson-2018-citation-2.2`, `comparison-s3-z3-routing`, `conclusion-aut-bullet-s3`, `kirmse-paragraph-addition`); all four passed verbatim, with the one stylistic note that the brief's em-dashes appear as the LaTeX en-dash convention (`--`) in the source. An independent diff of [paper/main_2026-05-25.tex](main_2026-05-25.tex) against [paper/main.tex](main.tex) returned six hunks, all six classified as user-requested; no unintended differences were detected. Change 5 (Gresnigt–Gourlay–Varma 2023 plus the "$S_3$ symmetry" token in the §8 list) is a subsequent same-day addition, web-verified against arXiv:2306.13098 (Eur. Phys. J. C **83** (2023), 747) before insertion.

---

*For background: the previous revision and the brief driving the present one are recorded in [paper/v3_to_v4_summary.md](v3_to_v4_summary.md) and [paper/update_brief_2026-05-25_v4.pdf](update_brief_2026-05-25_v4.pdf).*

— Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger, 2026-06-07
