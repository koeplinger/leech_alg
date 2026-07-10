# v5 → v6 revision notes

**v5** — dated 7 June 2026, 20 pages. Frozen at [paper/main_2026-06-07.tex](main_2026-06-07.tex). The v5 freeze includes the reviewer-response additions recorded as Change 9 of [paper/v4_to_v5_changelog.md](v4_to_v5_changelog.md) (span-defect remarks, exhaustive minimal-shell facts, φ comparison, §1 scope paragraph, §7 stabiliser clause).

**v6** — dated 10 July 2026, 20 pages. Source at [paper/main.tex](main.tex); frozen snapshot [paper/main_2026-07-10.tex](main_2026-07-10.tex).

v6 is a clarity revision of the span-defect material introduced in v5. No mathematical content was added or changed; the computed indices ($2^{16}$ for $\star$, $2^{25}$ for $\varphi$) and all other v5 content are unchanged.

---

## Changes

### Change 1 — §6 φ span comparison: restructured

- **What changed.** The single dash-linked sentence comparing the span defects of $\star$ and $\varphi$ was split into three shorter sentences (dash construction removed). No clarifying footnote is attached: the possible misreading of the comparison (as if the Baez–Egan construction failed to close) is instead resolved at the source, by the index explanation added to Remark 5.4 (Change 3) together with the anatomy paragraph's own opening, which already states that the Baez–Egan closure lives on the 27-dimensional ambient lattice and that $\varphi$ is the projection read-off.
- **Where.** §6, end of the "Anatomy of the Baez–Egan closure" paragraph.
- **Rationale.** As written in v5, the comparison could be misread as claiming the Baez–Egan construction fails to close; the fix grounds the reader in what a span defect is before the comparison arrives.

### Change 2 — (folded into Changes 1 and 3)

An intermediate revision this same day added a complementarity sentence to Remark 5.4 and a clarifying footnote to §6; both were superseded within the v6 cycle by the sharper index explanation of Change 3 and removed as redundant. They do not appear in the v6 freeze.

### Change 3 — Remark 5.4: the index explained

- **What changed.** The redundant closure-versus-reach sentence (added earlier the same day) was removed. In its place, after the computed index, a plain explanation of what the index means: the index counts the cosets of $S$ in Λ; Λ decomposes into $[\Lambda : S]$ disjoint translates of $S$; every product and every integer combination of products lies in $S$ alone; a span-surjective product would have index 1; here $2^{16} - 1$ of every $2^{16}$ cosets are unreachable.
- **Where.** §5, Remark 5.4.
- **Rationale.** "Index" was used without definition; the surjective baseline (index 1) makes the size of the defect legible to a non-specialist.

### Change 4 — Dash-aside sweep (paper-wide)

- **What changed.** All spaced en-dash parenthetical asides (` -- ... -- ` and ` -- ...`) were removed from the prose, roughly 26 constructions across §1, §2.3, §4, Remark 4.9, §6, §7, §8, Appendix A.1, Appendix B, and Appendix C. Each was rewritten as a sentence split, a colon, parentheses, or (where already present) left to the enclosing footnote. Unspaced en-dashes for name joins (Baez--Egan, Cayley--Dickson) and numeric ranges (1923--1946, page ranges) are retained as standard typography.
- **Where.** Paper-wide.
- **Rationale.** Author's style directive: dash-as-punctuation interrupts the sentence flow. The rule is now recorded in `evidence_and_reasoning/editorial_standards.md` (standard 4, adopted 2026-07-10).

### Change 5 — Version and date

- Subtitle "7 June 2026 (v5)" → "10 July 2026 (v6)".

---

## Page count and structural impact

Page count: v5 = 20, v6 = 20. The added footnote and sentence are absorbed without page-level reflow. No sections, remarks, theorem numbering, or cross-references changed; the bibliography is unchanged.

---

*The revision chain: [v3 → v4](v3_to_v4_summary.md), [v4 → v5](v4_to_v5_changelog.md), v5 → v6 (this file).*

— Claude (Anthropic), at the direction of Jens Köplinger, 2026-07-10
