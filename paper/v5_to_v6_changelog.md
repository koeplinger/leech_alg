# v5 → v6 revision notes

**v5** — dated 7 June 2026, 20 pages. Frozen at [paper/main_2026-06-07.tex](main_2026-06-07.tex). The v5 freeze includes the reviewer-response additions recorded as Change 9 of [paper/v4_to_v5_changelog.md](v4_to_v5_changelog.md) (span-defect remarks, exhaustive minimal-shell facts, φ comparison, §1 scope paragraph, §7 stabiliser clause).

**v6** — dated 10 July 2026, 20 pages. Source at [paper/main.tex](main.tex); frozen snapshot [paper/main_2026-07-10.tex](main_2026-07-10.tex).

v6 is a clarity revision of the span-defect material introduced in v5. No mathematical content was added or changed; the computed indices ($2^{16}$ for $\star$, $2^{25}$ for $\varphi$) and all other v5 content are unchanged.

---

## Changes

### Change 1 — §6 φ span comparison: restructured, clarifying footnote added

- **What changed.** The single dash-linked sentence comparing the span defects of $\star$ and $\varphi$ was split into three shorter sentences (dash construction removed). A footnote was added stating explicitly that a span defect is an attribute of the induced product, not a qualification of closure: closure places every product inside Λ, while the span measures how much of Λ the products reach; both $\star$ and $\varphi$ are closed on Λ and differ in how large a sublattice their images span; the Baez–Egan result itself (the doubled Jordan product on the 27-dimensional lattice $\widetilde\Lambda$) is untouched, $\varphi$ being the projection read-off described in §6.
- **Where.** §6, end of the "Anatomy of the Baez–Egan closure" paragraph.
- **Rationale.** As written in v5, the comparison could be misread as claiming the Baez–Egan construction fails to close. The footnote pre-empts that misreading.

### Change 2 — Remark 5.4 (*Span of the image*): complementarity sentence

- **What changed.** After the opening "Closure does not mean surjectivity.", one sentence was added: "Theorem 1.1 places every product $u \star v$ inside Λ; the question addressed here is the complementary one, namely how much of Λ the products reach."
- **Where.** §5, Remark 5.4.
- **Rationale.** Anchors the closure-versus-span distinction at the point where the span defect is first introduced, so the §6 comparison inherits the framing.

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
