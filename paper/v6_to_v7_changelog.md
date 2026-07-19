# v6 → v7 revision notes

**v6** — dated 18 July 2026, 23 pages. Frozen at [paper/main_2026-07-18.tex](main_2026-07-18.tex); the last released version.

**v7** — working source [paper/main.tex](main.tex), dated 19 July 2026. **In progress, not frozen;** until it is frozen it carries the calendar date of its most recent update. This changelog is self-contained: it records every v6 → v7 change on its own, and the frozen v6 artifacts (including [v5_to_v6_changelog.md](v5_to_v6_changelog.md) and everything earlier) are not touched.

v7 is a publication-preparation revision aimed at an experimental-mathematics venue. The emphasis shifts toward the experimentation and methodology, the experimental discovery of closure and its exact-arithmetic certification, presented as a path that may lead to deeper structural understanding without claiming that understanding here.

---

## 1. Abstract shortened to under 150 words

- **What changed.** The abstract was rewritten from 309 words to 149, within a 150-word cap. It keeps the existence statement, the construction (nine octonion-multiplication blocks under cyclic $\mathbb{Z}_3$ routing), and Wilson's transposition $\sigma$, and it foregrounds the experimental arc: closure "first observed experimentally, on over 12 million random pairs of minimal vectors," followed by "the proof reduces to four lemmas certified in exact integer arithmetic." The mod-$2L$ Witt decomposition and the open polarization question are compressed to one sentence ("we ask which polarizations lift"). The v6 second paragraph surveying the Section 5 structural facts (idempotents, square-zero elements, automorphisms, span and index, the permutation-cycle census) is dropped from the abstract; those results stay in Section 5. The historical paragraph is compressed to a single sentence naming Bruck's repair; the full Dickson–Kirmse–Mahler–Bruck–Coxeter account stays in Appendix B.
- **Where.** Abstract only.
- **Rationale.** The target venue caps the abstract length, and the experimental-methodology framing matches the audience. The abstract remains standalone (no `\cite`, no `\ref`); the v6 inline locator "[J. Algebra 322 (2009), 2186–2190]" for Wilson is dropped, the full citation remaining in the bibliography.
- **Accuracy.** "exact integer arithmetic" is exact: the four-lemma certification runs in doubled-integer coordinates throughout (`symbolic_proof_checks.py`, which doubles all coordinates "to avoid fractions entirely"), and the Appendix A certificate tables are integer 8-tuples. "over 12 million random pairs of minimal vectors" matches the >12,000,000-pair experimental run recorded in Appendix C. "$\sigma$ ... is an octonion-algebra isomorphism" is the same statement as v6's "an algebra isomorphism of the octonions"; the precise form, $\sigma\colon (\mathbb{O}, \cdot_\sigma) \to (\mathbb{O}, \cdot)$, is Definition 3.2. "$\sigma$ leaves the $E_8$ lattice $L$ invariant ... but moves Wilson's sublattices $L\bar s$ and $Ls$" is Lemma 4.1 ($\sigma(L) = L$) together with $\sigma(Ls) \ne Ls$, $\sigma(L\bar s) \ne L\bar s$ (Remark 4.5).

---

## 2. Section 6: recent computer-assisted experimental work in adjacent fields

- **What changed.** Three groups working experimentally in nearby fields are added to Section 6 (Related work), before the Corradetti entry, introduced by a one-line framing sentence and ordered chronologically by their *Experimental Mathematics* contribution:
  - **Kirschmer–Nebe (2022):** binary Hermitian lattices over number fields, classified by computer (*Exp. Math.* 31); together with Nebe's extremal even unimodular 72-dimensional lattice of minimum 8, built from the Leech lattice (*Crelle* 673, 2012).
  - **Höhn–Mason (2023):** the 290 fixed-point sublattices of the Leech lattice (*J. Algebra* 448, 2016) and the $N=4$ superconformal structure on most vertex superalgebras of odd unimodular rank-24 lattices (*Exp. Math.* 32); together with the Höhn–Seysen computer-assisted determination of the order of the Monster (arXiv:2508.01037, 2025).
  - **Dotsenko (2025):** identities in nonassociative algebras by computer experiment, with conjectures on nilpotence and the nil property (*Exp. Math.* 34).
- **Where.** Section 6; six new `\bibitem` entries (Dotsenko2025, HohnMason2016, HohnMason2023, HohnSeysen2025, KirschmerNebe2022, Nebe2012). Page count 23 → 24.
- **Rationale.** Placing computer-assisted experimentalists in the body rather than the historical appendix keeps the experimentation-and-methodology approach of nearby fields prominent, matching the venue.
- **Verification.** Every field of every reference was checked against live sources (journal pages, arXiv, author publication lists) before insertion. Two discrepancies in the initial list were corrected: Höhn's *Exp. Math.* paper is the odd-unimodular / $N=4$ vertex-superalgebra paper (not an "even unimodular rank-24" paper), and the Nebe entry was split into the two distinct works it had conflated — the 2012 *Crelle* 72-dimensional lattice and the 2022 Kirschmer–Nebe *Exp. Math.* paper on binary Hermitian lattices. The Dotsenko article is volume 34, the 2025 volume (online-first 2024); the citation and its label use 2025 for volume-year consistency.

---

*Revision chain: [v3 → v4](v3_to_v4_summary.md), [v4 → v5](v4_to_v5_changelog.md), [v5 → v6](v5_to_v6_changelog.md), v6 → v7 (this file).*

— Claude (Anthropic), at the direction of Jens Köplinger; opened 2026-07-19, in progress
