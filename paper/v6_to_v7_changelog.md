# v6 → v7 revision notes

**v6** — dated 18 July 2026, 23 pages. Frozen at [paper/main_2026-07-18.tex](main_2026-07-18.tex); the last released version.

**v7** — working source [paper/main.tex](main.tex). **In progress, not frozen.** This changelog is self-contained: it records every v6 → v7 change on its own, and the frozen v6 artifacts (including [v5_to_v6_changelog.md](v5_to_v6_changelog.md) and everything earlier) are not touched.

v7 is a publication-preparation revision aimed at an experimental-mathematics venue. The emphasis shifts toward the experimentation and methodology, the experimental discovery of closure and its exact-arithmetic certification, presented as a path that may lead to deeper structural understanding without claiming that understanding here.

---

## 1. Abstract shortened to under 150 words

- **What changed.** The abstract was rewritten from 309 words to 149, within a 150-word cap. It keeps the existence statement, the construction (nine octonion-multiplication blocks under cyclic $\mathbb{Z}_3$ routing), and Wilson's transposition $\sigma$, and it foregrounds the experimental arc: closure "first observed experimentally, on over 12 million random pairs of minimal vectors," followed by "the proof reduces to four lemmas certified in exact integer arithmetic." The mod-$2L$ Witt decomposition and the open polarization question are compressed to one sentence ("we ask which polarizations lift"). The v6 second paragraph surveying the Section 5 structural facts (idempotents, square-zero elements, automorphisms, span and index, the permutation-cycle census) is dropped from the abstract; those results stay in Section 5. The historical paragraph is compressed to a single sentence naming Bruck's repair; the full Dickson–Kirmse–Mahler–Bruck–Coxeter account stays in Appendix B.
- **Where.** Abstract only.
- **Rationale.** The target venue caps the abstract length, and the experimental-methodology framing matches the audience. The abstract remains standalone (no `\cite`, no `\ref`); the v6 inline locator "[J. Algebra 322 (2009), 2186–2190]" for Wilson is dropped, the full citation remaining in the bibliography.
- **Accuracy.** "exact integer arithmetic" is exact: the four-lemma certification runs in doubled-integer coordinates throughout (`symbolic_proof_checks.py`, which doubles all coordinates "to avoid fractions entirely"), and the Appendix A certificate tables are integer 8-tuples. "over 12 million random pairs of minimal vectors" matches the >12,000,000-pair experimental run recorded in Appendix C. "$\sigma$ ... is an octonion-algebra isomorphism" is the same statement as v6's "an algebra isomorphism of the octonions"; the precise form, $\sigma\colon (\mathbb{O}, \cdot_\sigma) \to (\mathbb{O}, \cdot)$, is Definition 3.2. "$\sigma$ leaves the $E_8$ lattice $L$ invariant ... but moves Wilson's sublattices $L\bar s$ and $Ls$" is Lemma 4.1 ($\sigma(L) = L$) together with $\sigma(Ls) \ne Ls$, $\sigma(L\bar s) \ne L\bar s$ (Remark 4.5).

---

*Revision chain: [v3 → v4](v3_to_v4_summary.md), [v4 → v5](v4_to_v5_changelog.md), [v5 → v6](v5_to_v6_changelog.md), v6 → v7 (this file).*

— Claude (Anthropic), at the direction of Jens Köplinger; opened 2026-07-18, in progress
