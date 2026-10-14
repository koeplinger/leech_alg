# Paper

The write-up and everything that supports it.  Genres per
[`DOCUMENT_GENRES.md`](../DOCUMENT_GENRES.md): every released version, note,
brief and review is frozen; the working source `main.tex` (v8, in progress)
is edited in place, and each release is frozen as a dated snapshot.

The paper is published: J. Köplinger, *An order on the Leech lattice from a
Z₃-symmetric triple-octonion product*,
[doi:10.13140/RG.2.2.22093.19686](https://doi.org/10.13140/RG.2.2.22093.19686)
(ResearchGate), corresponding to the frozen v7 source `main_2026-07-19.tex`.

## The paper

| File | |
|---|---|
| [`main.tex`](main.tex), [`main.pdf`](main.pdf) | **The paper.** Working source, v8 in progress (opened 12 October 2026). Bibliography inline via `\bibitem`. Build: `pdflatex main.tex` twice. |
| [`main_2026-04-29.tex`](main_2026-04-29.tex) | v3, frozen |
| [`main_2026-05-25.tex`](main_2026-05-25.tex) | v4, frozen |
| [`main_2026-06-07.tex`](main_2026-06-07.tex) | v5, frozen |
| [`main_2026-07-18.tex`](main_2026-07-18.tex) | v6, frozen |
| [`main_2026-07-19.tex`](main_2026-07-19.tex), [`main_2026-07-19.pdf`](main_2026-07-19.pdf) | v7, frozen; the last released version |
| [`v3_to_v4_summary.md`](v3_to_v4_summary.md), [`v4_to_v5_changelog.md`](v4_to_v5_changelog.md), [`v5_to_v6_changelog.md`](v5_to_v6_changelog.md), [`v6_to_v7_changelog.md`](v6_to_v7_changelog.md), [`v7_to_v8_changelog.md`](v7_to_v8_changelog.md) (in progress) | what changed between versions |

## Companion and standalone notes

Each is self-contained and stands apart from the paper; none is in its referee
scope.

| File | |
|---|---|
| [`companion.tex`](companion.tex) | Companion to the paper: worked examples and extended tables |
| [`automorphism_group_2026-07-12.tex`](automorphism_group_2026-07-12.tex) | The automorphism group of the triple-octonion product. Aut(Λ, +, ⋆) ≅ C₆ × S₃ of order 36, the ambient group G₂ × G₂ × S₃, the containment in Co₀, and what failed verification |
| [`hermes_mod2_structure_supplement.tex`](hermes_mod2_structure_supplement.tex) | The subalgebra and ideal structure of the mod-2 quotient L/2L: two externally proposed discoveries, verified and corrected |
| [`kirmse_1924_exposition.tex`](kirmse_1924_exposition.tex) | Kirmse 1924, re-derived from the original, with the forensic reconstruction of his count |
| [`historical_appendix.tex`](historical_appendix.tex), [`historical_appendix_expanded.tex`](historical_appendix_expanded.tex) | Drafts of the 1923–1946 historical synthesis, now Appendix B of the paper |
| [`Baez-Egan_Leech_Verification.pdf`](Baez-Egan_Leech_Verification.pdf) + [`_addendum_2026-05-24`](Baez-Egan_Leech_Verification_addendum_2026-05-24.tex) | Verification of the Baez–Egan construction, and the addendum comparing their L/2L treatment with ours |

## Update briefs

Drafted against a frozen version, then folded into the next one.

| File | Folded into |
|---|---|
| [`update_brief_2026-04-29_v3.tex`](update_brief_2026-04-29_v3.tex), [`update_brief_2026-04-29_v3_kirmse.tex`](update_brief_2026-04-29_v3_kirmse.tex) | v4 |
| [`update_brief_2026-05-25_v4.tex`](update_brief_2026-05-25_v4.tex) | v5 |

## Reviews

[`reviews/`](reviews/) holds every review of the paper, dated, with the model
that produced it named in the filename.  They are frozen on writing: a review
records what was found on the day it was run.

A review may come with its evidence.  [`reviews/2026-10-13_claude_opus_5_5_review_v8/`](reviews/2026-10-13_claude_opus_5_5_review_v8/)
is the external review of v8 by Claude Opus 5.5 (13 October 2026), incorporated
verbatim: its write-up (also at
[`reviews/2026-10-13_claude_opus_5_5_review_v8.md`](reviews/2026-10-13_claude_opus_5_5_review_v8.md)),
the reviewer's working scripts, and four independent reproductions in Python,
GAP and PARI/GP with their outputs.  Its `reproduction/final/run_all.sh` reads
the paper's source through the environment variable `MAIN_TEX`; from that
folder the working source is `MAIN_TEX=../../../../main.tex`.  Its GAP
extractor reads Appendix A by line range; for an exact rerun use the paper as
reviewed, `git show 6f7b37c:paper/main.tex`.  The folder is frozen as a whole.
