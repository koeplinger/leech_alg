# Paper

Created: 2026-07-12.  Updated: 2026-10-24.

The write-up and everything that supports it.  Genres per
[`DOCUMENT_GENRES.md`](../DOCUMENT_GENRES.md): each released version is kept
as a dated snapshot, a public record that is never edited; the external
review of 13 October 2026 is kept as received; every other file here, the
working source `main.tex` (v8, in progress) among them, describes its
current state and is corrected in place.

The paper is published: J. Köplinger, *An order on the Leech lattice from a
Z₃-symmetric triple-octonion product*,
[doi:10.13140/RG.2.2.22093.19686](https://doi.org/10.13140/RG.2.2.22093.19686)
(ResearchGate), corresponding to the released v7 source `main_2026-07-19.tex`.

## The paper

| File | |
|---|---|
| [`main.tex`](main.tex), [`main.pdf`](main.pdf) | **The paper.** Working source, v8 in progress. Bibliography inline via `\bibitem`. Build: `pdflatex main.tex` three times (from a fresh checkout, two runs leave a rerun warning and an unsettled PDF outline). |
| [`main_2026-04-29.tex`](main_2026-04-29.tex) | v3, released |
| [`main_2026-05-25.tex`](main_2026-05-25.tex) | v4, released |
| [`main_2026-06-07.tex`](main_2026-06-07.tex) | v5, released |
| [`main_2026-07-18.tex`](main_2026-07-18.tex) | v6, released |
| [`main_2026-07-19.tex`](main_2026-07-19.tex), [`main_2026-07-19.pdf`](main_2026-07-19.pdf) | v7, the last released version |
| [`v3_to_v4_summary.md`](v3_to_v4_summary.md), [`v4_to_v5_changelog.md`](v4_to_v5_changelog.md), [`v5_to_v6_changelog.md`](v5_to_v6_changelog.md), [`v6_to_v7_changelog.md`](v6_to_v7_changelog.md), [`v7_to_v8_changelog.md`](v7_to_v8_changelog.md) (in progress) | what changed between versions |

## Companion and standalone notes

Each is self-contained and stands apart from the paper; none is in its referee
scope.

| File | |
|---|---|
| [`companion.tex`](companion.tex) | Companion to the paper: an expository note on the σ-twist, from the Kirmse–Coxeter example at the E₈ level to the twisted triple product on Λ, with worked examples (AI-authored training material) |
| [`automorphism_group_2026-07-12.tex`](automorphism_group_2026-07-12.tex) | The automorphism group of the triple-octonion product. Aut(Λ, +, ⋆) ≅ C₆ × S₃ of order 36, the ambient group G₂ × G₂ × S₃, the containment in Co₀, and what failed verification |
| [`hermes_mod2_structure_supplement.tex`](hermes_mod2_structure_supplement.tex) | The subalgebra and ideal structure of the mod-2 quotient L/2L: two externally proposed discoveries, verified and corrected |
| [`kirmse_1924_exposition.tex`](kirmse_1924_exposition.tex) | Kirmse 1924, re-derived from the original, with the forensic reconstruction of his count |
| [`historical_appendix.tex`](historical_appendix.tex), [`historical_appendix_expanded.tex`](historical_appendix_expanded.tex) | Drafts of the 1923–1946 historical synthesis behind Appendix B of the paper |
| [`Baez-Egan_Leech_Verification.pdf`](Baez-Egan_Leech_Verification.pdf) + [`_addendum_2026-05-24`](Baez-Egan_Leech_Verification_addendum_2026-05-24.tex) | Verification of the Baez–Egan construction, and a note on the L/2L parallel and how our paper builds on it |

## Update briefs

Drafted against a released version, then folded into the next one.

| File | Folded into |
|---|---|
| [`update_brief_2026-04-29_v3.tex`](update_brief_2026-04-29_v3.tex), [`update_brief_2026-04-29_v3_kirmse.tex`](update_brief_2026-04-29_v3_kirmse.tex) | v4 |
| [`update_brief_2026-05-25_v4.tex`](update_brief_2026-05-25_v4.tex) | v5 |

## Reviews

[`reviews/`](reviews/) holds every review of the paper, dated, with the model
that produced it named in the filename.  Each records an assessment and, where
it is wrong about a fact, is corrected in place; the external review below is
kept as received.

A review may come with its evidence.  [`reviews/2026-10-13_claude_opus_5_5_review_v8/`](reviews/2026-10-13_claude_opus_5_5_review_v8/)
is the external review of v8 by Claude Opus 5.5 (13 October 2026), incorporated
verbatim: its write-up (also at
[`reviews/2026-10-13_claude_opus_5_5_review_v8.md`](reviews/2026-10-13_claude_opus_5_5_review_v8.md)),
the reviewer's working scripts, and four independent reproductions in Python,
GAP and PARI/GP with their outputs.  The second pass on that review, with the
first-hand re-derivation of its index claims, is
[`reviews/2026-10-14_claude_fable_5_1_index_one_verification.md`](reviews/2026-10-14_claude_fable_5_1_index_one_verification.md);
ledger claims 017 to 021 cite it.  The package's
`reproduction/final/run_all.sh` runs `f1_setup.py`, which reads the paper's
source from the environment variable `MAIN_TEX` and parses its three Appendix A
coefficient tables; the working source has none, and on it `f1_setup.py` stops
with a `ValueError`.  The tables are those of the paper as reviewed,
`git show 6f7b37c:paper/main.tex`: save that source outside the repository and
pass its path as `MAIN_TEX`.  The GAP extractor
(`reproduction/independent/gap/extract_paper_tables.sh`) reads the same tables
by fixed line ranges of that source (on the working source it finds no rows):
pass the saved source as the first argument and an output path outside the
repository as the second (by default it writes `paper_tables.g` into the
current directory).  The folder is
kept as received, as a whole: each `run_all.sh` regenerates its outputs in place, so a
full rerun is made from a copy outside the repository.
