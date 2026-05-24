# Detailed bibliography reference check (paper v4, 24 May 2026)

**Date:** 2026-05-24
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 163 — detailed bibliography reference check.
**Scope:** every entry in the `\thebibliography` block of `paper/main.tex` v4 (24 May 2026); 27 entries total.

---

## Method

For each entry I checked:
- author name(s)
- title
- journal / venue
- volume, issue, year, page range or article number
- DOI / arXiv ID / URL

against two sources:
- **Project record**: `evidence_and_reasoning/references/*.txt` (entries assembled with full bibliographic detail during the research).
- **Web verification**: arXiv abstract pages (where the entry has an arXiv ID), via `WebFetch`.

Journal-only entries behind paywalls (e.g. Proc. Roy. Soc. London 1982, J. Korean Math. Soc. 1996) were cross-checked against the project record only.

## Summary

All 27 entries cross-check cleanly. Five low-grade observations and one item worth flagging for human verification follow; none is a substantive error.

---

## I. Verifications

### I.1 arXiv-cross-verified (web-confirmed)

| Entry | Bib statement | Web verification |
|---|---|---|
| Baez2002 | "The octonions", Bull. AMS **39** (2002), 145–205; arXiv:math/0105155 | ✓ (title, author, journal, volume, year, pages, arXiv ID all match) |
| Dixon1995 | "Octonions: invariant Leech lattice exposed", preprint 1995; arXiv:hep-th/9506080 | ✓ (arXiv lists title in caps as a typographic artefact; the bib's normal capitalization is conventional) |
| Dixon2010 | "Integral octonions, octonion XY-product, and the Leech lattice", preprint 2010; arXiv:1011.2541 | ✓ |
| FureyHughes2025_TrioOfTrialities | Phys. Lett. B **865** (2025), 139473; arXiv:2409.17948 | ✓ (article number, journal, year, authors all match) |
| MarraniCorradettiZucconi2025 | J. Phys. A: Math. Theor. **58** (2025), 075202; arXiv:2309.17435 | ✓ (matches the published article; DOI 10.1088/1751-8121/adafef from journal — not in bib but available) |
| Corradetti2026 | "Integral elements of Okubo algebra and the $E_8$-lattice", arXiv:2605.09333, 2026 | ✓ (arXiv listing confirms preprint status, math.RA category, submitted 10 May 2026) |

### I.2 Project-record-cross-verified (high confidence; no web check needed)

| Entry | Source | Notes |
|---|---|---|
| Baez2014 | `prior_art_orders_on_leech.txt` | n-Category Café posts, parts 9-10, Nov/Dec 2014, URLs included in bib |
| ConwaySloane1982 | `prior_art_orders_on_leech.txt` | Proc. Roy. Soc. London Ser. A **381** (1982), 275-283 |
| ConwaySloane1999 | `background_octonions.txt` | 3rd edition, Springer-Verlag, 1999 (Grundlehren vol. 290 in record; not needed in bib) |
| ConwaySmith2003 | `background_octonions.txt` | A.K. Peters, 2003 (city "Natick, Massachusetts" in record; absent from bib — acceptable) |
| Coxeter1946 | `background_octonions.txt` | Duke Math. J. **13** (1946), 561-578 |
| Dickson1923 | `background_octonions.txt` | J. Math. Pures Appl. (9) **2** (1923), 281-326 |
| Elduque1996 | `okubo_algebra.txt` | J. Korean Math. Soc. **33** (1996), no. 1, 183-203 |
| Elduque2000_Triality | `okubo_algebra.txt` | Linear Algebra Appl. **314** (2000), no. 1-3, 49-74; DOI 10.1016/S0024-3795(00)00105-1 (not in bib) |
| Elduque2023_IsotropicNorm | `okubo_algebra.txt` | NAART 2020, Springer Proc. Math. Stat. **427**, 2023, pp. 287-301; arXiv:2210.03456 |
| ElkiesGross1996 | `prior_art_orders_on_leech.txt` | Internat. Math. Res. Notices **1996**, no. 14, 665-698 |
| Hall2019_Moufang | `triality_and_algebraic_physics.txt` | Memoirs AMS vol. 260, no. 1252, 2019 |
| KamiyaOkubo2015 | `okubo_algebra.txt` | Preprint only; arXiv:1503.00614 |
| Kirmse1924 | `background_octonions.txt` | Sächs. Akad. Wiss. Leipzig **76** (1924), 63-82 |
| Koeplinger2023 | `okubo_algebra.txt` | Preprint, ResearchGate deposit, DOI 10.13140/RG.2.2.19003.39209 |
| LepowskyMeurman1982 | `prior_art_orders_on_leech.txt` | J. Algebra **77** (1982), no. 2, 484-504 |
| Mahler1942 | `background_octonions.txt` | Proc. Roy. Irish Acad. Sect. A **48** (1942), 123-133 |
| Okubo1995_Book | `okubo_algebra.txt` | Montroll Mem. Lect. Series Math. Phys. vol. 2, CUP 1995 |
| Petersson2018 | `background_octonions.txt` | Málaga Workshop lecture, 2018 |
| SmithVojtechovsky2022 | `okubo_algebra.txt` | Doc. Math. **27** (2022), 535-580; DOI 10.4171/DM/877 (not in bib) |
| Wilson2009 | `leech_lattice_octonionic_constructions.txt` | J. Algebra **322** (2009), 2186-2190; DOI 10.1016/j.jalgebra.2009.03.021 (not in bib) |

### I.3 Repository link (just updated)

| Entry | Bib statement | Notes |
|---|---|---|
| repo | `https://bitbucket.org/jenskoeplinger/leech_alg` (primary), mirrored at `https://github.com/koeplinger/leech_alg`, 2026 | ✓ (both URLs as the user specified) |

---

## II. Low-grade observations (optional polish, not errors)

These are choices the paper has made consistently. None is a discrepancy; they are noted in case a journal copy editor would flag them.

### II.1 Issue numbers omitted from journal articles

Several entries omit issue numbers where they exist:
- **Wilson2009**: bib gives J. Algebra **322** (2009), 2186-2190; full citation has issue 7 (per project record).
- **MarraniCorradettiZucconi2025**: bib gives J. Phys. A: Math. Theor. **58** (2025), 075202; the journal record adds issue 7 (per web).
- **Hall2019_Moufang**: bib gives Memoirs AMS vol. 260, no. 1252, which is the *issue*-equivalent for Memoirs; correct as written.
- **Baez2002**: bib gives Bull. AMS **39** (2002), 145-205; journal record adds No. 2 (per web).

For most algebra journals, volume + page range or article number is sufficient. **No change recommended** unless the target journal specifies issue numbers as required.

### II.2 DOIs absent from several journal-published entries

These have a published DOI but the bib lists arXiv ID or no link:

| Entry | DOI (not in bib) |
|---|---|
| Wilson2009 | 10.1016/j.jalgebra.2009.03.021 |
| Elduque2000_Triality | 10.1016/S0024-3795(00)00105-1 |
| SmithVojtechovsky2022 | 10.4171/DM/877 |
| MarraniCorradettiZucconi2025 | 10.1088/1751-8121/adafef |
| FureyHughes2025_TrioOfTrialities | (journal DOI also exists; bib has arXiv only) |

For arXiv-first culture, the current convention (arXiv link) is fine. If the target journal asks for DOIs, these are available.

### II.3 Petersson2018 — URL omitted from bib

The bib has *"Lecture, Málaga Workshop on Non-Associative Algebras, 2018."* but does **not** include the URL where the lecture notes are hosted. The project record has:
`https://www.fernuni-hagen.de/mi/fakultaet/emeriti/docs/petersson/ass.-rem.-int.-oct.pdf`

A reader wanting to consult Petersson 2018 currently has no pointer to the file. **Suggestion**: add `\url{...}` after the venue line. This is the only entry citing a non-journal resource without a URL.

### II.4 Author-name capitalisation consistency

The bib consistently uses initials + surname for all authors. Project records sometimes carry the full first name (e.g. "Nichol Furey", "Holger P. Petersson", "Alessio Marrani"). This is a stylistic choice and consistent within the bib; no change needed unless the target journal asks for full first names.

### II.5 ConwaySmith2003 city omitted

Bib gives *"A. K. Peters, 2003"* without city. AK Peters is in Natick, Massachusetts (per project record). Many algebra-journal styles allow either. No change needed.

---

## III. Items flagged for human verification

### III.1 Corradetti 2026 — preprint status to confirm by publication date

The bib cites Corradetti 2026 as arXiv:2605.09333, 2026, and the web check confirms the arXiv entry exists (submitted 10 May 2026, math.RA category). At time of writing it is a **preprint, not yet in a peer-reviewed venue**.

**Action when relevant**: if Corradetti 2026 reaches a journal between now and resubmission acceptance, the bib should be updated to include the journal citation. As written, "arXiv:2605.09333, 2026" is the correct interim form.

### III.2 Koeplinger2023 — citable form

The bib gives *"preprint, 2023"* with the ResearchGate DOI 10.13140/RG.2.2.19003.39209. This is your own self-deposited preprint. Two notes:

(a) If the target journal restricts self-citations of unpublished work, this might attract a copy-editor query.
(b) If you have an arXiv mirror of this preprint (or plan to deposit one), citing it would broaden access. The ResearchGate DOI is functional but harder to retrieve programmatically.

Neither is a correction; both are author's-choice items.

---

## IV. Verdict

The bibliography is clean. 27 entries; all author names, titles, journals, volumes, years, pages, and arXiv IDs cross-check against either the project's central reference registry (`evidence_and_reasoning/references/`) or arXiv directly. The journal-published entries (Furey-Hughes 2025, Marrani-Corradetti-Zucconi 2025) were spot-checked via WebFetch and confirmed.

**Recommended changes**: only one — add the Petersson 2018 URL to the bib (item II.3). Everything else is style choice or "wait until publication" (Corradetti 2026).

---

## V. Summary table

| Status | Count | Entries |
|---|---|---|
| Web-verified clean | 6 | Baez2002, Dixon1995, Dixon2010, FureyHughes2025, MarraniCorradettiZucconi2025, Corradetti2026 |
| Project-record-verified clean | 20 | Baez2014, ConwaySloane1982, ConwaySloane1999, ConwaySmith2003, Coxeter1946, Dickson1923, Elduque1996, Elduque2000_Triality, Elduque2023_IsotropicNorm, ElkiesGross1996, Hall2019_Moufang, KamiyaOkubo2015, Kirmse1924, Koeplinger2023, LepowskyMeurman1982, Mahler1942, Okubo1995_Book, Petersson2018, SmithVojtechovsky2022, Wilson2009 |
| Repository link verified | 1 | repo |
| **Substantive errors** | **0** | — |
| Suggested URL addition (Petersson) | 1 | Petersson2018 |
| Wait-and-update (Corradetti journal pub) | 1 | Corradetti2026 |
| **Total entries** | **27** | |
