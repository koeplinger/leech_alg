# Documentation consistency sweep — post-v5 (2026-06-07)

Reviewer: Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
Scope: post-v5 cross-document consistency audit across five clusters — changelogs, reference files, frozen artifacts, READMEs and key claims, and live-docs cross-check.

## Summary

- Total findings: 36
- Required: 11
- Recommended: 11
- Note: 14
- Clusters checked: 5

The v5 release itself (date 7 June 2026, page count 19, frozen snapshots, changelog paths) is internally consistent across the seven live docs and the v3-to-v4 summary is accurate against the v4 freeze. Most issues are localised one- or two-line edits in auxiliary documents (changelog, two reference files, two READMEs, two frozen artifacts, and two stale lines in CURRENT_STATE.md). One required addition (the Gresnigt-Gourlay-Varma 2023 reference) and two required wording fixes are needed to align the reference files with the corrected `σ` framing established in the May commits and the standing user-memory note.

---

## Cluster 1 — Changelogs

### Required

- **`paper/v4_to_v5_changelog.md`, Verification section (around line 77)** — The claim "An independent diff of `paper/main_2026-05-25.tex` against `paper/main.tex` returned six hunks, all six classified as user-requested" is factually wrong. Running `diff main_2026-05-25.tex main.tex` returns 18 hunks, not 6. Several hunks are not covered by the 7 enumerated changes, so the "all six classified as user-requested" assurance also fails.
  *Fix:* Rerun the diff and report the correct hunk count (18), or rewrite the paragraph to list each hunk explicitly. At minimum, replace "six hunks" with the actual number and acknowledge the additional unenumerated edits.

- **`paper/v4_to_v5_changelog.md`, Summary of changes list (lines 13-19)** — The 7-item change list is incomplete. The diff against the v4 freeze shows additional substantive edits not recorded:
  (a) §2.3 footnote at lines 248-251 (v4) / 250-254 (v5) was rewritten — the claim changed from "closure under left-multiplication by `L` holds for `Lσ̄` but fails for `Lσ`" to "`Lσ` is not closed under the octonion product itself (`Lσ·Lσ ⊄ Lσ`)". This is a meaningful one-sided → two-sided closure clarification.
  (b) §5 properties table at line 584 (v4) / 587 (v5) was relabeled "Power-associativity" → "Quartic power-associativity".
  (c) §6 Comparison closing sentence at lines 701-703 (v4) / 705-708 (v5) was rewritten — `Σ(Λ)` notation became `σ(Λ) ⊂ O³` with explicit block-wise wording.
  (d) The bibliography was re-alphabetized: `Mahler1942` moved (was after `Koeplinger`, now between `LepowskyMeurman` and `MarraniCorradettiZucconi`); `SmithVojtechovsky2022` moved to its correct alphabetical position after `Petersson2018`.
  (e) `J.~Koeplinger` → `J.~K\"oplinger` in the repo bibitem (line 1624 / 1645).
  *Fix:* Add these five edits to the summary list (changes 8-12) with Per-change entries, or fold them into a "minor mechanical/editorial" bucket. The footnote rewording (a) is substantive enough to deserve its own numbered change record.

### Recommended

- **`paper/v4_to_v5_changelog.md`, Change 4 (line 46)** — Line-number citation "lines 1357-1362 of `paper/main.tex`" is off. The Kirmse new content actually spans lines 1359-1364 in the current `main.tex`.
  *Fix:* Update line range to 1359-1364, or remove explicit line numbers since they drift with subsequent edits.

- **`paper/v4_to_v5_changelog.md`, Change 2 (line 34)** — Line-number citation "lines 698-701" is off by one and incomplete. The §6 Comparison change actually spans lines 699-702 of `main.tex`, and there is also a separate hunk at lines 705-708 in the same paragraph (the `Σ(Λ)` → `σ(Λ) ⊂ O³` rewording) that Change 2 does not mention but that travels with the §6 edit.
  *Fix:* Update line range to 699-702 and either extend Change 2 to cover the related `σ(Λ)`-notation rewording at 705-708, or record that rewording as a separate enumerated change.

- **`paper/v4_to_v5_changelog.md`, opening sentence (line 7)** — The phrasing "v5 is a minimal revision implementing the four issues recorded in `update_brief_2026-05-25_v4.pdf` ... the four edits clarify three statements that were already in v4 and restore one citation" counts the brief-driven edits as four, but the per-change record lists seven changes and the diff contains more. The summary opening understates the scope of the revision.
  *Fix:* Adjust the opening to read along the lines of "v5 implements the four issues from the update brief, adds the Gresnigt-Gourlay-Varma 2023 citation and the Petersson acknowledgment, makes mechanical preamble updates, and includes a small number of editorial clarifications (see Summary of changes below)."

### Note

- **`paper/v4_to_v5_changelog.md`, Page count and structural impact section (line 71)** — Page count claim (v4 = 19, v5 = 19) is verified by `pdfinfo`; absorption of new content via removal of the `\clearpage` at line 1460 of the v4 freeze is correct (confirmed by diff). Frozen-source path claims verified.
- **`paper/v3_to_v4_summary.md`** — Spot-checked against the v4 freeze (`main_2026-05-25.tex`): bibliography additions (`Dickson1923`, `Mahler1942`, `Corradetti2026`, `MarraniCorradettiZucconi2025`), new Section A.1 mod-2-quotient subsection, new Appendix B "Historical note", and §5 table expansion (Cube/Quartic power-associativity split, Symmetric composition row) are all present as described. No drift from v5 is being flagged.
- **Cross-document terminology consistency** — Section name capitalization, path-reference style (`[paper/<file>](<file>)`), date format (DD Month YYYY), and author-signature format are consistent across both files.

**Cluster assessment:** `v3_to_v4_summary.md` is accurate and needs no changes. `v4_to_v5_changelog.md` has one factual error (the "six hunks" claim) and is materially incomplete: at least five edits visible in the v4-to-v5 diff are not enumerated. Two line-number citations are off by one to two lines. Terminology and path conventions are consistent.

---

## Cluster 2 — Reference files

### Required

- **`evidence_and_reasoning/references/triality_and_algebraic_physics.txt`** — Gresnigt-Gourlay-Varma 2023 (Eur. Phys. J. C 83, 747; arXiv:2306.13098) is not present anywhere in the `references/` directory. It was added to the paper bibliography in prompt 194 and is now cited alongside `Koeplinger2023`, `FureyHughes2025_TrioOfTrialities`, and `MarraniCorradettiZucconi2025` in the Outlook physics paragraph (`\cite` at line 849 of `main.tex`; `\bibitem` at lines 1569-1574). This file is the natural home — it already houses `FureyHughes2025_TrioOfTrialities`, the construction exposes an `S_3` cross-block symmetry (`S_3 ⊆ Aut(Λ,+,⋆)`, line 809), and the new paper passage explicitly groups `S_3` symmetry with triality and symmetric-composition-algebra structures.
  *Fix:* Add a `[GresnigtGourlayVarma2023]` entry. Authors: Niels G. Gresnigt, Liam Gourlay, Abhinav Varma. Title: "Three generations of colored fermions with S_3 family symmetry from Cayley-Dickson sedenions". Eur. Phys. J. C 83 (2023) 747. DOI 10.1140/epjc/s10052-023-11923-y. arXiv:2306.13098 [hep-ph]. Notes: cited in Outlook as a contemporary physics application invoking `S_3` family symmetry, paralleling the `S_3` cross-block symmetry of the present construction. Status: bibliographic entry only.

- **`evidence_and_reasoning/references/prior_art_orders_on_leech.txt:17`** — Summary describes the construction as "three copies of a transposition-twisted octonion product with Z_3 cross-block routing". This "transposition-twisted octonion product" framing is stale and contradicts the corrected framing established in subsequent commits (32efdcc, 168bdd0, 4af176f, 99aeb3a) and the standing user-memory note that `σ` acts on `Lσ`/`Lσ̄`, not on `L`. The paper's current abstract calls it a `Z_3`-symmetric triple-octonion product on `R^24` where `σ` is a transposition of two imaginary basis units fixing the octonion algebra and `L` but moving `Lσ`/`Lσ̄`.
  *Fix:* Rewrite the parenthetical clause to read along the lines of: "The specific construction in this project — a `Z_3`-symmetric triple-octonion product on `R^24`, with `σ` a coordinate transposition of two imaginary basis units that is an octonion-algebra isomorphism fixing `L` but exchanging the Wilson sublattices `Lσ` and `Lσ̄` — does not appear in the literature."

- **`evidence_and_reasoning/references/computer_algebra_systems.txt:37`** — Describes the GAP verification target as "closure of the Leech lattice under the transposition-twisted triple product". Same stale framing as `prior_art`: the product is the `Z_3`-symmetric triple-octonion product. Also, this file is the only `references/*.txt` without a `Last updated:` line.
  *Fix:* Replace "closure of the Leech lattice under the transposition-twisted triple product" with "closure of the Leech lattice under the `Z_3`-symmetric triple-octonion product" (or "the triple-octonion product of Theorem 1.1"). Add a `Last updated:` line.

### Recommended

- **`evidence_and_reasoning/references/background_octonions.txt:262`** — "Last updated: 2026-04-11" lags ~2 months behind today. The Kirmse and Dickson entries were materially edited after that date (the May commits 995bbfc, 99aeb3a touched Petersson/Kirmse credit).
  *Fix:* Bump to 2026-06-07 with a brief one-line note describing the most recent change (Kirmse positive contributions / Petersson1).

- **`evidence_and_reasoning/references/triality_and_algebraic_physics.txt:69`** — "Last updated: 2026-04-12" will be ~2 months stale once the Gresnigt-Gourlay-Varma entry above is added.
  *Fix:* Bump to 2026-06-07 with note "(Gresnigt-Gourlay-Varma 2023 added to accompany v5 paper bibliography update; `S_3` family symmetry context)".

### Note

- **`evidence_and_reasoning/references/leech_lattice_octonionic_constructions.txt:69`** — "Last updated: 2026-04-12" is ~2 months old; content (Wilson 2009, Dixon 2010) has not changed. Refresh on next touch.
- **`evidence_and_reasoning/references/okubo_algebra.txt:189`** — "Last updated: 2026-04-12" is ~2 months old; content appears consistent with current paper citations. Refresh on next touch.
- **`evidence_and_reasoning/references/prior_art_orders_on_leech.txt:261`** — "Last updated: 2026-05-23" notes only the `Corradetti2026` addition. Once the summary rewording is applied, advance the date and append "summary reworded to align with current `Z_3`-symmetric triple-octonion framing".

**Cluster assessment:** One required addition (Gresnigt-Gourlay-Varma 2023) and two required wording fixes ("transposition-twisted octonion product" / "transposition-twisted triple product") that conflict with the paper's current framing and the user-memory note. All five `Last updated:` stamps lag the current paper state; `computer_algebra_systems.txt` has no `Last updated:` line at all. Otherwise the reference files match the current paper bibliography.

---

## Cluster 3 — Frozen artifacts

### Required

- **`paper/Baez-Egan_Leech_Verification_addendum_2026-05-24.tex`, Section 5 ("What it would take to close the gap further"), line 205** — Factually wrong dimension: the text reads "our mod-2 picture works within `24`-dimensional `L/2L`". The same document defines `\bar L := L/2L ≅ F_2^8` in Step (2) of Section 3 (line 123), i.e. `L/2L` is 8-dimensional over `F_2` (since `L = D_8^+` has Z-rank 8). The 24 must be 8 (or, if the intended contrast was with the 27-dimensional Jordan ambient that hosts the Leech lattice, the sentence should be rephrased). This is a self-contradiction inside the artifact, not a date-driven framing issue.
  *Fix:* Replace "`24`-dimensional `L/2L`" with "`8`-dimensional `L/2L` (the `F_2`-algebra of a single octonion block)". If the contrast intended was with the 27-dim Jordan ambient, recast as "whereas our mod-2 picture works within the 8-dimensional `F_2`-algebra `L/2L` inherited from one octonion block of `L` alone".

### Recommended

- **`paper/hermes_mod2_structure_supplement.tex`, Abstract (line 50), Setting Section 1 (lines 87-88), and recurring throughout** — Section A.1 title is quoted as "The mod-2-quotient reading of Tables A.2 and A.3", but in `main.tex` (both the v4 frozen snapshot `main_2026-05-25.tex` and the current v5 `main.tex`) the appendix tables are not letter-prefixed: they render as Tables 1, 2, 3. The actual rendered subsection title is "The mod-2-quotient reading of Tables 2 and 3". This was already inconsistent at the artifact's freeze date (24 May 2026), so it is not a v4→v5 drift; but it remains a factually wrong cross-reference now.
  *Fix:* Replace every occurrence of "Tables A.2 and A.3" with "Tables 2 and 3" (lines 50, 87-88, 201, and any other instances). Alternatively, change `main.tex` to reset the table counter inside Appendix A — but the simpler fix is to update the supplement to match the actual numbering.

**Cluster assessment:** Across the seven frozen artifacts, two issues meet the bar for flagging. All v4 self-references are honest to artifact dates and should not be touched. All cited verification scripts exist at the paths given. All `\cite` keys resolve inside each artifact's own bibliography. The `companion.tex`, `historical_appendix.tex`, `historical_appendix_expanded.tex`, `kirmse_1924_exposition.tex`, and `update_brief_2026-05-25_v4.tex` artifacts pass: their factual claims (Dickson's three vs seven, Kirmse's eight vs seven, the `J_1` closure witness, the 21 transpositions / 7 orders count, Mahler's Euclidean bound 1/2 vs Dickson's 5/4, and the Coxeter-Bruck attribution) all agree with the current `main.tex` Appendix B and Section 5 wording. The lemma numbering (4.3, 4.4) referenced in the `hermes_mod2` supplement and the Baez-Egan addendum still matches v5, and Section A.1 / `app:mod2quotient` still exists with the same label.

---

## Cluster 4 — READMEs and key claims

### Required

- **`gap_project/README.md`, line 84** — References to "Definition 3.1, Prop. 3.2" are out of date. In current paper v5 (and v4), there is no Proposition 3.2 — the three Wilson-condition propositions are now `prop:cond1`, `prop:cond2`, `prop:cond3` in Section 4. The v3 paper had a `prop:iso` as Prop. 3.2 (the algebra-isomorphism statement), which was removed when transitioning to v4. Definition 3.1 currently refers to `def:twist` (Transposition); Definition 3.2 is `def:twisted-product` (the σ-twisted octonion product).
  *Fix:* Replace "paper Section 3 (Definition 3.1, Prop. 3.2)" with "paper Section 3 (Definitions 3.1, 3.2 — Transposition and σ-twisted product)". Optionally reference the algebra-isomorphism identity that was formerly Prop. 3.2 but is now stated in the displayed prose just after Def. 3.2 (lines 309-311 of `main.tex`).

- **`gap_project/README.md`, lines 85-86, 88, 151** — "Theorem 1.2" should be "Theorem 1.1". In v5 the only Theorem in Section 1 is `thm:main` = Theorem 1.1 (followed by Corollary 1.2 = `cor:order`, then Remark 1.3 = `rem:order-closure`). The earlier v3/v4 numbering may have had the main theorem as 1.2; v5 has it as 1.1.
  *Fix:* Change every occurrence of "Theorem 1.2" to "Theorem 1.1" in the gap_project README (three places).

- **`gap_project/README.md`, line 87** — "paper Definition 2.4" is wrong. In v5 the only Definition in Section 2 is `def:wilson-leech` = Definition 2.1.
  *Fix:* Change "paper Definition 2.4 / Wilson Section 3" to "paper Definition 2.1 / Wilson Section 3".

- **`gap_project/README.md`, lines 88-89, 157-158** — "Section 4 fix on multiplicative identity" is wrong. The non-existence of a multiplicative identity is treated in Section 5 (`sec:properties`) of `main.tex` (line 164: `(Section~\ref{sec:properties})`, lines 616-626). Section 4 is the proof of closure and contains no statement about a multiplicative identity.
  *Fix:* Change "Section 4 fix on multiplicative identity" to "Section 5 (Algebraic properties): no multiplicative identity" in both the layout block (line 88) and the "What each test file verifies" paragraph (line 157).

- **`prompt_logs/README.md`, lines 33-102** — The "Contents" table stops at prompt 069 (companion paper authorship & abstract, ~mid-April 2026) but the directory now contains prompts through 199 (`post_sweep_docs_consistency.txt`, 7 June 2026). Roughly 130 prompt-log entries (070-199) are not indexed. The most recent indexed entry pre-dates v4 freeze (prompt 189), the v5 work (prompts 193-197), and the editorial/docs sweep prompts (198, 199).
  *Fix:* Either (a) extend the table to cover prompts 070-199, or (b) replace the per-prompt table with a brief explanation that the canonical list is the directory listing itself (e.g. "For the full enumerated record see the directory listing; entries are immutable per Manifesto §12") and keep the table only for a curated set of milestone prompts.

### Recommended

- **`gap_project/README.md`, line 130** — "companion Section 8 (`L·_σ L ⊆ L`, so the twist does not undo Coxeter)" references the wrong companion section. Section 8 of `companion.tex` is "The standard triple product on `R^24`"; the section that addresses whether the σ-twist undoes Coxeter's fix is Section 10 ("Could the σ-twist undo Coxeter's fix?"). Section 6 ("The twist at the E_8 level changes the story") is also a candidate for the closure claim `L·_σ L ⊆ L` at the E_8 stratum.
  *Fix:* Replace "companion Section 8" with "companion Section 10" (or, if the intent is the E8-level closure statement, "companion Section 6").

- **`gap_project/README.md`, line 155** — "companion Example 4.7" does not match current companion numbering. Companion uses `\newtheorem{example}[section]`, so the standard-triple-product failure example (label `ex:triple-fails`, line 580) is Example 8.1 (section 8, first example). The (a,b) pair in Section 4 — referenced as "Example 4.5 pair" on line 164 — is `ex:Ls-fails` at line 333 of `companion.tex`, which is Example 4.2, not 4.5.
  *Fix:* Change "Example 4.5 pair" (line 164) to "Example 4.2 pair" (the `ex:Ls-fails` instance), and change "companion Example 4.7" (line 155) to Example 8.1 (the `ex:triple-fails` instance).

- **`gap_project/README.md`, line 119** — "Remark 4.5" matches v5 numbering for `rem:Ls-not-closed`. Consistent — but the `test_sublattices.g` claim that it verifies `Lσ·Lσ ⊄ Lσ` should also note the witness lives in the companion (Example 4.2 = `ex:Ls-fails`), which the README already cross-references elsewhere.
  *Fix:* No edit strictly required; optionally clarify that the companion's hand-checkable witness for `Lσ·Lσ ⊄ Lσ` is Example 4.2 (`ex:Ls-fails`), not Example 4.5.

- **`gap_project/README.md`, line 194** — "paper Remark 3.5 all transposition twists are equivalent up to relabelling" does not match v5. Remark 3.5 does not exist in v5; Section 3 contains only Definitions 3.1, 3.2, 3.3 followed by no remarks. The "21 transpositions yield seven maximal orders in seven sets of three" claim now lives in Section 7 (Related work / historical, ~line 1377 of `main.tex`). The independence-up-to-relabelling claim itself is no longer set off as a paper Remark; it is empirically verified in `python_project/src/consistency_checks.py`.
  *Fix:* Change the conventions note to: "The transposition σ is `(e_1 e_2)` throughout (the canonical choice; the 21 transpositions are equivalent up to `GL(3, F_2)` relabelling — verified in `python_project/src/consistency_checks.py`)". Drop the spurious "paper Remark 3.5" citation.

### Note

- **`gap_project/README.md`, line 126** — "Definition 3.2 of the paper" is correct: `def:twisted-product` is indeed Definition 3.2 in v5. Flagging that this reference is internally consistent with the file even though line 84 says "Definition 3.1, Prop. 3.2". Ensure the line-84 fix harmonises with line 126.
- **`evidence_and_reasoning/README.md`** — No v4/v5 references appear; the file lists key claims 001-008 and trials 001-007 with neutral statuses. All path references resolve to files that exist on disk. Optionally, the trial-007 row could update its references file name to match the canonical Python file (`trial_007_kirmse_twist.py`), but that's a stylistic note only.
- **`evidence_and_reasoning/key_claims/008_transposition_twist_order.txt`** — The forward-correction note (dated 2026-05-23) is internally consistent with the v4/v5 paper structure: it correctly states `σ(L)=L` (Lemma 4.1), that Lemma 4.5 was demoted to Remark `rem:Ls-not-closed`, and that v4 §4 has four lemmas plus three condition-propositions. This matches v5 exactly. The original-text references in the file (the file name `trial_007_triple_octonion_swap.py` versus current `trial_007_kirmse_twist.py`) are explicitly preserved per Manifesto §12. No update required.

**Cluster assessment:** `evidence_and_reasoning/README.md` and `key_claims/008` are in good shape. The two READMEs needing significant updates are `gap_project/README.md` and `prompt_logs/README.md`. The gap_project README has multiple stale paper cross-references from earlier paper versions (v3/v4): "Theorem 1.2" should be "Theorem 1.1", "Definition 2.4" should be "Definition 2.1", "Prop. 3.2" no longer exists, "Section 4 fix on multiplicative identity" should be "Section 5", and the companion-section/example numbers (Section 8, Example 4.5, Example 4.7, Remark 3.5) don't match current `companion.tex` structure. The prompt_logs README's table stops at prompt 069 but the directory now contains prompts through 199 — a ~130-entry gap. None of the issues affect mathematical correctness; they are documentation-consistency problems introduced by paper revisions outpacing the auxiliary READMEs.

---

## Cluster 5 — Live-docs cross-check

### Required

- **`CURRENT_STATE.md:54`** — The verification table for trial 007 references `trial_007_triple_octonion_swap.py` as the source file, but no such file exists in `python_project/src/`. The actual file is `trial_007_kirmse_twist.py` (which `README.md` correctly cites at line 90). The other rows of the same table (`trial_007_scaled_test.py`, `trial_007_fast.py`, `trial_007_exhaust.py`) all match existing files.
  *Fix:* Replace `trial_007_triple_octonion_swap.py` with `trial_007_kirmse_twist.py` in the "Initial (trial 007 base)" row of the verification-status table in `CURRENT_STATE.md`.

### Recommended

- **`CURRENT_STATE.md:472`** — Repository-map block says "# 157 tests verifying foundations", but the prose at line 161 of the same file says "197 tests verify the foundations" and `pytest --collect-only` confirms 197 tests. The "157" looks like a stale count that was updated in prose but not in the ASCII tree.
  *Fix:* Change "# 157 tests verifying foundations" to "# 197 tests verifying foundations" on line 472 of `CURRENT_STATE.md` to match the prose at line 161 and the actual pytest count.

### Note

- **All seven docs (page count)** — Only `CURRENT_STATE.md` (line 360) and `v4_to_v5_changelog.md` (lines 3, 5, 71) explicitly state "19 pages" for v5. `README.md`, `editorial_standards.md`, `2026-05-22_plan.md`, `terminology.md`, `research_result.md`, and `research_statement.md` do not give an explicit page count. This is internally consistent (no doc disagrees with 19). Optionally add "(19 pages)" to the `README.md` v5 line for parity.
- **All seven docs (date 7 June 2026)** — Appears consistently in `CURRENT_STATE.md` (lines 3, 359), `README.md` (line 16), `editorial_standards.md` (line 70), `2026-05-22_plan.md` (lines 3, 5, 25), `research_result.md` (line 97), and `terminology.md` (line 347). `research_statement.md` uses "Last updated: 2026-05-23" — correct because that file is a research-program statement, not a paper-status doc, and its content has not changed since v3. No inconsistency.
- **All seven docs (frozen snapshot paths)** — `paper/main_2026-05-25.tex` (v4) and `paper/main_2026-06-07.tex` (v5) consistently referenced in `CURRENT_STATE.md` (lines 369, 372), `README.md` (lines 16, 19), `editorial_standards.md` (lines 69-70), and `2026-05-22_plan.md` (lines 9-10). Both files exist on disk.
- **All seven docs (changelog paths)** — `paper/v3_to_v4_summary.md` and `paper/v4_to_v5_changelog.md` consistently referenced in `CURRENT_STATE.md` (lines 373-374), `README.md` (lines 17-18), `editorial_standards.md` (line 73), and `2026-05-22_plan.md` (lines 10-11). Both files exist on disk.

**Cluster assessment:** The seven live docs are largely consistent on the v5 metadata (page count = 19, date = 7 June 2026, frozen snapshots, changelog paths). Two pre-existing internal inconsistencies in `CURRENT_STATE.md` were not caught by prompt 198's sweep: (1) the ASCII repository map says "157 tests" while the prose and actual `pytest` count both say 197 — likely a stale number from before the GAP/LOOPS work or before a later batch of tests was added; (2) the trial-007 verification table cites `trial_007_triple_octonion_swap.py`, a filename that does not exist on disk (the actual base trial-007 script is `trial_007_kirmse_twist.py`, which `README.md` cites correctly). Both fixes are localised one-line edits. No issues found in the other six docs; `research_statement.md` keeping its 2026-05-23 "Last updated" stamp is correct since its mathematical content didn't change between v4 and v5.

---

## Verdict

The project documentation is largely consistent for the v5 release on the high-stakes metadata that matters most: v5 date (7 June 2026), page count (19), frozen-snapshot paths, and changelog paths all agree across the seven live docs, and the v3-to-v4 summary is accurate against the v4 freeze. Eleven required items need attention before this sweep can be considered closed: one factually wrong claim and at least five unenumerated edits in `v4_to_v5_changelog.md`, one missing reference entry (Gresnigt-Gourlay-Varma 2023) and two stale "transposition-twisted" framings in the `references/` directory, one self-contradictory dimension claim in the Baez-Egan addendum, four out-of-date paper cross-references in `gap_project/README.md`, a 130-entry gap in `prompt_logs/README.md`, and one wrong-filename row in `CURRENT_STATE.md`'s trial-007 verification table. None of these affect the mathematical correctness of the v5 paper or its companion artifacts — they are documentation-consistency problems where auxiliary docs lag behind the paper revisions and the corrected σ-on-sublattices framing established in the May commits. The fixes are localised (mostly one- to three-line edits) and can be applied without disturbing any of the frozen v4/v5 PDFs or the underlying mathematical record.

---

*Reviewer: Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger, 7 June 2026.*
