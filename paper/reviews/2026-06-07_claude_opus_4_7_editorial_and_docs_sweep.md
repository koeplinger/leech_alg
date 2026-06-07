# Editorial sweep + documentation sweep (v5, 2026-06-07)

## Summary
- Total findings: 45
- Required (must-fix): 7
- Files checked: 10

## Paper sweep

### Paper editorial standards

**Recommended**
- **main.tex:141-148 vs 299-334 vs 817 (Definition/Construction sections and Conclusion's open-question bullet)** — Vocabulary inconsistency: "octonion product" (dominant, ~20+ uses, used in Theorem statement, lemmas, Definition 3.2, and Related Work) and "octonion multiplication" (4 uses: lines 142, 148, 817, 1407) are used interchangeably for the same object. The compound noun "octonion-multiplication blocks" in the abstract (line 86) further conflates the two. Per standard 4, fix the convention.
  *Suggested fix:* Pick "octonion product" as canonical (it is already dominant and appears in the Theorem statement, the corollary, and Definition 3.2). Rewrite line 142 ("using the same octonion multiplication in all three blocks" -> "...same octonion product..."), line 148 ("twisting the octonion multiplication on each block" -> "...the octonion product..."), line 817 ("conjugate the octonion multiplication into a triple product" -> "...octonion product..."), and line 1407 ("octonion multiplication satisfies the composition property" -> "octonion product satisfies..."). Also reconsider the abstract phrase "nine octonion-multiplication blocks" (line 86) -> "nine octonion-product blocks" for consistency.
- **main.tex:650 vs 1381, 1383 vs 759** — Three different orderings of the same set of historical names appear: "Kirmse--Bruck--Coxeter relationship" (line 650), "Coxeter--Bruck transposition" (lines 1381, 1383), and "Bruck and Coxeter (1946)" (line 759). The first two refer to the same transposition device — once via the chronological triad, once via the credit pair. The mixed orderings on the same object are mild inconsistencies. Per standard 4, fix.
  *Suggested fix:* Pick one convention for the transposition-credit pair. Recommend "Coxeter--Bruck transposition" as canonical (matches the historical appendix's own attribution that Coxeter records Bruck's repair, so Coxeter's paper is the published source and Bruck is the originator of the device). Leave "Kirmse--Bruck--Coxeter" in its one occurrence as the chronological triad — but consider rewording to "Kirmse/Bruck/Coxeter" or "as in the Kirmse--Bruck--Coxeter trajectory" so the ordering is clearly chronological-historical, not credit-style.
- **main.tex:627-630** — Redundancy. The two-sentence paragraph "The order does not satisfy any of the classical algebraic identities (associativity, alternativity, flexibility, power-associativity, commutativity). The defining property is \emph{closure}: the Leech lattice is stable under the product." restates what the property-by-property table directly above (lines 568-593) already exhibits with numerical fail rates, and the second sentence re-states the order/closure claim already made in Corollary 1.2 and Remark 1.3 (lines 162-187). Per standard 3, this is a redundant restatement after its proper place.
  *Suggested fix:* Cut the paragraph, or compress to a single linking sentence if a closing flourish for the section is wanted (e.g., "None of the classical algebraic identities holds; the defining property is closure."). Consider deleting entirely — the table already speaks for itself.

**Note**
- **main.tex:41, 43, 54, 56** — Four occurrences of `---` (em-dash). All are inside LaTeX comments (`% ---------`) used as visual separators between section headers in the source, not in prose. These are exempt under the explicit comment-exemption in the task.
  *Suggested fix:* No change needed. (Optional: if strict zero-em-dash policy is desired even in comments, replace with `% ====` or `% ----` style; not recommended.)
- **main.tex:334** — "All nine block-pair products use the same multiplication~$\cdot_\sigma$" — here "multiplication" refers to the operation symbol $\cdot_\sigma$. Reads naturally, but inconsistent with the surrounding "octonion product" vocabulary. Minor.
  *Suggested fix:* Optional: change to "...use the same product~$\cdot_\sigma$" for term consistency. Low priority.
- **main.tex:813-814 vs 815-818 (Conclusion open-questions bullets)** — Bullet 4 ("Classification ... non-associative, non-alternative, non-flexible, and non-unital ... outside the classical families") restates the table on lines 568-593 a second time in close succession (the third restatement of the same fact, after the table and the 627-630 paragraph). The Outlook section (lines 825-843) then makes the same observation a third time when motivating the ternary direction.
  *Suggested fix:* If the 627-630 paragraph is cut per the previous finding, this bullet can remain as the proper home for the classification-as-open-question. If 627-630 stays, consider trimming this bullet to just "Classification of $(\RR^{24}, +, \star)$ outside the classical families" without re-enumerating the negations.
- **main.tex:88-92 (abstract) vs 147-152 (introduction)** — Near-verbatim restatement of the $\sigma$-acts-on-$E_8$-but-moves-sublattices claim between abstract and introduction. This is explicitly allowed by the Corollaries ("An introduction should read as a compressed version of the document") and the abstract is by convention a compressed version too, so this is not a violation. Flagged only for awareness.
  *Suggested fix:* No change needed. The abstract/intro relationship is permitted by the editorial standards.

*Overall:* main.tex (v5) is largely compliant: no em-dashes in prose (all `---` are in LaTeX comments), σ is uniformly typeset as $\sigma$, "twist/transposition" vocabulary is consistent. Two real cleanup items remain — the "octonion product" vs "octonion multiplication" interchange (4 spots) and the redundant restating paragraph at lines 627-630 immediately after the properties table; a smaller issue is the mixed "Coxeter--Bruck" vs "Kirmse--Bruck--Coxeter" ordering for the same historical reference.

### Paper cross-references (/home/jens/Projects/leech_alg/paper/main.tex)

**Note**
- **bibliography (lines 1561, 1567, 1595, 1607)** — Four \bibitem entries are defined but never cited in the body: FureyHughes2025_TrioOfTrialities (line 1561), GresnigtGourlayVarma2023 (line 1567), Koeplinger2023 (line 1595), MarraniCorradettiZucconi2025 (line 1607). These produce no LaTeX warnings (unlike undefined cites), but they will not appear in the rendered references when using a numeric/sorted scheme that prunes uncited entries; with the manual thebibliography environment they will appear regardless. Decide whether to cite them somewhere or remove them.
  *Suggested fix:* Either add \cite{...} for each in an appropriate spot (e.g., FureyHughes2025 / GresnigtGourlayVarma2023 in the related-work / outlook section on triality and physics; MarraniCorradettiZucconi2025 near the Corradetti2026 discussion; Koeplinger2023 in the historical or methodology context) or remove the unused \bibitem lines.
- **labels never referenced (intro/construction/properties/proof sections)** — Fifteen \label{...} declarations are not targeted by any \ref/\eqref in the body: cor:order, def:twist, def:twisted-product, eq:P1, eq:P2, eq:P3, eq:triple, eq:twist-formula, prop:cond1, prop:cond2, prop:cond3, rem:consecutive-twists, rem:order-closure, sec:intro, sec:outlook. These are not dangling (they cause no warnings), but they are dead anchors. Some (sec:intro, sec:outlook, cor:order, prop:cond1-3) are reasonable to keep as stable hyperlink targets; others (eq:P1/P2/P3 superseded by eq:P1-std/P2-std/P3-std, def:twist, def:twisted-product, eq:triple, eq:twist-formula, rem:consecutive-twists, rem:order-closure) appear orphaned.
  *Suggested fix:* Keep section-level labels (sec:intro, sec:outlook) and named-statement labels (cor:order, prop:cond1, prop:cond2, prop:cond3) as stable anchors. Consider deleting clearly orphaned labels: eq:P1, eq:P2, eq:P3 (replaced by *-std variants), def:twist, def:twisted-product, eq:triple, eq:twist-formula, rem:consecutive-twists, rem:order-closure -- or add explicit references to them if they were intended as forward pointers.

*Overall:* Every \ref/\eqref target resolves to an existing \label and every \cite key resolves to an existing \bibitem -- no dangling references or undefined citations; the only issues are four uncited bib entries and a handful of orphan labels.

### Paper bibliography

**Required**
- **main.tex:1614 (Mahler1942)** — Bibitems are not in strict alphabetical order by cite key. Mahler1942 appears AFTER MarraniCorradettiZucconi2025 (line 1607), but alphabetically 'Mahler' < 'Marrani'.
  *Suggested fix:* Move the Mahler1942 bibitem above MarraniCorradettiZucconi2025 (i.e., place Mahler1942 between LepowskyMeurman1982 and MarraniCorradettiZucconi2025).
- **main.tex:1631 (Petersson2018)** — Petersson2018 appears AFTER SmithVojtechovsky2022 (line 1626), breaking alphabetical-by-cite-key order ('Petersson' < 'Smith').
  *Suggested fix:* Move the Petersson2018 bibitem above SmithVojtechovsky2022 (i.e., place Petersson2018 between Okubo1995_Book and SmithVojtechovsky2022).
- **main.tex:1643 (repo)** — Author name in the 'repo' bibitem is written 'J.~Koeplinger' (no umlaut), while the Koeplinger2023 bibitem on line 1596 writes the same author as 'J.~K\"oplinger' (with umlaut). The author should be spelled consistently across both bibitems.
  *Suggested fix:* Change 'J.~Koeplinger' to 'J.~K\"oplinger' on line 1643 to match line 1596 (and the \author{} field on line 76).

**Recommended**
- **main.tex:1611 (MarraniCorradettiZucconi2025)** — Journal abbreviation uses '\ ' instead of '~' between 'J.' and 'Phys.': 'J.\ Phys.\ A: Math.\ Theor.\'. Every other journal name in the bibliography uses 'J.~' (e.g., 'J.~Algebra' line 1605/1640, 'J.~Math.\ Pures Appl.' line 1523, 'J.~Korean Math.\ Soc.' line 1540).
  *Suggested fix:* Change 'J.\ Phys.\ A: Math.\ Theor.\' to 'J.~Phys.\ A: Math.\ Theor.\' for consistency with other J.* journal entries.
- **main.tex:1509-1513 (Corradetti2026)** — The arXiv identifier is written twice: once as plain text ('arXiv:2605.09333, 2026.') and once inside an \href ('\href{...}{arXiv:2605.09333}'). Other arXiv-cited entries (Baez2002, Dixon1995/2010, Elduque2023, FureyHughes2025, GresnigtGourlayVarma2023, KamiyaOkubo2015) give the arXiv ID only once, inside the \href.
  *Suggested fix:* Remove the bare 'arXiv:2605.09333, ' from the citation text, leaving only the year ('2026.') followed by the \href{...}{arXiv:2605.09333} link, matching the format of the other arXiv entries.
- **main.tex:1635 (Petersson2018)** — Link uses \url{...} whereas all other web/arXiv/DOI references in the bibliography (Baez2002, Dixon1995, Dixon2010, Elduque2023, FureyHughes2025, GresnigtGourlayVarma2023, KamiyaOkubo2015, Koeplinger2023) use \href{URL}{display-text}. The 'repo' bibitem and Baez2014 also use \url, but those are bare URLs without a citation label, so they are arguably a separate idiom; Petersson2018 is a labeled external lecture and would normally be rendered with \href like the rest.
  *Suggested fix:* Either (a) convert to \href{https://www.fernuni-hagen.de/mi/fakultaet/emeriti/docs/petersson/ass.-rem.-int.-oct.pdf}{lecture notes (PDF)} for consistency with other linked references, or (b) accept that bare-URL \url is used for non-arXiv/non-DOI web resources (Baez2014, repo) and leave as-is. Pick one convention and apply it uniformly.

**Note**
- **main.tex:1509-1513 (Corradetti2026)** — The arXiv identifier '2605.09333' parses as a May-2026 submission. The paper is also dated 2026 in the entry, so this is internally consistent, but the ID should be double-checked against the actual arXiv listing since it is unusual to cite a future-month preprint.
  *Suggested fix:* Verify the arXiv number against arxiv.org/abs/2605.09333 (or the author's published version) and correct if it was a typo for a different ID.
- **main.tex:1477 \begin{thebibliography}{XXXXXXXX}** — The widest label argument 'XXXXXXXX' (8 chars) is narrower than several actual labels used (e.g., 'Elduque2023_IsotropicNorm', 'FureyHughes2025_TrioOfTrialities', 'GresnigtGourlayVarma2023', 'MarraniCorradettiZucconi2025'). This may not affect rendering if the document uses author-year-style citations, but the conventional argument should be at least as wide as the longest label.
  *Suggested fix:* If labels are rendered (default LaTeX \thebibliography), widen the argument to e.g. 'MarraniCorradettiZucconi2025' (the longest label) or to a placeholder of equal width. If labels are suppressed by the citation style, this is cosmetic only.

*Overall:* Author-name initial formatting (\,/~ convention) is consistent throughout; the bibliography has two clear ordering errors (Mahler1942 and Petersson2018 out of alphabetical position), a name-spelling inconsistency for the author K\"oplinger between the Koeplinger2023 and repo bibitems, one journal-abbreviation spacing inconsistency in MarraniCorradettiZucconi2025, and a duplicated arXiv-ID rendering in Corradetti2026.

### Paper prose

**Required**
- **main.tex:704 (§6 Comparison paragraph)** — Undefined notation: "Whether $\Sigma(\LL)$ overlaps with one of the embeddings catalogued by Egan, or sits outside that family, is a question this paper does not address." The symbol $\Sigma(\LL)$ appears exactly once in the paper and is never defined. A first-time reader has no way to parse it. This same issue was already present in v4.
  *Suggested fix:* Replace $\Sigma(\LL)$ with the intended object — e.g. "Whether the embedding $\LL \subset \OO^3$ underlying our construction overlaps with one of the embeddings catalogued by Egan…" or "Whether the order $(\LL, +, \star)$ corresponds to a Leech embedding…". (If $\Sigma$ was meant to denote a specific Leech-embedding map, define it at the point of first use.)
- **main.tex:242-253 (§2.3 footnote on $L\bar s$ and $Ls$ as sublattices)** — Technical inconsistency between the footnote and what the proof actually uses. The footnote claims: "Closure under left-multiplication by $L$, by contrast, holds for $L\bar s$ but fails for $Ls$ — and that failure is precisely the obstruction resolved in this paper (Remark 4.6)." But the obstruction actually resolved in §4 is $Ls \cdot Ls \not\subseteq Ls$ (Prop 4.10 introduction line 508–510: "the block sum would require $Ls \cdot Ls \subseteq Ls$, which does not hold"), not $L \cdot Ls \not\subseteq Ls$. These are different statements (the former is the bilinear closure needed by condition (3); the latter is an $L$-module closure that is not what the proof's load-bearing lemma asserts). Cross-reference to Remark 4.6 reinforces the mismatch.
  *Suggested fix:* Rewrite the relevant sentence in the footnote to align with what is actually proved, e.g.: "$Ls$ is moreover not closed under the octonion product itself — $Ls \cdot Ls \not\subseteq Ls$ — and this is the obstruction resolved in the present paper (Remark 4.6)."

**Recommended**
- **main.tex:1383-1387 (Appendix B closing paragraph)** — Slightly inverts the v5 framing recorded in the project memory ("σ acts on Ls / L̄s, not on L"). The current text reads: "the present $\sigma$ acts on the coordinate-symmetric placement $L = D_8^+$ (where it is inert at the $E_8$ stratum) and does its work on the sublattice $Ls$". Saying σ "acts on $L$ … where it is inert" is a self-contradicting locution; readers internalize the first half (σ acts on L) before the qualifier lands. The cleaner v5 phrasing is that σ acts on $Ls$ and $L\bar s$, leaving $L$ (and the octonion algebra) invariant.
  *Suggested fix:* Reorder, e.g.: "the present $\sigma$ leaves the coordinate-symmetric placement $L = D_8^+$ invariant (it is inert at the $E_8$ stratum) and does its work on Wilson's sublattices $Ls$ and $L\bar s$".
- **main.tex:147-152 (§1 statement of main result)** — "… but it moves Wilson's sublattices such that inclusion holds." Reads awkwardly; "inclusion" has no referent at this point in the prose (the three inclusion-like conditions are only described informally in the previous sentence). A first-time reader has to guess what "inclusion" refers to.
  *Suggested fix:* Replace with explicit naming, e.g.: "… but it moves Wilson's sublattices $L\bar s$ and $Ls$ onto closed sublattices, restoring the third Wilson condition." (Compare the abstract on line 92, which is concrete; the introduction is vaguer.)
- **main.tex:144-145 (§1 sentence before the main-result paragraph)** — "The third, however, is generally not closed under octonion products, including Wilson's." Two issues: (a) the antecedent of "the third" is grammatically the third sublattice condition, but a condition is not the kind of thing that is or isn't "closed under octonion products" — the underlying sublattice is; (b) "including Wilson's" is ambiguous ("Wilson's octonion product"? "Wilson's third condition"? Wilson does not pick a distinguished octonion product). The intended content is that the sublattice $Ls$ is not closed under the standard octonion product in any of the seven $E_8$-aligned multiplication conventions.
  *Suggested fix:* Rephrase, e.g.: "The third condition, however, asks that the block sum land in $Ls$, and $Ls$ is not closed under the standard octonion product — neither Wilson's choice nor any of the other six $E_8$-aligned conventions."
- **main.tex:82-96 (abstract) vs §1 (introduction)** — The abstract previews a Witt decomposition of $L/2L$ into a complementary pair of Lagrangians as one of two headline structural results, but §1 makes no mention of this material at all. A reader who reads the abstract and then the introduction will not see this thread previewed before it appears (deep inside Appendix A.1). Either the abstract is over-promising or §1 is under-previewing.
  *Suggested fix:* Add a one-sentence pointer in §1 (after the main-result paragraph, or as a short paragraph before §2) along the lines of: "A complementary mod-2 reading of the proof identifies the construction as producing a Witt decomposition of $L/2L$ into a pair of Lagrangians; see Appendix A.1 (§A.1)." Or, conversely, tone down the abstract claim.
- **main.tex:83-86 (abstract, opening sentences)** — "… assembled from nine octonion-multiplication blocks under cyclic $\ZZ_3$ permutation of three octonion factors." The phrase "three octonion factors" is potentially confusing — readers may parse "factor" as "input of a binary product" (in which case $\star$ has two factors). The intended meaning is "three copies of $\OO$ in the decomposition $\RR^{24} = \OO_1 \oplus \OO_2 \oplus \OO_3$".
  *Suggested fix:* Replace "three octonion factors" with "three octonion blocks" or "three copies of $\OO$" (consistent with §3's terminology, which uses "blocks").
- **main.tex:586-587 (§5 properties table)** — The table reports a row labelled simply "Power-associativity" with 12.75% pass rate, distinct from "Cube power-associativity $x^2 \star x = x \star x^2$" (30.70%). The surrounding prose discusses "cube and quartic power-associativity" — implying the unlabelled "Power-associativity" row is the quartic identity — but the row itself does not state which identity is being tested. A reader scanning the table cannot tell which form of power-associativity is meant.
  *Suggested fix:* Rename the row to "Quartic power-associativity" and include the explicit identity, paralleling the cube row (e.g., $(x \star x) \star (x \star x) = ((x \star x) \star x) \star x$, or whichever specific quartic identity is tested).

**Note**
- **main.tex:573 (§5 properties table, first row)** — "Multiplicative identity | No | n/a (structural)" — putting "Multiplicative identity" in a table of identities tested by sampling is a category error (it's not an identity in $u, v$, and is established structurally in the remark below the table, not by sampling). Readers may briefly be puzzled by the "n/a (structural)" entry.
  *Suggested fix:* Either drop this row from the table (the structural remark below already covers it) or move it to a separate "structural facts" mini-list above the table.
- **main.tex:698-700 (§6 Comparison, S_3 phrasing)** — "… the cross-block $S_3$ extending the $\ZZ_3$ routing is supplied by the routing of Definition~3.3 rather than by the octonion product itself…" The grammar parses awkwardly because "the cross-block $S_3$ extending the $\ZZ_3$ routing" reads as one heavy noun phrase, and then "is supplied by the routing" repeats "routing". The v4-to-v5 changelog flags this as a deliberate rewrite from the v4 phrasing, and the intended fix is mathematically right, but the prose result is mildly clunky.
  *Suggested fix:* Consider e.g.: "… the cross-block $S_3$ symmetry that extends the $\ZZ_3$ routing is built into Definition~3.3, not into the octonion product itself…"
- **main.tex:815-818 (§7 Conclusion, final open-question bullet)** — "Whether other linear coordinate permutations of $\RR^8$ (beyond simple transpositions of imaginary axes) conjugate the octonion multiplication into a triple product that closes on $\LL$…" — this open question is partially addressed by Remark 4.7 (Consecutive twists), which already enumerates 2-fold composites. The conclusion bullet does not acknowledge Remark 4.7, making the open question look broader than it is.
  *Suggested fix:* Add a parenthetical pointer, e.g.: "… (cf. Remark 4.7, which enumerates two-fold composites of transpositions; the open question concerns the remaining linear automorphisms of $\RR^8$ stabilising $L$)."

*Overall:* The paper is internally consistent with the v5 σ-framing and free of v3/v4 residue overall, but two genuine issues need fixing — an undefined symbol $\Sigma(\LL)$ in §6, and a footnote in §2.3 that misstates the obstruction the paper resolves — plus several recommended prose tweaks around the §1 introduction's vagueness, the abstract/§1 preview asymmetry on the Witt-decomposition thread, and a §5 properties-table labelling ambiguity.

### Paper numerical claims (paper/main.tex)

**Note**
- **/home/jens/Projects/leech_alg/paper/main.tex:75** — Subtitle confirmed: '7 June 2026 (v5)'.
  *Suggested fix:* No change.
- **/home/jens/Projects/leech_alg/paper/main.tex:269, 276-280** — Minimal-vector count 196,560 and decomposition 720 + 11,520 + 184,320 confirmed. Arithmetic checks out: 3*240=720; 3*240*16=11,520; 3*240*256=184,320; sum=196,560.
  *Suggested fix:* No change.
- **/home/jens/Projects/leech_alg/paper/main.tex:563** — Section 5 sample size N = 10^6 per property confirmed.
  *Suggested fix:* No change.
- **/home/jens/Projects/leech_alg/paper/main.tex:685-686, 1273-1275** — Polarization count 270 x 64 = 17,280 stated correctly in two places; both attribute the lattice-level Lagrangian/complementary-pair picture to Baez/Egan (\cite{Baez2014}).
  *Suggested fix:* No change.
- **/home/jens/Projects/leech_alg/paper/main.tex:352, 385, 397, 408, 419, 457, 476, 513, 541** — Theorem-environment counter shared (lemma/proposition/remark) with [section] reset. In Section 4: lemma 4.1 (sigma-L), 4.2 (L-order), 4.3 (L-sigmaLsbar), 4.4 (sigmaLs-closed), remark 4.5 (Ls-not-closed), proposition 4.6 (cond1), 4.7 (cond2), 4.8 (cond3), remark 4.9 (consecutive twists). Matches the requested numbering exactly.
  *Suggested fix:* No change.
- **/home/jens/Projects/leech_alg/paper/main.tex:805-810** — S_3 subseteq Aut(Lambda,+,star) subsetneq Co_0 statement found with -I_{24} witness, mathematically consistent with the prompt. Two minor presentation notes: (a) it lives in Section 7 ('Conclusion'), not Section 8 ('Outlook'); the user said 'in section 7' which matches the section number. (b) It is rendered as a bullet item, not a footnote.
  *Suggested fix:* No change required to math; if the user expected a literal \footnote, none is present. Otherwise no change.

*Overall:* All six numerical and structural claims verified in paper/main.tex; no factual mismatches, only the minor presentation observation that the S_3/Co_0 statement is a bullet item in Section 7 (Conclusion) rather than a footnote.

## Documentation sweep

### CURRENT_STATE.md

**Recommended**
- **CURRENT_STATE.md:3** — The 'Last updated' header still reads '2026-05-23', but the document body references v5 of the paper dated 2026-06-07 and the v4->v5 changelog. The header date is stale relative to the v5 content described below.
  *Suggested fix:* Update line 3 to 'Last updated: 2026-06-07' to match the v5 paper date now described in the document.

**Note**
- **CURRENT_STATE.md:161 vs CURRENT_STATE.md:472** — Inconsistent Python test count: line 161 states '197 tests verify the foundations', while the repository map at line 472 says '157 tests verifying foundations'. (Outside the five points requested, but flagged as it appeared while verifying lemma/appendix references.)
  *Suggested fix:* Reconcile the two numbers (verify the current pytest count and update both lines to agree).

*Overall:* The five requested items (19 pages, v5 / 7 June 2026, both frozen snapshot paths, both changelog references, and the Lemma 4.1-4.4 + Appendix A.1 structure) are all consistent with the paper; the only blemish in scope is the stale 'Last updated' header date, plus an out-of-scope test-count mismatch worth noting.

### README.md

No findings.

*Overall:* README.md is fully consistent with the v5 / 7 June 2026 release: all four required elements (version/date, frozen snapshot path, both changelogs, preserved v4 snapshot) are correctly stated at lines 15-19, and all four referenced files exist on disk.

### editorial_standards.md

No findings.

*Overall:* The "Note on paper/main.tex" section is consistent: v5 is identified as current, both frozen snapshots (main_2026-05-25.tex and main_2026-06-07.tex) are named, and both changelogs (v3_to_v4_summary.md and v4_to_v5_changelog.md) are pointed to; all four referenced files exist on disk.

### 2026-05-22_plan.md

**Recommended**
- **/home/jens/Projects/leech_alg/evidence_and_reasoning/2026-05-22_plan.md:17-18** — The second paragraph still reads as the original v3-era framing: "`paper/main.tex` is frozen at **v3 (29 April 2026)** and under review. A window has opened for a substantial revision. **This update produces v4 — a real new released version**". With v4 frozen and v5 issued, this paragraph is stale relative to the top-of-file status and could mislead a reader skimming past the status block.
  *Suggested fix:* Either delete this paragraph (the status block now carries the same information accurately) or prepend a phrase such as "Original framing (v3-era):" so it is clearly marked as historical context rather than current state.
- **/home/jens/Projects/leech_alg/evidence_and_reasoning/2026-05-22_plan.md:20-24** — The revision-notes list stops at "Revision 4 (2026-05-23): Phases A–D complete; v4 released." There is no Revision 5 line marking the 2026-06-07 v5 minimal revision, even though the top-of-file status documents it. The revision log is therefore one entry behind the status block.
  *Suggested fix:* Add a line: "*Revision 5 (2026-06-07): v5 minimal revision released; see `paper/v4_to_v5_changelog.md` and the update brief `paper/update_brief_2026-05-25_v4.pdf`.*"

**Note**
- **/home/jens/Projects/leech_alg/evidence_and_reasoning/2026-05-22_plan.md:3-14** — Status block at top correctly reflects v5 dated 7 June 2026, and names both frozen snapshots (main_2026-05-25.tex for v4, main_2026-06-07.tex for v5), both changelogs (v3_to_v4_summary.md and v4_to_v5_changelog.md), and the v5 update brief (update_brief_2026-05-25_v4.pdf). All three required references are present and the referenced files exist on disk.
  *Suggested fix:* No change needed — the three required items are all satisfied in the status block.
- **/home/jens/Projects/leech_alg/evidence_and_reasoning/2026-05-22_plan.md:12-14** — The v5 update brief is named as `paper/update_brief_2026-05-25_v4.pdf` — i.e. the filename still carries the v4 tag even though the brief drove the v5 minimal revision. The file does exist at that path (verified), so the reference is correct, but the name itself encodes a small versioning mismatch (a brief named "v4" that produced v5). Not a defect in this plan file, just a thing to be aware of.
  *Suggested fix:* No change to the plan file required. If desired for future clarity, the brief could be renamed (e.g. `update_brief_2026-05-25_v4_to_v5.pdf`) and the reference here updated to match — but only if the file is renamed in lockstep.

*Overall:* The three required elements (v5 status, both frozen snapshots, both changelogs, v5 update brief) are all present and accurate in the top status block; the only inconsistencies are a stale v3-era second paragraph and a missing Revision 5 line in the revision-notes list.

### terminology / research_result / research_statement

**Required**
- **/home/jens/Projects/leech_alg/evidence_and_reasoning/terminology.md:347-348** — Stale 'Last updated' line dated 2026-05-23 and explicitly refers to 'main paper v4 Lemma 4.1', whereas the paper is now v5 (2026-06-07). The σ-on-L statement still matches v5's Lemma 4.1 (lem:sigma-L) verbatim, so only the dating/version attribution is stale.
  *Suggested fix:* Update to '2026-06-07' and replace 'main paper v4 Lemma 4.1' with 'main paper v5 Lemma 4.1' (or simply 'main paper §4'). The technical content is fine; only the trailing tag is out of date.
- **/home/jens/Projects/leech_alg/evidence_and_reasoning/research_result.md:97-98** — Stale 'Last updated: 2026-05-23' line that parenthetically notes the additions 'both completed for v4'. The §4 lemma proof and §5 N=10^6 tests still match v5 (no changes from v4 to v5 in those sections per v4_to_v5_changelog.md), but the v4 attribution lags two weeks behind v5 (2026-06-07).
  *Suggested fix:* Bump 'Last updated' to 2026-06-07 and either drop 'completed for v4' or change to 'unchanged from v4 in v5' to acknowledge the v5 freeze.

**Recommended**
- **/home/jens/Projects/leech_alg/evidence_and_reasoning/research_statement.md:76-78** — 'Last updated: 2026-05-23' lags v5 by ~2 weeks. The statement text itself is the open-question framing and contains nothing v5-specific to revise; only the dating is stale.
  *Suggested fix:* Bump 'Last updated' to 2026-06-07 with a brief note like 'no substantive changes for v5 freeze' (or leave the parenthetical as-is and only update the date).

**Note**
- **/home/jens/Projects/leech_alg/evidence_and_reasoning/research_result.md:21-27 and §'What makes this simple'** — Describes the cross-block symmetry as 'Z₃ cross-block routing' only. v5 §6 and §7 now record that the resulting product is fixed by S₃ extending the Z₃ routing (update brief item 2; changelog Change 2 and Change 3, with the §7 bullet stating S₃ ⊆ Aut(Λ,+,⋆) ⊊ Co₀). This is not an outright inconsistency — the paper itself still names the construction 'Z₃-symmetric triple-octonion product' in the theorem statement — but research_result.md never mentions the S₃ refinement that v5 emphasises.
  *Suggested fix:* Optional: add one sentence noting that the product is in fact fixed by S₃ on the three blocks (extending the Z₃ routing built into the definition), matching v5 §6/§7.
- **/home/jens/Projects/leech_alg/evidence_and_reasoning/terminology.md:347 (file-level)** — No reference to 'Proposition 3.2' or '19 pages'/'20 pages' anywhere in the three files (checked by grep). The construction-level description of σ as a 'transposition of two imaginary basis units' / 'linear involution of R⁸' matches v5 throughout; the phrase 'transposition-twisted algebra swap(s,t)' on line 35 refers to the algebra (O, ·_σ) qua algebra-with-twisted-product, which is consistent with v5 Definition 3.2.
  *Suggested fix:* No fix needed for these — recording as a clean check.

*Overall:* The three files are technically consistent with v5 — the σ description, the formula x ·_σ y = σ(σ(x)·σ(y)), the lemma structure (lem:sigma-L, lem:L-order, lem:L-sigmaLsbar, lem:sigmaLs-closed = Lemmas 4.1–4.4 ≡ A–D), and the 192-basis-products count all match; the only real issues are stale 'Last updated' dates and a v4 attribution in terminology.md that should be bumped to v5, plus an optional S₃-vs-Z₃ refinement opportunity in research_result.md.

## Verdict

The paper and its documentation are substantively ready to ship as v5: every numerical claim verified, all cross-references resolve, the σ-on-sublattices framing is consistent across the paper body, README.md, and editorial_standards.md, and the abstract/intro/lemma/proof structure matches the project memory. Seven required items remain — three bibliography ordering/spelling fixes (Mahler1942 and Petersson2018 alphabetisation, K\"oplinger umlaut consistency in the `repo` bibitem), two paper-prose corrections (the undefined $\Sigma(\LL)$ in §6 and the §2.3 footnote that misstates the obstruction as $L \cdot Ls$ rather than $Ls \cdot Ls$), and two stale-date/version-tag fixes in terminology.md and research_result.md. None block the v5 freeze in principle, but the two §-level prose issues (undefined symbol, mis-stated obstruction) are the only items that could mislead a careful reader and should be fixed before any external circulation; the rest are housekeeping that can be batched into a routine cleanup pass.
