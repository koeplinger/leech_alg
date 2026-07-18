# Forensic pre-freeze review of v6 (18 July 2026)

**Date:** 2026-07-18
**Reviewer:** Claude Fable 5, at the direction of Jens Köplinger
**Manuscript:** `paper/main.tex`, v6 (18 July 2026), 23 pp., not frozen
**Prompt:** 264.
**Method:** eight read-only review dimensions (mathematical correctness;
every-number-versus-its-script audit; reproducibility of the computer proofs;
bibliography forensics; house-style editorial sweep; argument-flow logic;
Experimental-Mathematics-referee pass; LaTeX build audit on an isolated copy),
followed by deduplication and two independent adversarial refuters per
finding. 85 agents, ~3.5M tokens. 49 raw findings, 38 after deduplication:
31 confirmed, 5 split, 2 refuted. **Every finding rated blocker or
should-fix, and every script- or build-related claim, was then re-verified
first-hand by the orchestrator** (spot-check scripts in the session
scratchpad; live fetches for the two URL claims; both discovery scripts
executed; the PDF bookmark file decoded; `main.log` re-scanned in Python
after `grep` silently treated it as binary).

Line numbers refer to the v6 state at review time (the state committed as of
18 July 2026, after prompt 263).

---

## Headline

The mathematics of the main theorem chain (Sections 2–4, Appendix A
certificates) is sound, and no number in the census, automorphism,
idempotent, square-zero, or span results is wrong. The build is clean of
errors and overfull boxes at 23 pages.

The review found two blockers, both of the same species as the last round's:
the text asserting something other than what the computation behind it
computes. **B1 is the serious one**: the displayed formula defining the
permutation-twisted product in §5.6 is not the operation the census table
reports, and under the formula as printed the table is false for every
non-involution cycle type. The fix is one inverse sign; no number changes.
**B2** is a normalization clash in the two sentences that open the
introduction. Beyond these: one bibliography URL whose content has moved out
from under the entry, a cluster of stale-prose issues in cited discovery
scripts, and the usual editorial residue.

---

## Blockers

### B1. §5.6 defines a different product than the one the census computes (line 813)

Line 813 prints:

> "For a general such permutation~$\pi$, the product
> $x \cdot_\pi y := \pi(\pi(x) \cdot \pi(y))$ may or may not close on~$\LL$."

The operation the census actually computes is
$x \cdot_\pi y := \pi^{-1}(\pi(x) \cdot \pi(y))$. The cited script chain is
explicit: `verify_consecutive_twists_exact.py` line 125 has
`def dot(a,b): return ap(pinv, omul(ap(perm,a), ap(perm,b)))`, where
`ap(p,v)` is the map $e_j \mapsto e_{p(j)}$, i.e. $\pi$, and `ap(pinv,·)` is
$\pi^{-1}$; `verify_all_permutations_exact.py` (the table's footnote) reuses
that `closes()` unchanged.

The two conventions coincide exactly when $\pi$ is an involution
($\pi^{-1} = \pi$), so Definition 3.2, the whole of Section 4, and
Theorem 1.1 are untouched: for $\sigma = (1\,2)$ the printed and computed
products are identical. But for general $\pi$ the printed formula is the
computed product composed with $\pi^2$ on the output, and the table is then
false for every non-involution row. Verified three independent ways:

1. Two review-workflow refuters each recomputed the census exactly from the
   paper's own Fano triples and Wilson conditions. Under the printed
   formula: 3-cycles close 0/70 (table says 35), (3,2) 0/420 (says 126),
   5-cycles 0/504 (says 252), (3,3) 0/280 (says 112), (4,2) 0/630 (says
   294), 6-cycles 42/840 (says 210), (5,2) 0/504 (says 42), 7-cycles 0/720
   (says 336). Under the $\pi^{-1}$ form, every entry, the totals
   1,764/5,039 (35.0%), the identity failure, and the 3-cycle
   one-orientation-per-subset description all reproduce exactly.
2. The orchestrator's independent spot-check (scratchpad,
   `spotcheck_pi_convention.py`, reusing the repo's own `closes()` and
   basis): 3-cycle (123) closes under $\pi^{-1}$, fails under printed
   $\pi$; (132) fails under both; the transposition (12) closes under both.
3. Code reading of `ap()`/`pinv` as above.

**Fix:** print $x \cdot_\pi y := \pi^{-1}(\pi(x) \cdot \pi(y))$, optionally
noting that $\pi$ is then an algebra isomorphism
$(\OO, \cdot_\pi) \to (\OO, \cdot)$ and that for the transposition $\sigma$
this agrees with Definition 3.2 since $\sigma^{-1} = \sigma$. No table
number changes.

A same-sentence secondary (from the logic dimension): what closes or fails
on $\LL$ is not the $\RR^8$-valued product $x \cdot_\pi y$ itself but the
triple product assembled from it. Suggested combined wording: "the triple
product assembled blockwise from $x \cdot_\pi y := \pi^{-1}(\pi(x) \cdot
\pi(y))$ may or may not close on~$\LL$."

### B2. Introduction asserts unimodularity and a norm-8 minimal shell of the same lattice (lines 145–148; echo at §6 lines 856–857)

Lines 145–148:

> "The Leech lattice~$\LL$ is the unique even unimodular lattice of
> rank~24 with no vectors of squared norm~2 \cite{ConwaySloane1999}. Its
> minimal shell $\operatorname{Min}(\LL)$ consists of 196,560 vectors of
> squared norm~8..."

No lattice satisfies both sentences. The unique even unimodular rank-24
lattice with no roots has minimal norm 4 (Conway–Sloane, cited in the same
sentence); the lattice with minimal shell of norm 8 is its $\sqrt{2}$-scaled
copy, of Gram determinant $2^{24}$, which is the paper's $\LL$ (the plain
Euclidean norm is used throughout: $N(u) = 3N(x)$ in §5.4, $q = N/2 \bmod 2$
in A.1). The paper states the exactly analogous scaling for $L$ with care
(lines 273–278, "a $\sqrt{2}$-scaled copy of a maximal order") but never
states it for $\LL$. Section 6 lines 856–857 repeat the clash: "and proves
that~$\LL$ is even and self-dual" is true only in Wilson's halved-norm
convention, which the paper's own `evidence_and_reasoning/terminology.md`
says must be translated when cited.

**Fix (one clause in each place):** in §1, e.g. "...with no vectors of
squared norm~2. In the normalization used throughout this paper ($E_8$
roots of norm~2), $\LL$ is realized as the $\sqrt{2}$-scaled copy of this
unimodular lattice; its minimal shell $\operatorname{Min}(\LL)$ consists of
196,560 vectors of squared norm~8..." In §6: "proves that $\LL$ is even and
self-dual (in his halved-norm convention; the $\sqrt{2}$-scaled copy used
here has determinant $2^{24}$)" or equivalent.

---

## Should-fix

### S1. The Dixon2005 URL now serves a different document (bibliography, lines 1804–1809)

`www.7stones.com/Homepage/8Lattice.pdf` is live (HTTP 200) but today serves
"Integral Octonions, Octonion XY-Product, and the Leech Lattice," dated
November 21, 2010 — essentially the Dixon2010 paper — not the 2005 preprint
"Integral octonions, octonion multiplication and the Leech lattice" that
Wilson's reference [12] cited at that URL in 2009. Verified by downloading
and reading the PDF's title page this session. Hot-linking a URL as the
address of the 2005 preprint when it serves the 2010 paper will confuse any
reader who follows it.

**Fix options:** (a) drop the hot link and cite the preprint by its locator
only: "preprint, 10 March 2005; reference [12] of [Wilson2009]"; or (b)
keep the URL but add "(the URL now serves a later, 2010, version)".
Recommendation: (a), as the locator through Wilson's published paper is the
stable one.

### S2. Appendix C claims the tests verify "the $E_8$ lattice is indeed a maximal order" (lines 1671–1672)

Contradicts the paper's own §2.2 scoping: at this normalization $L$ is a
non-unital order ($e_0 \notin L$), a $\sqrt{2}$-scaled copy of a maximal
order of the unit-scale octavians. No test verifies maximality at the
paper's normalization. **Fix:** "...and that the $E_8$ lattice is closed
under the octonion product, a $\sqrt{2}$-scaled copy of a maximal order of
the integral octonions" or similar.

### S3. Remark 5.1: "the order has no idempotent at all" (lines 677–678)

$0$ is an idempotent and lies in $\LL$; §5.2's own classification lists $0$
among the eight. **Fix:** "no nonzero idempotent."

### S4. The two discovery scripts cited in Appendix A.1 print stale mid-discovery prose that contradicts the paper (footnote at lines 1545–1546)

Both scripts were executed this session. Their computed results support the
paper; their prose does not:

- `verify_discovery1_W_subalgebra.py` prints "MISMATCH at (w0, w0)...
  Total mismatches: 2 of 16" and "(But L/2L itself is commutative from (a)
  -- contradiction.)" — it compares a superseded proposed multiplication
  table against the actual one, and the "contradiction" aside reflects
  mid-discovery confusion that the paper's non-commutativity footnote
  (lines 1551–1559) later resolved correctly.
- `verify_discovery2_V_isotropic.py`'s docstring asserts "W has mixed
  isotropy," the opposite of the paper's Lagrangian claim (lines
  1567–1571). Its own computation contradicts its docstring: the actual
  distribution on $W$ is $q=0 \to 16$, $q=1 \to 0$ (totally isotropic),
  against its printed "(proposed: q=0 -> 8, q=1 -> 8)".

A referee reproducing this footnote sees output reading as a failed or
self-contradictory verification of the very facts it is cited for.
**Fix recommendation:** write one clean verification script (e.g.
`verify_polarization_lagrangians.py`) asserting exactly the paper's three
claims — $V$ two-sided ideal; $V$ and $W$ both totally isotropic, hence a
complementary Lagrangian pair; $\bar\mu$ non-commutative on $\bar L$ — with
ALL PASS semantics, and cite it in the footnote (keeping or dropping the
discovery scripts as the author prefers). The discovery scripts themselves
are frozen-genre referenced programs and should not be rewritten in place.

### S5. Appendix B: "Each was obtained in the original and re-derived in exact arithmetic in this project" (lines 1593–1596)

Overstates the Mahler coverage: `verify_mahler_1942.py`'s covering-radius
check maximizes over 6,000 random sample points (exact rational arithmetic,
but sampled, so it certifies a lower bound on the covering radius, not the
theorem), and the ideal-principality statement is not re-derived at all.
Under the house rule that non-exhaustive runs are flagged, this blanket
sentence needs a qualifier, e.g. "Each was obtained in the original, and its
key computational content re-derived in exact arithmetic in this project
(for Mahler's covering-radius theorem, by exact-arithmetic sampling)."

### S6. Remark 4.5 asserts $Ls$ is not closed with no witness or citation (lines 481–490)

The remark's explicit witnesses establish $\sigma(Ls) \ne Ls$ and
$\sigma(L\bar s) \ne L\bar s$, but nothing in the paper witnesses the
lead claim "$Ls$ is \emph{not} closed under the standard octonion product"
— the only comparable claim in the paper carrying no evidence. The claim is
true; the orchestrator found an explicit witness this session (scratchpad,
using the repo's own `omul`/`inLs`): $u_1 = (-1,0,0,0,1,0,1,1)$ and
$u_2 = (-1,1,0,0,0,1,0,1)$ both lie in $Ls$ (norm 4), and
$u_1 \cdot u_2 = (0,-2,2,1,-2,-1,-1,-1) \in L\bar s \setminus Ls$.
**Fix:** add the witness pair to the remark (one sentence), optionally with
a small verification script.

### S7. §1: "the sublattice it involves is not closed under octonion products, including Wilson's" (lines 172–174)

The bare plural overclaims: the paper's central mechanism is that this
sublattice IS closed under the $\sigma$-twisted octonion product.
**Fix:** e.g. "is in general not closed under the standard octonion
product, Wilson's included" — leaving room for the twisted product that the
paper then constructs.

### S8. Two computational-reproduction claims cite only bare \cite{repo} (lines 878–882 and 911–914)

"We have independently reproduced this construction in the supplemental
code~\cite{repo}" (Dixon XY-product) and "Details and code are
in~\cite{repo}" (Baez–Egan reading) are the only two reproduction claims in
the paper without a file path, against the paper's own script-citation
convention. Matching scripts exist: `leech_dixon.py` (with `e8_dixon.py`,
`octonions.py`) for the former; `verify_phi_span_index.py` and
`egan_baez_count.py` are among those for the latter. **Fix:** add
`\scriptref` footnotes naming the entry-point scripts.

### S9. "catalogued" (line 926)

US-English house rule: "cataloged."

### S10. Keywords hyphenation clashes with the body (line 133–135 vs. 1050, 1072)

The keyword list has "nonassociative algebra; ... non-unital order" — mixing
closed-up and hyphenated forms side by side — while the body consistently
writes "non-associative" (and "non-alternative, non-flexible, non-unital").
Note MSC 17 is officially "Nonassociative rings and algebras," which may be
why the keyword took that form. **Fix:** either "non-associative algebra"
in the keywords (local, one word), or adopt "nonassociative" body-wide
(larger diff). Recommendation: the former.

### S11. Abstract wording residue (lines 106–121)

"We further provide a historical account from the original sources
(1923--1946) in a historical appendix" — "historical ... historical" in one
sentence; and the passage stacks "Furthermore" (106), "We further provide"
(116), "a further narrowing" (120). **Fix:** e.g. "A historical appendix
gives an account from the original sources (1923--1946)." — which removes
both the repetition and one "further."

### S12. The PDF bookmark for §4.4 drops $\bar{s}$ and asserts something false (line 536)

The heading "Condition~(2): pairwise sums lie in~$L\bar{s}$" produces the
PDF outline entry "Condition (2): pairwise sums lie in L" (verified by
decoding `main.out`): hyperref drops the math, losing $\bar s$ entirely, so
the bookmark asserts membership in $L$ — false — and the build emits
"Token not allowed in a PDF string" warnings. §4.5's bookmark renders as
"Condition (3): the block sum lies in Ls" (correct, since $Ls$ survives).
**Fix:** `\texorpdfstring{$L\bar{s}$}{L s-bar}` in the §4.4 heading (and
optionally `\texorpdfstring{$Ls$}{Ls}` in §4.5 for symmetry, which also
silences the warnings).

### S13. The three [H]-anchored coefficient tables produce badness-10000 pages (lines 1274, 1370; PDF pages 14 and 16)

`main.log` (re-scanned in Python; the file contains control bytes that make
`grep` treat it as binary and print nothing) shows Underfull \vbox
(badness 10000) while \output is active: page 16 is roughly 45% blank below
Table 2 and page 14 has a stretched bottom, because \begin{table}[H] can
neither float nor break. **Fix:** change [H] to [htbp] (or [p]) for
Tables 2 and 3 and let them float, or accept the white space deliberately.

---

## Nits

### N1. `verify_egan_lower_bound.py` docstring dates part 10 as "12 December 2014"

The live page (fetched this session) is posted December 1, 2014; 12 December
is part 11's date. The paper's bibliography (lines 1752–1753) is correct;
the script docstring is not. One-line docstring fix (referenced program;
fix before freeze with a prompt-log record, or leave with this review as
the record).

### N2. `verify_aut_lambda_star.py` carries no runtime note

It runs well over 10 minutes; the comparable long-running census scripts
document their runtimes in their docstrings and the src README. One-line
docstring/README addition.

### N3. Appendix B footnote run-on (lines 1596–1600)

"...\scriptref{verify\_coxeter\_1946.py} Verification notes accompany
each..." — missing period after the last script citation.

### N4. Five bibliography labels expose internal key suffixes (lines 1822, 1828, 1841, 1854, 1899)

Printed labels "[Elduque2000\_Triality]", "[Elduque2023\_IsotropicNorm]",
"[FureyHughes2025\_TrioOfTrialities]", "[Hall2019\_Moufang]",
"[Okubo1995\_Book]" against plain Author-Year labels on the other 25
entries. Plain labels (Elduque2000, Elduque2023, FureyHughes2025, Hall2019,
Okubo1995) collide with nothing (checked against the full bibitem list).
Optional-label-only change; \cite keys stay as they are.

### N5. "Further structural readings ... for further investigation" (lines 221–222)

Repeats "further" within one sentence. (The sentence is the author's own
15 July dictation; flagged only for the repetition.) E.g. "as starting
points for investigation."

### N6. Tense mismatch (lines 895–897)

"Baez left open the question ... and conjectures that the answer is
negative." Either "left ... and conjectured" or "leaves ... and
conjectures."

### N7. Two missing serial commas (lines 775, 842)

"have equal norms and sum to zero" (§5.4) and "The $3$-cycles, $5$-cycles
and $(3,2,2)$ types" (§5.6); the paper otherwise uses the serial comma.

### N8. "come in useful" (line 1642)

Chiefly British idiom; e.g. "would prove useful here."

### N9. $\FF_2^8$ vs. $\FF_2^{\,8}$ (line 907 vs. 1485, 1527, 1567)

Same object typeset two ways; three thin-space occurrences vs. one plain.

### N10. Symbol-noun compound hyphenation (lines 1066, 1068 vs. 919, 1085)

"$\ZZ_3$-symmetry," "$\ZZ_3$-action" vs. "$\ZZ_3$ routing," "$S_3$
symmetry." Pick one convention (unhyphenated after a math symbol is the
more common journal style).

### N11. "All entries are integer" (line 1133)

Standard: "All entries are integers."

### N12. "squared norm" vs. "norm" for the same $N$ (lines 146, 270 vs. 284 on)

The early occurrences precede the definition of $N$ and may be deliberate
disambiguation; if so, a parenthetical at first use ("squared norm; simply
'norm' hereafter") would make the convention explicit.

### N13. §7 census bullet vs. §5.6 (lines 1052–1057)

"Open: a reason for ... the exact one-half rates" does not exclude the
3-cycles, whose one-half rate §5.6 already accounts for exhaustively by the
orientation description; and the bullet's "explained only for the
$3$-cycles" upgrades §5.6's more careful "admits a description." E.g.
"...for the exact one-half rates beyond the $3$-cycles."

### N14. "Type~1/2/3" vs. "type-3$\times$type-3" (line 1712 vs. 327–331)

Capitalization of the family names differs between §2.3 and Appendix C.

### N15. Duplicated opener (lines 879–880, 898)

"We have independently reproduced..." opens two nearby passages in §6.
(Carried over from the previous review's nit list, still unruled.)

### N16. Outlook justification breakdown (lines 1064–1083, PDF page 13)

One line is set at badness 10000 (visible inter-word gaps around "and
others"). Rephrasing or a \sloppypar would fix it; cosmetic only.

### N17. Script-path footnote lines at badness 10000 (footnotes at 1545–1546 and 1596–1600)

The unbreakable \nolinkurl paths leave visible gaps in the 8pt footnote
lines. Cosmetic; could be accepted, or the footnotes' script lists broken
manually.

### N18. Loose bibliography lines (Baez2014, Elduque2023, Baez2002 entries)

Ten moderately underfull lines (badness 1496–5388), driven by the long
golem.ph.utexas.edu URLs. Cosmetic; \raggedright inside those entries or
acceptance.

### N19. `\geometry{verbose}` marginal-notes warning (line 5)

Every build logs "The marginal notes overrun the paper" from the class's
121pt \marginparwidth on letter paper; the document uses no \marginpar, so
there is no output effect. Silence with `\geometry{marginparwidth=90pt}` or
ignore.

---

## Refuted findings (recorded for completeness)

- *Abstract "hinted ... at a further narrowing" vs. Appendix B "unstated":*
  the two phrases target different objects (the narrowing's existence,
  which Kirmse hinted at, vs. its criterion, which he did not state) and
  are consistent; the abstract wording is the author's deliberate 11 July
  ruling.
- *Global build health:* aside from the items above, the build is clean;
  final page count 23.

---

## Recommendation

Fix B1 and B2 before freezing; both are one-to-two-clause repairs that
change no result. S1–S13 are each small and worth ruling on individually
before the freeze, since v6 freezes exactly once. The nit list is optional
polish; N1 and N2 touch only script docstrings, not the paper.
