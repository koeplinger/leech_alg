# Pre-resubmission journal-style review of v6 (15 July 2026)

**Date:** 2026-07-15
**Reviewer:** Claude Fable 5, at the direction of Jens Köplinger
**Manuscript:** `paper/main.tex`, v6 (15 July 2026), 23 pp., not frozen
**Prompt:** 232.
**Method:** eight read-only review dimensions matching the author's checklist
(references; formal math; computer-proof reproducibility; logic;
proof-vs-computation separation; editorial standards; prior art; submission
readiness), followed by independent adversarial verification of every finding
rated blocker or should-fix. 40 agents, ~2.2M tokens. 59 findings: 4 blockers,
28 should-fix, 27 nits. **31 of 32 verifications CONFIRMED, 0 refuted**; the
one unverified blocker (its verifier was lost to an infrastructure outage) was
then **confirmed by the orchestrator by direct construction** (see B1).

Line numbers refer to the v6 state at review time; three small user-directed
edits landed mid-review (Conclusion sentence, §6 half-sentence, §6 heading
full names, §5.6 opening) and are reconciled in the text below.

---

## Headline

The mathematics is sound and the references are in unusually good shape: all
27 bibliography entries were verified against live sources field by field, and
every hand-checkable derivation in the paper was independently re-derived,
including all 192 appendix certificate rows recomputed from the LaTeX source.
No number in the paper is wrong.

What the review did find is a small cluster of **scope overclaims** — places
where the text asserts more than the computation behind it covers. These matter
precisely because the resubmission's posture is "computer proof, honestly
labeled." Four are blockers in that sense. All are cheap to fix.

---

## Blockers

### B1. The square-zero census is a stratum, not a census (§5.1 Remark, §5.4)

Λ contains square-zero vectors of norm **24**: the sum of two blockwise
orthogonal hexagonal triples is purely imaginary, has equal block norms 8, sums
to zero, lies in Λ, and squares to zero. **Verified by explicit construction**
with the repo's own Λ-membership test and ⋆ implementation (orchestrator, this
review). The paper's §5.4 sentence is literally true (the *root* triples give
4,032 vectors, all of norm 12 — that stratum count is exhaustive), but Remark
5.1's "it does, however, contain square-zero vectors, at norm 12" and the
surrounding framing invite the census reading "the square-zero vectors of Λ are
the 4,032 at norm 12," which is false.

**Fix:** state the scope: square-zero norms lie in 12ℤ; the norm-12 stratum
consists of exactly the 4,032 hexagonal root triples; higher strata exist
(norm 24 witness). One sentence in §5.4.

**Repo ripple:** CURRENT_STATE.md, the claims ledger (012), research_result.md
and terminology.md all say "Λ contains 4,032 square-zero vectors, all of norm
12" as a census. Those need the same scoping.

### B2. "L is a maximal order" at the paper's normalization (§2.2) [CONFIRMED]

At the paper's root-norm-2 normalization, e₀ ∉ L (L is even; N(e₀) = 1), so L
is not literally a maximal order in the unital sense of the cited classical
works — it is a √2-scaled copy of one. The paper elsewhere handles exactly this
distinction carefully (the non-unital gloss of "order" in the abstract), so
this is an internal-consistency gap, not a conceptual error.

**Fix:** say L is closed under multiplication (Table A.1) and is a √2-scaled
copy of a maximal order of the integral octonions; reserve unital maximality
for the unit-scale octavians of the cited literature.

### B3. "Every cycle type of S₇" covers only the single-cycle types (§5.6, §7) [CONFIRMED]

S₇ has 14 non-identity cycle types (partitions of 7); the census covers the six
single-cycle types — 2,365 of the 5,039 non-identity permutations. The mixed
types ((2,2), (3,2), (2,2,2), …) are not in the table; the (2,2) data was
deliberately removed from the paper this week. "We find for every cycle type of
S₇" (§5.6) and "Every cycle type of S₇ has been enumerated exactly" (§7)
overstate the coverage.

**Fix:** "for every single-cycle type of S₇" / "every *n*-cycle, n = 2,…,7", in
both places. The §7 open-questions bullet then correctly names the mixed types
as unenumerated.

### B4. The product-norm remark is the only §5 claim with no script and no evidence class (§5.1 Remark) [CONFIRMED]

The value set {16,…,128}, "all multiples of 16," mode 64, "neither 0 nor 8
occurs" are **sampled** facts (off-diagonal pairs are sampled; ~3.9 × 10¹⁰
pairs exist), yet the remark states them flatly with no footnote, no N, no SE —
uniquely in Section 5. What *is* exhaustive is the diagonal: all 196,560
squares u ⋆ u (script exists: `verify_idempotents_min_shell.py`). Remark 1.3
inherits the same issue ("never," "always land on strictly higher shells").

**Fix:** split by evidence class: exhaustive facts about squares (cited), then
sampled facts about pair products (cited, with N), and soften "always"/"never"
to the sampled framing. This is the paper's own separation standard applied to
its one non-compliant remark.

---

## Should-fix, by the author's checklist

### (1) References — all 27 entries verified live; one attribution error

Every `\bibitem` checked against the live web: authors, titles, journals,
volumes, pages, years, DOIs, URLs — **all correct** (details in the clean
report; includes confirming the Baez part dates, Calderbank's 2 April 2023
comment, Dickson's Theorem XV and the p. 321 Hurwitz assumption, Wilson 2009 in
full text, and that both repo mirrors are live). Two exceptions:

- **R1 [CONFIRMED].** The §6 footnote says Wilson's correction targets the
  formula in *Dixon 2010*; Wilson's bracketed note actually corrects Dixon's
  **2005 preprint** (Wilson's ref. [12], "Integral octonions, octonion
  multiplication and the Leech lattice"), and the locator "Section 5" is
  inexact in the archived version. Retarget the sentence.
- **R2 (nit).** [Elduque2023] "(NAART 2020)" — the printed book says NAART II,
  Coimbra 2022. Cosmetic; matches Elduque's own citation as is.

### (2) Formal math — clean

Everything re-derived independently: the σ-twist identity, the parity argument,
the basis-reduction logic, **all 192 appendix certificate rows recomputed from
the LaTeX source**, the Gram/index facts ([L:Ls] = [L:L s̄] = 16, sum L, meet
2L), the Remark 4.5 witnesses, the blockwise identities behind Propositions
4.6–4.8, the coset arithmetic of §5.5, the D ⊕ T ideal proofs, the
compactness-to-Co₀ chain, and both directions of the hexagonal-triple
characterization. Two wording-level items:

- **M1 [CONFIRMED].** The §4 WLOG footnote does not discharge the WLOG:
  Theorem 1.1 quantifies over any transposition, the proof covers (1 2), and
  the orbit fact as stated does not transport the proof (transpositions are
  conjugate by permutations that need not preserve the Fano orientation; the
  correct transport is by the order-21 Fano-preserving subgroup, under which
  the 21 transpositions form a single orbit *and* the construction is
  equivariant). One corrected footnote sentence fixes it.
- **M2 [CONFIRMED].** §5.3 calls the C₆ generator "a blockwise octonion
  automorphism"; it is an automorphism of the *twisted* product (O, ·_σ), not
  of the standard one. Specify the algebra once.

### (3) Computer-proof reproducibility — ran everything; five gaps

Every cited script exists at its cited path, runs, and its output matches the
claim it is attached to, with these exceptions/qualifications [all CONFIRMED]:

- **C1.** The Appendix A.1 polarization paragraph makes three computational
  claims (V two-sided ideal; both Lagrangians; μ̄ non-commutative) with **no
  script footnote** — the only such cluster in the paper. Cite
  `verify_discovery1_W_subalgebra.py` / `verify_discovery2_V_isotropic.py`.
- **C2.** The Kirmse eight-via-unit-stability reconstruction is verified by
  `verify_kirmse_1924_forensic.py`, which the paper never cites (it cites only
  `verify_kirmse_1924.py`, which does the 30/7/J₁ facts). Add the forensic
  script to the Appendix B footnote.
- **C3.** Appendix C says both implementations "reproduce the same integer
  coefficients displayed in the appendix" — neither prints coefficients; both
  verify memberships (and the GAP basis differs by sign on 7 of 8 vectors).
  Soften to "verify the same basis-pair identities."
- **C4.** §5.1's table is not reproduced by the script's default invocation
  (default reproduces an earlier ~802k run); the exact N = 10⁶ table came from
  two seeded runs summed. Record the two invocations in the repo (script
  docstring or README), or in the footnote.
- **C5.** `verify_all_cycles_exact.py` runs ≳15 min single-threaded; a runtime
  note in the docstring would save a referee a surprise. `aut_lambda_star.g`
  silently needs GAP's SmallGrp library for `StructureDescription` — add a
  requirements line.

### (4) Logic — sound; three seams

The dependency graph from definitions through Theorem 1.1 to §5 is acyclic and
complete; Propositions 4.6–4.8 discharge the theorem exactly [CONFIRMED
clean]. Seams: the σ = (1 2) scope is declared "throughout" in §4 but §5's
machine facts silently reuse it (one sentence at the head of §5); "quartic
power-associativity" is never defined in the paper (footnote: all five degree-4
parenthesizations agree); the φ span-index 2²⁵ is a new quantitative result of
the paper but is announced nowhere outside §6.

### (5) Proof vs. computation — the posture holds, with the B1/B4 caveats

The separation is real and mostly visible: Section 4 is proof; the §5.1 table
is labeled sampled with N and SE; the exhaustive computations cite their
scripts. Beyond the blockers: "strictly stronger" (quartic vs. cube, §5.1) is
an observed implication on samples, not a theorem — the algebra has zero
divisors, so no cancellation argument is available [CONFIRMED]; conversely two
*symbolic* classifications (the eight idempotents, the square-zero iff) read as
if they were computations — they deserve the word "symbolic," since no
enumeration could establish them over ℝ²⁴ [CONFIRMED]. The abstract's "mostly
exact arithmetic" hedges the wrong axis: everything load-bearing is exact; the
qualifier belongs on sampled-vs-exhaustive, not exact-vs-inexact [CONFIRMED].

### (6) Editorial standards — good compliance; residue

"stabiliser" in the Egan-bound footnote (the paper's one British spelling,
introduced 15 July); the §4 footnote's "They" has no antecedent after the
lead-in trim; "type-2 root" (§2.3) collides with the Type 1/2/3 minimal-vector
vocabulary defined 35 lines later; minor: "block-wise" once vs. "blockwise"
elsewhere, "cyclic perms" vs. "cyclic permutations" within one table, arXiv id
printed twice in [Corradetti2026], footnote punctuation drift. [heading
first-names finding obsolete — resolved by prompt 235 during the review.]

### (7) Prior art — fair and accurate, with two precision items

The Wilson/Dixon/Baez–Egan/Silverman/Calderbank/Corradetti credits check out
against the sources, including the part-by-part dates and the comment thread
[CONFIRMED clean]. Two corrections:

- **P1 [CONFIRMED].** The Corradetti paragraph credits an "order on E₈ via the
  Okubo algebra"; Corradetti's abstract explicitly disclaims that (the Okubo
  product forces ℚ(√3) coefficients and does not preserve the same ℤ-order).
  Rephrase to his actual claim (closed ℤ[√3]-structure / 2-adic saturation),
  keeping the 2-primary contrast that §6 already draws.
- **P2 [CONFIRMED].** "Egan's sharper count of 17,280" — the 244M and 17,280
  count *different things* (all Leech copies vs. Jordan-closed ones); "sharper"
  implies one refines the other. Say "Egan's count of at least 17,280 Jordan
  rings." Also, part 10 normalizes the Jordan product as ½(xy+yx), so our
  "doubled" is Baez's "quadrupled" — a parenthetical would spare a
  cross-reading referee [nit].

### (8) Submission readiness

- **No MSC 2020 codes, no keywords** — both venues expect them. Suggest MSC
  11H56, 17A75, 17D99, 11R52, 68V05 (author to confirm) and 4–6 keywords.
- **The abstract does not stand alone**: it contains `\cite{repo}` and a
  `\ref` to the appendix; abstracts circulate without bibliographies. Replace
  with plain text ("a historical appendix"; a spelled-out repository URL).
- **Repo permanence**: the bibliography names Bitbucket primary, every hotlink
  targets GitHub, and neither is archival. Recommend a Zenodo DOI snapshot at
  the submission commit, added to [repo].
- **The companion note as sole completeness reference (§5.3)**: the order-36
  completeness rests on a repo PDF authored by an AI agent and marked not
  peer-reviewed. For a journal, either state in one main-text sentence what
  completeness rests on (the ambient classification with its one
  non-enumerative step), or expect a referee to ask. [CONFIRMED]
- **Declarations block** (funding / competing interests / data availability)
  missing; AACA requires it. The data-availability statement is nearly free
  here given the public-record posture.
- Subtitle carries the working date/version — drop in the submitted copy.
- "(‘Donald’)" in the abstract: charming, but abstracts are indexed verbatim;
  consider moving the nickname to Appendix B.
- The §1 four-application list: flagged by two dimensions as the one sentence
  reading as overclaim. **The author has already ruled** (prompt 230): it stays,
  grounded in real prior interest. Recorded here as reviewed-and-decided, with
  one available softening if ever wanted: ground it as "fields from which
  interest in earlier versions of this work came."

---

## Clean-report coverage (what was checked and found in order)

All 27 bibliography entries field-by-field against live sources; the full
symbolic chain §2–§4 re-derived including all 192 certificate rows; every
cited script run with output matched to its claim (two runtime caveats); the
census, span, idempotent, square-zero, automorphism and Egan-bound numbers
reproduced; the argument's dependency graph; the Kirmse hedging ("our
reconstruction") intact; both repo mirrors live; build health (verified
separately) clean.

## Verdict

Ready for resubmission after a focused pass: the four blockers are all
scope-of-claim wording (B1–B4), fixable in under a page of edits total, plus
the reference retarget (R1), the Corradetti and Egan precision items (P1, P2),
and the submission apparatus (MSC/keywords, standalone abstract, declarations,
DOI snapshot). Nothing touches the mathematics, which held up under full
independent re-derivation. The paper's posture — computer proof for an applied
audience, honestly bounded — survives review *because* of these fixes, not
despite them: they are exactly the places where the text currently claims a
hair more than the machine checked.
