# Fine-comb review of v6 (15 July 2026)

**Date:** 2026-07-15
**Reviewer:** Claude Opus 4.8, at the direction of Jens Köplinger
**Manuscript:** `paper/main.tex`, v6 (15 July 2026), 23 pp.
**Prompt:** 226.
**Method:** mechanical pre-check (build, cross-references, overfull, stale
patterns), then five read-only review dimensions (technical correctness,
today's-edits regression, terminology, editorial, applied-audience fit and
honesty boundary), with adversarial verification of correctness findings.
Every number was re-checked by running the relevant script. Findings are
assessed against the state before this session's small post-launch edits
(date, Petersson de-naming, Bruck); those are reconciled below.

## Headline

No math error and no wrong number survived verification: the technical
dimension re-ran the scripts and reproduced the eight idempotents, the
126 × 32 = 4,032 square-zero count, the σ(Ls) 32/64 and 21/64 ratios, the
census table, and the 2¹⁶ / 2²⁵ span indices. The mechanical pass found the
paper clean: 0 undefined references, 0 overfull boxes, no dangling "Remark 4.9",
no orphaned labels from the un-remarking. 17 findings, all editorial or
boundary: 7 should-fix, 10 nits.

## Should-fix

### S1. §6: the 244,035,421 figure listed under "independently reproduced" (honesty boundary)

§6 opens "We have independently reproduced the full Baez–Egan verification
chain: … the derivation of the lower bound 244,035,421 …, and Egan's sharper
count of 17,280 Jordan rings." No runnable artifact in `python_project/src/`
reproduces 244,035,421; `egan_baez_count.py` reproduces only 270 and
17,280 = 270 × 64. Prompt log 089 records the intent precisely: *"Egan's count
of 17,280 … cannot be reverified without re-running the enumeration; we verified
the* derivation *of the lower bound 244,035,421 … Be honest about what was
verified."*

**Nuance the reviewer missed:** the paper already says "the *derivation of* the
lower bound 244,035,421", which matches prompt 089 exactly — the derivation was
checked, not the enumeration. So the wording is defensible. The residual risk is
only that it sits in a list under the umbrella "independently reproduced" where
the neighbouring items have runnable scripts, so a reader may assume one exists
for 244M too. **Judgment call for the author:** leave as-is (it is accurate), or
attribute the figure to Egan ("Egan's reported lower bound of 244,035,421") to
remove any ambiguity.

### S2. §1 scope paragraph: the four-application list (honesty boundary)

"…pointers for a reader interested in applications such as medical imaging,
cryptography, quantum computing, or phenomenology in the physics of fundamental
particles." Of the four, only particle-physics phenomenology has any pointer in
the paper (the Outlook references). The concrete list reads more promotional
than the restrained voice used elsewhere. **This is the author's own wording
(prompt 223), chosen deliberately** to signal applied relevance, so it is his
call. Options if he wants to soften: trim to the one application the paper
actually supports, or keep the list but frame it explicitly as prospective
("domains where such structures have been proposed to be useful").

### S3. §8 Outlook: "Programmes" → "Programs" (US spelling)

Line 1041, the lone British spelling against the US-English convention and the
paper's own "research program" (lines 112, 996, 1522). Clean fix.

### S4. §5.6: duplicated title (regression from the un-remarking pass)

The subsection heading "5.6 Closure of permutation cycles" is immediately
followed by "Remark 5.x (Closure of permutation cycles)" — word for word. This
is a leftover seam: 5.6 kept both a heading and a titled Remark wrapper, while
5.2–5.5 became plain prose. Fix: drop the Remark's bracketed title (keeping the
Remark environment for the table is fine), or unwrap it to prose like 5.2–5.5.

### S5. Abstract line 110: "also close" / "also produces" echo

Adjacent sentences at the seam between the applied-exposition list and the
"Reduced modulo 2L" continuation. Drop one "also".

### S6. §1 scope paragraph: "Further … further investigation" echo

"Further structural readings are offered as starting points for further
investigation…". (Note: this session's prompt-227 edit already removed the
companion "verified facts"/"machine-obtained facts" echo the reviewer also
flagged; only the "Further…further" remains.) Reword the opener, e.g. "The
structural readings that follow are offered as starting points for further
investigation."

### S7. Abstract topic order no longer matches §5 order

The abstract lists idempotents → square-zero → automorphisms → span; §5 now runs
idempotents (5.2) → automorphisms (5.3) → square-zero (5.4) → span (5.5), after
today's move of square-zero. An abstract list is not bound to section order, so
this is borderline; if alignment is wanted, reorder to "…its idempotents, its
automorphisms, its square-zero elements, the span of its image…".

## Nits

- **N1.** rem:order-closure says products of minimal vectors "always" land on
  strictly higher shells, and the value set {16,…,128} / "all multiples of 16" /
  "never 8" is stated flatly — but for the full Min × Min (≈ 3.9 × 10¹⁰ pairs)
  these are *sampled*, not proven. Since Λ does contain norm-12 vectors, "a
  multiple of 16" is a genuine non-trivial universal, not a corollary of
  closure. Either prove it or soften "always" to an observed-on-samples framing.
- **N2.** Abstract uses `L` (in `Ls̄`, `Ls`, `2L`, `L/2L`) without naming it;
  add "leaves the E₈ lattice `L` invariant" at first use.
- **N3.** Abstract "mostly exact arithmetic" is an unexplained hedge (the reason,
  §5.1 sampling, is invisible to the abstract reader).
- **N4.** "an undercount of the total" (abstract, Dickson) is elliptical; the
  total (seven) is never stated in the abstract.
- **N5.** Abstract closing sentence ≈ the introduction's closing sentence
  (near-verbatim "developed through a systematic human–AI collaboration").
- **N6.** §5.3 footnote hand-codes the PDF path in `\texttt` rather than a
  `\scriptref`-style `\nolinkurl` macro; the monospace treatment differs from
  every other repo footnote.
- **N7.** §6 two consecutive related-work paragraphs both open "We have
  independently reproduced…".
- **N8.** "catalogued" (§6 line 882) is the British-preferred form; "cataloged"
  for strict US consistency (both are accepted US variants).
- **N9.** §5.3 could add a two-word signpost that the exact order 36 is the
  enumerated/computed result (the containment argument is proof; the count is
  computation).

## Verdict

The paper is in strong shape and technically sound. The only findings with any
substance are the two honesty-boundary items (S1, S2), and both are the author's
call rather than errors — the paper's wording on S1 is already careful, and S2
is his deliberate framing. Everything else is polish. Recommended before freeze:
S3, S4, S5 (mechanical); consider S1/S2 phrasing and N1's "always".
