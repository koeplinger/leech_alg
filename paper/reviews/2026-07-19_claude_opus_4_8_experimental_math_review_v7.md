# Journal-quality review of v7 for an Experimental Mathematics audience (19 July 2026)

**Date:** 2026-07-19
**Reviewer:** Claude Opus 4.8, at the direction of Jens Köplinger
**Manuscript:** `paper/main.tex`, v7 (19 July 2026), 24 pp., in progress
**Prompt:** 276.
**Method:** seven read-only referee dimensions chosen for the *Experimental
Mathematics* readership (reproducibility; the experiment↔proof relationship;
mathematical correctness; engagement with the experimental-math literature;
presentation; significance/verdict; a LaTeX build audit), followed by
deduplication and two independent adversarial refuters per finding. 36
agents. 17 raw findings, 14 after deduplication: **2 confirmed, 1 split, 11
refuted**. The safety classifier was unavailable for 7 of the verify agents;
the orchestrator therefore re-verified **both** confirmed findings first-hand
against the paper and the scripts before writing this review.

Line numbers refer to the v7 state at review time.

---

## Headline

The paper is in strong shape for this venue. The two dimensions *Experimental
Mathematics* weighs most heavily — reproducibility and the honest separation
of experiment from proof — came back essentially clean: every computational
claim is traced to a named, runnable script, exhaustive computation is the
default, and sampled figures are labeled as sampled. The build is clean at 24
pages, the mathematics (heavily reviewed through v6) surfaced no new error,
and the newly cited experimental-math works are described accurately. The
v7 additions (the 149-word abstract and the Section 6 experimental circle)
are internally consistent with the body.

Two **should-fix** items remain, both about epistemic/citation hygiene and
both cheap. Neither is a mathematical error.

---

## Should-fix

### SF1. Remark 1.3 states a sampled distinct-pair bound as a categorical fact (§1, lines 190–193)

> "The product of two minimal vectors has norm at least~16
> (Remark~\ref{rem:product-norms}), never the minimal value $\Nrm = 8$: such
> products always land on strictly higher shells."

The quantifier ("two minimal vectors", "never", "always") ranges over all
ordered pairs $(u,v)$, but the cited `rem:product-norms` splits the evidence
by epistemic status: the **self-product** case is exhaustive over all 196,560
minimal vectors ("Exhaustively over the minimal shell, the self-product norms
$\Nrm(u\star u)$ ... neither $0$ nor $8$ occurs"), whereas the **distinct-pair**
case is sampled ("On $10^6$ sampled pairs..."). So for $u \ne v$ the
"never/always" bound rests on $10^6$ draws out of $\sim\!3.9\times10^{10}$
ordered pairs, not on proof. This is exactly the distinction an *Experimental
Mathematics* referee weighs most closely, and it is the one place in the paper
where the house rule "sampled results are flagged as sampled" is dropped at the
point of a categorical claim. The gap is not merely pedantic: norm 12 is a
legitimate $\LL$-norm (the smallest square-zero vector has norm 12, and every
$\LL$-norm lies in $4\mathbb{Z}$), so "never 8, always higher shells" for
distinct pairs is an empirical observation the proven lattice geometry does not
force.

An exhaustive spot check (8 representative $u$, each ranged over **all**
196,560 $v$ — about 1.57M exact star-products) found the minimum always 16, so
the statement is almost certainly true; but that is still sampling in $u$, not
a proof.

*Note:* the 15 July fine-comb review already raised this as its item N1
("soften 'always' to an observed-on-samples framing"); the wording still
stands in v7, so this is a corroborated, not-yet-applied point rather than a
new one.

**Fix options (author's call — this remark's wording was deliberately set):**
(a) restrict the categorical claim to the exhaustively verified self-product
case, which alone serves the remark's purpose of severing order-closure from
shell-closure: "The self-product of a minimal vector has norm at least 16
(Remark~4.5), never the minimal value $\Nrm = 8$..."; or (b) keep the
general-pair statement and flag its basis: "...and on $10^6$ sampled distinct
pairs the same bound held."

### SF2. The §6 Kirschmer–Nebe header reverses the published byline (§6, line 987)

> "**Gabriele Nebe and Markus Kirschmer (2022).** In computational lattice
> theory, Kirschmer and Nebe classify binary Hermitian lattices over number
> fields by computer [KirschmerNebe2022]; ..."

The header names the authors in the order *Nebe, Kirschmer*, but the sentence
one line below, the bibitem ("M.~Kirschmer and G.~Nebe"), and the paper as
published in the target journal (*Exp. Math.* 31(1), byline "Markus Kirschmer
and Gabriele Nebe") all give *Kirschmer, Nebe*. The two sibling headers added
in v7 ("Gerald Höhn and Geoffrey Mason (2023)", "Vladimir Dotsenko (2025)")
each follow their paper's byline order; only this one is reversed, and the
"(2022)" pins it to that specific paper. An editor or referee at the target
venue recognizes a paper from their own journal and reads a flipped byline —
two words above the correct order — as a careless slip.

**Fix options:** (a) set the header to "**Markus Kirschmer and Gabriele Nebe
(2022).**" to agree with the byline, the sentence, the bibliography, and the
sibling-header convention; or (b) if the intent is to feature Nebe as the
research-circle member (she also authors the solo Nebe2012), name her alone:
"**Gabriele Nebe (2022).**", keeping "Kirschmer and Nebe" in the sentence.

---

## Split verdict

### Nebe's 72-dimensional lattice: "built from the Leech lattice" (§6, lines 990–991)

One refuter accepted the phrasing, one asked for precision: the lattice is the
tensor product of the **Barnes** lattice and the Leech lattice over the ring
of integers of $\mathbb{Q}(\sqrt{-7})$, so "built from the Leech lattice"
is accurate but partial. Optional: "...is a tensor product of the Barnes
lattice and the Leech lattice." For a one-line related-work mention the current
phrasing is defensible; flagged for the author's taste.

---

## Considered and set aside (refuted), two worth a second look

- **Bremner not cited alongside Dotsenko** (§6). Murray Bremner's computational
  work on identities in nonassociative algebras is the canonical body for this
  audience; the reviewers refuted this as not required (the paper cites Dotsenko
  as one example, not a survey), but a single Bremner citation would resonate
  with the readership. *Optional.*
- **The two in-body results tables are unnumbered/uncaptioned** (§5.1 identities,
  §5.6 census) while the three appendix certificate tables are numbered and
  captioned; the Conclusion refers back to one as "the table" across a section
  boundary. Refuted as stylistic, but numbering and captioning the two headline
  tables would help a referee navigate. *Optional.*
- **§1's four applied-fields sentence** (medical imaging, cryptography, quantum
  computing, particle phenomenology) carries no citation and a serious-journal
  referee may query it. This is the author's deliberate, already-softened choice
  ("fields from which interest ... came"); recorded for awareness, not as a
  defect.
- **Sampled figures cite their scripts rather than stating seeds in-text.** The
  cited scripts carry the two-seed convention and exact invocations; judged
  adequate for reproducibility. Optional to surface a seed in the caption.
- **Pre-existing cosmetics:** the three `[H]`-anchored appendix tables leave
  their pages ~40–50% blank (forensic S13) and ten bibliography/footnote lines
  are underfull from long URLs (forensic N18); both were ruled cosmetic in v6.
  The Outlook list's "non-power-associative" among closed-up "nonassociative /
  noncommutative / nonalternative" is the **deliberate** hyphenation exception,
  not an inconsistency.

---

## Verdict, as an Experimental Mathematics referee

**Minor revision.** The contribution — a genuinely new bilinear order on the
Leech lattice under a $\mathbb{Z}_3$-symmetric twisted-octonion product, with a
machine-obtained structural anatomy — clears the venue's bar on novelty and,
crucially, on reproducibility: a reader can re-run every number from the public
repository, and the exhaustive-vs-sampled discipline is exactly what this
journal expects. The one thing a skeptical referee will probe is significance:
the algebra satisfies none of the classical identities, so its interest rests
on existence plus the experimental structural facts. The paper's honest,
bounded framing — existence and anatomy now, deeper theory left open — is the
right defense, and v7's move to foreground the experimentation and methodology
is well aimed at this readership. Resolve SF1 (epistemic hygiene, the venue's
sharpest lens) and SF2 (a byline slip in a paper from the target journal), and
the manuscript reads as venue-ready.
