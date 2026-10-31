# Historical-accuracy note — the Kirmse / Coxeter characterization

Created: 2026-05-17.  Updated: 2026-10-28.

**Date:** 2026-05-17
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 098. Concerns the description of the "Kirmse twist" / Kirmse–Coxeter
correction in `main.tex` (Remark *Kirmse's convention versus ours*; Remark
*$\sigma$ vs. a Fano-line permutation*; §6 *Kirmse (1924)* paragraph) and in
`companion.tex`. **No manuscript was changed** — v3 is frozen and several
points below need a primary source before any revision.

## Where the paper's current account came from

The paper's Kirmse account derives from **Petersson (2018)**, *Integral
octonions* (Málaga lecture), digested in prompt 037 and recorded in
`evidence_and_reasoning/terminology.md` and `background_octonions.txt`.
Petersson is a reliable source here — he obtained and read hard-copies of
Kirmse's and Mahler's original papers.

However, the specific phrase the manuscript uses — Coxeter corrected Kirmse
by **"permuting / swapping two basis indices in the Fano-plane labelling"** —
is **not Petersson's**. It was an AI-written gloss introduced in the prompt-037
summary and never checked against Coxeter's 1946 paper. That gloss is the
source of the sloppiness; its substance is right, since Coxeter (1946, §4)
records Bruck's correction as "transposing two of the $i$'s" (see below).

## What the sources actually say

Verbatim, Petersson (2018):

> "[Kirmse] exhibited a lattice $M \subseteq \mathbb{O}$ containing $1_\mathbb{O}$
> … Thus $(M, n_\mathbb{O}|_M)$ is a positive definite unimodular integral
> quadratic lattice of rank 8 … He then proceeded to claim **without proof**
> that … $M$ is closed under multiplication. Coxeter … was able to show that
> **Kirmse's assertion is false** … Coxeter then turned to **Bruck**, one of
> the leading experts in non-associative systems at the time, **who remedied
> the defect by a modification of Kirmse's original construction.**"

> "Mahler relied on a paper by **Dickson [5], dating back to 1923** … In this
> paper, Dickson constructs and investigates … an algebra of integral
> octonions which is **easily seen to be isomorphic to Coxeter's**. Thus the
> Coxeter octonions are more than twenty years older than is often assumed."

Baez (*Integral Octonions*, Part 6) and Wikipedia (*Octonion*) describe the
correction concretely: take the Kirmse integers and **switch the coordinate
$a_0$ (the identity, "$\infty$") with one imaginary coordinate $a_j$**; the
seven choices $j=1,\dots,7$ give the seven maximal orders.

### Confirmed facts

1. **Bruck's correction transposes two imaginary units.** Coxeter (1946,
   §4) states it for Kirmse's lattice $J_1$: "Bruck's domain $J$ can be
   derived from $J_1$ by transposing two of the $i$'s"; the 21 transpositions
   fall into seven sets of three and give the seven maximal orders.
   Interchanging the identity element $e_0$ with an imaginary unit $e_j$, as
   Baez and Wikipedia describe it and as Köplinger recalled ("replace the
   unit element … then factor through the affected triples"), also carries
   $J_1$ onto a maximal order, and the seven choices of $j$ give the seven
   maximal orders. The manuscript's "swap two basis indices in the Fano-plane
   labelling" matches Coxeter's description.
2. **Bruck supplied the correction.** Coxeter found the error; Bruck fixed it.
   The manuscript credits only Coxeter — an attribution gap.
3. **Dickson (1923)** independently gave a *correct* integral octonion algebra,
   isomorphic to Coxeter's, **23 years before Coxeter and before Kirmse's 1924
   paper.** Dickson did **not** repeat Kirmse's error; his work predates it.
   Coxeter acknowledged Dickson's priority in a postscript.
4. Dates and papers, 1920s–1940s:
   - L. E. Dickson, *A new simple theory of hypercomplex integers*,
     J. Math. Pures Appl. (9) **2** (1923), 281–326.
   - J. Kirmse, *Über die Darstellbarkeit …*, Ber. Verh. Sächs. Akad. Wiss.
     Leipzig, Math.-Phys. Kl. **76** (1924), 63–82.
   - K. Mahler, *On ideals in the Cayley–Dickson algebra*, Proc. Roy. Irish
     Acad. Sect. A **48** (1942), 123–133.
   - H. S. M. Coxeter, *Integral Cayley numbers*, Duke Math. J. **13** (1946),
     561–578.
   - **Bruck:** no separate publication is identified; his correction appears
     within Coxeter (1946).

### Settled by the primary sources

- **Did Kirmse himself propose the seven correct orders?** Wikipedia states
  the seven "were constructed by Kirmse (1924), Dickson and Bruck" and that
  Kirmse "thought there were eight … rather than seven." But Petersson — who
  read Kirmse's actual paper — describes Kirmse exhibiting **one** lattice and
  wrongly claiming it closed, and credits the correction to Bruck; Baez's
  account likewise does not credit Kirmse with the seven. **Kirmse's own text
  settles it** (`2026-05-20_claude_opus_4_7_kirmse_1924_verification.md`,
  `2026-05-20_claude_opus_4_7_kirmse_forensic.md`): he states eight maximal
  orders (p. 70) and names one, $J_1$, which is not closed; in the forensic
  reconstruction his eight are $J_1$ and the seven maximal orders.
- **Is the manuscript's Remark *$\sigma$ vs. a Fano-line permutation* correct?**
  It asserts Coxeter's correction "rewrites the multiplication table without
  being induced by a linear map on $\mathbb{R}^8$." **Coxeter (1946, §4)
  settles it:** the correction transposes two imaginary units, a linear
  coordinate transposition like $\sigma$, so the claimed dichotomy collapses;
  the difference is one of role, not kind
  (`2026-05-22_claude_opus_4_7_coxeter_1946_verification.md`).

## Recommendation (for the v4 revision; not for frozen v3)

1. Obtain **Coxeter (1946)** directly and consult **Conway–Smith (2003),
   ch. 12, "The octavian integers"** before rewording anything. These settle
   the precise form of the correction and the seven-orders attribution.
2. In the revision: (a) credit **Bruck**; (b) describe the correction
   accurately: a transposition of two imaginary units (Coxeter 1946, §4), the
   21 transpositions giving the seven maximal orders in seven sets of three;
   (c) re-examine Remark *$\sigma$ vs. a Fano-line permutation*: Coxeter's
   correction and $\sigma$ are both linear coordinate transpositions, and the
   contrast is one of role, not "linear map" vs. "not a linear map"; (d) consider citing
   Dickson (1923) and Mahler (1942), and noting Dickson's priority.
3. Köplinger's two instincts are both borne out: the wording *is* sloppy, and
   there *is* a substantive worry — the uncredited Bruck and the
   incorrect "not a linear map" claim. Neither is load-bearing for
   Theorem 1.1, but both belong in the revision.

## Primary sources to obtain (acquisition list)

The *historical* record — who constructed what, when, and what they claimed —
cannot be recomputed; it requires the original documents. (By contrast the
*mathematical* facts in this area — $E_8$ is a maximal order, there are seven,
the swap yields a closed order — should be verified computationally in-project,
as `L·L ⊆ L` already was, not taken on anyone's authority.)

1. **Kirmse (1924)** — *required.* J. Kirmse, "Über die Darstellbarkeit
   natürlicher ganzer Zahlen als Summen von acht Quadraten und über ein mit
   diesem Problem zusammenhängendes nichtkommutatives und nichtassoziatives
   Zahlensystem," Berichte über die Verhandlungen der Sächsischen Akademie der
   Wissenschaften zu Leipzig, Math.-Phys. Klasse **76** (1924), 63–82. German.
   Settles the contested question (one lattice, or seven?) and the journal
   name (below). Hardest to find — needs a research library holding
   the Saxon Academy's *Berichte*, or a German digitisation archive.
2. **Coxeter (1946)** — *required.* H. S. M. Coxeter, "Integral Cayley
   numbers," Duke Math. J. **13** (1946), 561–578.
   DOI 10.1215/S0012-7094-46-01347-6. English. Paywalled on Project Euclid
   ($30, "not available for individual sale"); free via a university library
   with DMJ access. Reportedly reprinted in Coxeter, *The Beauty of Geometry:
   Twelve Essays* (Dover, 1999), pp. 21–39 — worth checking (cheap, common),
   but verify the reprint actually contains this essay. Settles the precise
   form of the correction and how Bruck is credited.
3. **Dickson (1923)** — *recommended; already free online, no library needed.*
   L. E. Dickson, "A new simple theory of hypercomplex integers," J. Math.
   Pures Appl. (9) **2** (1923), 281–326. English (in a French journal). Open PDF on Numdam:
   `numdam.org/item/JMPA_1923_9_2__281_0.pdf`. Verifies Dickson's correct,
   pre-Coxeter construction.
4. **Mahler (1942)** — *optional, not load-bearing.* K. Mahler, "On ideals in
   the Cayley–Dickson algebra," Proc. Roy. Irish Acad. Sect. A **48** (1942),
   123–133. English. No claim in the paper needs it; it is the historical
   link (how Coxeter learned of Dickson) and concerns one-sided ideals of the
   integral octonions — which could become relevant if the ideal-theoretic
   structural questions (Conclusion; the table-reduction note) are pursued.

**Bruck:** there is *no separate Bruck paper* to obtain — per Petersson his
correction is recorded within Coxeter (1946). Reading Coxeter 1946 covers it.

**Journal name:** Kirmse's journal is the *Berichte über die Verhandlungen
der Sächsischen Akademie der Wissenschaften zu Leipzig*, as `main.tex`,
Petersson and `background_octonions.txt` give it
(`2026-05-20_claude_opus_4_7_kirmse_1924_verification.md`).

## Secondary sources consulted

- H. P. Petersson, *Integral octonions*, Málaga Workshop lecture, 2018
  (`fernuni-hagen.de/.../ass.-rem.-int.-oct.pdf`) — read in full.
- J. C. Baez, *Integral Octonions (Part 6)*, math.ucr.edu/home/baez/octonions/.
- *Octonion*, Wikipedia, "Integral octonions" section.
- J. H. Conway, D. A. Smith, *On Quaternions and Octonions* (2003), ch. 12
  — useful cross-check, but secondary; not where the history lives.
