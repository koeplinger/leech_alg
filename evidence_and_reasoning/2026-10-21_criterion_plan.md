# 2026-10-21 — Plan: the closure criterion in Section 5.6

Created: 2026-10-21.  Updated: 2026-10-24.

**Status: executed.**  This plan is item 2 of the closed-form plan,
`2026-10-20_closed_form_plan.md`.  The author's decisions of §1: (d), (e)
and (f) as recommended, and "is given here" for the tension around F7;
(g) the paper carries no credit note for Proposition 5.3, since the work
is a human–AI collaboration (Section 1, the declaration on the use of AI,
Appendix C) and its results are joint, while the records trace where each
result came from.  The four edits of §3 are in `paper/main.tex` as
written here: Proposition 5.3 (`prop:criterion`) and the reading of the
census table in Section 5.6 (`sec:census`), Section 1's overview, and
Section 7's third question (`sec:conclusion`).  `paper/v7_to_v8_changelog.md` states
the change under Section 5.6; the records of §5 are part of item 7 of
the parent plan.  Drafted by Claude Opus 5.5 at the direction of Jens
Köplinger.

---

## 1. The decisions of this item

The parent plan asks decisions (d) to (g) again before this item.  The edits
of §3 follow the recommendations.

- **(d) Twist.**  An invertible linear map $g$ of $\mathbb R^8$ with $g(L) =
  L$ (*recommended*).  Every twist the paper uses already preserves $L$ (the
  permutations of $e_1, \ldots, e_7$, $\sigma$ among them, and the
  reflections in the roots of $L$, $\rho$ among them), and the proposition
  then decides closure for every twist.  Maps with $g(L) \ne L$ are set
  aside by definition, though some close: $g = 2\,\mathrm{id}$ gives
  $2\star_{\mathrm{id}}$, which maps $\Lambda$ into $(2L)^3 \subseteq
  \Lambda$.  Otherwise the definition stays, the proposition carries the
  hypothesis $L \cdot_g L \subseteq L$, and Section 7 keeps a question for
  the twists outside it.
- **(e) The $\sigma$ instance.**  One sentence after the proposition: for
  $g = \sigma$ it gives Theorem 1.1 again.  No remark on Lemma 4.3; Section 4
  unchanged (*recommended*).
- **(f) Section 7, third question.**  Narrowed to what stays open: its last
  clause, "and which twists other than those of Section 5.6 close on
  $\Lambda$", which the proposition answers, is deleted (*recommended*); or
  the question is deleted.
- **Tension (F7).**  Since July, the sentence after the census table has said
  that no description like that of the 3-cycles, and no reason for the
  failure of every 4-cycle, "is known" (in its present wording since prompt
  314).  F7 of the parent plan (new results, pointed out only) nearly gives
  both in terms of $\pi$, through the Fano plane of (2.1): the failure of
  every 4-cycle and the number of closing permutations of every even cycle
  type, hence the exact halves, proved up to one lattice identity and one
  nonclosure witness.  The sentence therefore says "is given here" in place
  of "is known" (*recommended*: true whatever becomes of F7), or keeps "is
  known".  Section 7 keeps the question either way, as one the paper leaves
  open.
- **(g) Provenance.**  The criterion is Theorem A of the external review of
  13 October 2026 (ledger claim 017; the parent plan, §0, says what prompts
  315 and 317 added).  The paper names no review for its results so far,
  also not for $\rho$'s index 1 from the same review.  Besides the AI agent's
  general role (Section 1, the declaration on the use of AI, Appendix C), it
  names one specific AI contribution, in Appendix C: the Hermes proposals
  behind the Lagrangian framing of Appendix A, credited at the author's
  request (prompt 157) and in the author's wording (prompts 281 and 282).
  Options: none in the paper, the records carrying it (*recommended*), the
  question kept for the phase of prompt 319 in which the author absorbs and
  owns the work; or a footnote at Proposition 5.3 naming the review of 13
  October 2026.

---

## 2. What changes

- **Section 5.6, the twist.**  "A twist is an invertible linear map $g$ of
  $\mathbb R^8$." gains "with $g(L) = L$", and the sentence on $\sigma$ cites
  Lemma 4.1 for $\sigma(L) = L$.
- **Section 5.6, the criterion.**  After the definition of $\star_g$, a new
  Proposition 5.3: $\star_g$ closes on $\Lambda$ if and only if $L\bar s$ and
  $Ls$ are closed under $\cdot_g$, that is, if and only if $g(L\bar s)$ and
  $g(Ls)$ are closed under the standard octonion product.  Its proof takes
  about a dozen lines: the substitution $a = x+y$, $b = x+z$, $w = x+y+z$
  describes $\Lambda$ exactly; two identities of the routing give the block
  sum and a pairwise sum of a product; the test vectors $(a, 0, -a)$ and
  $(w, w, -w)$ give the converse.  One sentence then reads $g = \sigma$:
  $\sigma(L\bar s) = 2\mathcal O$ and $\sigma(Ls) = p\mathcal O$ are closed,
  by the closure of $\mathcal O$ and Lemma 4.4, which is Theorem 1.1 again.
- **Section 5.6, after the census table.**  The identity permutation's
  failure is read off the proposition ($Ls$ is not closed, Remark 4.5), in
  place of its separate argument.  A new sentence: each permutation $\pi$
  preserves $L$, as $\sigma$ does (Lemma 4.1), so by the proposition the
  column "Close" counts the $\pi$ for which $\pi(L\bar s)$ and $\pi(Ls)$ are
  closed.  The last sentence says that no description in terms of $\pi$ like
  that of the 3-cycles "is given here" for the other types that close in
  part, nor a reason for the failure of every 4-cycle (the tension of §1).
  The data sentences (the 4-cycles, the $(5,2)$ type, the exact halves, the
  3-cycle orientations) stay.
- **Section 7, third question.**  Its last clause is deleted (decision (f)).
- **Section 1.**  "It also tests twists other than $\sigma$ for closure on
  $\Lambda$." becomes "It also characterizes the twists other than $\sigma$
  that close on $\Lambda$."

Untouched: Section 4; the opening paragraph of 5.6; the census lead-in (the
author's wording, prompt 298) and the nonisomorphism paragraph (prompts 299
to 301); the paragraphs on the span, on deciding by finitely many products,
on the parity and on $\rho$; the census table and its footnotes.  The script
`verify_twist_span_index_census.py`, already cited at the table, checks the
criterion on all 5,040 permutations; no new script.

---

## 3. The edits, old → new

Each block is replaced exactly once (the edit script asserts it).

### 3.1. Section 1, overview of Section 5

Old:

```latex
index of the span of its products.  It also tests twists other
than~$\sigma$ for closure on~$\LL$.
```

New:

```latex
index of the span of its products.  It also characterizes the twists
other than~$\sigma$ that close on~$\LL$.
```

### 3.2. Section 5.6, the twist and the closure criterion

Old:

```latex
A \emph{twist} is an invertible linear map~$g$ of~$\RR^8$.  It defines
the twisted product $x \cdot_g y := g^{-1}\bigl(g(x) \cdot g(y)\bigr)$
on~$\OO$.  The triple product~$\star_g$ is that of
Definition~\ref{def:triple}, with $\cdot_g$ in all nine block
products.  Since $\sigma^{-1} = \sigma$, the twist $g = \sigma$ gives
the product of Definition~\ref{def:twisted-product}, and
$\star_\sigma = \star$.
```

New:

```latex
A \emph{twist} is an invertible linear map~$g$ of~$\RR^8$ with
$g(L) = L$.  It defines the twisted product $x \cdot_g y :=
g^{-1}\bigl(g(x) \cdot g(y)\bigr)$ on~$\OO$.  The triple
product~$\star_g$ is that of Definition~\ref{def:triple}, with
$\cdot_g$ in all nine block products.  Since $\sigma^{-1} = \sigma$,
the twist $g = \sigma$ (Lemma~\ref{lem:sigma-L}) gives the product of
Definition~\ref{def:twisted-product}, and $\star_\sigma = \star$.

\begin{proposition}\label{prop:criterion}
For every twist~$g$, the product~$\star_g$ closes on~$\LL$ if and only
if $L\bar s$ and~$Ls$ are closed under~$\cdot_g$, that is, if and only
if $g(L\bar s)$ and~$g(Ls)$ are closed under the standard octonion
product.
\end{proposition}

\begin{proof}
Write $xy$ for $x \cdot_g y$, while $L\bar s$ and~$Ls$ remain the
lattices of~\eqref{eq:disp1}.  Since $g(L) = L$,
Lemma~\ref{lem:L-order} gives $L \cdot_g L = g^{-1}(L \cdot L)
\subseteq L$, so every block of a product of vectors of~$L^3$ lies
in~$L$.  For $u = (x, y, z) \in L^3$ put $a = x + y$, $b = x + z$ and
$w = x + y + z$.  As $y + z = a + b - 2x$ and $2L \subseteq L\bar s$
(Section~\ref{subsec:wilson-leech}), $u \in \LL$ exactly when $a, b
\in L\bar s$ and $w \in Ls$, and then $u = (a + b - w,\ w - b,\ w -
a)$.  For two such vectors $u$ and~$u'$, Definition~\ref{def:triple}
gives
\[
  P_1 + P_2 + P_3 = ww', \qquad
  P_1 + P_2 = ab' + ba' - aa' + 2(bb' - wb' + ww' - bw'),
\]
and the other two pairwise sums follow by the $\ZZ_3$ symmetry of the
routing.  If $L\bar s$ and~$Ls$ are closed, the block sum lies in~$Ls$
and, since $2L \subseteq L\bar s$, the pairwise sums lie in~$L\bar s$.
Conversely, $(a, 0, -a)$ and $(w, w, -w)$ lie in~$\LL$ for $a \in L\bar
s$ and $w \in Ls$; the product of two vectors of the first kind has $P_1
+ P_2 = -aa'$, and that of two of the second kind has block sum $ww'$,
so closure forces both inclusions.
\end{proof}

For $g = \sigma$, the lattices $\sigma(L\bar s) = 2\MO$ and $\sigma(Ls) =
p\MO$ of~\eqref{eq:multiples} are closed under the standard octonion
product, by the closure of~$\MO$ and Lemma~\ref{lem:sigmaLs-closed}, so
the proposition gives Theorem~\ref{thm:main} again.
```

### 3.3. Section 5.6, the reading of the census table

Old:

```latex
The $2$-cycles, the transpositions of Definition~\ref{def:twist},
close by Theorem~\ref{thm:main}.  The identity permutation does not
close: for $a, b \in Ls$, the vectors $(a, a, -a)$ and $(b, b, -b)$
lie in~$\LL$, and the block sum of their product is $a \cdot b$, so
closure would require $Ls \cdot Ls \subseteq Ls$, which fails
(Remark~\ref{rem:Ls-not-closed}).  Among the cycle types in the table,
the $4$-cycles are the only one failing without exception; the $(5,2)$ type has the lowest
nonzero rate.  The $3$-cycles, $5$-cycles, and $(3,2,2)$ types
close at exactly half.  For the $3$-cycles the closing half admits a
description: for each three-element subset of imaginary indices, one
of the two cyclic orientations closes.  No such description is known
for the other types that close in part, and no reason is known for the
failure of every $4$-cycle.
```

New:

```latex
The $2$-cycles, the transpositions of Definition~\ref{def:twist},
close by Theorem~\ref{thm:main}.  By Proposition~\ref{prop:criterion},
the identity permutation does not close, since $Ls$ is not closed
under the standard octonion product (Remark~\ref{rem:Ls-not-closed}).
Among the cycle types in the table,
the $4$-cycles are the only one failing without exception; the $(5,2)$ type has the lowest
nonzero rate.  The $3$-cycles, $5$-cycles, and $(3,2,2)$ types
close at exactly half.  For the $3$-cycles the closing half admits a
description: for each three-element subset of imaginary indices, one
of the two cyclic orientations closes.  Each permutation~$\pi$
preserves~$L$, as $\sigma$ does (Lemma~\ref{lem:sigma-L}), so by the
same proposition the column ``Close'' counts the permutations~$\pi$
for which $\pi(L\bar s)$ and~$\pi(Ls)$ are closed under the standard
octonion product.  No description in terms of~$\pi$ like that of the
$3$-cycles is given here for the other types that close in part, nor a
reason for the failure of every $4$-cycle.
```

### 3.4. Section 7, the third open question

Old:

```latex
\item A reason for the closure pattern of Section~\ref{sec:census}
  beyond the $2$-cycles, in particular for the failure of every
  $4$-cycle and for the exact one-half rates; and which twists other
  than those of Section~\ref{sec:census} close on~$\LL$.
```

New:

```latex
\item A reason for the closure pattern of Section~\ref{sec:census}
  beyond the $2$-cycles, in particular for the failure of every
  $4$-cycle and for the exact one-half rates.
```

---

## 4. Checks done before this hand-over

- A scratch copy with the four edits builds with three `pdflatex` runs: no
  errors, no undefined references; 21 pages (20 before); the known 1.3 pt
  overfull vbox, and two pages stretched by the reflow (a gap after the
  proof on page 10, space at the foot of page 11), layout only.  Section 5.6
  was read in the PDF (pages 10 and 11).
- The proof's two identities, the membership of the test vectors and of the
  parametrized vectors, and the criterion against direct closure (on 22
  permutations, the identity and $\sigma$ among them) were checked in exact
  arithmetic; the criterion was checked on all 5,040 permutations, all 120
  root reflections, all 17,280 ordered pairs of the parent plan's F4 and
  several hundred twists in $\mathrm{GL}(L)$ at prompt 317.
- An independent review in three lenses (the mathematics, the paper as a
  reader sees it, fidelity to the parent plan and the records), each
  finding checked by a second agent that tried to refute it, recomputed with
  its own exact code: the parametrization (index $2^{12}$ in $L^3$), both
  identities for an arbitrary bilinear product, the $\mathbb Z_3$ step, both
  directions of the proof, and the criterion against direct closure on all
  5,040 permutations (every Close count of the table), on the 120 root
  reflections (56 close, $\rho$ among them) and on 600 twists in
  $\mathrm{GL}(L)$ that are not isometries.  It found no mathematical error.
  Applied from it: "the standard octonion product" in four places where the
  draft shortened it; in the proof, the clause that keeps $L\bar s$ and $Ls$
  in their meaning, $L \cdot_g L = g^{-1}(L \cdot L)$ made explicit, and the
  step $y + z = a + b - 2x$; Lemma 4.1 cited where $\sigma$ and the
  permutations are used as twists; Proposition 5.3 cited for the identity's
  failure; "the column Close" for "each row"; the qualifier "in terms of
  $\pi$" dropped from the 4-cycle clause and from Section 7, where it
  implied a reason in other terms; the tension of §1; the alternative of (d);
  the options of (g); the records of §5.  Refuted as taste: a signal in the
  untouched paragraph on finitely many products, and a second notation point
  in the proof.

---

## 5. On execution, and what stays for later items

The four edits are applied to `paper/main.tex`; the paper is built and the
changed pages read; the changelog records the change; the linter runs.
The parent plan's addendum notes decision (g), the provenance of
Proposition 5.3 in the paper, for the phase in which the author absorbs and
owns the work.  Appendix A's research-program sentence and Section 7's fourth
question (on polarizations) stay until item 3, which answers them.  The
records follow in item 7, as the parent plan provides:

- ledger claim 017, extended to every twist with $g(L) = L$ and every
  bilinear $\mu$ with $\mu(L, L) \subseteq L$, with the exact proof of F1;
- the glossary entry "Twist";
- CURRENT_STATE: in open item 1, the paragraph "Not in the paper, for the
  author's decision ..."; in open item 5, the clause "the closure of twists
  that are not coordinate permutations, beyond the reflection ρ", which the
  proposition answers for every twist with $g(L) = L$; and, in the paragraph
  "The span index depends on the twist", the sentence that gives the
  criterion as checked exhaustively, with its proof in the review;
- README's open bullet "A structural characterization of which permutation
  cycles close";
- the index line of `verify_twist_span_index_census.py` in
  `python_project/src/README.md`, which names only the span-index column and
  the reflection as what the paper cites it for.

Parked, outside this item: the paragraph after the census table restates
three facts of the table in prose (the 4-cycles, the $(5,2)$ type, the exact
halves), the kind of restating that prompt 320 removed from Section 5.1.
They stay in this item: the later sentences and Section 7's question refer
to the 4-cycles and the exact halves, and the $(5,2)$ clause shares the
4-cycles' sentence.
