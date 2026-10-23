# 2026-10-22 — Plan: Section 6, the Baez–Egan block

**Status (2026-10-23, prompt 330): engaged and executed.**  The author:
"Your recommendation appear good. Continue as planned."  Decisions (l),
(m) and (n′) as recommended, the rulings of prompt 250 kept, Appendix A
unchanged.  The author committed this plan and `verify_egan_construction.py`
at 4854e73, before the paper cites the script.  The five edits of §3 were
applied to `paper/main.tex` from this file, each old block found exactly
once; the result equals the reviewed scratch build.  Build: 21 pages, no
undefined references or citations, no hyperref warnings.  The script's
index line now names the paper's §6.  Changelog item 25.  The records of
§5 follow in item 7 of the parent plan.

**Status (2026-10-22, prompt 329): proposed, not executed.**  Item 4 of
`2026-10-20_closed_form_plan.md` (engaged at prompt 318), as its own plan.
Drafted by Claude Opus 5.5 at the direction of Jens Köplinger, and reviewed
before hand-over (§4).  Nothing in the paper has changed; the edits below are
applied only after the author has read them.  The new script of decision
(n′) is in the repository, not yet cited.

---

## 1. The decisions of this item

The parent plan asks decisions (l) and (m) again before this item.  The
research and the review for this plan raised three more points for the
author.  The edits of §3 follow the recommendations.

- **(l) The Comparison paragraph.**  The author asked for it at prompt 080
  ("It does not require a larger 27-dim construction ..."; "I'll review this
  paragraph, and may still decide not to have it at all (and let the reader
  make their own opinion)").  Step C (prompt 315) reduced it to two
  sentences, without its contrast with Dixon's XY-product.  The first
  restates Sections 3 and 5.4; its point on the 27-dimensional lattice stays
  in the Anatomy ("not on $L_L$ alone"; "distinct bilinear products on
  $\mathbb R^{24}$").  The second says that whether $\sigma(\Lambda)$ is one
  of the lattices "cataloged by Egan" is "a question this paper does not
  address"; the paper answers it: by (4.1) and Lemma 4.1, $\sigma(\Lambda)$
  is Egan's lattice for the pair $(2\mathcal O, p\mathcal O)$.
  *Recommended:* delete the paragraph and state the answer in the paragraph
  on Egan's construction.  *Or:* keep one sentence that gives the answer.
- **(m) The count $17{,}280 = 270 \times 64$.**  *Recommended:* in lattice
  terms (270 sublattices of $L$ isometric to $\sqrt2\,L$; each meets 64 of
  them in $2L$), with a pointer to Appendix A for the reading modulo 2, so
  that no term is used before Appendix A defines it.  *Or:* the reading
  modulo 2 in place (Lagrangians and their complements), before its
  definition.
- **(n′) The script behind the reproduction chain.**  The paper says that
  the full Baez–Egan chain was reproduced independently.  At prompt 082 the
  author asked to "mention that we've done a computer proof"; at prompt 089
  the author took the report of 16 April 2026, already in the paper folder,
  to be enough for that mention ("That's enough for us to briefly mention in
  the paper that we've independently verified their construction as
  correct, right?") and asked to "just state that we've done so and defer
  to the online repo for the code".  At prompt 230 the author had the lower
  bound reproduced by a script rather than left as it stood; the review that
  prompt answered (15 July 2026) said that the other links already had
  runnable scripts.  Two links have scripts: the lower bound
  (`verify_egan_lower_bound.py`) and the count (`egan_baez_count.py`).  The
  other four (Egan's simple roots and the two sublattices built from them;
  the lemma that these meet in twice that $E_8$ lattice; the Leech-lattice
  theorem for that pair; the doubled-Jordan-product closure) rest only on
  the report `paper/Baez-Egan_Leech_Verification.pdf` of 16 April 2026,
  which checks some steps on examples ("checked representative values",
  "spot-checked").  *Recommended:* cite the new exact script
  `python_project/src/verify_egan_construction.py` at the chain.  It checks
  the four links in exact arithmetic: Egan's rows from part 9 (sums, dot
  products, Gram matrices, $L_1 \cap L_2 = 2L_0$), the lattice they give
  (index $2^{12}$, half its Gram matrix even of determinant 1, none of the
  179,280 vectors of norm 4 of $L_0^3$ in it), and, at the standard scale in
  the paper's conventions, the closure of $\mathbb Z^3 \oplus L_L$ under
  $2(xy+yx)$ on all 378 pairs of basis elements and its failure under
  $xy+yx$.  It runs in about ten seconds, ALL PASS; it is indexed in
  `python_project/src/README.md` as not yet cited, and is committed before
  the paper cites it.  *Or:* no footnote at those four links (the chain
  keeps its new labels and place), with the April report as their record
  (the script then stays as a record only, or goes).
- **Appendix A and prompt 309.**  The draft of this plan proposed restoring
  in Appendix A the sentence "Egan's construction of $\Lambda$ (part 9 of
  [Baez2014]) simplifies Wilson's by eliminating the octonions", on the
  reading that item 3 had dropped a remark the author asked to keep at
  prompt 309 ("make sure not to drop the remark that this polarization
  directly relates to the Baez2014 'simplification', and the acknowledging
  references in that discussion").  The review showed that reading wrong.
  At prompt 309, Appendix A's "Relation to prior work" read "The reading of
  the Leech construction in terms of the preimages in $L$ of a
  complementary pair of such Lagrangians is David Calderbank's ...", with
  Silverman credited; the "simplifies" sentence was written later, at Step
  B (prompt 314).  Item 3, approved at prompt 327, removed it (the author
  disagreed with the "simplification" wording at prompt 184) and kept the
  relation ("construct a lattice isometric to $\Lambda$ from any
  polarization ... Egan's construction (Section 6) is of this form") and
  the acknowledging references (Silverman, Calderbank).  *Recommended:* no
  change; Appendix A stays as item 3 left it, for the author's exact read
  (prompt 329).  *Or:* Baez's word, attributed to Baez: "Egan's construction
  (Section 6), which Baez describes as simplifying Wilson's (part 9 of
  [Baez2014]), is of this form."
- **Two rulings of prompt 250, against the parent plan.**  The parent
  plan's item 4 asks for "both names of the product" and "the exact count"
  in the Jordan paragraph.  At prompt 250 the author declined the
  "doubled"/"quadrupled" distinction ("would distract from the fact that
  we're talking about the same thing") and approved "Egan's count of at
  least 17,280 Jordan rings" in the reproduction chain; that phrase left the
  chain at Step C (prompt 315), and the Jordan paragraph's "at least 17,280"
  dates from prompt 060.  The change of prompt 250 followed a review's
  reading (prompt log 236) that 17,280 counts a more constrained set and is
  not an exact count.  The parent plan finds otherwise: 17,280 is exact for
  Egan's construction (part 10's second Result 1; part 11), and for every
  sublattice $M \subset \mathbb O^3$ of the integral octonions isometric to
  the Leech lattice, $\mathbb Z^3 \oplus M$ is closed under $2(xy+yx)$ (F3,
  F7).  *Recommended:* keep the rulings of prompt 250: one name,
  "doubled", and "at least 17,280", the lower bound of part 10, for the
  Jordan-closed lattices, linked to the exact 17,280 of Egan's
  construction.  *Or:* the parent plan's version (both names; the exact
  count in the Jordan paragraph).

---

## 2. What changes

- **Section 6, the Baez–Egan block, in a new reading order.**
  - *"Baez–Egan (2014)"*: Egan's construction, credited to part 9, stated
    exactly in the form of part 11 as a display, with the condition $L_1
    \cap L_2 = 2L$, which the present text omits; $\Lambda = \Lambda(L\bar
    s, Ls)$ (Definition 2.1) and, with $\sigma$ applied in each block,
    $\sigma(\Lambda) = \Lambda(2\mathcal O, p\mathcal O)$; through an
    isometry of $L$ that carries Wilson's pair to $(L_1, L_2)$ (Appendix
    A), the untwisted product closes on $\Lambda(L_1, L_2)$ exactly when
    $L_1$ and $L_2$ are closed (Proposition 5.3): on $\sigma(\Lambda)$ it
    closes, which is Theorem 1.1, and on $\Lambda$ it does not; the count in
    lattice terms (decision (m)), all 17,280 lattices distinct, with Egan's
    comments of 11 and 13 December 2014; the pointer to Appendix A, which
    the author asked for at prompt 281, now without restating "ideal and
    subalgebra" (Appendix A says it once).
  - *The Jordan paragraph*: its lead sentence, "Baez–Egan (2014) [Baez2014]
    study closure properties of integer sublattices of the exceptional
    Jordan algebra $\mathfrak h_3(\mathbb O)$" (prompt 060), goes, and the
    heading now leads the paragraph on Egan's construction.  The Leech
    lattice $L_L$ of the Baez series lies in the integral octonions cubed
    and at the standard scale: "an off-diagonal Leech lattice $L_L$ at the
    standard scale, with entries in a maximal order of integral octonions
    (here $\mathcal O$), so that $L_L \cong \Lambda/\sqrt2$" (the present
    text writes $L_L \cong \Lambda$, inside the real octonions).  The closure under $2(xy+yx)$ and its failure under
    $xy+yx$ are credited to Egan, who found them (comment of 26 November
    2014).  "Off-diagonal Leech embeddings $L_L \subset \mathbb O^3$" becomes
    "such lattices $L_L$", tied to the lattices of Egan's construction with
    $\mathcal O$ in place of $L$.  "Doubled" and "at least 17,280" stay
    (prompt 250), and so does the sentence on Baez's open question (the
    author's cut of prompt 234), with "any such lattice" made "the
    27-dimensional lattice of any such $L_L$", since $L_L$ is now its
    nearest antecedent; F7 of the parent plan notes that the question has an
    elementary negative answer (pointed out only).  The reproduction chain
    moves to the end of this paragraph, before the Anatomy, where it stood
    before prompt 315 (the parent plan's §2 lists it after the Anatomy),
    with "the $E_8$ simple-roots data" and "$E_8^3$-orbit Leech embeddings"
    replaced by the wording verified at prompt 316, "Egan's simple roots of
    $E_8$ and the two sublattices built from them" (now without its symbols
    $L_1, L_2$, which would clash with the display's) and "sublattices of
    $L^3$ isometric to $\Lambda$", and "the lemma $L_1 \cap L_2 = 2L$"
    replaced by the review's "the lemma that these meet in twice that $E_8$
    lattice" (Egan's pair lies in Egan's $E_8$ lattice, not in $L$), and
    with the footnote of decision (n′); "likewise" (the author's, prompt 267)
    stays.
  - *Comparison*: deleted (decision (l)).
  - *Anatomy*: its opening at the right scale ("The Baez–Egan closure is on
    the 27-dimensional lattice $\mathbb Z^3 \oplus L_L$, not on $L_L$ alone.
    At the paper's scale, ..."); the rest unchanged.
- **Section 1.**  "explored Jordan-algebraic structures on $\Lambda$"
  becomes "explored Jordan-algebraic structures containing the Leech
  lattice" (the scale, and the closed object is the 27-dimensional
  lattice).
- **Bibliography.**  [Baez2014] links parts 11 and 12 to their Café posts,
  as for parts 9 and 10, and dates part 11 by its Café post (11 December
  2014; the old link, to the math.ucr.edu copy over http, was dated 12
  December).  New: [Egan2014], Egan's comments of 26 November 2014 on part 9
  and of 11 and 13 December 2014 on part 11, with their anchors, like
  [Silverman2016] and [Calderbank2023].

The credits of F4 are in place since item 3: Appendix A credits Lepowsky
and Meurman (with their note on Tits), Silverman and Calderbank, and Section
6 no longer names Silverman.  Untouched: the Wilson, Dixon, precursors,
history and Corradetti paragraphs of Section 6 (the precursors stay, as the
parent plan's §2 provides; the pointer to Appendix A leads the reader to
the 1982 construction); the Anatomy after its opening; Appendix A.

---

## 3. The edits, old → new

Each block is replaced exactly once (the edit script asserts it).

### 3.1. Section 1, the sentence on Baez and Egan

Old:

```latex
John C.~Baez and Greg Egan \cite{Baez2014} explored
Jordan-algebraic structures on~$\LL$ within the exceptional Jordan
```

New:

```latex
John C.~Baez and Greg Egan \cite{Baez2014} explored
Jordan-algebraic structures containing the Leech lattice within the exceptional Jordan
```

### 3.2. Section 6, Egan's construction, the Jordan setting and the reproduction chain (the Comparison paragraph deleted)

Old:

```latex
\textbf{Baez--Egan (2014)} \cite{Baez2014} study closure
properties of integer sublattices of the exceptional Jordan algebra
$\mathfrak{h}_3(\OO)$.  In parts~9--10 of the series, Baez considers
the 27-dimensional lattice $\ZZ^3 \oplus L_L \subset \mathfrak{h}_3(\OO)$,
with integers on the diagonal and an off-diagonal Leech sublattice
$L_L \cong \LL$ embedded as an octonionic triple.  This lattice is
\emph{not} closed under the Jordan product $x \circ y = xy + yx$,
but \emph{is} closed under the doubled product $2(xy + yx)$.
Egan further enumerates at least $17{,}280$ distinct off-diagonal
Leech embeddings $L_L \subset \OO^3$ for which the associated
27-dimensional lattice is closed under this doubled product
(see part~10 of~\cite{Baez2014} for Egan's computational counts).
Baez leaves open the question whether any such lattice is closed under
the undoubled Jordan product, and conjectures that the answer is
negative.

Egan's construction of~$\LL$ uses a pair of sublattices
$L_1, L_2 \subset L$, each isometric to~$\sqrt{2}\,L$.  We have
likewise independently reproduced the full Baez--Egan verification
chain: the $E_8$ simple-roots data, the lemma $L_1 \cap L_2 = 2L$, the
Leech-lattice theorem for that pair, the doubled-Jordan-product
closure, the
derivation of the lower bound $244{,}035{,}421$ on the number of
$E_8^3$-orbit Leech embeddings,\footnote{An orbit-stabilizer bound, $\lceil 3!\,\lvert W(E_8)\rvert^3 / \lvert\Co_0\rvert\rceil$; \scriptref{verify\_egan\_lower\_bound.py}.} and Egan's
count.\footnote{\scriptref{egan\_baez\_count.py}.}  The count factors as
$17{,}280 = 270 \times 64$: the $270$ is the number of sublattices $L_1 \subset L$
isometric to $\sqrt{2}\,L$ (part~11 of \cite{Baez2014}), equivalently
the number of $4$-dimensional maximal totally isotropic subspaces of
$L/2L \cong \FF_2^8$ under the quadratic form induced by the
$E_8$ Gram form (Appendix~\ref{app:mod2quotient}).  The $64$ is
the number of $L_2$ in this collection that meet a given $L_1$ only
in $2L$, equivalently the number of $4$-dimensional totally isotropic
subspaces complementary to a given one; Appendix~\ref{app:mod2quotient} reads one such
pair, $\sigma(L\bar s)/2L$ and $\sigma(Ls)/2L$, as ideal and
subalgebra under the product that $L/2L$ inherits from the standard
octonion product.

\textbf{Comparison.}  The product~$\star$ is bilinear on~$\RR^{24}$
itself, and its block-permutation symmetry~$S_3$ comes from the
routing of Definition~\ref{def:triple}.  Whether the Leech embedding $\sigma(\LL) \subset
\OO^3$ obtained by applying $\sigma$ blockwise overlaps with one of
the embeddings cataloged by Egan, or sits outside that family, is
a question this paper does not address.
```

New:

```latex
\textbf{Baez--Egan (2014).}  In part~9 of the series \emph{Integral
octonions}~\cite{Baez2014}, Greg Egan constructs the Leech lattice from
the $E_8$ lattice alone; in the form of part~11, for sublattices
$L_1, L_2 \subset L$, each isometric to~$\sqrt{2}\,L$, with
$L_1 \cap L_2 = 2L$, the lattice
\begin{equation}\label{eq:egan}
  \LL(L_1, L_2) \;:=\; \bigl\{ (x,y,z) \in L^3 :
    x{+}y,\, x{+}z,\, y{+}z \in L_1;\; x{+}y{+}z \in L_2 \bigr\}
\end{equation}
is isometric to~$\LL$.  Wilson's lattice is $\LL = \LL(L\bar s, Ls)$
(Definition~\ref{def:wilson-leech}), and, with~$\sigma$ applied in each
block, $\sigma(\LL) = \LL(2\MO, p\MO)$ by~\eqref{eq:multiples} and
Lemma~\ref{lem:sigma-L}.  Applied in each block in the same way, an
isometry~$g$ of~$L$ with $g(L\bar s) = L_1$ and $g(Ls) = L_2$, which
exists by Appendix~\ref{app:mod2quotient}, carries $\LL$ onto
$\LL(L_1, L_2)$ and~$\star_g$ onto the untwisted product~$\star_{\mathrm{id}}$
of Section~\ref{sec:census}.  By Proposition~\ref{prop:criterion}, then,
$\star_{\mathrm{id}}$ closes on $\LL(L_1, L_2)$ exactly when $L_1$
and~$L_2$ are closed under the standard octonion product: it closes
on~$\sigma(\LL)$, which is Theorem~\ref{thm:main}, and not on~$\LL$
(Remark~\ref{rem:Ls-not-closed}).  There are $270$ sublattices of~$L$
isometric to~$\sqrt{2}\,L$, and each meets exactly $64$ of them in~$2L$,
so Egan's construction gives $17{,}280 = 270 \times 64$ lattices, all
distinct (part~11, and Egan's comments there~\cite{Egan2014}).
Appendix~\ref{app:mod2quotient} reads these pairs modulo~$2$, with the
product that $L/2L$ inherits from the standard octonion product.

In parts~9 and~10, Baez considers the 27-dimensional lattice
$\ZZ^3 \oplus L_L \subset \mathfrak{h}_3(\OO)$ of the exceptional Jordan
algebra, with integers on the diagonal and an off-diagonal Leech
lattice~$L_L$ at the standard scale, with entries in a maximal order of
integral octonions (here~$\MO$), so that $L_L \cong \LL/\sqrt{2}$.  Egan found~\cite{Egan2014} that this lattice
is \emph{not} closed under the Jordan product $x \circ y = xy + yx$ but
\emph{is} closed under the doubled product $2(xy + yx)$.  Egan further
enumerates at least $17{,}280$ distinct such lattices~$L_L$ for which
the associated 27-dimensional lattice is closed under this doubled
product: those of Egan's construction with~$\MO$ in place of~$L$ (see
part~10 for Egan's computational counts).  Baez leaves open the question
whether the 27-dimensional lattice of any such~$L_L$ is closed under the
undoubled Jordan product, and conjectures that the answer is negative.  We have likewise independently reproduced the full
Baez--Egan verification chain: Egan's simple roots of~$E_8$ and the two
sublattices built from them, the lemma that these meet in twice that
$E_8$ lattice, the Leech-lattice theorem for that pair, the
doubled-Jordan-product closure,\footnote{\scriptref{verify\_egan\_construction.py}.}
the derivation of the lower bound $244{,}035{,}421$ on the number of
sublattices of~$L^3$ isometric to~$\LL$,\footnote{An orbit-stabilizer bound, $\lceil 3!\,\lvert W(E_8)\rvert^3 / \lvert\Co_0\rvert\rceil$; \scriptref{verify\_egan\_lower\_bound.py}.} and Egan's
count.\footnote{\scriptref{egan\_baez\_count.py}.}
```

### 3.3. Section 6, the opening of the Anatomy

Old:

```latex
\textbf{Anatomy of the Baez--Egan closure on $\LL$.}  The
Baez--Egan closure is on the 27-dimensional lattice
$\ZZ^3 \oplus \LL \subset \mathfrak{h}_3(\OO)$, not on $\LL$ alone.
For off-diagonal Hermitian matrices $M_u, M_v \in \mathfrak{h}_3(\OO)$
```

New:

```latex
\textbf{Anatomy of the Baez--Egan closure on $\LL$.}  The
Baez--Egan closure is on the 27-dimensional lattice
$\ZZ^3 \oplus L_L$, not on~$L_L$ alone.  At the paper's scale, for
off-diagonal Hermitian matrices $M_u, M_v \in \mathfrak{h}_3(\OO)$
```

### 3.4. Bibliography, [Baez2014]

Old:

```latex
\bibitem[Baez2014]{Baez2014}
J.\,C.~Baez (with contributions by G.~Egan, T.~Silverman and others),
\emph{Integral octonions} (parts 1--12),
The $n$-Category Caf\'e, 2013--2016.
The material discussed here appears in parts 9--12: Egan's Leech
construction (part~9, 23 November 2014) and the Jordan-ring count
(part~10, 1 December 2014),
\url{https://golem.ph.utexas.edu/category/2014/11/integral_octonions_part_9.html},
\url{https://golem.ph.utexas.edu/category/2014/12/integral_octonions_part_10.html};
the 270 rescaled copies of the first $E_8$ shell (part~11, 12 December
2014) and Tim Silverman's account of the integral octonions modulo~2
(part~12, 31 January 2016),
\url{http://math.ucr.edu/home/baez/octonions/integers/integers_11.html},
\url{http://math.ucr.edu/home/baez/octonions/integers/integers_12.html}.
```

New:

```latex
\bibitem[Baez2014]{Baez2014}
J.\,C.~Baez (with contributions by G.~Egan, T.~Silverman and others),
\emph{Integral octonions} (parts 1--12),
The $n$-Category Caf\'e, 2013--2016.
The material discussed here appears in parts 9--12: Egan's Leech
construction (part~9, 23 November 2014) and the Jordan-ring count
(part~10, 1 December 2014),
\url{https://golem.ph.utexas.edu/category/2014/11/integral_octonions_part_9.html},
\url{https://golem.ph.utexas.edu/category/2014/12/integral_octonions_part_10.html};
the 270 rescaled copies of the first $E_8$ shell (part~11, 11 December
2014) and Tim Silverman's account of the integral octonions modulo~2
(part~12, 31 January 2016),
\url{https://golem.ph.utexas.edu/category/2014/12/integral_octonions_part_11.html},
\url{https://golem.ph.utexas.edu/category/2016/01/integral_octonions_part_12.html}.
```

### 3.5. Bibliography, Egan's comments (new)

Old:

```latex
\bibitem[Elduque1996]{Elduque1996}
```

New:

```latex
\bibitem[Egan2014]{Egan2014}
G.~Egan,
comments of 26 November 2014 on Baez~\cite{Baez2014}, part~9, and of
11 and 13 December 2014 on part~11,
\url{https://golem.ph.utexas.edu/category/2014/11/integral_octonions_part_9.html#c047749},
\url{https://golem.ph.utexas.edu/category/2014/12/integral_octonions_part_11.html#c047862}
and
\url{https://golem.ph.utexas.edu/category/2014/12/integral_octonions_part_11.html#c047873}.

\bibitem[Elduque1996]{Elduque1996}
```

---

## 4. Checks done before this hand-over

- A scratch copy with the five edits builds with three `pdflatex` runs: no
  errors, no undefined references or citations, no hyperref warnings; 21
  pages; the known 1.3 pt overfull vbox and the space at the foot of page
  11.  Section 6 was read in the PDF (pages 12 and 13).
- `verify_egan_construction.py` (decision (n′)) runs ALL PASS in about ten
  seconds.
- The research (prompt 329): the Baez series read and quoted from the pages
  downloaded at prompt 325 and three more (part 8's thread, the
  math.ucr.edu copies of parts 9 and 10): Egan's construction and Egan's
  conditions (part 9, Egan's text; Egan's comment of 23 November 2014 on
  part 8); Baez's comment of 25 November 2014 that introduces the
  27-dimensional lattices and Egan's of 26 November 2014 that finds the
  closure under $2(xy+yx)$ and its failure under $xy+yx$; part 10's Theorem
  1′, its two statements labeled "Result 1", its Theorem 4 and the names
  of the products; the scale of $L_L$; the URLs and dates of parts 9 to 12.
  Every mathematical statement was recomputed by two independent programs
  in exact arithmetic: $\Lambda = \Lambda(L\bar s, Ls)$; $\sigma(\Lambda) =
  \Lambda(2\mathcal O, p\mathcal O)$; the transport $g^3(u \star_g v) =
  g^3u \star_{\mathrm{id}} g^3v$ on twelve isometries; closure of
  $\star_{\mathrm{id}}$ on $\Lambda(L_1, L_2)$ against both halves closed,
  over all 17,280 pairs (4,704 close); $\Lambda(L_1, L_1) = L_1^3$, not a
  Leech lattice; theta series of complementary pairs; $\Lambda/\sqrt2$ the
  Leech lattice; the Anatomy's $\Delta$, $\operatorname{tr}\Delta = 8\langle
  u, v\rangle$ and $R \in \Lambda$ for all placements of the conjugates;
  and the footnote's bound, 244,035,421.
- An independent review in three lenses (mathematics and sources; the paper
  as a reader sees it; fidelity, provenance and the author's rulings), each
  finding checked by a second agent that tried to refute it.  It found no
  mathematical error.  Applied from it: $\sigma$ and $g$ said to act in each
  block (the deleted Comparison paragraph was the only place that said so);
  $L_L$ placed in the integral octonions, where the count of 17,280 means
  something, and tied to Egan's construction; the general form of the
  construction credited to part 11; the chain's labels without symbols that
  clashed with the display; Egan's comments cited with anchors
  ([Egan2014]); "containing" in Section 1; the script of decision (n′); the
  draft's restoration of the "simplifies" sentence in Appendix A withdrawn,
  its premise being wrong (§1); the history and options of (l) and of the
  prompt-250 rulings; the records of §5.  Refuted: one point of wording.
  A second check of the revision, in three parts (the paper text, this
  plan's claims, and the new script, which a separate agent tried to make
  pass on false data), found nothing false.  Applied from it: "the
  27-dimensional lattice of any such $L_L$" in Baez's question; $\mathcal O$
  named once; the script now checks that its basis of $\Lambda$ is one,
  computes the index $2^{12}$ instead of arguing it, tests $2\bar a(ax) = x$,
  and gives the right reason for $aL = \mathcal O$ (definition (2.2)); the
  history of (n′) (prompts 082, 089 and 230) and its option; §2's account
  of the chain's labels.

---

## 5. On execution, and what stays for later items

`verify_egan_construction.py` is committed first; then the five edits are
applied to `paper/main.tex`, the paper is built and the changed pages read,
the script's index line in `python_project/src/README.md` names the paper's
§6 in place of "not yet cited", changelog item 25 records the change, and
the linter runs.  The records follow in
item 7, as the parent plan provides: a ledger row for $\sigma(\Lambda)$ as an
Egan lattice; the two registry entries [Baez2014] (`background_octonions.txt`
and `prior_art_orders_on_leech.txt`), for the Café links and part 11's date
(11 December 2014; `prior_art_orders_on_leech.txt` says 12 December), and in
`prior_art_orders_on_leech.txt` "a Leech sublattice L_L ≅ Λ off-diagonal" (at
the standard scale, $L_L \cong \Lambda/\sqrt2$); a registry entry for
[Egan2014]; a dated note in `verify_egan_lower_bound.py`, whose correction
note says "12 December is part 11's date"; and CURRENT_STATE's open item 1,
which says that (4.1) "bears on §6's sentence on Egan's catalog".  The
glossary names neither the Comparison nor the reproduction chain's old
labels (checked).
