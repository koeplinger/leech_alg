# 2026-10-22 — Plan: Appendix A around the closure criterion modulo 2

Created: 2026-10-22.  Updated: 2026-10-24.

**Status: executed.**  This plan is item 3 of the closed-form plan,
`2026-10-20_closed_form_plan.md`.  The author's decisions of §1: (i), (j)
and (k) as recommended, with one addition to (j) at the author's
direction: Silverman's comment that reads Egan's construction modulo 2
is cited as well.  So the last paragraph of Appendix A credits Silverman
with the identification of these sublattices with the Lagrangians (part
12), the reading of Egan's construction modulo 2, and the finding that
Wilson's pair is a polarization, and credits Calderbank in a sentence of
its own; the entry [Silverman2016] names both comments, with their
anchors.  With that addition, the edits of §3 are in `paper/main.tex`:
Appendix A (`app:mod2quotient`), Section 1's overview, Section 7's open
questions (`sec:conclusion`) and the bibliography.  The sentence of §3.2
is not in the paper: Section 6 gives the count in the form of item 4's
plan, `2026-10-22_section6_plan.md`.  `paper/v7_to_v8_changelog.md` states the
change under Appendix A; the records of §5 are part of item 7 of
the parent plan.  Drafted by Claude Opus 5.5 at the direction of Jens
Köplinger.

---

## 1. The decisions of this item

The parent plan asks decisions (i) to (k) again before this item.  Reading
Lepowsky and Meurman's paper in the original (the author's copy, prompt 326)
changed the facts behind (j).  The edits of §3 follow the recommendations.

- **(i) The words for the subspaces.**  On $\bar L = L/2L$ the paper uses
  the quadratic form $q$, half the norm modulo 2, and calls a subspace on
  which $q$ vanishes "totally isotropic", as Wilson's book does.
  Lepowsky–Meurman and Taylor call the same subspaces "totally singular",
  and in Taylor's book "totally isotropic" means instead that the associated
  bilinear form vanishes: read that way, there would be 2,295
  four-dimensional subspaces with 1,024 complements each, not 270 and 64.
  The paper's definition ("$q$ vanishes on it") settles which is meant.
  *Recommended:* keep "totally isotropic", with that definition, and no
  second name, per the standing rule against redundant terms (prompt 291);
  keep "Lagrangian" (a 4-dimensional such subspace) and "polarization" (an
  ordered pair of complementary Lagrangians), both defined in Appendix A;
  replace "plus type", "Witt index" and "Witt decomposition", which the paper
  never defines, by plain words, also in the one Section 6 sentence that
  edit 3.2 touches.  *Or,* as an exception to that rule for the author to
  rule on: add "or totally singular" once at the definition, for readers of
  Lepowsky–Meurman.
- **(j) Who is credited for what.**  The paper now credits the picture of
  $\bar L$ as an orthogonal space to Tim Silverman's part 12 of the Baez
  series (2016), in Section 6 and in Appendix A.  The sources say otherwise:
  - The form $q$ on $\bar L$ and the action of the Weyl group of $E_8$ on it
    as the full orthogonal group, with kernel $\{\pm 1\}$, are classical;
    Wilson's book states both (§3.12.4, p. 104).
  - Lepowsky and Meurman (1982, §2) prove that the lattices between $2L$ and
    $L$ isometric to $\sqrt2\,L$ are exactly the preimages of the
    Lagrangians (Proposition 2.1, with the uniqueness of $E_8$ invoked before
    Proposition 2.4), and build a lattice isometric to $\Lambda$ from any
    complementary pair of them (Theorem 2.2, Corollary 2.3), in the form
    $\{(x_1 + y, x_2 + y, x_3 + y)\}$ that the proof of Proposition 5.3
    also uses.  For a pair $(L_1, L_2)$ it is the lattice of Egan's
    construction, which Egan states by conditions (pairwise sums in $L_1$,
    the total in $L_2$).  That every sublattice of $L$ isometric to
    $\sqrt2\,L$ contains $2L$ is not in their §2; Appendix A adds the
    one-line reason, so that Section 6's "equivalently" stays covered.
  - Their introduction (p. 484) says that Tits's work (their refs. 18 and
    19, both 1980) "includes what amounts to the same construction of the
    Leech lattice from a different viewpoint".  The project has not read
    Tits's papers, so the paper does not cite them.
  - Silverman, in part 12 (31 January 2016), identified the 270 sublattices
    with these subspaces, saying he assumed it and had not proved it, and
    he calls the reduction of $E_8$ modulo 2 "a well-known thing that people
    do".  In a comment of 18 February 2016 he did not know who first built
    the Leech lattice from two copies of $E_8$ ("somebody (dunno who)").
    In a comment of 19 February 2016 he found that Wilson's $L\bar s$ and $Ls$ reduce to two disjoint
    such subspaces, that is, that Wilson's pair is a polarization.  The
    parent plan named his comment of 16 February 2016 instead, which reads
    Egan's construction modulo 2 with the parametrization $(x+p, x+q,
    x+r)$; Lepowsky and Meurman's Theorem 2.2 already has that form, so the
    plan credits it to them.
  - Calderbank's comment of 2 April 2023 relates Egan's construction to
    Tits and to Lepowsky–Meurman.  It names no corollary.

  *Recommended:* each credited once, in Appendix A's last paragraph, which
  keeps its heading "Relation to prior work":
  "The orthogonal space $(\bar L, q)$ is classical [Wilson's book, §3.12.4].
  Every sublattice $M$ of $L$ isometric to $\sqrt2\,L$ contains $2L$, since
  its dual lattice, $\frac12 M$, contains $L$.  Lepowsky and Meurman [§2]
  identify the lattices between $2L$ and $L$ that are isometric to
  $\sqrt2\,L$ with the preimages of the Lagrangians, and construct a lattice
  isometric to $\Lambda$ from any polarization; in their introduction they
  note that work of Tits contains what amounts to the same construction.  Egan's construction
  (Section 6) is of this form.  In the Baez series, Tim Silverman arrived at
  the identification of these sublattices with the Lagrangians independently
  (part 12) and found that Wilson's pair is a polarization [his comment of
  19 February 2016], and David Calderbank related Egan's construction to that
  of Lepowsky and Meurman [his comment]."  Wilson's book and Silverman's
  comment are new bibliography entries; Calderbank's gets the comment's
  anchor.  Section 6's clause "that mod-2 orthogonal-space picture is Tim
  Silverman's" becomes a pointer to Appendix A now, so that the paper does
  not contradict itself until item 4 rewrites Section 6.  *Or:* also cite
  Silverman's comment of 16 February 2016, a revised repost of his comment
  of 28 January 2016 in the part 11 thread, where he first posted his
  reading of Egan's construction modulo 2 with the parametrization $(x+p,
  x+q, x+r)$; or leave out the clause on Tits, so that the paper does not
  name him; or leave Silverman and Calderbank out.  The last would leave both
  uncredited in the paper's text: their credits replaced, at the v6 review,
  the credit to the Baez series that the author asked for at prompts 183 and
  184 (Baez and Egan stay credited in Section 6).
- **(k) Section 7's fourth question.**  It asks which polarizations "lift to
  a closed bilinear product on $\Lambda$".  For the twisted products of
  Section 5.6, the closure criterion answers it: every polarization is
  $(g(L\bar s)/2L, g(Ls)/2L)$ for an isometry $g$ of $L$, and $\star_g$
  closes exactly when both halves are subalgebras of $(\bar L, \bar\mu)$.
  Beyond twisted products the question is not well posed: "lift" is defined
  nowhere in the paper, and closed bilinear products exist trivially (the
  zero product), so nothing open remains to narrow it to.  *Recommended:*
  the answer in Appendix A, with the short reason why every polarization
  arises, and the question deleted.  *Or:* in addition, the number of
  polarizations that close (4,704 of the 17,280), which needs a new exact
  script in `python_project/src/`: it has been computed several times,
  among them by the review of 13 October 2026 (its §3.3, with the scripts of
  its frozen reproduction package), but no script in `python_project/` or
  `gap_project/` computes it.  This option adds a computed classification
  of subalgebra pairs, the kind of structure the author's ruling SF10 of
  prompt 284 leaves until a referee asks for it.

---

## 2. What changes

- **Appendix A, rebuilt** under the title "The closure criterion modulo 2",
  in five parts: the descent of the product to $\bar L$ and the criterion
  read modulo 2; Lagrangians and polarizations, defined; Wilson's pair and
  $\sigma$'s pair (Wilson's pair is already a polarization, and $\sigma$
  moves it to one whose halves are an ideal and a subalgebra); which
  polarizations close (decision (k)); the relation to prior work (decision
  (j)).  Gone: the lead on Lemmas 4.3 and 4.4; the term "mod-2 quotient";
  the sentence on structure constants (left over from the retired tables);
  display (A.3), now one sentence; "plus type", "Witt index" and "Witt
  decomposition"; the heading "Polarization of $\bar L$ by $\sigma$" (it
  suggested that $\sigma$ creates the polarization); the sentence that
  Egan's construction "simplifies Wilson's by eliminating the octonions"
  (the author disagreed with the "simplification" wording at prompt 184);
  Calderbank's credit for the preimage reading (it is Lepowsky and
  Meurman's); and the research-program sentence, approved at prompts 157
  and 291/292, whose last clause is the author's wording of prompt 281: it
  is replaced by its answer, as the parent plan's item 3 provides.  Kept:
  the descent, the definition of $q$, $V$ and $W$, the two-sided ideal with
  its footnote, the two script footnotes, and the heading "Relation to prior
  work" with its closing sentence on what the appendix adds, whose role the
  author set at prompts 184 and 281.  New words: "associated bilinear form",
  Wilson's and Lepowsky–Meurman's term, for the inner product modulo 2.
- **Section 1.**  "Appendix A reads two of the lemmas modulo 2" becomes
  "reads the closure criterion modulo 2".
- **Section 6.**  In the sentence on the 270, "the plus-type quadratic form
  ...; that mod-2 orthogonal-space picture is Tim Silverman's (part 12)"
  becomes "the quadratic form ... (Appendix A)" (decisions (i) and (j)).
- **Section 7.**  The fourth question is deleted (decision (k)).
- **Bibliography.**  Wilson's *The Finite Simple Groups* and Silverman's
  comment are added; Calderbank's comment gets its anchor.

Untouched: the rest of Section 6 (item 4); Appendix C, whose "Lagrangian
framing of Appendix A" stays accurate; the scripts.  Left out on purpose,
as new results beyond item 3: that no complement of $L\bar s/2L$ is a
subalgebra, so that $\sigma$ has to move both halves; which of the 270
Lagrangians are subalgebras or ideals (142, 16 right, 8 two-sided); the
count of decision (k).

---

## 3. The edits, old → new

Each block is replaced exactly once (the edit script asserts it).

### 3.1. Section 1, overview of the appendices

Old:

```latex
Section~\ref{sec:outlook} is an outlook.  Appendix~\ref{sec:appendix-tables}
reads two of the lemmas modulo~$2$, Appendix~\ref{app:history} traces
```

New:

```latex
Section~\ref{sec:outlook} is an outlook.  Appendix~\ref{sec:appendix-tables}
reads the closure criterion modulo~$2$, Appendix~\ref{app:history} traces
```

### 3.2. Section 6, the count read modulo 2

Old:

```latex
$L/2L \cong \FF_2^8$ under the plus-type quadratic form induced by the
$E_8$ Gram form; that mod-2 orthogonal-space picture is Tim
Silverman's (part~12 of \cite{Baez2014}).  The $64$ is
```

New:

```latex
$L/2L \cong \FF_2^8$ under the quadratic form induced by the
$E_8$ Gram form (Appendix~\ref{app:mod2quotient}).  The $64$ is
```

### 3.3. Section 7, the fourth open question

Old:

```latex
  $4$-cycle and for the exact one-half rates.
\item Which polarizations of~$L/2L$, beyond the one of
  Appendix~\ref{app:mod2quotient}, lift to a closed bilinear product
  on~$\LL$.
\end{itemize}
```

New:

```latex
  $4$-cycle and for the exact one-half rates.
\end{itemize}
```

### 3.4. Appendix A

Old:

```latex
\section{The mod-2-quotient reading of Lemmas~\ref{lem:L-sigmaLsbar}
and~\ref{lem:sigmaLs-closed}}\label{sec:appendix-tables}\label{app:mod2quotient}
% ===============================================================

Lemmas~\ref{lem:L-sigmaLsbar} and~\ref{lem:sigmaLs-closed} are
statements about the sublattices $\sigma(L\bar s)$ and~$\sigma(Ls)$
of~$L$; we now record what they say after reduction modulo~$2L$.

Because $L \cdot L \subseteq L$ (Lemma~\ref{lem:L-order}) and the
standard octonion product is $\ZZ$-bilinear, $(a + 2c) \cdot b \equiv a \cdot
b \pmod{2L}$ for $c \in L$, and likewise in the right factor.  The
product therefore descends to a well-defined $\FF_2$-bilinear product
$\bar\mu$ on the \emph{mod-2 quotient}
\begin{equation}\label{eq:disp14}
  \bar L := L/2L \cong \FF_2^8, \qquad
  \bar\mu \colon \bar L \times \bar L \to \bar L,
\end{equation}
the reduction of the order~$L$ modulo~$2$.  The
structure constants of $\bar\mu$ are the coefficients, in a $\ZZ$-basis
of~$L$, of the products of that basis, reduced modulo~$2$.

The sublattices $\sigma(L\bar s)$ and~$\sigma(Ls)$ contain~$2L$, since
$L\bar s \cap Ls = 2L$ (Section~\ref{subsec:wilson-leech}) and
$\sigma(2L) = 2L$ (Lemma~\ref{lem:sigma-L}).  They are therefore the \emph{full preimages}, under
$L \to \bar L$, of the $\FF_2$-subspaces
\begin{equation}\label{eq:disp15}
  V := \sigma(L\bar s)/2L = 2\MO/2L, \qquad
  W := \sigma(Ls)/2L = p\MO/2L,
\end{equation}
where the second equalities are~\eqref{eq:multiples}.  From
$\sigma(L\bar s) + \sigma(Ls) = L$
and $\sigma(L\bar s) \cap \sigma(Ls) = 2L$ one obtains
$\bar L = V \oplus W$, with $\dim V = \dim W = 4$ since both
sublattices have index~$16$ in~$L$.  Reducing the closure statements
modulo~$2L$ then gives
\begin{equation}\label{eq:disp16}
  \text{Lemma~\ref{lem:L-sigmaLsbar}} \iff
    \bar\mu(\bar L, V) \subseteq V, \qquad
  \text{Lemma~\ref{lem:sigmaLs-closed}} \iff
    \bar\mu(W, W) \subseteq W.
\end{equation}
That is: Lemma~\ref{lem:L-sigmaLsbar} says that $V$ is a
\emph{left ideal} of $\bar L$, and Lemma~\ref{lem:sigmaLs-closed}
that $W$ is a \emph{subalgebra} of it.  The descent and both
equivalences were checked in exact arithmetic.\footnote{\scriptref{verify\_mod2\_quotient.py}.}

\paragraph{Polarization of $\bar L$ by $\sigma$.}
Two further facts follow
from~\eqref{eq:disp15}.\footnote{\scriptref{verify\_polarization\_lagrangians.py}.}
First, $V$ is a \emph{two-sided} ideal of~$\bar L$: since
$L \subseteq \MO$ and $\MO$ is closed, $2\MO$ is a two-sided ideal
of~$L$ (Remark~\ref{rem:Ls-not-closed}).\footnote{Hence
$(\bar L, \bar\mu)$ is not an octonion algebra in the
composition-algebra sense: those are simple.}

Second, equip $\bar L$ with the quadratic form
\begin{equation}\label{eq:disp17}
  q \colon \bar L \to \FF_2, \qquad q(v + 2L) \;:=\; \Nrm(v)/2 \pmod 2,
\end{equation}
well-defined because $L$ is even
($\Nrm(v + 2u) \equiv \Nrm(v) \pmod 4$).  A subspace of~$\bar L$ is
\emph{totally isotropic} when $q$ vanishes on it.  For $L = D_8^+$, $q$
is the standard plus-type form on~$\FF_2^8$, of Witt index~$4$: its
maximal totally isotropic subspaces, the \emph{Lagrangians}, have
dimension~$4$.  Both $V$ and~$W$ are totally isotropic: every vector
of $2\MO$ and of~$p\MO$ has norm divisible by~$4$, since
$\Nrm(2x) = \Nrm(px) = 4\Nrm(x)$ with $\Nrm(x) \in \ZZ$ for
$x \in \MO$.  Being $4$-dimensional, each is a Lagrangian, and the
decomposition $\bar L = V \oplus W$ is a Witt decomposition of
$(\bar L, q)$ into a complementary pair of Lagrangians, a
\emph{polarization} of the orthogonal $\FF_2$-space.

\paragraph{Relation to prior work.}
The picture of $\bar L$ as an orthogonal space under~$q$, with its
Lagrangians, is Silverman's (part~12 of~\cite{Baez2014}); see
Section~\ref{sec:related}.  Egan's construction of~$\LL$ (part~9
of~\cite{Baez2014}) simplifies Wilson's by eliminating the octonions.
David Calderbank reads Egan's construction through a polarization
of~$\bar L$: its two sublattices of~$L$ are the preimages of a
complementary pair of Lagrangians (discussion thread of part~9
\cite{Calderbank2023}).  The Lagrangians $V$ and~$W$ are one such
pair; the present appendix adds the product~$\bar\mu$ inherited from
the standard octonion product, under which they are a two-sided ideal
and a subalgebra.  This opens a research
program: which polarizations of~$\bar L$ arise from the construction
of Section~\ref{sec:construction}, and which polarizations, beyond
the one constructed here, lift to a closed bilinear product on~$\LL$
itself.
```

New:

```latex
\section{The closure criterion modulo~\texorpdfstring{$2$}{2}}\label{sec:appendix-tables}\label{app:mod2quotient}
% ===============================================================

Because $L \cdot L \subseteq L$ (Lemma~\ref{lem:L-order}) and the
standard octonion product is $\ZZ$-bilinear, $(a + 2c) \cdot b \equiv a \cdot
b \pmod{2L}$ for $c \in L$, and likewise in the right factor.  The
product therefore descends to an $\FF_2$-bilinear product~$\bar\mu$ on
\begin{equation}\label{eq:disp14}
  \bar L := L/2L \cong \FF_2^8,
\end{equation}
the reduction of the order~$L$ modulo~$2$.  A lattice~$M$ with
$2L \subseteq M \subseteq L$ is the full preimage of~$M/2L$, so it is
closed under the standard octonion product exactly when $M/2L$ is a
subalgebra of~$(\bar L, \bar\mu)$.  For a twist~$g$
(Section~\ref{sec:census}), the lattices $g(L\bar s)$ and~$g(Ls)$
contain~$2L$, since $L\bar s \cap Ls = 2L$
(Section~\ref{subsec:wilson-leech}) and $g(2L) = 2L$.
Proposition~\ref{prop:criterion} therefore reads modulo~$2$: $\star_g$
closes on~$\LL$ if and only if $g(L\bar s)/2L$ and~$g(Ls)/2L$ are
subalgebras of~$(\bar L, \bar\mu)$.

\paragraph{Lagrangians and polarizations.}
The quadratic form
\begin{equation}\label{eq:disp17}
  q \colon \bar L \to \FF_2, \qquad q(v + 2L) \;:=\; \Nrm(v)/2 \pmod 2,
\end{equation}
is well defined because $L$ is even ($\Nrm(v + 2u) \equiv \Nrm(v)
\pmod 4$), and its associated bilinear form, the inner product
modulo~$2$, is nondegenerate because $L$ is unimodular.  A subspace
of~$\bar L$ is \emph{totally isotropic} when $q$ vanishes on it; it
then lies in its orthogonal complement, so its dimension is at
most~$4$.  The $4$-dimensional ones are the \emph{Lagrangians}, and a
\emph{polarization} of~$\bar L$ is an ordered pair of complementary
Lagrangians.

\paragraph{Wilson's pair and $\sigma$'s pair.}
The subspaces $L\bar s/2L$ and~$Ls/2L$ form a polarization.  They are
complementary, since $L\bar s + Ls = L$ and $L\bar s \cap Ls = 2L$;
each has dimension~$4$, since $[L : L\bar s] = [L : Ls] = 16$; and $q$
vanishes on both, since $\Nrm(\ell\bar s) = \Nrm(\ell s) =
2\Nrm(\ell) \in 4\ZZ$ for $\ell \in L$.  The isometry~$\sigma$ of~$L$
maps this polarization to
\begin{equation}\label{eq:disp15}
  V := \sigma(L\bar s)/2L = 2\MO/2L, \qquad
  W := \sigma(Ls)/2L = p\MO/2L,
\end{equation}
where the second equalities are~\eqref{eq:multiples}.\footnote{\scriptref{verify\_polarization\_lagrangians.py}.}
The identity twist does not close, since $Ls/2L$ is not a subalgebra
(Remark~\ref{rem:Ls-not-closed}).  The twist~$\sigma$ closes:
modulo~$2$, Lemma~\ref{lem:L-sigmaLsbar} says that $V$ is a left ideal
of~$(\bar L, \bar\mu)$, and Lemma~\ref{lem:sigmaLs-closed} that $W$ is
a subalgebra.\footnote{\scriptref{verify\_mod2\_quotient.py}.}
Moreover, $V$ is a two-sided ideal, since $2\MO$ is one in~$L$
(Remark~\ref{rem:Ls-not-closed}).\footnote{Hence
$(\bar L, \bar\mu)$ is not an octonion algebra in the
composition-algebra sense: those are simple.}

\paragraph{Which polarizations close.}
The Weyl group of~$E_8$, the isometry group of~$L$, acts on~$\bar L$
as the full orthogonal group of~$q$, with kernel~$\{\pm 1\}$
\cite[\S3.12.4]{WilsonBook2009}.  This group permutes the polarizations
transitively.  Indeed, the associated bilinear form pairs the two
halves $U_1$ and~$U_2$ of a polarization, so a basis $u_1, \ldots, u_4$ of~$U_1$
has a dual basis $u'_1, \ldots, u'_4$ in~$U_2$, and then
$q\bigl(\sum_i (a_i u_i + b_i u'_i)\bigr) = \sum_i a_i b_i$; the linear
map between two such bases preserves~$q$.  Hence every polarization is
$(g(L\bar s)/2L, g(Ls)/2L)$ for an isometry~$g$ of~$L$, and $\star_g$
closes on~$\LL$ exactly when both halves are subalgebras
of~$(\bar L, \bar\mu)$.

\paragraph{Relation to prior work.}
The orthogonal space~$(\bar L, q)$ is classical
\cite[\S3.12.4]{WilsonBook2009}.  Every sublattice~$M$ of~$L$
isometric to~$\sqrt2\,L$ contains~$2L$, since its dual lattice,
$\tfrac12 M$, contains~$L$.  Lepowsky and Meurman
\cite[\S2]{LepowskyMeurman1982} identify the lattices between~$2L$
and~$L$ that are isometric to~$\sqrt2\,L$ with the preimages of the
Lagrangians, and construct a lattice isometric to~$\LL$ from any
polarization; in their introduction they note that work of Tits
contains what amounts to the same construction.  Egan's construction (Section~\ref{sec:related}) is
of this form.  In the Baez series, Tim Silverman arrived at the
identification of these sublattices with the Lagrangians independently
(part~12 of~\cite{Baez2014}) and found that Wilson's pair is a
polarization \cite{Silverman2016}, and David Calderbank related Egan's
construction to that of Lepowsky and Meurman \cite{Calderbank2023}.
The present appendix adds the product~$\bar\mu$ inherited from the
standard octonion product, under which $V$ and~$W$ are a two-sided
ideal and a subalgebra.
```

### 3.5. Bibliography, Calderbank's comment

Old:

```latex
comment of 2 April 2023 on Baez~\cite{Baez2014}, part~9,
\url{https://golem.ph.utexas.edu/category/2014/11/integral_octonions_part_9.html}.
```

New:

```latex
comment of 2 April 2023 on Baez~\cite{Baez2014}, part~9,
\url{https://golem.ph.utexas.edu/category/2014/11/integral_octonions_part_9.html#c062182}.
```

### 3.6. Bibliography, Silverman's comment (new)

Old:

```latex
\bibitem[SmithVojtechovsky2022]{SmithVojtechovsky2022}
```

New:

```latex
\bibitem[Silverman2016]{Silverman2016}
T.~Silverman,
comment of 19 February 2016 on Baez~\cite{Baez2014}, part~12,
\url{https://golem.ph.utexas.edu/category/2016/01/integral_octonions_part_12.html#c050282}.

\bibitem[SmithVojtechovsky2022]{SmithVojtechovsky2022}
```

### 3.7. Bibliography, Wilson's book (new)

Old:

```latex
J.~Algebra \textbf{322} (2009), 2186--2190.

\bibitem[repo]{repo}
```

New:

```latex
J.~Algebra \textbf{322} (2009), 2186--2190.

\bibitem[WilsonBook2009]{WilsonBook2009}
R.\,A.~Wilson,
\emph{The Finite Simple Groups},
Graduate Texts in Mathematics, vol.~251, Springer, 2009.

\bibitem[repo]{repo}
```

---

## 4. Checks done before this hand-over

- A scratch copy with the seven edits builds with three `pdflatex` runs: no
  errors, no undefined references or citations, no hyperref warnings; 21
  pages; the known 1.3 pt overfull vbox and the stretched page 11 of item 2,
  layout only.  Appendix A was read in the PDF (pages 15 and 16).
- The sources, read and quoted (prompts 325 and 326): Lepowsky–Meurman §2
  and their reference list, in the author's copy; Wilson's book, §3.12.4,
  p. 104, in the publisher's preview and in the 2007 draft; Silverman's part
  12 and his comments of 28 January, 16, 18 and 19 February 2016; Calderbank's
  comment, with its anchor.
- Every mathematical statement of the new Appendix A except one was
  recomputed by two independent programs in exact arithmetic (prompt 325):
  the descent; Wilson's
  pair and $\sigma$'s pair (complementary, $q$ zero on each half, $V$ a
  two-sided ideal, $W$ a subalgebra, $Ls/2L$ not one); the Weyl group of
  $E_8$ onto the orthogonal group of $q$ (order 348,364,800) with kernel
  $\{\pm 1\}$; transitivity on the 17,280 ordered polarizations; the
  criterion modulo 2 against direct closure on all 5,040 permutations, the
  120 root reflections, 2,000 random elements of the Weyl group and all
  17,280 polarizations (4,704 close).  The exception, the containment step
  added after the review (every sublattice of $L$ isometric to $\sqrt2\,L$
  contains $2L$), was checked by two programs of the review: exhaustively,
  through the 2,160 vectors of norm 4 of $L$ and the 270 preimages, and on
  400 sampled sublattices.
- An independent review in three lenses (mathematics and sources; the paper
  as a reader sees it; fidelity and credits), each finding checked by a
  second agent that tried to refute it.  It found no mathematical error and
  confirmed that every source supports its credit once the credits were
  narrowed.  Applied from it: Lepowsky–Meurman's credit restricted to
  lattices between $2L$ and $L$, with the one-line containment step;
  Silverman's credit stated as the identification and Wilson's pair;
  "a lattice isometric to $\Lambda$"; Wilson's book cited once for the
  action, once for the space; the heading "Relation to prior work" and its
  closing sentence kept; "associated bilinear form" for "polar form", next
  to "polarization"; "the Weyl group of $E_8$" in words, since $W$ names a
  subspace there; "plus-type" dropped from Section 6; the title's math in
  `\texorpdfstring`; "vol.~251"; decisions (i) to (k) completed with the
  facts above; the records of §5.  Refuted: two points, as taste or as
  already provided by the parent plan.  A second check of the revision found
  one slip in the paper text (the note on Tits is in Lepowsky and Meurman's
  introduction, not their §2) and six points in this plan's prose, all
  applied.

---

## 5. On execution, and what stays for later items

The seven edits are applied to `paper/main.tex`; the paper is built and the
changed pages read; the changelog records the change; the linter runs.
The records follow in item 7, as the parent plan provides:

- the glossary: the entry "Polarization of L/2L" and the words of decision
  (i);
- the reference registry: Lepowsky–Meurman read (Proposition 2.1, Theorem
  2.2, Corollary 2.3, and their note on Tits); Wilson's book; Silverman's
  comment; Calderbank's entry, its heading "Complementary-Lagrangians
  reading" and its sentence on what Appendix A credits him for; and the
  [Baez2014] sentence of `prior_art_orders_on_leech.txt` on where the paper
  uses parts 11 and 12 (part 12 is then cited in Appendix A only);
- dated notes in `verify_mod2_quotient.py`,
  `verify_polarization_lagrangians.py`, `verify_product_span_structure.py`
  and `verify_discovery1_W_subalgebra.py`, which name the old paragraphs or
  "the mod-2-quotient reading", and the
  index line of `verify_polarization_lagrangians.py` in
  `python_project/src/README.md` ("Polarization paragraph", "Witt
  decomposition");
- README's open-question bullet on polarizations, removed with Section 7's
  fourth question;
- CURRENT_STATE: the description of Appendix A in the paper's structure (its
  title and "Lagrangian-polarization framing"), and the lead and item 1 of
  "What remains open", which say the general question stays open;
- the ledger, on the polarization answer, with claim 017.
