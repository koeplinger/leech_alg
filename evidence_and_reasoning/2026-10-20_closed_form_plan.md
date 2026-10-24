# 2026-10-20 — Plan: the closed-form v8 (the closure criterion, Egan's construction, Section 5.1)

Created: 2026-10-20.  Updated: 2026-10-24.

**Status: executed.**  The seven items of §3 are done, items 1 to 5 each
through its own plan; the body is the plan as engaged.  In this plan "§"
refers to the plan's own sections and "Section" to the paper's.  The
author's decisions of §5, as taken: (a) item 1 first; (b) the witnesses
as a short table of explicit vectors, so that Section 5.1 defines no
symbols and Section 5.2 keeps its definition of $\varepsilon_i$ and
$\omega$, and the exhaustive value set of Remark 5.1 in the new
appendix; (c) the new appendix is Appendix D; (d) to (f) as
recommended, and "is given here" for the tension around F7; (g) no credit
note in the paper for results found on the AI side of the work: the work
is a human–AI collaboration and its results are joint, while the records
trace where each result came from; (h) Appendix A rebuilt around the
criterion modulo 2; (i) to (k) as recommended, with the credits that
Lepowsky and Meurman's paper, read in the original, supports, and with
Silverman's comment that reads Egan's construction modulo 2 cited as
well; (l) and (m) as recommended, with the exact script
`verify_egan_construction.py` cited at the reproduction chain and two
earlier rulings of the author kept against item 4 (one name of the
product, and "at least 17,280"); (n) item 5, by hand and more simply than
§3 has it, without F6(a)'s lemma, whose literature check is therefore not
run, and without the X-product, which does not carry the step; (o) F6(b)
to (e) stay out of v8; (p) the abstract and Section 7 each name the
criterion in one sentence; (q) as recommended, with the corrections made
in place: the Hermes supplement and the W-concrete review state $w_0
\equiv s^2 \pmod{2L}$ (F8), and the stabilizer statements of the ledger,
CURRENT_STATE and research_result name the group $\mathrm{Aut}(\mathbb O,
\cdot_\sigma)$ in which they hold.  The items, with their own plans and the
parts of the paper they changed (`paper/v7_to_v8_changelog.md` states the
changes part by part):

- Item 1: `2026-10-21_plan.md`; Section 5.1
  (`sec:identities`) and Appendix D (`app:samples`).
- Item 2: `2026-10-21_criterion_plan.md`; Proposition
  5.3 (`prop:criterion`) in Section 5.6 (`sec:census`).
- Item 3: `2026-10-22_appendix_a_plan.md`; Appendix A
  (`app:mod2quotient`), and Section 7 (`sec:conclusion`) with three open
  questions.
- Item 4: `2026-10-22_section6_plan.md`; Section 6
  (`sec:related`), with Egan's construction as display (6.2) (`eq:egan`).
- Item 5: `2026-10-23_span_plan.md`; Proposition 5.4
  (`prop:span`).
- Item 6: the abstract and Section 7.
- Item 7, the records: CURRENT_STATE, the indexing READMEs,
  research_result, the ledger (the paper's labels in its "Recorded in"
  column; claim 017 for every twist; rows 023, $\sigma(\Lambda)$ as one of
  Egan's lattices, and 024, the witnesses), the glossary, the reference
  registry (F8; Lepowsky and Meurman as read in the original; the entries
  [Egan2014], [Silverman2016] and [WilsonBook2009], and, for the
  attribution of the X-product, [Dixon1994] and
  [CederwallPreitschopf1995]), and the docstrings of the cited scripts, in
  the paper's numbering.

Open, outside the items, for the author: the keyword "computer-assisted
proof" and the MSC code 68V05, to be reconsidered for the closed-form
proofs; the paragraph after the census table, which restates three facts
of the table.  Drafted by Claude Opus 5.5 at the direction of Jens
Köplinger, and reviewed before hand-over by four independent readers.

---

## 0. Goal, method, provenance

**Goal** (prompt 317): a compact v8 that contains all of the new closed-form
proofs, worked one finding at a time, with an item that is not trivial split
into its own plan, so that the author stays ahead of the changes, controls
them, and sees what remains to be learned.  **Roles of the sections**
(prompt 316): Section 5 holds the structural facts; Section 6 the tally of
prior work; Appendix A what is interesting but too far from the main line for
Section 5.  Section 7 is a summary followed by the open questions (prompt
309).  **Ruling (5) of prompt 317:** Section 5.1 keeps only general
properties, each failure easily verifiable by the reader, and the sampled
table moves to an appendix.

**Evidence.**  Two read-only workflows at prompt 317.  The first made an
inventory of the records and side documents and ran five verifications
(Egan's construction; the criterion and the polarizations; the census; terms
and attribution; the span index), each followed by an independent attempt to
refute it.  The second searched for hand-checkable witnesses for Section 5.1,
with an independent check.  All arithmetic was exact, and the core statements
were also re-derived first-hand.  Prompt log 317 names the workflows, the
checks behind each count, and the scratch files.

**Provenance.**  The closure criterion (F1) is Theorem A of the
external review of 13 October 2026
(`paper/reviews/2026-10-13_claude_opus_5_5_review_v8.md`, lines 125–169).  The
review states it in lattice form for twists $g \in W(E_8)$: $\Lambda \star_g
\Lambda \subseteq \Lambda$ iff $L\bar s \cdot_g L\bar s \subseteq L\bar s$ and
$Ls \cdot_g Ls \subseteq Ls$.  It proves it by reduction modulo $2L$, using
only $g(L) = L$, and reads it modulo 2 in lines 148–150.  Prompt 287 verified
it, and ledger claim 017 records it in lattice form.  The same review answered
Appendix A's question on which polarizations lift, for these twists (lines
167–169; its §4, item 3; its §6.1(d)), and computed in its §3.1 and §3.3 the
transitivity on the 17,280 ordered pairs, the 4,704 closing pairs, the 142
subalgebras, 16 right ideals and 8 two-sided ideals among the 270 subspaces,
and the span of $\star$ as the diagonal of $Ls$ plus $2\Lambda$, shared by all
588 closing odd permutations.  Prompt 315 stated the criterion modulo 2 for
every twist $g$ with $g(L) = L$, isometric or not.  Prompt 317 found that the
routing identities behind it hold exactly over $\mathbb Z$, not only modulo
$2L$, so the criterion holds for every bilinear $\mu$ with $\mu(L, L)
\subseteq L$.

---

## 1. The findings

### F1. The closure criterion, in lattice form

Let $\mu$ be a $\mathbb Z$-bilinear product on $\mathbb R^8$ with $\mu(L, L)
\subseteq L$, and $\star_\mu$ the triple product of Definition 3.3 built from
$\mu$.  Then $\star_\mu$ closes on $\Lambda$ if and only if $\mu(L\bar s,
L\bar s) \subseteq L\bar s$ and $\mu(Ls, Ls) \subseteq Ls$.  More finely,
Wilson's condition (1) always holds, condition (2) holds for all products iff
the first inclusion holds, and condition (3) iff the second.

*Proof* (about ten lines).  Put $a = x+y$, $b = x+z$, $w = x+y+z$.  Since $2L
\subseteq L\bar s$ (Section 2.3), this inverts exactly:
$\Lambda = \{(a+b-w,\ w-b,\ w-a) : a, b \in L\bar s,\ w \in Ls\}$.  For two
such vectors the routing gives, for any bilinear $\mu$,
$P_1+P_2+P_3 = \mu(w,w')$ and
$P_1+P_2 = \mu(a,b') + \mu(b,a') - \mu(a,a') + 2\bigl(\mu(b,b') - \mu(w,b') +
\mu(w,w') - \mu(b,w')\bigr)$, with the $\mathbb Z_3$-images for the other two
pairs.  The doubled terms lie in $2L \subseteq L\bar s$, which gives
sufficiency.  Necessity: $(a, 0, -a)$ and $(w, w, -w)$ lie in $\Lambda$, and
their products have $P_1 + P_2 = -\mu(a, a')$ and block sum $\mu(w, w')$.

For a twist $g$ with $g(L) = L$ and $\mu = \cdot_g$: $\star_g$ closes on
$\Lambda$ iff $L\bar s$ and $Ls$ are closed under $\cdot_g$, that is, iff
$g(L\bar s)$ and $g(Ls)$ are closed under the standard product.  Both contain
$2L$, so the statement reads modulo 2 as: iff $g(L\bar s)/2L$ and $g(Ls)/2L$
are subalgebras of $(L/2L, \bar\mu)$, which is the review's reading in its
lines 148–150 and the form of prompt 315.

*Instances.*  $g = \sigma$: $\sigma(L\bar s) = 2\mathcal O$ is closed (closure
of $\mathcal O$) and $\sigma(Ls) = p\mathcal O$ is closed (Lemma 4.4), which is
Theorem 1.1 again; Lemma 4.3, the stronger left-ideal property, is not needed
on this route.  $g = \mathrm{id}$: $Ls$ is not closed (Remark 4.5), so
condition (3) fails, as Section 5.6 already shows.  $L\bar s$ is closed under
the standard product (indeed $L\bar s \cdot L\bar s \subseteq 2L$), which the
paper does not state; so condition (3) is the only failure, as Appendix C's
failure analysis reports.

*Checks.*  Direct closure agrees on all 5,040 permutations of $e_1, \ldots,
e_7$ (the census of Section 5.6, row by row), all 120 root reflections (56
close), all 17,280 ordered pairs of F4, several hundred twists in
$\mathrm{GL}(L)$ that are not isometries, and the 35,280 permutations that
move $e_0$ (12,348 close).  The repository script
`verify_twist_span_index_census.py`, which Section 5.6 cites, has checked the
criterion on the 5,040 permutations since prompt 287.

*Scope.*  Outside it: twists with $L \cdot_g L \not\subseteq L$, all of which
have $g(L) \ne L$, and which Section 5.6's definition of a twist allows.  A
twist with $g(L) \ne L$ but $L \cdot_g L \subseteq L$ is covered by the first
form: $g = 2\cdot\mathrm{id}$ gives $\star_g = 2\star_{\mathrm{id}}$, which
closes, as the criterion says.

*In the paper today.*  Section 5.6 says "No such description is known for the
other types that close in part" (main.tex 916–918).  The criterion gives a
test for every type, not a description like the 3-cycle rule, and no uniform
reason for the 4-cycles (of the 210, $L\bar s$ fails for 168, $Ls$ for 147,
both for 105), so that sentence is narrowed rather than refuted.  Section 7's
third question asks "which twists other than those of Section 5.6 close on
$\Lambda$" (1110–1113), which the criterion answers for every $g$ with $g(L) =
L$.

### F2. Egan's construction, exactly; $\Lambda$ and $\sigma(\Lambda)$ are two of its lattices

Egan's construction (Baez, part 9, quoting Egan; parts 10 and 11; Egan's
comments): for an ordered pair $(L_1, L_2)$ of sublattices of $L$, each
isometric to $\sqrt2\,L$, with $L_1 \cap L_2 = 2L$, the lattice
$\Lambda(L_1, L_2) := \{(x,y,z) \in L^3 : x+y,\ x+z,\ y+z \in L_1;\ x+y+z \in
L_2\}$ is a Leech lattice at the paper's scale.  There are 270 such
sublattices and $17{,}280 = 270 \times 64$ ordered pairs, giving 17,280
distinct lattices (Egan's comment of 13 December 2014; part 11).  The
condition $L_1 \cap L_2 = 2L$ carries the construction: with $L_1 = L_2$ the
set contains $(v, 0, 0)$ for every $v \in L_1$ of norm 4.  Section 6's
statement of the construction (main.tex 993–994) omits it; the condition
appears there only as a lemma of the reproduction chain (996) and in the
count of the 64 (1007–1009).  Baez's one-paragraph summary theorem at the top
of part 9 omits Egan's condition that the two images of each simple root
$r_i$ sum to $r_i$, and is false as literally stated, so the paper states
Egan's condition, not that summary.

Definition 2.1 is $\Lambda(L\bar s, Ls)$, as Egan himself notes.  By (4.1) and
Lemma 4.1, $\sigma(\Lambda) = \Lambda(2\mathcal O, p\mathcal O)$.  So
$\sigma(\Lambda)$ is one of Egan's 17,280 lattices: the Comparison paragraph's
open question (main.tex 1017–1020) is answered by the paper's own (4.1).  The
blockwise map $\sigma \oplus \sigma \oplus \sigma$ carries $\star$ to the
untwisted product, so by F1 the untwisted product closes on $\Lambda(L_1,
L_2)$ iff $L_1$ and $L_2$ are closed; it closes on $\sigma(\Lambda)$ (Theorem
1.1, through Section 4.2) and not on Wilson's $\Lambda$ (Remark 4.5).

### F3. Section 6's Jordan and Anatomy paragraphs, and Section 1's sentence: precision

- Baez and Egan's $L_L \subset \mathbb O^3$ is the Leech lattice at its
  standard scale (minimum 4); the paper's $\Lambda$ is the $\sqrt2$-scaled
  copy, so "$L_L \cong \Lambda$" (main.tex 982) holds only up to $\sqrt2$.
- The Anatomy paragraph (main.tex 1022–1024) calls $\mathbb Z^3 \oplus
  \Lambda$ the Baez–Egan lattice; it computes the analogue at the paper's
  scale ($\operatorname{tr}\Delta = 64$, span index $2^{25}$).  The same scale
  point.
- 17,280 is the exact number of Egan's lattices, and all of them are closed
  under $2(xy + yx)$ (part 10, Theorem 4).  "At least 17,280" (main.tex
  985–988) is the wording of part 10's first "Result 1", a lower bound on
  Jordan-closed sublattices.
- Part 9 calls $2(xy+yx)$ the "doubled" product; part 10, which the sentence
  cites, calls $xy + yx$ "doubled" and $2(xy+yx)$ "quadrupled".
- Baez introduced the 27-dimensional lattice (comment of 25 November 2014);
  the closure, and its failure for $xy + yx$, are Egan's findings (26 November
  2014).
- Section 1 (main.tex 144–147), "Jordan-algebraic structures on $\Lambda$ ...
  a doubled Jordan product", has the same scale and naming points.
- Accurate as it stands: the report of Baez's open question (main.tex
  989–991).

### F4. The reading modulo 2: terms, credits, and the polarization question

*Terms.*  Two conventions coexist in characteristic 2.  "Totally singular
subspace" ($q$ vanishes) is the usage of Taylor, Cameron and Shimakura.
Silverman's part 12 defines "singular" subspaces ($q$ vanishes) against
"isotropic" ones (the polar form vanishes), but his tables call the $q$-zero
points "isotropic", his subspaces are "maximal flats", and he never writes
"totally singular".  Elman–Karpenko–Merkurjev and Wilson's book say "totally
isotropic subspace" for $q$, the paper's usage; Baez and Calderbank say
"maximal isotropic".  In Elman–Karpenko–Merkurjev a *form* is totally singular
when its polar form vanishes.  Read with the polar form, the counts would be
2,295 four-dimensional subspaces with 1,024 complements each, not 270 and 64
(computed).  The paper uses "totally isotropic" in Section 6 (main.tex 1004,
1009) before Appendix A defines it (1242–1243); "plus type", "Witt
decomposition" and "lift" are not defined.

*Credits.*  The orthogonal space $(L/2L, q)$ of plus type, on which $W(E_8)$
acts as the full orthogonal group with kernel $\{\pm1\}$, is older than the
blog series: Wilson, *The Finite Simple Groups*, §3.12.4 (2009; draft of 23
July 2007).  Silverman's part 12 identifies the 270 sublattices with its
maximal flats, and his comment of 16 February 2016 reads Egan's construction
modulo 2, with the parametrization $(x+p, x+q, x+r)$ that the criterion uses.
The paper credits Silverman twice (main.tex 1006–1007 and 1255–1256) and, for
the picture itself, the wrong source.  The construction of $\Lambda$ from a
complementary pair is also older: Calderbank's comment of 2 April 2023, which
Appendix A cites, traces it to Tits (1980), Lepowsky and Meurman (1982; their
Corollary 2.3, per Shimakura), Conway and Sloane (1982) and Borcherds' thesis
(1984).  Section 6's precursors paragraph cites Lepowsky–Meurman ($\Lambda$ as
a sublattice of $E_8^3$ containing $2E_8^3$) and Conway–Sloane (a Turyn-style
construction from three copies of $E_8$) without linking them to Egan's
construction (main.tex 1061–1066); Tits and Borcherds are not cited.  Neither
Tits 1980 nor Lepowsky–Meurman's Corollary 2.3 has been read in the project.
A product modulo 2 on a unital maximal order, an octonion algebra over
$\mathbb F_2$, comes up in the part-12 thread (Baez, Silverman; Wilson's book,
§4.4); the paper's $\bar\mu$ is the reduction of the nonunital order $L$.

*The polarization question.*  Every ordered complementary pair of
four-dimensional totally singular subspaces of $L/2L$ is $(g(L\bar s)/2L,
g(Ls)/2L)$ for an isometry $g$ of $L$: $W(E_8)$ maps onto the full orthogonal
group of $q$, and Witt's theorem gives transitivity.  By F1, $\star_g$ closes
exactly when both halves are subalgebras: 4,704 of the 17,280 pairs.  The
review of 13 October computed both (its §3.3); prompt 317 recomputed them
exhaustively, three times.  This answers Section 7's fourth question and
Appendix A's research program within twisted triple products.  For bilinear
products in general the question is not well posed, since trivially closed
bilinear products exist, and "lift" is undefined.

*Wilson's pair and $\sigma$'s pair.*  Appendix A's paragraph "Polarization of
$\bar L$ by $\sigma$" suggests that $\sigma$ creates the polarization.
Wilson's own pair $(L\bar s/2L, Ls/2L)$ is already a complementary pair of
totally singular subspaces.  $\sigma$ moves both halves (Remark 4.5; modulo 2,
$L\bar s/2L$ and $\sigma(L\bar s)/2L$ share only a 3-dimensional subspace, as
do $Ls/2L$ and $\sigma(Ls)/2L$).  The first half keeps its property:
$L\bar s/2L$ and $V = 2\mathcal O/2L$ are both two-sided ideals on which
$\bar\mu$ vanishes.  The second gains one: $Ls/2L$ is not a subalgebra, $W =
p\mathcal O/2L$ is.  Moving the first half is forced: no complement of
$L\bar s/2L$ is a subalgebra (all computed), so by F1 no twist $g$ with $g(L)
= L$ and $g(L\bar s) = L\bar s$ closes.  The sentence on structure constants
(main.tex 1199–1201) is left over from the retired tables.

### F5. Section 5.1: every failure has a hand-checkable witness (ruling 5)

Two observations make the failures checkable by hand.  (i) For vectors with
one nonzero block, Definition 3.3 gives $[a]_\alpha \star [b]_\beta = [a
\cdot_\sigma b]_{\alpha\circ\beta}$, where $\alpha\circ\beta = \alpha$ for
$\alpha = \beta$ and otherwise the third index.  This routing is commutative
but not associative: $(1\circ1)\circ2 = 3$, $1\circ(1\circ2) = 2$.  (ii) With
$\varepsilon_1 = (e_0, 0, 0)$, $\varepsilon_2$, $\varepsilon_3$ and $\omega =
\varepsilon_1 + \varepsilon_2 + \varepsilon_3$ (today defined in Section 5.2,
after 5.1), vectors whose blocks are real multiples of $e_0$ multiply by
$\varepsilon_\alpha \star \varepsilon_\beta = \varepsilon_{\alpha\circ\beta}$,
so the arithmetic is on three numbers.  Witnesses, each verified three times
independently in exact arithmetic:

| Property | Witness | The two sides |
|---|---|---|
| associativity | $\varepsilon_1, \varepsilon_2, \varepsilon_3$ | $(\varepsilon_1\star\varepsilon_2)\star\varepsilon_3 = \varepsilon_3$, $\varepsilon_1\star(\varepsilon_2\star\varepsilon_3) = \varepsilon_1$ |
| left alternativity | $\varepsilon_1, \varepsilon_2$ | $(\varepsilon_1\star\varepsilon_1)\star\varepsilon_2 = \varepsilon_3$, $\varepsilon_1\star(\varepsilon_1\star\varepsilon_2) = \varepsilon_2$ |
| right alternativity | $\varepsilon_2, \varepsilon_1$ | $(\varepsilon_2\star\varepsilon_1)\star\varepsilon_1 = \varepsilon_2$, $\varepsilon_2\star(\varepsilon_1\star\varepsilon_1) = \varepsilon_3$ |
| quartic power-associativity | $x = \varepsilon_1 - \varepsilon_2$ | $(x\star x)\star(x\star x) = -3\,x\star x$, $((x\star x)\star x)\star x = 3\,x\star x$ (the cube identity holds at $x$) |
| norm multiplicativity | $\omega$, $\varepsilon_1 - \varepsilon_2$ | $\omega \star (\varepsilon_1 - \varepsilon_2) = 0$: zero divisors |
| commutativity | $(e_3,0,0)$, $(e_4,0,0)$ | the two products are $(e_6,0,0)$ and $(-e_6,0,0)$ |
| symmetric composition | $\varepsilon_1$, $a = (e_3,0,0)$ | $\langle \varepsilon_1\star a, a\rangle = 1$, $\langle \varepsilon_1, a\star a\rangle = -1$ |
| flexibility | $u = (e_3, e_4, 0)$, $v = \varepsilon_1$ | the two sides differ by $(0, \mp 2e_6, \pm 2e_6)$ |
| cube power-associativity | $u = (e_0, e_3, e_4)$ | the two sides differ by $(0, \mp 4e_6, \pm 4e_6)$ |
| multiplicative identity | Remark 5.2, unchanged | |

$\sigma$ fixes $e_0, e_3, e_4, e_6$, so the reader uses table (2.1) directly.
Multiplying each witness by 4 puts it in $(2L)^3 \subseteq \Lambda$ ($4e_i \in
2L$); the identities are homogeneous, so the failures hold on the order.
Remark 5.1's statement that $\star$ does not preserve the minimal shell has the
one-line witness $u = (2\lambda, 0, 0)$ with $\lambda = e_3 + e_4$, a Type 1
minimal vector (Section 2.3): $u \star u = (-8e_0, 0, 0)$, of norm 64.

*In Section 5.1 today.*  Sampled: the sampling paragraph with the standard
error (main.tex 635–644), the frequencies of the table (648–673), "On the
$10^6$ samples ... cube identity" (682–683), and Remark 5.1's second sentence
on sampled pair-product norms (689–694).  Exhaustive, but checkable only by
machine: Remark 5.1's first sentence (all 196,560 self-product norms).
Definitions that the witnesses need and that stay in 5.1: each identity as
written in the table's first column (654–670), and the quartic identity, that
all five parenthesizations of $u \star u \star u \star u$ coincide (644–646).
Symbolic, and staying: Remark 5.2 and the closing sentence on closure
(approved at prompt 284).  Pointed out, not proposed: over the whole minimal
shell the cube identity holds on exactly $4/13$ of $\mathrm{Min}(\Lambda)$ and
the quartic identity on exactly $58/455$, and the quartic identity implies the
cube identity on all of $\mathrm{Min}(\Lambda)$; the sampled values agree.

### F6. Closed-form proofs of claims the paper now supports by computation

- **(a) The span of $\star$.**  $S = (\Lambda \cap D) + 2\Lambda$, with $D$ the
  diagonal ideal of Section 5.4 and $\Lambda \cap D = \{(a,a,a) : a \in Ls\}$;
  $\Lambda/S \cong (\mathbb Z/2)^{16}$, so $[\Lambda : S] = 2^{16}$ and
  $2\Lambda \subseteq S$.  Section 5.6 states the index and $2\Lambda
  \subseteq S$ by Hermite normal form; the identification of $S$ is the
  review's (its §3.1) and ledger claim 015's, by Hermite normal form only.
  New at prompt 317: a proof by hand, from the coordinates of F1, (4.1) to
  (4.3), the closure of $\mathcal O$, and one lemma that appears to be new
  (literature not searched): in an octonion algebra $C$ over a field, the
  products $(ay)b$ span $C$ for every $y \ne 0$.  Applied in $\mathcal
  O/2\mathcal O$, the split octonion algebra over $\mathbb F_2$, to $y = e_1 +
  e_2$ and $y = p$, it shows that $\langle(\mathcal O y)\mathcal O\rangle$ has
  odd index in $\mathcal O$; since it also contains $(yy)\mathcal O =
  N(y)\,\mathcal O$, with $N(y) = 2$ or $4$, the index is a power of 2, hence
  1.  With (4.2) this gives $\langle L\cdot L\rangle = L$ and $\langle
  p\mathcal O \cdot p\mathcal O\rangle = p\mathcal O$.  Verified exactly in
  two passes.
- **(b) The odd column of the census.**  All 588 closing odd permutations
  share the span of $\star$ (the review's §3.1).  New: the reason.  Each maps
  $L\bar s$ onto $2\mathcal O^{(n)}$ for one of the seven maximal orders that
  contain the basis units, and then (a) applies.  The even column ($2^{10}$)
  stays computational.
- **(c) $\rho$'s index 1.**  Closed form apart from three eight-vector lattice
  identities, one of them through the image of $\mathcal O^{(6)}$ under
  conjugation by a unit of order 3, an order that is not among the seven.
- **(d) The eight idempotents (Section 5.2) and the ambient square-zero
  classification (Section 5.5)** through $D \oplus T$, with $T \cong
  (\mathbb O \otimes \mathbb C,\ x * y = \tau(xy))$, $\tau$ the complex
  conjugation of the factor $\mathbb C$, via $(x,y,z) \mapsto \sigma(x + \zeta
  y + \zeta^2 z)$, $\zeta = e^{2\pi i/3}$.  This is the tensor description the
  author kept out of the paper at prompt 284 (ruling SF10 there); prompt 314
  also excluded the $D \oplus T$ reading of the idempotents as a new result.
- **(e) The review's index formula** (its Theorem B: for every closing $g$,
  $[\Lambda : S_g] = [Ls : \langle Ls \cdot_g Ls\rangle]\cdot[L\bar s :
  \langle L\bar s \cdot_g L\bar s\rangle]^2$) has an elementary proof through
  $D \oplus T$.  Excluded at prompt 314 as a new result.

### F7. New results beyond the paper's claims (pointed out only)

- **The census, nearly in closed form**, through the Fano plane $P$ of (2.1)
  and its mirror $P'$ (the lines of $P$ under $i \mapsto -i \bmod 7$).
  $\pi(L\bar s)/2L$ is a subalgebra iff $\pi(P)$ shares a line with $P$
  (proved).  The rule for the $Ls$-half is proved up to one lattice identity
  and one nonclosure witness.  Consequences: no 4-cycle closes, and an even
  cycle type $K$ closes for $(|K| - |K \cap G|)/2$ of its permutations, $G$
  the 168 permutations that preserve the lines of $P'$, which gives the exact
  halves of the 3-cycles, 5-cycles and $(3,2,2)$.  It needs about a page of
  Fano-plane and code machinery.
- **Baez's puzzles (part 10).**  Every Leech sublattice of $\mathbb O^3$
  contains $2\mathbb O^3 = (\mathbb O^3)^*$, so all of them, at least
  244,035,421, are closed under $2(xy+yx)$ (Puzzle 3); none is closed under
  $xy + yx$ (Puzzle 1).
- **$(L/2L, \bar\mu)$.**  Of the 270 totally singular subspaces, 142 are
  subalgebras, 16 are right ideals and 8 are left ideals (the review's §3.3).
  The 8 left ideals are exactly the two-sided ones, and new at prompt 317 is
  their identification: the seven $2\mathcal O^{(n)}/2L$ and $L\bar s/2L =
  2J_1/2L$ (Kirmse's set, doubled).  The other 8 right ideals are not left
  ideals.  Also: a cheaper subalgebra test by a trilinear form, and $\bar\mu$
  as the isotope $x * y = (x(e_1+e_2))y$ of $\mathcal O/2\mathcal O$.
- **Overlattices.**  $\Lambda + \mathbb Z\varepsilon_1 + \mathbb
  Z\varepsilon_2 + \mathbb Z\varepsilon_3$ is $\star$-closed of index 16; the
  projections to $D$ and $T$ give a closed overlattice of index $3^8$; the
  maximal order $\Omega$ above $L^3$ (ledger claim 020).
- **Conventions.**  Signed permutations; right multiples, where the parity
  law flips; permutations that move $e_0$.

### F8. Records outside the paper that are wrong or imprecise

- The provenance of the criterion and of the counts (§0); corrected in this
  hand-over.
- The Hermes supplement (frozen; cited in Appendix C's footnote; lines
  178–184) and the W-concrete review
  (`paper/reviews/2026-05-24_claude_opus_4_7_W_concrete.md`, frozen; lines
  55–57 and 139) state "$2e_0 \in 2L$" and conclude that $w_0$ represents $s$
  modulo $2L$; the review also concludes "$s^2 \equiv s \pmod{2L}$" (lines
  77–79).  These are false: $e_0 \notin L$, so $2e_0 \notin 2L$; $w_0 = s +
  2e_0$ lies in $\sigma(Ls)$ while $s$ does not; and $s^2 - s = -(e_0 +
  \cdots + e_7) \notin 2L$.  The multiplication table of $W$, the idempotent
  $w_0$ of $W$, and the absence of an identity in $L/2L$ are unaffected.
  Each document takes only a dated addendum.
- Ledger claim 010a, CURRENT_STATE and research_result give
  $|\mathrm{Stab}(Ls)| = 168$ and $|\mathrm{Stab}(L\bar s)| = 12096$ without
  naming the product.  Both hold in $\mathrm{Aut}(\mathbb O, \cdot_\sigma)$,
  the companion note's setting; under the standard product the orders are 21
  and 1344 ($|\mathrm{Stab}(L)| = 1344$ under both, since $\sigma(L) = L$).
- The registry entry [Baez2014] says that only the 17,280 lattices are
  compatible with Jordan closure; part 10 says "at least 17,280" (its first
  Result 1) and asks in Puzzle 3 whether there are more.  The entry
  [Calderbank2023] omits his credits.
- `egan_baez_count.py`, cited for Egan's count, builds its Gram matrix in
  floating point with integrality assertions; its counts are right.

---

## 2. Target structure of v8

- **Abstract.**  One optional clause naming the criterion (decision (p)).
- **Section 1.**  The overview clauses for Section 5 (the criterion) and for
  Appendices A and D; the Baez–Egan sentence at the right scale.
- **Section 4.**  Unchanged.
- **Section 5.1.**  General properties only: each identity, as now in the
  table's first column, with its witness (F5); $\varepsilon_\alpha$ and
  $\omega$ defined here; the consequences for composition and symmetric
  composition algebras (the Elduque and Okubo citations, first full names,
  stay here); Remark 5.1 with its witness; Remark 5.2; the closing sentence on
  closure.
- **Section 5.2.**  One clause: $\varepsilon_\alpha$ and $\omega$ as in
  Section 5.1.  **Sections 5.3 to 5.5.**  Unchanged.
- **Section 5.6.**  A twist is an invertible linear map $g$ with $g(L) = L$.
  The criterion as a proposition with its proof; $\sigma$ and the identity as
  instances; the census read through the proposition; the span paragraphs,
  with F6(a) if the author takes it.
- **Section 6, the Baez–Egan block, in reading order.**  Egan's construction,
  exactly, with the count $270 \times 64$ in lattice terms (pointing to
  Appendix A for the reading modulo 2), its two instances $\Lambda$ and
  $\sigma(\Lambda)$, and the consequence of the criterion; the Jordan setting
  at the right scale, followed by the Anatomy with its opening at the right
  scale; the reproduction chain with defined labels.  The Comparison
  paragraph goes.  The precursors, the history pointer and Corradetti stay.
- **Section 7.**  The summary, with one optional clause on the criterion
  (decision (p)); the third question narrowed to what stays open; the fourth
  question deleted.
- **Appendix A**, retitled for the criterion modulo 2: the descent of the
  product to $L/2L$; the criterion modulo 2 in one sentence; Wilson's pair and
  $\sigma$'s pair; $q$ and its subspaces in the terms of decision (i), each
  defined before use; the answer to the polarization question (decision
  (k)); prior work, credited once.
- **Appendix B.**  Unchanged.  **Appendix C.**  The Hermes sentence keeps its
  pointer to Appendix A; its "Lagrangian framing" follows decision (i).
- **Appendix D (new).**  The sampled frequencies: the table, the sampling
  method with the standard error, the quartic-implies-cube observation, the
  sampled pair-product norms, and (decision (b)) the exhaustive self-product
  norms.
- **Bibliography.**  The links and dates of parts 11 and 12 in [Baez2014];
  the anchor of Calderbank's comment; any citation decision (j) adds.

Estimated length: about even, with Section 6 and Appendix A shorter and
Sections 5.1 and 5.6 longer; 20 pages should hold.

---

## 3. Items, in order (one hand-over each)

*Own plan.*  An item marked *own plan* gets its own short plan file in
`evidence_and_reasoning/` (status line, frozen once engaged, indexed), with
the exact old → new text, as the span-index findings got
`2026-10-14_plan.md`; it is applied only after the author has read it.
*Every item*, before its hand-over: its changelog entry, an index line in
`python_project/src/README.md` for any new script, a build, and the linter.

**Item 1.  Section 5.1 by witnesses; the sampled material to Appendix D**
(F5; ruling 5).  *Own plan.*  $\varepsilon_\alpha$ and $\omega$ defined in 5.1,
with 5.2's defining clause edited; the identities' definitions kept in 5.1;
the witnesses; Remark 5.1 with its witness; Appendix D; Section 1's overview
clause for Appendix D; a new exact script for the witnesses, committed before
the paper cites it; the scripts behind the table cited from Appendix D.
Decisions (a), (b), (c), with (h) asked together with them (it decides the
appendix letters).

**Item 2.  The closure criterion in Section 5.6** (F1).  *Own plan.*  The
twist restricted to $g(L) = L$; the proposition, its proof and its two
instances; after the census table, "No such description is known ..."
narrowed to the criterion's reading of the table; Section 7's third question
narrowed; Section 1's "It also tests twists other than $\sigma$" updated.  The
census sentences on the 4-cycles, the $(5,2)$ type, the exact halves and the
3-cycle orientations stay as data.  Untouched: Section 4, the approved lead of
5.6, and the parity, $\rho$ and nonisomorphism paragraphs.  Decisions (d) to
(g).

**Item 3.  Appendix A around the criterion modulo 2** (F4).  *Own plan*, which
first reads Lepowsky–Meurman's Corollary 2.3 in the original.  The title, and
Section 1's overview clause for Appendix A to match; the descent of $\bar\mu$
without the structure-constants sentence; the criterion modulo 2; Wilson's
pair and $\sigma$'s pair, with (A.3) condensed; the terms of decision (i),
including Section 6's two uses of "totally isotropic" (main.tex 1004, 1009) if
the terms change, so that the paper stays consistent; the polarization answer
of decision (k); prior work credited once (decision (j)), with Silverman's
full name, since Section 6 no longer names him after item 4; the
research-program sentence replaced by the answer; Section 7's fourth question
deleted; Appendix C's Hermes wording per decision (i).  Decisions (i) to (k).

**Item 4.  Section 6, the Baez–Egan block** (F2, F3, the credits of F4).
*Own plan.*  The reading order of §2; Egan's construction with $L_1 \cap L_2
= 2L$, the count in lattice terms, the instances $\Lambda$ and
$\sigma(\Lambda)$, and the consequence of the criterion; the Jordan paragraph
at the right scale, with the exact count, both names of the product, and
Egan's finding credited to Egan; the Anatomy's opening at the right scale;
the reproduction chain with the wording verified at prompt 316 (item 4c of
its log); the Comparison paragraph deleted; Section 1's Baez–Egan sentence;
the bibliography entries.  Decisions (l), (m).

**Item 5.  The span of $\star$ by hand** (F6(a)).  *Own plan*, with a
literature check of the lemma first.  Section 5.6's sentence "In particular,
$\star$ has span index $2^{16} = 65{,}536$, and its span contains $2\Lambda$"
(main.tex 922–924) becomes the statement $S = (\Lambda \cap D) + 2\Lambda$
with its proof; the Hermite-normal-form scripts stay as cross-checks.
Decision (n).

**Item 6.  The abstract and Section 7's summary.**  Trivial: one clause each
naming the criterion, in the author's words, or no change.  Decision (p).

**Item 7.  Records, the changed sections read, the linter.**  The
cross-cutting records: CURRENT_STATE, README, `python_project/src/README.md`,
this plan's status in `evidence_and_reasoning/README.md` and CURRENT_STATE,
research_result, the ledger (claim 017, already in lattice form, extended to
every $g$ with $g(L) = L$ and every bilinear $\mu$ with $\mu(L, L) \subseteq
L$, with the exact proof of F1; rows for $\sigma(\Lambda)$ as an Egan lattice
and for the witnesses; F6(a) if taken), the glossary (the entries "Twist" and
"Polarization of L/2L", and the terms of decision (i)), the registry (F8;
Lepowsky–Meurman's Corollary 2.3 if read), dated notes in the cited scripts
whose docstrings name the old Appendix A or the Section 5.1 table, the
corrections of decision (q), and memory.  Then the changed sections are read
in the PDF and the linter is run once more.

*Order.*  Item 1 is independent of the mod-2 thread and is the author's most
recent ruling; with (h) decided, it fixes the appendix letters first.  Items
2 to 4 are the mod-2 thread in dependency order: the criterion is the
statement the other two cite; Appendix A restates it modulo 2 and answers the
polarization question, so Section 7's fourth question goes there; Section 6
then points to both.  The paper is consistent after each item.  After item 2,
Appendix A's research program and Section 7's fourth question are still open
in the paper, since the answer is not yet stated.  After item 3, Section 6
keeps its remaining mod-2 wording until item 4.  Item 5 is optional and
independent of items 3 and 4.  Items 6 and 7 close.

*Where the items parked at prompt 316 go* (its log, Response, item 4): (a)
the census sentence and the third question: item 2, decisions (d) and (f);
(b) the Comparison paragraph: item 4, decision (l); (c) the Egan paragraph's
three subjects and undefined labels: item 4; (d) the credits: items 3 and 4,
decision (j); (e) the terms: item 3, decision (i); (f) the fourth question and
"lift": item 3, decision (k); (g) "Polarization of $\bar L$ by $\sigma$" and
(h) the structure-constants sentence: item 3; (i) the repetitions between
Section 6 and Appendix A: Egan's construction stated once, in Section 6,
which Appendix A cites (item 4), and "ideal and subalgebra" once, in Appendix
A (item 3); (j) the [Baez2014] links: item 4.  The nits of that log's item 5
are not in this plan, except where an item rewrites the sentence anyway.

---

## 4. What to learn, item by item

| Item | Concepts | Where to read |
|---|---|---|
| 1 | The identities: associative, alternative, flexible, power-associative; why a polynomial identity that fails on $\mathbb R^{24}$ fails on a lattice of full rank (scale by 4) | R. D. Schafer, *An Introduction to Nonassociative Algebras* (1966), Ch. I and III |
| 2 | The substitution $a = x+y$, $b = x+z$, $w = x+y+z$, which inverts Wilson's conditions; an "if and only if" proved by sufficiency plus test vectors; automorphisms of $L$ versus isometries | Wilson (2009), §3; Conway–Sloane, *Sphere Packings, Lattices and Groups*, Ch. 4, §8 |
| 3 | Quadratic forms over $\mathbb F_2$ and their polar forms; the two conventions for "totally isotropic" and "totally singular"; Witt index and plus type; $W(E_8)/\{\pm1\}$ as the full orthogonal group of $q$ ($O_8^+(2){:}2$ in ATLAS notation) and Witt's theorem; the complementary-pair construction | R. A. Wilson, *The Finite Simple Groups* (2009), §3.12.4; D. E. Taylor, *The Geometry of the Classical Groups* (1992), Ch. 11; Lepowsky–Meurman (1982), §2 |
| 4 | Egan's construction and its count; the exceptional Jordan algebra $\mathfrak h_3(\mathbb O)$ and its Jordan product; the standard and $\sqrt2$ scales of the Leech lattice | Baez, *Integral octonions*, parts 9 to 11 with Egan's comments; Baez, *The octonions* (2002), §3.4; Conway–Sloane, Ch. 4, §11 |
| 5 | Octonion algebras over a finite field; $\mathcal O/2\mathcal O$ as the split octonions over $\mathbb F_2$; the span of a set of products and its index | R. A. Wilson, *The Finite Simple Groups* (2009), §4.4 |

---

## 5. Decisions for the author

Each decision is asked again before its item; the recommendation is marked.

- **(a) Order.**  Item 1 first (*recommended*), or after the mod-2 thread.
- **(b) Section 5.1.**  The witnesses as a short table (*recommended*) or in
  sentences.  $\varepsilon_\alpha$ and $\omega$ defined in 5.1, with 5.2's
  clause edited (*recommended*), or the witnesses written as explicit vectors
  and 5.2 unchanged.  Remark 5.1: its statement, with the witness, stays in
  5.1; the exhaustive value set of the 196,560 self-product norms moves to
  Appendix D with the sampled pair-product norms (*recommended*: the letter of
  ruling (5), each failure checkable by the reader), or stays in 5.1 (exact,
  but checkable only by machine).
- **(c) The new appendix.**  Appendix D, after C (*recommended*: no
  relettering, provided (h) keeps Appendix A), or elsewhere with relettering.
- **(d) Twist.**  An invertible linear map with $g(L) = L$ (*recommended*);
  otherwise the proposition carries the hypothesis $L \cdot_g L \subseteq L$,
  and Section 7 keeps a question for the twists outside it.
- **(e) The $\sigma$ instance.**  One sentence in 5.6 stating that the
  proposition gives Theorem 1.1 again, without a remark on Lemma 4.3; Section
  4 unchanged (*recommended*).
- **(f) Section 7, third question.**  Narrowed to what stays open
  (*recommended*), or deleted.
- **(g) Provenance in the paper.**  None, as for $\rho$'s index 1 from the same
  review; the records carry it (*recommended*).  Or a footnote naming the
  review of 13 October 2026.
- **(h) Appendix A.**  Rebuilt around the criterion modulo 2
  (*recommended*), or deleted, with the count's reading moved to Section 6.
  Asked together with (a) to (c), before item 1.
- **(i) Terms.**  The paper's present terms ("totally isotropic", with
  "Lagrangian" and "polarization"), each used only after its definition, and
  "plus type" and "Witt decomposition" defined or replaced by plain words in
  item 3's own plan (*recommended*: the standing rule "established terms
  only; never rename"; the paper's usage is that of Elman–Karpenko–Merkurjev
  and of Wilson's book).  Or, as an exception to that rule for the author to
  rule on: "totally singular" for $q = 0$, the usage of Taylor, Cameron and
  Shimakura, defined before use, with plain words ("a complementary pair")
  instead of "Lagrangian", "polarization" and "Witt decomposition".
- **(j) Credits.**  The orthogonal space $(L/2L, q)$ as classical, with the
  citation chosen in item 3's own plan (Wilson's book, §3.12.4, is the source
  found).  Silverman once, in Appendix A, for the identification of the 270
  sublattices with the maximal subspaces on which $q$ vanishes (part 12) and
  for the reading of Egan's construction modulo 2 (his comment of 16 February
  2016).  The complementary-pair construction credited through Calderbank's
  comment; Lepowsky–Meurman's Corollary 2.3 cited directly only after it is
  read; Tits 1980 not cited, since unread (*recommended*).
- **(k) The polarization answer.**  In Appendix A, with Section 7's fourth
  question deleted.  The characterization: every ordered complementary pair
  is $(g(L\bar s)/2L, g(Ls)/2L)$ for an isometry $g$ of $L$, and $\star_g$
  closes exactly when both halves are subalgebras (*recommended*).  Or, in
  addition, the count 4,704 of 17,280 with an exact repository script: a
  computed classification of subalgebra pairs, the kind of structure the
  author's ruling SF10 of prompt 284 leaves until a referee asks for it.
- **(l) The Comparison paragraph.**  Deleted (*recommended*), or reduced to
  one sentence with the answer.
- **(m) Section 6's count.**  In lattice terms with a pointer to Appendix A
  (*recommended*), or with the reading modulo 2 in place.
- **(n) Item 5.**  The span of $\star$ by hand into v8, after a literature
  check of the lemma (*recommended*: it replaces the computational support of
  a headline property), or the Hermite-normal-form support kept.
- **(o) The rest of F6, and F7.**  This narrows ruling (2), "all of the new
  closed-form proofs", for these reasons: F6(b) explains only the odd column,
  and the even column stays computational; F6(c) still rests on three
  computed lattice identities; F6(d) uses the tensor description ruled out at
  prompt 284 (SF10); F6(e) and the $D \oplus T$ reading of the idempotents
  were excluded at prompt 314 as new results.  So F6(b) to (e) stay out of
  v8, each a possible plan of its own afterwards (*recommended*); or (b) and
  (c) join item 5.  F7: new results, pointed out only.
- **(p) The abstract and Section 7's summary.**  One clause each naming the
  criterion, in the author's words (*recommended*), or unchanged.
- **(q) Records.**  Dated addenda to the Hermes supplement and to the
  W-concrete review, for the $2e_0$ slip and what was drawn from it; the
  product named in the stabilizer statements of the ledger, CURRENT_STATE and
  research_result (*recommended*, in item 7).
