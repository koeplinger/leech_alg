# The span index of twisted triple products: verification of the external v8 review, and what index 1 means (14 October 2026)

**Date:** 2026-10-14
**Author of this note:** Claude Fable 5.1 (Anthropic), at the direction of Jens Köplinger
**Subject:** the external review of v8 by Claude Opus 5.5 (13 October 2026),
incorporated verbatim with its evidence as
[`2026-10-13_claude_opus_5_5_review_v8/`](2026-10-13_claude_opus_5_5_review_v8/)
(write-up also at [`2026-10-13_claude_opus_5_5_review_v8.md`](2026-10-13_claude_opus_5_5_review_v8.md)).
**Prompt:** 287.
**Status:** machine-obtained; every computational statement below was obtained
in at least two independent implementations; the structural identification of
Section 4 in two; the maximal order of Section 3.3 in two.  Not
peer-reviewed.

**Addendum (14 October 2026, prompt 288).**  Three corrections; the text
below stands as written.

1. Section 3.1's parenthesis "(the minimal vectors included, which are
   never products, since product norms are at least 16)" gives the wrong
   reason.  The bound $N \geq 16$ holds only for products of two minimal
   vectors under $\sigma$ (paper, Remarks 1.3 and 5.1); under $s_r$ such
   products have norms 20, 24, 28 and others.  The valid argument: the block
   sum of $u \star_g v$ is $(x+y+z) \cdot_g (x'+y'+z')$, a product of two
   vectors of $Ls$, whose norm lies in $16\mathbb Z$.  So the 189,360 minimal
   vectors whose block sum has norm 4, 8 or 12 are never single products,
   under any twist by an isometry of $L$ (count computed at prompt 288); the
   remaining 7,200 were not settled.  Index 1 says that every vector of
   $\Lambda$ is an integer combination of products, not that it is a
   product.
2. Section 3.2's "chain $\Lambda \subset \Lambda + \mathbb Z(1_r, 0, 0)
   \subset L^3$ ... with indices 4 and $2^{10}$" is wrong: $1_r \notin L$,
   so $\Lambda + \mathbb Z(1_r, 0, 0)$ is not contained in $L^3$.  It is
   $\star_r$-closed of index 4 over $\Lambda$ and, like
   $\Lambda + \mathbb Z\varepsilon_1$ for $\sigma$, not contained in $L^3$
   (checked exactly at prompt 288).
3. Section 2's "Not re-verified here: the review's spectrum ..., its rule
   via $2 \cdot 1_g$" was overtaken after this note was written.  The
   adversarial pass on the assessment (prompt 287) re-derived the spectrum
   and, at the level of the 142 subalgebra Lagrangians, the rule
   ($h = 1$ for all 120 not containing the class of $2e_0$; $2^5$ or $2^8$
   otherwise).  At prompt 288 two further agents found the rule, Theorem A
   and Theorem B to hold on random samples of $W(E_8)$ (1,115 and 439
   closing twists).  The rule's direction "index $> 1$ when $2 \cdot 1_g$
   lies in a Wilson sublattice" has a short proof; the converse and the
   exact values are computed, not proved.

**Addendum (14 October 2026, prompt 290).**  Two corrections; the text
below stands as written.

1. Section 4, item 4: "this remains a finite check ($p$ does not
   normalize the order, so the Moufang identity does not give it)" is
   wrong.  The premise holds ($pO^{(4)} \neq O^{(4)}p$), but closure needs
   no normalizing.  For purely imaginary $p$, $p^2 = -N(p)$, and the left
   Moufang identity gives $(pa)(pb) = p((ap)b)$ for all octonions $a, b$.
   Since $p = e_2 - e_4 + e_5 - e_6$ is purely imaginary and lies in
   $O^{(4)}$, $p\,O^{(4)}$ is closed, and Lemma 4.4 follows from the closure
   of $O^{(4)}$ as Lemma 4.3 does (ledger claim 021; symbolic and exact
   checks at prompt 290).
2. Section 4, item 1: $L\bar s = 2J_1$ and $L = J_1 \cdot s$ are Wilson's,
   and this note does not cite him: in his coordinates "$L\bar s = 2B$",
   with $B$ the set whose assumed closure "is known as Kirmse's mistake"
   (Wilson 2009, p. 2187), and "$Bs = L$" (p. 2189).  New beyond Wilson:
   $\sigma(L\bar s) = 2O^{(4)}$ and $\sigma(Ls) = p\,O^{(4)}$, and the closure
   argument of item 1.  The closing remark of Section 4 on "twice Bruck's
   order" was not checked: $O^{(4)}$ is one of the seven, and which of them
   is Coxeter's $J$ was not determined.

Notation as in the paper (v8): $L = D_8^+$, Wilson's left-multiple sublattices
$L\bar s$ and $Ls$, $\Lambda \subset L^3$; for an isometry $g$ of $\mathbb R^8$
with $g(L) = L$, $x \cdot_g y := g^{-1}(g(x) \cdot g(y))$, $\star_g$ the triple
product of Definition 3.3 built from $\cdot_g$, $S_g := \mathbb Z\text{-span}\{u
\star_g v : u, v \in \Lambda\}$; the paper's $\sigma = (1\,2)$ is the reflection
in the root $e_1 - e_2$.

---

## 1. Method

The review's computational claims were re-derived without its code:

1. **Repository script** `python_project/src/verify_twist_span_index_census.py`
   (exact arithmetic; own integer Hermite normal form; the lattices and
   membership tests of `verify_consecutive_twists_exact.py`): the reflection
   twist, then all 5,040 permutation twists with the review's Theorems A and B
   tested on each.  Output ALL PASS (1,169 s on 12 cores).
2. **An independent Python implementation** written from the paper's
   definitions by a second-pass agent, with its own 24-product certificate for
   index 1 (determinant $+4096$ in $L^3$ coordinates, Gram determinant
   $2^{24}$), cross-checked against sympy and PARI on the witnesses; its own
   census of all 5,040 permutations and 120 reflections.
3. **An independent GAP 4.12.1 and PARI/GP 2.15.4 implementation** written from
   the paper's definitions by a third agent (elementary divisors instead of
   Hermite normal forms), including the full reflection table; adversarially
   re-checked by a fourth agent.

The orchestrator checked Theorem A's congruence identities modulo $2L$ on
random lattice pairs for $\sigma$ and for the reflection, the necessity
witnesses $(v, 0, v)$ and $(w, w, w)$, closure of $L$ under both twisted
products, and the identifications of Section 4.

## 2. Results

**The index-1 witness.**  $r = \tfrac12(e_0 + e_1 + \dots + e_6 - e_7)$ is a
root of $L$, $s_r(L) = L$, and $s_r(e_0) \notin L$.  All 576 basis-pair
products $u \star_r v$ lie in $\Lambda$, and $[\Lambda : S_r] = 1$: the
products $\mathbb Z$-span $\Lambda$.  Obtained by own Hermite normal form (index
$[L^3 : S_r] = 4096 = [L^3 : \Lambda]$), by elementary divisors $1^{24}$ in GAP
and PARI, and by a normal-form-free certificate (24 single products forming a
$\mathbb Z$-basis of $\Lambda$).  $L \cdot_r L \subseteq L$ on all 64 basis
products.

**The permutation census** (all 5,040 permutations $\pi$, the paper's
Section 5.6 convention $x \cdot_\pi y = \pi^{-1}(\pi(x)\cdot\pi(y))$): 1,764
close, in the paper's cycle-type counts; every closing even permutation has
$[\Lambda : S_\pi] = 2^{10}$ and every closing odd one (the 21 transpositions,
$\sigma$ among them) $2^{16}$.  The 3-cycle $5 \to 6 \to 7 \to 5$ closes with
$2^{10}$; its inverse does not close.

**The root reflections** (120, one per $\pm r$): 56 close, 28 with index
$2^{16}$ (the integral roots $e_i - e_j$ and $e_0 + e_j$) and 28 with index 1
(the half-integral roots with one or five minus signs, exactly those with
$\langle r, s\rangle = \pm 1$); the others do not close.  Reproduced row by row
in GAP and PARI.

**Theorem A of the review** ($\star_g$ closes on $\Lambda$ iff $L\bar s$ and
$Ls$ are both closed under $\cdot_g$): the proof's three identities modulo
$2L$,
$P_1 + P_2 \equiv aa' + ab' + ba'$, $P_1 + P_3 \equiv ab' + ba' + bb'$,
$P_1 + P_2 + P_3 \equiv ww'$ for $a = x+y$, $b = x+z$, $w = x+y+z$, hold on
random pairs for both twists; the proof was read line by line by two agents and
the orchestrator and is correct as written; the criterion holds as a numerical
equivalence on all 5,040 permutations and 120 reflections.  **Theorem B**
($[\Lambda : S_g] = h_W \, h_V^2$ with $h_W = [Ls : \langle Ls \cdot_g Ls\rangle]$,
$h_V = [L\bar s : \langle L\bar s \cdot_g L\bar s\rangle]$) holds on every
closing twist computed, with $(h_W, h_V) = (1, 2^8)$ for the $2^{16}$ class,
$(1, 2^5)$ for $2^{10}$, $(1, 1)$ for index 1; its proof is a sketch whose
steps were machine-checked, and the assessment agent wrote out the four
remaining lemmas.

Not re-verified here: the review's spectrum over all 17,280 Lagrangian pairs
(the values $2^5$, $2^8$ and the counts), its rule via $2\cdot 1_g$, the
signed-permutation census, and the right-multiple convention variant.

## 3. What index 1 means

**3.1 Two different properties.**  The span index $[\Lambda : S_g]$ measures
how far the products fall short of generating $\Lambda$ as a
$\mathbb Z$-module: index 1 says $S_g = \Lambda$, every vector of $\Lambda$
(the minimal vectors included, which are never products, since product norms
are at least 16) is an integer combination of products; the nonunital ring
$(\Lambda, \star_g)$ is idempotent, $\Lambda^2 = \Lambda$, and the mod-2
defect of Section 5.5 vanishes.  Maximality of the order is a different
property: no lattice $\Lambda' \supsetneq \Lambda$ is $\star_g$-closed.  The
first bounds the products from below, the second bounds the order from above;
neither implies the other.

**3.2 $\Lambda$ is not maximal under any twist by an isometry of $L$.**
For $g \in \operatorname{Aut}(L)$, $L \cdot_g L = g^{-1}(g(L) \cdot g(L))
\subseteq g^{-1}(L) = L$ by Lemma 4.2, so every block of $u \star_g v$ for
$u, v \in L^3$ is a sum of three elements of $L$: $L^3 \star_g L^3 \subseteq
L^3$, and $[L^3 : \Lambda] = 2^{12}$.  Hence $(\Lambda, \star_g)$ is not a
maximal order for any such $g$, the 28 index-1 reflections included (checked
exactly for $\sigma$, the 3-cycle and $s_r$).  For the index-1 reflection the
chain $\Lambda \subset \Lambda + \mathbb Z(1_r, 0, 0) \subset L^3$ of
$\star_r$-closed lattices has indices 4 and $2^{10}$, with $1_r = s_r(e_0)$
the identity of $\cdot_r$ (GAP, one agent; the analogous chain for $\sigma$ is
Section 7 of the paper).  Ledger claim 018.

The second pass added a reason that does not depend on the twist at all.
For any product assembled from a bilinear block product by the paper's
routing, $D$ and $T$ are mutually annihilating ideals, so the projections
$\operatorname{proj}_D$ and $\operatorname{proj}_T$ are multiplicative and
$P := \operatorname{proj}_D(\Lambda) \oplus \operatorname{proj}_T(\Lambda)$
is closed whenever $\Lambda$ is; and $P \ne \Lambda$, because $\Lambda \cap D
= \Delta(Ls)$, $\Lambda \cap T = L\bar s \otimes T_{\mathbb Z}$ and
$[\Lambda : (\Lambda \cap D) \oplus (\Lambda \cap T)] = 3^8$ (computed for
$\sigma$ and $s_r$).  Wilson's 3-adic gluing of the two parts is what rules
maximality out, for every closing twist of this kind, whether or not $g$
preserves $L$ (one agent; the ingredients are the paper's own ideal
decomposition and the review's Theorem B, step 3).

**3.3 An explicit maximal order above $L^3$** (assessment agent; closure,
chain, index and the maximality argument re-derived independently by the
adversarial second pass; ledger claim 020).  Let $O$ be a unital maximal order of $(\mathbb O, \cdot_g)$
containing $L$ (for $\sigma$: $O = \sigma(O^{(4)})$ with $O^{(4)}$ the
octavian order of Section 4).  With $\Delta(a) = (a,a,a)$ and $\omega$ a
primitive cube root of unity, put
$\Omega := \Delta(\tfrac13 O) \oplus \{(p,q,r) = \tfrac13(2x - y,\, -x + 2y,\,
-x - y) : x, y \in O\}$.  Then $\Omega$ is $\star_g$-closed, contains $O^3
\supseteq L^3 \supseteq \Lambda$, has $[\Omega : \Lambda] = 2^{24} \cdot 3^8$,
and is maximal: any closed $M \supsetneq \Omega$ contains the idempotent
$e_D = \tfrac13 \Delta(1_g)$, hence splits as $(M \cap D) \oplus (M \cap T)$,
and each part is forced to equal $\Omega$'s by maximality of $O$ and of
$O \otimes \mathbb Z[\omega]$.  $\Omega$ is decomposable and not an integral
lattice; nothing Leech-like survives at the top of the chain.

**3.4 Where the intuition comes from.**  For a unital order, $m = 1 \cdot m$
makes the span index 1 automatically, and the classical index that measures
distance from maximality is $[O_{\max} : O]$, the discriminant, never
$[O : O^2]$.  In the nonunital setting here, $[\Lambda : \Lambda^2]$ is a new
invariant, index 1 is not automatic ($\sigma$: $2^{16}$), and it says nothing
about maximality.  The review's remark that no element multiplies unimodularly
is correct ($\det$ of left multiplication by $u = (x,y,z)$ is
$N(x+y+z)^4$ times an integer, divisible by 16): index 1 is a collective
property of the products, not the effect of a unit-like element.  Two toy
counterexamples (two-dimensional, nonunital, exact) show that neither
implication holds in general.

**3.5 Verdict on significance.**  A necessary correction and a worthwhile
structural remark, not a breakthrough: it settles the paper's Section 7
question (a product built from octonion multiplication can have the minimum
index), it shows $2^{16}$ is a feature of $\sigma$ and not of the construction
(index 1 is the generic case over the closing twists, $\sigma$ the extreme
case), and it leaves the main theorem, maximality and the open question on
Lemma 4.4 untouched.  The review's own measured stance is right.

**3.6 Maximality or index 1 as a route to a table-free closure proof.**
Maximality presupposes closure and is false here; index 1 is a consequence of
the structure, not a cause.  Theorem A reduces closure to two subalgebra
conditions (it shrinks the certificate from three 64-entry integer tables to
closure of two 4-dimensional subspaces of $L/2L$) but does not remove finite
certificates.  What removes them for Lemma 4.3 is Section 4.

## 4. A structural reason for Lemma 4.3

Verified first-hand (exact; the assessment agent's computation, re-derived
independently by the orchestrator; ledger claim 019).  Let $C_K$ be Kirmse's
code, the $[8,4,4]$ code generated by $1^8$ and the words $\{0\} \cup \ell$
for the seven Fano lines $\ell$ of eq. (2.1), and $J_1 := \mathbb Z^8 +
\tfrac12 C_K$ Kirmse's set (Appendix B).  For $n = 1, \dots, 7$ let $O^{(n)}
:= \mathbb Z^8 + \tfrac12 (0\,n)C_K$ with $(0\,n)$ the coordinate transposition.
Then:

1. $L\bar s = 2 J_1$: Wilson's sublattice is twice Kirmse's set (and
   $L = J_1 \cdot s$, $Ls = J_1 \cdot s^2$).
2. $\sigma(L\bar s) = 2\, O^{(4)}$, and $O^{(4)}$ is the only one of the seven
   with this property.  All seven $O^{(n)}$ are closed under the standard
   product (Coxeter's theorem) and all seven contain $L$; $J_1$ is not closed.
3. Hence $L \cdot \sigma(L\bar s) = L \cdot 2 O^{(4)} \subseteq 2 \,(O^{(4)}
   \cdot O^{(4)}) \subseteq 2 O^{(4)} = \sigma(L\bar s)$, and likewise on the
   right: Lemma 4.3 follows from Coxeter's theorem, two-sided, with no table;
   and $h_V = [2O^{(4)} : \langle 2O^{(4)} \cdot 2O^{(4)}\rangle] = [2O : 4O]
   = 2^8$ explains the factor $2^8$ in $\sigma$'s index $2^{16}$.
4. $\sigma(Ls) = p \cdot O^{(4)}$ with $p = e_2 - e_4 + e_5 - e_6$ (norm 4).
   Lemma 4.4 becomes "the principal left-multiple set $p \cdot O^{(4)}$ is
   closed under multiplication"; this remains a finite check ($p$ does not
   normalize the order, so the Moufang identity does not give it).

Appendix B's remark that the Coxeter--Bruck transposition reappears "acting on
a different object" is therefore sharper than stated: $\sigma$ acts on $L\bar
s = 2J_1$, twice Kirmse's set itself, and carries it onto twice Bruck's order.

## 5. Other findings of the review, checked against the record

- **Section 6(a), the seven maximal orders:** correct.  The v8 sentence "For
  a fixed octonion multiplication, the integral octonions have seven maximal
  orders" is false as stated: the rational octonion algebra has infinitely
  many maximal orders (conjugation by $q = 1 + e_1 + e_2 + e_4$, $q^3 = -8$,
  is an automorphism carrying one of the seven to a maximal order outside the
  thirty code lattices); the classical count is of the maximal orders
  containing the Gravesian order $\mathbb Z\langle e_0, \dots, e_7\rangle$,
  which are the seven $O^{(n)}$ (Coxeter 1946; Conway--Smith 2003, Ch. 9).
  The same looseness is in Section 2.2 and in the glossary.  Paper wording is
  the author's decision; the glossary is corrected in this pass.
- **Section 6.2 (the record):** confirmed by two agents: prompt-log numbers
  073 and 112 do not exist (073 restarted as 074; 112 recorded inside 111);
  logs 025 and 026 carry elisions; log 281 records no provenance for its
  comments; until prompt 287 no file named the model of prompts 281 to 286
  (now recorded there: Claude Fable 5.1); the date change to 13 October is
  logged (286); the ResearchGate v7 is called "published" in several files;
  the docstring of `verify_discovery1_W_subalgebra.py` still calls a vector
  $e_0$ that is $s = L_1$ in $L$-coordinates (the octonion identity is not in
  $L$).  Disposition is the author's.

## 6. Recommendations for the paper (the author's decision)

1. Section 5.5: state that the span index depends on the twist ($2^{16}$ for
   $\sigma$ and all odd closing permutations, $2^{10}$ for the even ones, 1
   for the half-integral root reflections with $\langle r, s\rangle = \pm 1$);
   keep "a span-surjective product would have index 1" and add that such
   products exist.
2. Section 7: delete or replace the span-index bullet; keep the
   non-maximality paragraph; do not claim maximality or a table-free proof.
3. Optional, substantive: Theorem A with its short proof, and the Kirmse
   identification of Section 4, which replaces Table 2 by Coxeter's theorem
   and sharpens Appendix B.
4. Sections 1 and 2.2: the maximal-order sentences, as in Section 5 above.

Scripts of this note: `python_project/src/verify_twist_span_index_census.py`;
the agents' scratch scripts are not part of the record.
