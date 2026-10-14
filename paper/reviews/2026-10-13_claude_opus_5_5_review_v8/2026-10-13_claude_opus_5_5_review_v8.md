# Review of v8: the span index of twisted triple products, and other findings

- **Paper:** J. Köplinger, "An order on the Leech lattice from a
  $\mathbb Z_3$-symmetric triple-octonion product", v8 (13 October
  2026).
- **Reviewer:** Claude Opus 5.5 (Anthropic), in multi-agent review and
  reproduction workflows run on 13 October 2026 at the author's request,
  outside the public research repository.
- **Status:** machine-obtained and not peer-reviewed.  Section 5 says
  exactly what is proved, what is computed exhaustively, and what is
  open.

## Summary

1. **The paper's open question on the span index is already answered.**
   Section 7 of v8 asks whether a closed product built from octonion
   multiplication can have a smaller span index than $\star$, which has
   $2^{16}$.  It can:
   - in the paper's own Section 5.6 census, every one of the 1,176
     closing *even* permutations gives index $2^{10}$ (the 588 closing
     odd ones give $2^{16}$, as $\sigma$ does);
   - replacing $\sigma$ by the reflection in the root
     $r=\tfrac12(e_0+e_1+\dots+e_6-e_7)$ of $L$ gives a product that
     closes on $\Lambda$ with index **1**: its products span all of
     $\Lambda$.

   Note that $\sigma$ is itself a root reflection, in $e_1-e_2$.
2. **The results are robust.**  They were reproduced by six independent
   implementations in Python, GAP and PARI/GP.  Index 1 is certified
   without any normal-form computation, by an explicit
   $\mathbb Z$-basis of $\Lambda$ that consists of 24 single products
   (determinant $-1$).
3. **Index 1 is not trivial.**  No element of $\Lambda$ multiplies
   unimodularly on either side, for any twist: the determinant of
   $v\mapsto u\star_g v$ is $N(x+y+z)^4$ times an integer, hence
   divisible by 16.
4. **Index 1 is the generic case, and $\sigma$ is the extreme case.**
   Two theorems (the first with a short proof, the second with a proof
   sketch whose steps are machine-checked), and an exhaustive census
   over all twists $g$ in $W(E_8)=\operatorname{Aut}(L)$, give the full
   picture:
   - *closure criterion:* $\star_g$ closes on $\Lambda$ if and only if
     $L\bar s$ and $Ls$ are both closed under $\cdot_g$;
   - *index formula:* $[\Lambda:S_g]=[Ls:\langle Ls\cdot_gLs\rangle]\,
     [L\bar s:\langle L\bar s\cdot_gL\bar s\rangle]^2$;
   - only the values $1,2^5,2^8,2^{10},2^{16}$ occur, and 71% of the
     closing cases have index 1;
   - the index is decided by where twice the twisted identity
     $1_g=g^{-1}(e_0)$ lies.  Coordinate permutations fix $e_0$, and
     $2e_0\in L\bar s$, so they always land on the top values: $2^{10}$
     for even and $2^{16}$ for odd permutations.  $\sigma$ sits at the
     maximum.
5. **This is a substantive structural result, not an isolated
   anomaly.**  Taken together, these results reduce closure to two
   subalgebra conditions and locate the span defect in where
   $2\cdot1_g$ lies.  Two things remain open: why $\sigma(Ls)/2L$ is a
   subalgebra (the paper's own first open question), and why the
   factors are $2^5$ and $2^8$.  Whether any of it enters the paper is
   the author's decision; Section 4 gives options.  Other findings on v8
   are in Section 6.

## 1. Setting and conventions

Everything below uses the paper's definitions (v8, Sections 2 to 5).
The conventions were pinned to the paper's printed data: all 192 rows of
Tables 1, 2 and 3 (LL, LN, MM), the vectors $M_k$ and $N_k$, Remark 4.5,
and the 196,560 vectors of the three minimal-vector families are
reproduced exactly (`reproduction/final/f1_setup.out.txt`).

- **Octonions.**  Basis $e_0,\dots,e_7$, $e_0$ the identity, Fano triples
  $(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)$: for a triple
  $(a,b,c)$, $e_ae_b=+e_c$ with cyclic permutations positive and
  anticyclic ones negative, and $e_i^2=-e_0$.
- **The $E_8$ lattice** $L=D_8^+$: coordinates all integers with even
  sum, or all half-integers with odd sum; $L\cdot L\subseteq L$.
- **Wilson's sublattices** are formed with **left** multiples:
  $L\bar s=\{\ell\bar s:\ell\in L\}$ and $Ls=\{\ell s:\ell\in L\}$, with
  $s=\tfrac12(-e_0+e_1+\dots+e_7)$ and
  $\bar s=\tfrac12(-e_0-e_1-\dots-e_7)$.  Both have index 16 in $L$,
  $L\bar s+Ls=L$, $L\bar s\cap Ls=2L$, and $2e_0=s\bar s\in L\bar s$.
- **Leech lattice** $\Lambda=\{(x,y,z)\in L^3: x{+}y,\,x{+}z,\,y{+}z\in
  L\bar s;\ x{+}y{+}z\in Ls\}$.  It has $[L^3:\Lambda]=2^{12}$, Gram
  determinant $2^{24}$, minimum 8, and 196,560 minimal vectors.  These
  were verified in Python, GAP (`ShortestVectors`) and PARI/GP
  (`qfminim`).
- **Twists.**  For an isometry $g$ of $\mathbb R^8$ with $g(L)=L$, that
  is, $g\in W(E_8)=\operatorname{Aut}(L)$, put
  $x\cdot_g y:=g^{-1}(g(x)\cdot g(y))$.  This is an octonion algebra
  isomorphic to $\mathbb O$ via $g$, with identity element
  $1_g=g^{-1}(e_0)$.  Two families appear below.
  - *Permutation twists* (paper, Section 5.6): $\pi(e_a)=e_{\pi(a)}$,
    $\pi(e_0)=e_0$, and $x\cdot_\pi y=\pi^{-1}(\pi(x)\cdot\pi(y))$.
  - *Root reflections* $s_r(x)=x-\langle x,r\rangle r$ for roots $r$ of
    $L$; these are involutions.  **The paper's $\sigma=(1\,2)$ is $s_r$
    for $r=e_1-e_2$.**
- **Triple product** (Definition 3.3): $u\star_g v=(P_1,P_2,P_3)$ with
  $P_1=xx'+zy'+yz'$, $P_2=yy'+xz'+zx'$, $P_3=zz'+yx'+xy'$, every
  product being $\cdot_g$.
- **Closure** means $\Lambda\star_g\Lambda\subseteq\Lambda$; by
  bilinearity it suffices to check the 576 products of a
  $\mathbb Z$-basis.
- **Span index**: $S_g=\mathbb Z\text{-span}\{u\star_gv\}$.  The index
  $[\Lambda:S_g]$ and the structure of $\Lambda/S_g$ are computed from
  the 576 integer coordinate vectors of the basis-pair products.

**Two conventions matter.**
- The direction of a permutation twist matters for individual
  witnesses.  The 3-cycle $\pi(5)=6$, $\pi(6)=7$, $\pi(7)=5$ closes
  under $x\cdot_\pi y=\pi^{-1}(\pi(x)\cdot\pi(y))$, but its inverse does
  not.  The census counts and the parity law below do not depend on
  this.
- Left versus right multiples matters for most results.  With right
  multiples ($\bar s\ell$, $s\ell$) one gets a different Leech lattice
  $\Lambda_R$.  On it, $\sigma$ has index $2^{10}$, the parity law flips,
  and no root reflection reaches index 1 (`f7`).  The paper's printed
  tables fix the left-multiple convention used here.

## 2. Two theorems

All products in this section are $\cdot_g$.  Since $g(L)=L$, we have
$L\cdot_gL\subseteq L$, so by bilinearity $\cdot_g$ descends to $L/2L$.
As $2L\subseteq L\bar s\cap Ls$, membership in $L\bar s$ and in $Ls$
depends only on the class modulo $2L$.

**Theorem A (closure criterion).**  For $g\in W(E_8)$,
$$\Lambda\star_g\Lambda\subseteq\Lambda\iff L\bar s\cdot_gL\bar s\subseteq
L\bar s\ \text{ and }\ Ls\cdot_gLs\subseteq Ls .$$

*Proof.*  For $u=(x,y,z)\in\Lambda$ put $a=x+y$ and $b=x+z$, both in
$L\bar s$, and $w=x+y+z\in Ls$.  Modulo $2L$ we have $y+z\equiv a+b$,
hence $x\equiv w+a+b$, $y\equiv w+b$, $z\equiv w+a$.  The same holds for
$u'$ with primes.  Substituting into the blocks and reducing mod $2L$
(characteristic 2, so the $w$-terms cancel in pairs) gives the identities
$$P_1+P_2\equiv aa'+ab'+ba',\qquad P_1+P_3\equiv ab'+ba'+bb',\qquad
P_1+P_2+P_3\equiv ww' .$$
These were also checked symbolically (`f7`).
- *Sufficiency.*  If $L\bar s$ is closed, the first two sums lie in
  $L\bar s$, because $2L\subseteq L\bar s$.  The third pairwise sum is
  $(P_1+P_2)+(P_1+P_3)-2P_1$, which lies in $L\bar s$ too.  If $Ls$ is
  closed, the block sum lies in $Ls$.  Each block lies in $L$.
- *Necessity.*  For $v\in L\bar s$, the vector $(v,0,v)$ lies in
  $\Lambda$, with $a=v$ and $b\equiv w\equiv0$.  For two such vectors,
  $P_1+P_2\equiv vv'$, so $L\bar s\cdot_gL\bar s\subseteq L\bar s$.  For
  $w\in Ls$, the vector $(w,w,w)$ lies in $\Lambda$, with $a\equiv
  b\equiv0$ and block sum $\equiv w$.  Hence $Ls\cdot_gLs\subseteq Ls$.
  $\square$

In the paper's language, with $V_g=g(L\bar s)/2L$ and $W_g=g(Ls)/2L$
(a polarization of $\bar L=L/2L$): $\star_g$ closes exactly when both
Lagrangians are subalgebras of $(\bar L,\bar\mu)$.

**Consequences.**
- The paper's Lemma 4.3 ($\sigma(L\bar s)$ a left ideal) is stronger
  than closure needs.  Among closing twists, its analogue ($g(L\bar s)$
  a left ideal of $L$) holds exactly for those of index $2^{16}$:
  $\sigma$, all 588 closing odd permutations, and the 28 closing
  integral reflections.  It fails for all others.
- For the index-1 reflection, $V_g$ is a subalgebra and a right ideal
  but not a left ideal, so Lemma 4.3 itself fails.  The four-lemma
  proof still goes through in mirrored form.  Write
  $P_1+P_2=(x+z)x'+(y+z)y'+(x+y)z'$, and likewise for the other two
  pairs; then condition (2) needs only the right-ideal property.
- On 3,136 of the 3,360 index-1 pairs, $V_g$ is neither a left nor a
  right ideal; on the remaining 224 it is a right ideal only.  Neither
  form of the template applies to the 3,136, and Theorem A is needed
  (`f8_ideal_classes.out.txt`).
- The question at the end of Appendix A.1, which polarizations lift,
  is answered for twists of this kind: exactly the pairs of two
  subalgebras.

**Theorem B (index formula).**  For every closing twist $g\in W(E_8)$,
$$[\Lambda:S_g]=[Ls:\langle Ls\cdot_gLs\rangle]\cdot[L\bar s:\langle
L\bar s\cdot_gL\bar s\rangle]^2,$$
and $\Lambda/S_g\cong Ls/\langle Ls\cdot_gLs\rangle\oplus\bigl(L\bar s/\langle
L\bar s\cdot_gL\bar s\rangle\bigr)^2$.  On all 4,704 closing pairs
(Section 3.3) these quotients are elementary abelian 2-groups, so
$2\Lambda\subseteq S_g$ for every closing $g$.  This last part is
computed, not proved.

*Proof sketch* (steps 1 to 3 machine-checked in `f7`; the premise
$\langle L\cdot_gL\rangle=L$ of step 4 is checked for every twist in
`f2` and `f4`).
1. $\mathbb R^{24}=D\oplus T$ with $D=\{(a,a,a)\}$ and
   $T=\{p+q+r=0\}$.  For any bilinear block product,
   $D\star T=T\star D=0$, and $(a,a,a)\star(b,b,b)=3(ab,ab,ab)$.
2. Let $T_{\mathbb Z}=\{\alpha\in\mathbb Z^3:\alpha_1+\alpha_2+\alpha_3=0\}$,
   and let $\circ$ be the routing product on $\mathbb Z^3$,
   $(\alpha\circ\beta)_1=\alpha_1\beta_1+\alpha_3\beta_2+\alpha_2\beta_3$
   and cyclically.  On $T=\mathbb O\otimes T_{\mathbb Z}$ the product is
   $xy\otimes(\alpha\circ\beta)$, with $[T_{\mathbb Z}:\langle
   T_{\mathbb Z}\circ T_{\mathbb Z}\rangle]=3$.  Write
   $\Delta(a)=(a,a,a)$.
3. $\Lambda\cap D=\Delta(Ls)$ and $\Lambda\cap T=L\bar s\otimes
   T_{\mathbb Z}$, and $[\Lambda:(\Lambda\cap D)+(\Lambda\cap T)]=3^8$,
   a 2-adic unit.
4. $(2a,0,0)\star(2b,0,0)=(4\,a\cdot_g b,0,0)$ and $\langle
   L\cdot_gL\rangle=L$ give $(4L)^3\subseteq S_g$, so the index is a
   power of 2.

Because the gluing index $3^8$ is odd, $\Lambda$ splits 2-adically:
$\Lambda\otimes\mathbb Z_2=(\Lambda\cap D)\otimes\mathbb
Z_2\oplus(\Lambda\cap T)\otimes\mathbb Z_2$.  By step 1 the span splits
the same way, into two pieces:
$$\langle(\Lambda\cap D)\star(\Lambda\cap D)\rangle=3\Delta(\langle
Ls\cdot_gLs\rangle),\qquad \langle(\Lambda\cap T)\star(\Lambda\cap
T)\rangle=\langle L\bar s\cdot_gL\bar s\rangle\otimes\langle T_{\mathbb
Z}\circ T_{\mathbb Z}\rangle .$$
With $h_W=[Ls:\langle Ls\cdot_gLs\rangle]$ and $h_V=[L\bar s:\langle L\bar
s\cdot_gL\bar s\rangle]$, their indices are $3^8h_W$ and
$h_V^2\,3^8$:
- the exponent 2 is the rank of $T_{\mathbb Z}$;
- $3^8=[T_{\mathbb Z}:\langle T_{\mathbb Z}\circ T_{\mathbb
  Z}\rangle]^{\operatorname{rank}L\bar s}$.

As $2L\subseteq U$ gives $4L=\langle2L\cdot_g2L\rangle\subseteq\langle
U\cdot_gU\rangle$ for $U=Ls,L\bar s$, both $h_W$ and $h_V$ are powers of
2.  So the 2-parts of the two indices are $h_W$ and $h_V^2$.  By step 4,
$[\Lambda:S_g]$ is its own 2-part, which gives the formula.  It was
checked against direct normal-form computations on every closing twist
computed: 1,764 permutations, 56 reflections and 4,704 polarizations.
$\square$

## 3. Results

### 3.1 The paper's own family: permutation twists

All 5,040 permutations were computed exactly.  The closing counts match
the paper's Section 5.6 table in all 14 cycle-type rows (1,764 of 5,039;
the identity fails).  New is the span index:

| closing permutations | number | $[\Lambda:S]$ | $\Lambda/S$ | factors $(h_W,h_V)$ |
|---|---|---|---|---|
| odd (types 2, (3,2), (2,2,2), 6, (4,3), (5,2)) | 588 | $2^{16}$ | $(\mathbb Z/2)^{16}$ | $(1,2^8)$ |
| even (types 3, (2,2), 5, (3,3), (4,2), (3,2,2), 7) | 1,176 | $2^{10}$ | $(\mathbb Z/2)^{10}$ | $(1,2^5)$ |

Here $h_W=[Ls:\langle Ls\cdot_gLs\rangle]$ and
$h_V=[L\bar s:\langle L\bar s\cdot_gL\bar s\rangle]$.

- **Witness for $2^{10}$:** the 3-cycle $\pi(5)=6$, $\pi(6)=7$,
  $\pi(7)=5$ in the paper's convention.  Its exact Hermite normal form
  has pivots $1^{14}2^{10}$, and PARI's Smith form has ten 2s.  The
  reverse orientation does not close.
- **$\sigma$:** $2^{16}$, with sixteen 2s, as in the paper.  Its span is
  explicit: $S_\sigma=\Delta(Ls)+2\Lambda=\{(x,y,z):x,y,z\in Ls;\ x-y,\,
  x-z\in2L\bar s\}$.  All 588 closing odd permutations and the 28
  closing integral reflections below share this one span lattice; the 1,176 even ones
  give 14 distinct span lattices.
- The paper's 3-cycle remark (exactly one orientation per 3-subset
  closes) is confirmed.

### 3.2 Root reflections, and index 1

There are 120 reflections, one per $\pm r$.  The pattern is clean
(`f4`, independently in GAP and PARI):

| root $r$ | $\langle r,s\rangle$ | reflections | closes | $[\Lambda:S]$ |
|---|---|---|---|---|
| $e_i-e_j$ ($1\le i<j\le7$; includes $\sigma$) | 0 | 21 | yes | $2^{16}$ |
| $e_0+e_j$ | 0 | 7 | yes | $2^{16}$ |
| $e_0-e_j$ | $-1$ | 7 | no | – |
| $e_i+e_j$ | 1 | 21 | no | – |
| $\tfrac12(e_0+\sum\varepsilon_ie_i)$, one minus sign | 1 | 7 | yes | **1** |
| same, five minus signs | $-1$ | 21 | yes | **1** |
| same, three minus signs | 0 | 35 | no | – |
| same, seven minus signs ($r=-s$) | $-2$ | 1 | no | – |

So integral roots close exactly when $\langle r,s\rangle=0$, with index
$2^{16}$, and half-integral roots close exactly when $\langle
r,s\rangle=\pm1$, with index 1.

**The index-1 witness, in detail.**  $r=\tfrac12(e_0+e_1+\dots+e_6-e_7)$
is a root of $L$ (half-integral, coordinate sum 3), and $s_r(L)=L$.  The
twisted product $x\cdot_ry=s_r(s_r(x)\cdot s_r(y))$ is an octonion
product with identity $1_r=s_r(e_0)=\tfrac14(3e_0-e_1-\dots-e_6+e_7)$,
which is not in $L$.
- **Closure.**  All 576 basis-pair products lie in $\Lambda$, by three
  tests that agree product by product: integral $\Lambda$-coordinates,
  Wilson's conditions mod 2, and the definitional test "$t\in L\bar s$
  iff $(t\,s)/2\in L$".  Random lattice pairs were also checked.
- **Index 1, by five routes:**
  1. an exact integer Hermite normal form equal to the identity;
  2. $\mathbb F_2$-rank 24, which makes the index odd, combined with
     $(4L)^3\subseteq S_r$ (Theorem B, step 4), which makes it a power
     of 2;
  3. PARI `mathnf`/`matsnf`, with all elementary divisors 1;
  4. GAP elementary divisors $1^{24}$;
  5. the index formula, with factors $(1,1)$.
- **Certificates that use no normal form:**
  - *(a) A $\mathbb Z$-basis of $\Lambda$ made of 24 single products*
    $u_k\star_rv_k$, with $u_k,v_k\in\Lambda$.  All 48 factors and all
    24 products satisfy Wilson's definition, the coordinate determinant
    is $-1$, and the Gram determinant of the products is $2^{24}$
    (`reproduction/independent/refute/c3_certificate.gp(.out)`,
    re-verified in `f3`).
  - *(b) Two sets of 24 basis-pair products* with odd, coprime
    determinants: $79{,}194{,}800{,}595{,}163=11\cdot7{,}199{,}527{,}326{,}833$
    and $-81{,}930{,}166{,}777{,}915{,}407$.  The index divides their
    gcd, which is 1 (`cert_index1_final.gp`, re-checked from scratch in
    PARI by `f6`).
  - *(c)* GAP's two sets, with determinants $2^{44}$ and
    $3^4\cdot5\cdot7\cdot8{,}961{,}598{,}279$; gcd 1.

**Why index 1 is not trivial.**  For any twist and any
$u=(x,y,z)\in\Lambda$, left multiplication $L_u$ is block-diagonal for
$D\oplus T$.
- Its $D$-block is left multiplication by $x+y+z$ in
  $(\mathbb O,\cdot_g)$, with determinant $N(x+y+z)^4$, which is
  divisible by 16 because $L$ is even.
- Its $T$-block maps the lattice $L^3\cap T=L\otimes T_{\mathbb Z}$ into
  itself, because $u\in L^3$ and $L\cdot_gL\subseteq L$, so its
  determinant is an integer.

So $\det L_u$ is divisible by 16 (and is 0 when $x+y+z=0$), and the same
holds on the right.  This was checked exactly on 924 vectors, left and
right, for $\sigma$ and $s_r$: no unimodular case.  So $S=\Lambda$ does
not come from a unit-like element: no single left or right
multiplication maps $\Lambda$ onto $\Lambda$.

### 3.3 The whole group $W(E_8)$: spectrum and rule

The reflections act on $\bar L$ as transvections and reach all 17,280
ordered complementary pairs of Lagrangians from Wilson's pair.  So every
pair is $(V_g,W_g)$ for exactly $|W(E_8)|/17{,}280=40{,}320$ elements
$g$.

Applying $g$ blockwise is a ring isomorphism from $(\Lambda,\star_g)$
onto the lattice defined like $\Lambda$ with $g(L\bar s)$ and $g(Ls)$ in
place of $L\bar s$ and $Ls$, under the untwisted product.  So closure
and index of $\star_g$ depend only on $(V_g,W_g)$, and a computation
over the pairs covers all of $W(E_8)$.  Over the pairs (`f5`,
exhaustive, by direct lattice computation on each):

- there are 270 Lagrangians in two families of 135; 142 are
  subalgebras; 16 are right ideals, and 8 of these are two-sided (they
  are also the only left ideals);
- 4,704 of the 17,280 ordered pairs consist of two subalgebras, and
  exactly these close (Theorem A);
- **the index spectrum is** $1:3{,}360$, $2^5:448$, $2^8:224$,
  $2^{10}:448$ and $2^{16}:224$ closing pairs.  In terms of group
  elements, 189,665,280 of the 696,729,600 elements of $W(E_8)$ close,
  and 135,475,200 of them (19.4% of the group) give $S=\Lambda$.

**Rule (exhaustive over the 142 subalgebra Lagrangians).**  Let
$\varepsilon$ be the class of $2e_0$ in $\bar L$.  Since
$V_g=g(L\bar s)/2L$, $\varepsilon\in V_g$ if and only if
$2\cdot1_g=g^{-1}(2e_0)\in L\bar s$; likewise $\varepsilon\in W_g$ if and
only if $2\cdot1_g\in Ls$.  For a closing $g$:
$$[\Lambda:S_g]=\begin{cases}2^{2k}&\varepsilon\in V_g\ (2\cdot1_g\in L\bar s),\\
2^{k}&\varepsilon\in W_g\ (2\cdot1_g\in Ls),\\ 1&\text{otherwise,}\end{cases}
\qquad k=\begin{cases}5&\det g=+1,\\8&\det g=-1.\end{cases}$$

Every coordinate permutation fixes $e_0$, and $2e_0=s\bar s\in L\bar s$,
with $\det g=\operatorname{sign}\pi$.  This is the parity law of 3.1:
permutation twists can never go below $2^{10}$, and $\sigma$ (odd) sits
at the top of the spectrum.  The index-1 reflections move $e_0$ so that
$2\cdot1_g$ lies in neither sublattice.

### 3.4 Signed permutations

The 645,120 signed permutations of $e_1,\dots,e_7$ give 480 distinct
twisted tables, of which 168 close (112 at $2^{10}$, 56 at $2^{16}$).
This adds nothing new.  The 240 tables from maps preserving $L$ are
exactly the permutation tables.  The other 240 come from maps with
$g(L)\neq L$, and each is the opposite table $T'[a,b]=T[b,a]$ of a
permutation table, so $u\star'v=v\star u$ and closure and index are
identical.

## 4. Consequences for the paper (options for the author)

1. **Section 7, the span-index bullet.**  It cannot stay as an open
   question.  This is the second time the question has been settled,
   after the unrestricted form.
   - *Minimal, with no new content:* delete it.  Both indices,
     $2^{16}$ and $2^{25}$, already appear in Sections 5.5 and 6.
   - *Factual replacement* (needs a repository script under your
     protocol): "The span index depends on the twist.  Among the closing
     permutations of Section 5.6, every even one gives index $2^{10}$
     and every odd one $2^{16}$, as $\sigma$ does.  Replacing $\sigma$
     by the reflection in the root $\tfrac12(e_0+e_1+\dots+e_6-e_7)$ of
     $L$ gives a product that closes on $\Lambda$ and whose products
     span $\Lambda$ (index 1)."
2. **Section 5.5** should no longer read as if $2^{16}$ were a feature
   of the construction; it is a feature of $\sigma$.  The sentence "A
   span-surjective product would have index 1" can stay, since such
   products now exist.
3. **Optional, substantial: a short section on twists by isometries of
   $L$**, containing:
   - Theorem A, with its proof;
   - Theorem B, with a proof to be written out from the sketch;
   - the spectrum and the rule, as exhaustive computations;
   - the remark that the analogue of Lemma 4.3 ($g(L\bar s)$ a left
     ideal of $L$) holds, among closing twists, exactly for those of
     index $2^{16}$, and fails for all others.

   This answers the lifting question of Appendix A.1 for this family.
   It also gives part of the structural "why" that the Journal of
   Algebra referee asked for: a criterion for closure and a rule for
   the index, though not a reason why $\sigma(Ls)/2L$ is a subalgebra.  That cuts both ways with the paper's stated
   positioning (THAT, not WHY): the results are machine-obtained, with
   short proofs the author can check, and would have to be marked as
   such.  It is your decision.
4. **If any of this is used:** cite the convention (left multiples;
   $\pi(e_a)=e_{\pi(a)}$), add the scripts to the repository, and
   record the provenance (this review).

**Still open.**
- Why $\sigma(Ls)/2L$ is a subalgebra: the paper's own first open
  question, still certified only by Table 3.
- The exact values in the rule.  Let $U\subseteq\bar L$ be a subalgebra
  Lagrangian with full preimage $\tilde U\subseteq L$, and put
  $h(U)=[\tilde U:\langle\tilde U\tilde U\rangle]$ (untwisted product).
  A conceptual (non-exhaustive) proof is open that $h(U)=1$ when
  $\varepsilon\notin U$, and that otherwise $h(U)$ is exactly $2^5$ or
  $2^8$, according to the family of $U$.  Only the lower bound $h\ge2$
  for $\varepsilon\in U$ has a structural proof: $\tilde U$ then has
  integer coordinates and $\operatorname{Re}(xy)$ is even for all
  $x,y\in\tilde U$, while $2s\in2L\subseteq\tilde U$ has real part $-1$.
- What the index-1 orders $(\Lambda,\star_g)$ are as rings, and whether
  they are better objects for applications.  This was not examined.

## 5. What is proved, computed, or open

| statement | status |
|---|---|
| Closure criterion (Theorem A), both directions | proved (short proof above; identities also checked symbolically) |
| Index formula (Theorem B) | the formula: proof sketch with machine-checked steps; the elementary abelian quotient ($2\Lambda\subseteq S_g$): computed only; both hold on all 4,704 closing pairs, hence for every closing $g\in W(E_8)$ |
| Closure of all 5,040 permutation twists, 120 reflections and 17,280 polarizations, and the span indices of the closing ones (1,764, 56 and 4,704) | exhaustive exact computation, reproduced independently |
| Lemma 4.3's analogue holds exactly in the $2^{16}$ class; $V_g$ neither left nor right ideal on 3,136 index-1 pairs | exhaustive (`f8`) |
| Index 1 for $s_r$, $r=\tfrac12(e_0+\dots+e_6-e_7)$ | certified: an explicit $\mathbb Z$-basis of 24 products, plus coprime-determinant certificates |
| No unimodular multiplication | proved; checked on 924 vectors |
| Spectrum over $W(E_8)$ and the rule via $1_g$ | exhaustive over all pairs; a conceptual proof of the exact values is open |
| Convention dependence (right multiples flip the parity law) | computed |

## 6. Other findings on v8

These come from the review workflow of 13 October 2026, which had four
reviewers (mathematics, consistency, build and record, journal fit) and
a skeptic who tried to refute every finding.  **Confirmed** means real
as stated; **partly** means real but overstated, or the first proposed
fix was wrong (the fix given here is the corrected one).  Line numbers
refer to v8 (`main.tex`, 13 October 2026).

### 6.1 Mathematics and wording in the paper

**(a) Section 1, the seven maximal orders** (medium; confirmed).
"For a fixed octonion multiplication, the integral octonions have seven
maximal orders, each an $E_8$ lattice" is false as stated.  The rational
octonion algebra has infinitely many maximal orders: $G_2(\mathbb Q)$ is
infinite, and an automorphism such as conjugation by $q=1+e_1+e_2+e_4$
maps one of the seven to a maximal order that contains $e_0$ but not
$e_3$.  The classical count (Coxeter 1946; Conway–Smith 2003; Appendix
B's "the true count is seven") is of the maximal orders **containing the
Gravesian order** $\mathbb Z\langle e_0,\dots,e_7\rangle$.  This was
verified exactly: of the 30 $E_8$ lattices $\mathbb Z^8+C/2$ ($C$ a
doubly-even $[8,4,4]$ code), exactly 7 are closed under the paper's
multiplication.  Exactly 8 are stable under the units, the eighth being
Kirmse's $J_1$, and the $\sigma$-twisted multiplication has its own,
disjoint, 7.  Proposed wording (lines 136–138): "For a fixed octonion
multiplication, exactly seven maximal orders contain the Gravesian order
$\mathbb Z\langle e_0,\ldots,e_7\rangle$, each a copy of the $E_8$
lattice [Coxeter1946, ConwaySmith2003] (Appendix B); ...".  This is SF4,
which you marked for review in context.

**(b) Section 7, the non-maximality paragraph** (math correct; prose
low).  Every claim was re-checked exactly:
- $L^3\star L^3\subseteq L^3$, and $[L^3:\Lambda]=2^{12}$;
- $\varepsilon_1$ acts on either side as the swap of blocks 2 and 3,
  and $\varepsilon_1\star\varepsilon_1=\varepsilon_1$;
- $k\varepsilon_1\in\Lambda$ first at $k=4$;
- $\langle\varepsilon_1,b\rangle=-\tfrac12$ for a basis vector $b$.

Only the prose needs work.  "A second ⋆-closed overlattice is not
contained in $L^3$" asserts a property before naming the lattice, and
the next sentence runs on.  Proposed wording (lines 1051–1059): "A
second ⋆-closed overlattice is $\Lambda+\mathbb Z\varepsilon_1$.  The
idempotent $\varepsilon_1=(e_0,0,0)$ of Section 5.2 acts on either side
as the swap of blocks 2 and 3, $\varepsilon_1\star v=v\star\varepsilon_1=(x',z',y')$
for $v=(x',y',z')$.  This swap preserves $\Lambda$, because Wilson's
three conditions are symmetric in the blocks, and
$\varepsilon_1\star\varepsilon_1=\varepsilon_1$, so
$\Lambda+\mathbb Z\varepsilon_1$ is ⋆-closed by bilinearity."  The rest
stays.  This is SF3, which you also marked for review.

**(c) Abstract** (low; partly).  The headline result of Section 5.3, the
automorphism group $C_6\times S_3$ of the *order*, is not what the
abstract names; the clause reads as if the group were the algebra's,
which is $G_2\times G_2\times S_3$.  Proposed wording (lines 104–107):
"We determine the idempotents and the square-zero elements of the
algebra, the automorphism group $C_6\times S_3$ of the order, and the
index in $\Lambda$ of the span of its products."

**(d) Appendix A.1, the lifting question** (medium; partly).  It is
answered for twists of the paper's kind by Theorem A of Section 2.  At
minimum (lines 1638–1640): "... and which of them lift to a closed
bilinear product on $\Lambda$ itself; the closing permutations of
Section 5.6 already give 84, the one constructed here among them."

**(e) Section 5.1 lead-in** (low; confirmed).  It lists associativity,
which has no row in the table.  Proposed: "To give an initial idea of
the degree of violation of these and related properties, ...".  Also
add serial commas at line 604 ("power-associativity, and
commutativity") and line 715 ("$16k^2$, $48k^2$, or $24k^2$,
respectively").

**(f) Remark 1.3** (low; confirmed).  "on $10^6$ sampled distinct
pairs": the sampler does not enforce distinctness.  Write "sampled
pairs".

**(g) Hyperref warnings** (low; confirmed).  `\texorpdfstring` in the
title ($\mathbb Z_3$), in the subsection "The $E_8$ lattice" and in
"each block lies in $L$" removes the 8 hyperref warnings that v7 also
has (three from the PDF title, five from the two bookmarks), and fixes
the PDF title metadata.

### 6.2 The AI disclosure and the public record

**(h) "Every prompt ... end to end"** (medium; confirmed; your call).
The paper says the complete record is public: "every prompt,
computation, and verification", "the complete list of prompts" and
"all prompts".  The public log has some gaps:
- no files 073 and 112 (073 was interrupted and restarted as 074; 112
  survives as one line inside 111);
- elisions "[...]" under "User prompt (verbatim)" in 025 and 026;
- about 213 words of the user prompt in log 281 are verbatim passages
  of a non-public AI review (Claude Opus 5.5, 9 to 12 October 2026,
  `J_Alg_Comb/review_findings.md`) without that provenance.

Two routes:
- *Narrow the claims.*
  - Introduction: "the research record is public: the logged prompts,
    computations, and verifications, end to end".
  - Appendix C: "the prompt log", and "the prompt log, code, evidence,
    and git history".
  - Declarations: "the research record".
- *Disclose the gaps,* including a dated note in a log 285, e.g.: "The
  passages of prompt 281 on the ε₁ witness, the drafting sentence, the
  Hermes sentence, Python and GAP, the A.1 clause, and the Use-of-AI
  sentences were adopted from an AI-assisted review of the v7 journal
  package (Claude Opus 5.5, 9 to 12 October 2026, outside this
  repository); the $L^3$ observation first appears there too."

This review is itself a candidate for that record.

**(i) Models named for v8** (medium; confirmed).  The record names
Claude Fable 5.1, the reviewer at prompt 283.  It does not name the
model that edited the paper in sessions 281, 282 and 284.  The journal
AI statement needs both.  In Appendix C, "the AI agent (Claude Opus
4.7)" reads as the model for the whole project; "(at that stage Claude
Opus 4.7)" fixes that.

**(j) Who ran Hermes** (low; confirmed; your call, ideally with
Matthew Barley's consent).  Log 282 records that Matthew Barley ran the
Hermes checks.  Appendix C does not connect him to them, so a reader
may assume you ran Hermes.  Possible wording: "Structural proposals that
Matthew Barley obtained from a second AI system (...) and submitted to
the repository prompted ...".

**(k) The forwarded v8 is not the public v8** (medium as a
pre-submission step; confirmed).  The public `paper/main.tex` (head
59cd53a) says "12 October 2026 (v8)", and your copy says "13 October".
That is the only difference, and no log records it.  v8 is not frozen
(no `main_2026-10-13.tex`; the changelog says "in progress, not
frozen"), and neither mirror has a tag.  Before posting or submitting,
three steps:
- log the date change (285);
- freeze v8;
- tag the commit on both mirrors, and name the tag in the journal's
  data availability statement.

**(l) Public trackers** (low; partly).  README "What Remains Open" and
`CURRENT_STATE.md` open questions 1 and 3 repeat the span-index and
lifting questions; align them with whatever replaces the paper's text.
`CURRENT_STATE.md` still carries the date stamp 2026-10-12 although it
changed on 13 October.  The `CLOSING_NOTE.md` addendum leaves reason 1
("Journal submission is not the productive path right now")
unqualified.  Several files call the ResearchGate v7 "published".

**(m) Script docstring** (low; confirmed).
`verify_discovery1_W_subalgebra.py`, which v8 now links, still carries
the proposal-era expectations: that $L/2L$ is commutative, and that
$e_0$ is its identity.  The vector it labels $e_0$ is $L_1=s$ in
$L$-coordinates.  A dated docstring correction, as for
`symbolic_proof_checks.py` at prompt 283, fixes it.

### 6.3 Journal fit after v8

Not refutation-checked; this is advice.  v8's abstract dropped the
mod-2 Lagrangian sentence and contains no combinatorics.  Editors judge
scope from the cover letter, abstract, conclusion and references.  So
JACo becomes the weakest of the four candidates, and v8 strengthens the
case for Advances in Applied Clifford Algebras (AACA).

Its remaining risk is scope.  The paper never mentions Clifford
algebras and cites no AACA paper.  The cheapest honest de-risking is one
clause in the introduction citing Dixon's own AACA paper: G. Dixon,
"Division algebras, Clifford algebras, periodicity", Adv. Appl.
Clifford Algebras 28 (2018), art. 13, doi:10.1007/s00006-018-0820-8.
It proposes a Clifford-algebra periodicity of order 24 built over the
four division algebras.

The results of Sections 2 and 3, should they enter the paper, would
add structural content that also helps the case at Results in
Mathematics.

## 7. Reproducing the computations

Everything is in `reproduction/` next to this file, with relative paths
only.

- `reproduction/final/` is the consolidated pipeline: run `./run_all.sh`.
  It needs Python 3 with numpy and sympy, plus PARI/GP, and takes about
  2 minutes on 14 cores.  It reads `../../../main.tex` (override with
  the environment variable `MAIN_TEX`) to pin the conventions to the
  paper's printed tables.

  | script | what it does |
  |---|---|
  | `leech_final.py` | the library |
  | `f1_setup` | conventions and the lattices |
  | `f2_witnesses` | the witnesses, with exact HNF |
  | `f3_certificates` | normal-form-free certificates, including re-checks of the other agents' certificates |
  | `f4_census` | all 5,040 permutations in both conventions, 120 reflections, 645,120 signed permutations |
  | `f5_polarizations` | all 17,280 Lagrangian pairs, the spectrum, the rule |
  | `f6_pari_independent.gp` | an independent PARI/GP implementation |
  | `f7_structure_checks` | the proofs' identities, the D⊕T splitting, the no-unit test, the ideal properties, convention variants |
  | `f8_ideal_classes` | the ideal type of $V$ in every closing pair (added after the document check) |

  Each script's output is saved as `*.out.txt`.  After the folder was
  moved here, an unchanged rerun reproduced all ten generated files
  byte for byte: the seven outputs of f1 to f7, `f2_matrices.gp`,
  `cert_index1_final.gp` and `lambda_gram_final.gp`.  The run on 13
  October 2026 started at 01:03 UTC on 14 October.  Since then, f1's
  label for the paper's remark was corrected from 4.6 to 4.5, and f8
  was added.  The final complete rerun (14 October, 01:43 UTC) took about
  two minutes and regenerated all eleven files identically.
- `reproduction/independent/{py,gap,refute,structure}/` holds the four
  independent implementations, each written from the paper's
  definitions without seeing the others:
  - Python, with PARI/GP for `qfminim` and the cross-checks (`py`);
  - GAP 4 (`gap`; `run_*.g` with outputs; Leech minimum and the 196,560
    vectors by `ShortestVectors`);
  - PARI/GP and Python (`refute`; including the 24-product
    certificate);
  - Python with PARI (`structure`).

  Five large regenerable data files (about 12 MB) were not copied.
- `original_reviewer_scripts/` holds the first reviewer's working
  scripts, which found the 2¹⁰ and index-1 cases, kept for provenance.

Tool versions: Python 3.12.3, numpy 1.26.4, sympy 1.12, PARI/GP 2.15.4,
GAP 4.12.1.

## 8. Provenance and limits

- **Who did what.**  Claude Opus 5.5 orchestrated two workflows on 13
  October 2026:
  - a review of v8: four reviewers and one skeptic;
  - a reproduction: four independent implementations and one
    adjudicating skeptic, which reran the decisive checks itself.

  The first finding of index $2^{10}$ and index 1 came from the
  mathematics reviewer.  The headline numbers were then reproduced by
  at least three implementations in at least two systems: the census
  counts, the indices $2^{16}$, $2^{10}$ and 1, the reflection table,
  the spectrum, and the signed-permutation counts.

  Some results were reproduced in Python only:
  - the rule via $1_g$ (`final` and `independent/structure`);
  - the 14 distinct span lattices and the shared $2^{16}$ span lattice
    (`f4` and `independent/py`);
  - the right-multiple results (the permutation census in `f7` alone);
  - the ideal-type counts (`f8`).

  The numbers in Section 6 come from the review workflow, not from the
  reproduction.
- **Limits.**  Nothing here has been peer-reviewed.  The proof of
  Theorem B is a sketch with machine-checked steps, not a written
  proof.  The exact values in the rule are known only by exhaustion.
  As the paper's own methodology appendix says, errors are possible on
  both sides.  The scripts are included so that every statement can be
  checked independently.
