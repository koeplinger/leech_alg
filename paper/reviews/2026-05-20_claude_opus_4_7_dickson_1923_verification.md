# Primary-source verification — Dickson (1923), §19

**Date:** 2026-05-20
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 109. Second source paper of the 1920s–40s set (after Kirmse).
**Source read:** L. E. Dickson, *A new simple theory of hypercomplex
integers*, J. Math. Pures Appl. (9) **2** (1923), 281–326 — §19, *The
Integers of Cayley's Algebra*, pp. 319–324. The paper is in English (French
journal). Read from the Numdam scan.
**Reproducible check:** `python_project/src/verify_dickson_1923.py` — Dickson's
Cayley–Dickson product and his order (30) built solely from his paper.

---

## What Dickson did in §19

Dickson builds the octonions by the Cayley–Dickson doubling of the
quaternions: $x = q + Qe$ with $q,Q$ quaternions, and (his eq. 24)
$$(q+Qe)(r+Re) = (qr - \bar R Q) + (Rq + Q\bar r)\,e.$$
He then **determines the maximal orders** ("maximal sets of integers … having
properties C, U, N"). From the norm he derives that integer coordinates are
integers or half-integers, runs a case analysis, and arrives at an explicit
8-element basis (his (29): $i,j,k,p,e,W,Z,V$, with $p=\tfrac12(1{+}i{+}j{+}k)$)
spanning an order, his **system (30)** — call it $O_D$. He **asserts its
closure** ("It can be verified that this is true also of the product of any
two of the numbers (29)"), **proves it maximal**, and states:

> **Theorem XV.** *The only maximal systems of integers of Cayley's algebra
> having properties C, U, N are (30) and the two systems obtained from it by
> cyclic permutations of $i,j,k$.*

That is: Dickson claims there are **three** maximal orders.

His axioms (his §6): **C** = closure under multiplication, **U** = contains
the unit $1$, **N** = norm is a rational integer.

## Verification from scratch

Using only Dickson's product (eq. 24) and his basis (29):

1. **The product is a genuine octonion algebra** — $N(xy)=N(x)N(y)$ on every
   tested pair; basal units multiply to $\pm$ basal units.
2. **Dickson's order $O_D$ is genuinely closed** — all 64 products of the
   basis (29) lie in $O_D$. *Dickson's construction is correct:* his system
   (30) is a real maximal order. (This is the part Kirmse got wrong; it
   confirms Petersson's remark that Dickson's integral octonions are
   isomorphic to Coxeter's.)
3. **There are exactly seven maximal orders** for Dickson's product (30
   $E_8$-type candidate lattices; 7 closed) — not three. All seven satisfy
   C, U, N (each is closed, contains $J_0\ni 1$, and is an even lattice so has
   integral norm). Since Dickson's derivation also shows every octonion
   integer has half-integer coordinates, every "maximal system having C,U,N"
   is one of these seven rank-8 lattices.

**Theorem XV undercounts: it states three where there are seven.**

## Why three — and this time Dickson states the cause himself

In the §19 derivation (p. 321), having reached a case in which the order
contains $p=\tfrac12(1{+}i{+}j{+}k)$, Dickson writes:

> *"… so that the set contains (27) $p=\tfrac12(1+i+j+k)$, and hence all of
> Hurwitz's integral quaternions. **We shall assume that the latter occur in
> our set in all cases.**"*

This assumption — that every maximal order contains the Hurwitz quaternions,
i.e. contains $p$ — is **not** part of C, U, N; it is an extra restriction.
And it is false: of the seven maximal orders, **exactly three contain $p$**,
and those three are precisely $O_D$ and its two $(i\,j\,k)$-cyclic images —
exactly Dickson's Theorem XV triple. The four he missed are exactly the
maximal orders that do not contain $p$. (Under the $(i\,j\,k)$ 3-cycle Dickson
used to organise his case analysis, the seven split into orbits of sizes
$1,3,3$; his argument reached one orbit of three.)

So Theorem XV, **as printed, is false**; what Dickson's argument actually
establishes is the true *restricted* statement — *the only maximal orders
**containing the Hurwitz quaternions** are (30) and its two cyclic images*.
The qualifier, stated openly in the derivation, was dropped from the theorem.

## The two 1920s papers side by side

| | construction exhibited | claimed count | true count |
|---|---|---|---|
| **Dickson 1923** | $O_D$ — **a genuine maximal order** | **3** (Thm XV) | 7 |
| **Kirmse 1924** | $J_1$ — **not closed**, not an order | **8** (p. 70) | 7 |

Both miscount the maximal orders, in opposite directions, and the true number
seven lies between. Dickson's miscount comes from an *explicitly stated*
extra assumption (contains the Hurwitz quaternions) too strong by four;
Kirmse's, from a *silent* criterion (unit-invariance) too weak by one. The
constructions diverge sharply: Dickson's exhibited order is correct and
predates Kirmse's paper; Kirmse's exhibited order is the defective one.

## Status

Dickson (1923): **verified.** Kirmse (1924): verified (separate note).
Pending: Mahler (1942) — loaded; Coxeter (1946) — to be obtained. The
question of who first stated the count *seven* correctly remains open and is
deferred to Coxeter.

Facts vs. interpretation: the counts (7 maximal orders; 3 containing $p$; the
identity of those 3), the closure of $O_D$, and the orbit sizes are **verified
facts**. The identification of the *cause* of Theorem XV's undercount is, in
this case, not a reconstruction — Dickson states the offending assumption in
his own text; the verification only confirms that it accounts exactly for the
gap (3 = the orders satisfying it).
