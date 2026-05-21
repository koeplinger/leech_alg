# Forensic note — can Kirmse's "eight" maximal orders be constructed from his paper?

**Date:** 2026-05-20
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 103. Follows the verification in prompt 100
(`2026-05-20_..._kirmse_1924_verification.md`).
**Reproducible check:** `python_project/src/verify_kirmse_1924_forensic.py`.

## The question

Kirmse *states* there are eight maximal orders and then works with one
concrete lattice, $J_1$, that is not an order. That alone does not tell us
what his enumeration was, nor whether his paper contains enough to reconstruct
it. The question: **are there sufficient hints in Kirmse's paper to construct
the eight — and, done so, what do they turn out to be?**

## Part 1 — the explicit text does not construct them

A forensic re-read of Kirmse's §2 (pp. 68–71) shows the eight orders are
treated on **p. 70 only**, and there only as follows:

- the assertion that "a simple finite computation" finds *all* maximal orders,
  and that *"Man findet acht solche"* — one finds eight;
- a schematic form, *"deren jeder die Gestalt hat"*
  $\;J=[\,J_0,\ \tfrac12(i_{\nu_1}{+}i_{\nu_2}{+}i_{\nu_3}{+}i_{\nu_4}),\ \tfrac12\sum i_k\,]$;
- the structural remark that each order contains 14 "middle-kind" elements in
  7 complementary pairs;
- exactly **one** order given by an explicit generating set: $J_1$.

The finite computation itself is **not printed**, not even sketched, and the
eight are **not listed**. Moreover the schematic, read literally in Kirmse's
own bracket notation (defined p. 68 as the $\mathbb{Z}$-span), is
$J_0$ together with *two* half-integer elements — a module of index at most
$4$, not the index-$16$ of a maximal order. So the schematic cannot itself be
a construction; it states the *form* of the elements, not a recipe.

**Conclusion of Part 1:** Kirmse's explicit text does not construct the eight.
It asserts their number and form and exhibits one example — the defective
$J_1$.

## Part 2 — the eight are nonetheless reconstructible, and the reconstruction is decisive

Although Kirmse printed no computation, the following facts (table (1) only;
exact arithmetic) reconstruct it essentially uniquely.

A maximal order over $J_0$ is an even unimodular lattice, hence corresponds to
a doubly-even self-dual binary $[8,4]$ code. There are **30** such candidates.
Among them:

| condition on the candidate lattice | count |
|---|---|
| closed under multiplication (a genuine maximal order) | **7** |
| invariant under multiplication by the **units** | **8** |

and — the decisive point —

- the **8** unit-invariant lattices are *exactly* the **7** genuine maximal
  orders **together with Kirmse's $J_1$**;
- $J_1$ is the **unique** unit-invariant lattice that is not an order;
- $J_1$'s failure of closure is confined entirely to the products of two
  half-integer elements: of Kirmse's own generators $a_1,\dots,a_4$, ten of
  the sixteen products $a_i\cdot a_j$ leave $J_1$, while **all 128 products
  involving a unit close**.

So the number "eight" is not arbitrary. It is exactly the count of lattices
invariant under multiplication by the units — that is, lattices closed under
every product *except possibly* the product of two half-integer elements.
Full closure (a genuine order) additionally requires the half $\times$ half
products to land back; exactly seven of the eight pass that further test.

## Interpretation

The reconstruction is an inference — Kirmse printed no computation — but it is
supported by three independent matches: the unit-invariance test yields his
exact count (8); the eight it yields contain exactly the seven genuine orders;
and its unique defective member is exactly the lattice $J_1$ he chose to
exhibit. The natural reading of Kirmse's *"einfache endliche Rechnung"* is
therefore: he verified closure under the units — the products $e_i\cdot h$ and
$h\cdot e_i$ — and did not carry the test through to the products $h\cdot h'$
of two half-integer elements. That single missing condition is exactly where
$J_1$ fails (witness $a_1\cdot a_3=\tfrac12(i_1{+}i_2{+}i_4{+}i_7)\notin J_1$),
and it is what separates his eight from the true seven.

## Consequence for attribution

Kirmse's enumeration was an enumeration of **unit-invariant lattices**, not of
orders. Seven of his eight happen to be genuine maximal orders; one is not;
and the only one he wrote down explicitly is the defective one. His paper
therefore does **not** establish him as having identified, constructed, or
verified the seven maximal orders. What is securely his, and stands, is the
*framework* — the definition of orders and maximal orders, and the doubles
theorem (p. 69) that makes the search finite. The step from that framework to
"eight maximal orders" is a partial closure test mistaken for a complete one.

## Status

This sharpens the open point left in §6 of the preliminary exposition
(`paper/kirmse_1924_exposition.tex`): we can now say what Kirmse's eight most
likely were and where the computation fell short. Whether to fold this into
the exposition is deferred, with the rest, until Coxeter (1946) is read.

Facts vs. reconstruction are kept separate above: the counts 30 / 7 / 8, the
identity of the eight, and the location of $J_1$'s failure are **verified
facts**; the identification of Kirmse's unprinted computation with the
unit-invariance test is a **strongly-supported inference**.
