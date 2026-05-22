# Item (3) — analysis of Corradetti's Okubo / E8 / 2-adic work

**Date:** 2026-05-22
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 114 (Phase A of the v4 update plan).
**Primary source (read in full):** D. Corradetti, *Integral elements of Okubo
algebra and the E8-lattice*, arXiv:2605.09333, 10 May 2026 (17 pp.).
**Also surveyed (abstracts only):** arXiv:2605.09458, 2605.15075, 2605.18538.

---

## 1. What arXiv:2605.09333 does

Corradetti studies integral elements for the three eight-dimensional
*division* composition algebras — octonions $\mathbb{O}$, para-octonions
$p\mathbb{O}$, and the real Okubo algebra $\mathcal{O}$ — built from the one
maximal octonion order:

- **Foundation:** the **Coxeter–Dickson order** $\mathbb{O}_{E_8}$ — the
  maximal order of integral octonions, 240 units, underlying lattice $E_8$,
  shell formula $|\mathbb{O}_{E_8}(n)| = 240\sum_{d\mid n}d^3$. He cites
  Dickson 1923 and Coxeter 1946 — exactly the sources this project verified.
- **Para-octonions** ($x\bullet y = \bar x\cdot\bar y$): $\mathbb{O}_{E_8}$
  *is* closed; **Theorem 6** — a genuine $\mathbb{Z}$-order, underlying
  lattice still $E_8$. A clean positive result, parallel to the octonionic
  case.
- **Okubo product** ($x\ast y=\tau(\bar x)\cdot\tau^2(\bar y)$, $\tau$ the
  order-3 octonion automorphism, which carries a $\sqrt3$): **Theorem 7** —
  $\mathbb{O}_{E_8}$ is **not** closed under $\ast$ over $\mathbb{Z}$ (an
  explicit basis product has coefficient $-\tfrac32\sqrt3$). The natural ring
  is $\mathbb{Z}[\sqrt3]$. **Theorem 9** — a diagonal **2-adic scaling**
  $D=\operatorname{diag}(2,2,2,2,4,4,4,4)$ yields a closed
  $\mathbb{Z}[\sqrt3]$-order $\mathcal{O}_0$.
- **The conductor sublattice (Theorem 12):** the $\mathbb{Z}$-"metric shadow"
  $L_{\mathrm{Ok}} = D\cdot E_8$ is an even sublattice of $E_8$ of index
  $2^{12}$, determinant $2^{24}$, minimum 8 (no roots). $E_8$ itself is
  **not** the Okubo integral lattice; it is recovered (Theorem 14) as the
  **2-adic saturation** of $L_{\mathrm{Ok}}$ — gluing along a maximal
  isotropic $(\mathbb{Z}/2)^4\oplus(\mathbb{Z}/4)^4$. Corradetti is candid:
  this "does not make the title literal" — the recovery is metric-arithmetic,
  not multiplicative.

**The three newer papers** (abstracts only): *Integral Shell Polytopes*
(2605.09458) is a geometric companion — same Okubo / "two-adic hierarchy" /
Gosset-polytope gluing theme. *Non-crystallographic systems* (2605.15075) and
*Integral Planes and Unit-Norm Polytopes* (2605.18538) are adjacent
composition-algebra / golden-ring / polytope work, not on orders-on-$E_8$
with 2-adic reduction. **Only 2605.09333 (with 2605.09458 as companion) is
on point** for us.

## 2. Genuine parallels with our work

1. **Shared foundation.** Corradetti's whole construction sits on the
   Coxeter–Dickson integral octonions — the very order whose history this
   project verified (Kirmse 1924, Dickson 1923, Mahler 1942, Coxeter 1946
   pending). His Remark-4 shell formula $240\sum d^3$ is the $E_8$/Eisenstein
   identity Kirmse derived.
2. **The Okubo / symmetric-composition / ternary theme.** Our §8 Outlook
   explicitly invokes Okubo algebras, Elduque, symmetric composition
   algebras, and a "ternary reformulation" built on an order-3 automorphism.
   Corradetti's paper *is* the Okubo-arithmetic programme — built on exactly
   that order-3 automorphism $\tau$. It is contemporary independent work in
   the direction our Outlook gestures toward.
3. **2-primary structure.** Both works find octonion-derived arithmetic on
   $E_8$-type lattices governed by powers of 2: his $[E_8:L_{\mathrm{Ok}}]
   = 2^{12}$; our $[L^3:\Lambda]=2^{12}$ (verified in item 1c); our item-(2)
   reduction lives in $L/2L$.
4. **Order via a modified octonion product.** His para-octonion order on
   $E_8$ joins Dixon's XY-product and Baez–Egan as "an order from a modified
   octonion product" — the company our §6 Related Work already keeps.

## 3. Contrasts — to state precisely, not blur

- **Dimension/object.** Corradetti: 8-dimensional, the Okubo/para-octonion
  products *on* the octonion order. Ours: 24-dimensional, a $\mathbb{Z}_3$-routed
  *triple*-octonion product on the Leech lattice. **Corradetti does not touch
  the Leech lattice** — there is no overlap of results, no priority concern.
- **The two "2-adic"s are different operations.** His is a 2-adic *scaling*
  (the diagonal $D$), a conductor *sub*lattice, and 2-adic *saturation*
  (sublattice → overlattice $E_8$). Our item (2) is reduction *mod 2* — a
  *quotient* $L\to L/2L$ expressing the closure lemmas as $\mathbb{F}_2$-subspace
  conditions. Both are "2-primary in flavour," but they are not the same
  construction; **a technical identification of the two would be an
  overclaim** and must be avoided.
- **Closure over $\mathbb{Z}$.** Our $(\Lambda,\star)$ is a genuine
  $\mathbb{Z}$-order on $\Lambda$ itself. Corradetti's Okubo product does
  *not* close over $\mathbb{Z}$ on $E_8$ — it needs $\mathbb{Z}[\sqrt3]$, and
  even then the $\mathbb{Z}$-shadow is a proper sublattice. (His
  *para*-octonion order, by contrast, *is* a clean $\mathbb{Z}$-order on
  $E_8$.) So in the "closes over $\mathbb{Z}$ on the named lattice" sense our
  result is the unconditional one; his Okubo case is explicitly conditional —
  a fair point of contrast, provided his clean para-octonionic result is
  stated alongside it.
- **The twist.** His Okubo obstruction is a genuine $\sqrt3$ field extension
  (from $\tau$); our $\sigma$ is a $\mathbb{Z}$-linear coordinate transposition
  — no field extension.

## 4. Recommendation (decision DP5)

**Cite Corradetti (2605.09333). Yes.** It is recent, independent, on-topic
work; citing a preprint is consistent with the paper's existing practice
(Dixon, Furey–Hughes).

Where and how:
- **§8 Outlook** — the natural home: alongside the existing Okubo/Elduque/
  ternary sentences, as contemporary work realising the Okubo-arithmetic
  direction, *and* as an instructive contrast (its Okubo product requires
  $\mathbb{Z}[\sqrt3]$ and yields only a conductor sublattice of $E_8$,
  whereas our $\sigma$-twisted triple product closes on $\Lambda$ over
  $\mathbb{Z}$).
- **§6 Related Work** — optionally, one sentence: Corradetti's para-octonionic
  order on $E_8$ as a further "order from a modified octonion product."

What **not** to do: do not assert that our mod-2 reduction (item 2) and his
2-adic saturation are the same mechanism — they are thematically adjacent,
operationally distinct. The honest comparison is "both find 2-primary
structure governing octonion-derived orders," no more.

Bearing on item (2): the comparison sharpens our own exposition — when we
present the 2-adic nature of Tables A.2/A.3 we should name it precisely as a
*reduction mod 2 / $\mathbb{F}_2$-quotient*, so it is not conflated with a
2-adic scaling-and-saturation in the reader's mind.

## 5. Status

The analysis of 2605.09333 is from a full read; the three newer papers from
abstracts only (sufficient to confirm they are not more on-point than the
main one). Recommendation: cite Corradetti in §8 (and optionally §6) — final
wording to be settled when those sections are revised.
