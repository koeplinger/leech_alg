# Primary-source verification — Kirmse (1924)

**Date:** 2026-05-20
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 100. Step (1) of the Kirmse-twist plan: derive the history from
the primary sources, from scratch.
**Source read:** J. Kirmse, *Über die Darstellbarkeit natürlicher ganzer
Zahlen als Summen von acht Quadraten und über ein mit diesem Problem
zusammenhängendes nichtkommutatives und nichtassoziatives Zahlensystem*,
Ber. Verh. Sächs. Akad. Wiss. Leipzig, Math.-Phys. Kl. **76** (1924), 63–82.
Read in full from the scans in `LeechAlg_1920s-1940s/Kirmse_1924/`
(`seq_325`–`seq_344`), clear Antiqua type.
**Reproducible check:** `python_project/src/verify_kirmse_1924.py` — octonion
multiplication and $J_1$ built solely from Kirmse's paper; no external input.

---

## How far Kirmse got — what the paper actually contains

Kirmse got remarkably far. The 1924 paper is a near-complete arithmetic
theory of the integral octonions:

- **§1 (pp. 63–67).** Octonion arithmetic from an explicit multiplication
  table (his table (1), p. 64); norm $N(a)=a\bar a$, conjugate; norm-
  multiplicativity $N(ab)=N(a)N(b)$; the multiplication matrices and the
  "Grundnorm" $=N(a)^4$; the **Kirmse identities** $\bar a(ab)=N(a)b=(ab)\bar b$
  type relations (his (8)); coplanar numbers, their commutativity/
  associativity, and powers.
- **§2 (pp. 68–71).** Defines *module*, *finite module*, *Integritätsbereich*
  (an order: a finite module closed under multiplication) and *umfassendster
  Integritätsbereich* (a maximal order). Proves that doubles of order-elements
  have integer components. Asserts he has found **all** maximal orders
  over $J_0=[1,i_1,\dots,i_7]$, names one $J_1$, with $[J_1:J_0]=16$.
- **§3 (pp. 71–75).** Two-sided ideal theory; ideal norms; primitive ideals.
- **§4 (pp. 76–77).** Units: $J_0$ has 16, **$J_1$ has 240** — "a group in an
  extended sense" lacking only associativity (i.e. a Moufang loop); the
  $2160=240\cdot 9$ elements of norm 2.
- **§5 (pp. 77–82).** Prime decomposition; the ideal zeta function
  $Z(s)=\zeta(4s)\zeta(4s-3)$; all ideals of $J_1$ principal; and the
  **eight-squares theorem** $r_8(n)=16\,\sigma_3(n)$ for odd $n$. A postscript
  adds a Minkowski geometry-of-numbers proof that all ideals are principal.

The lattice/quadratic-form results (240 minimal vectors, the theta-series
counts, the eight-squares theorem) are correct — they concern the $E_8$
lattice as a quadratic form and do not depend on multiplicative closure.

## The error — and it is exactly the one Köplinger recalled

On **p. 70** (`seq_332`), verbatim:

> *"Dieser Satz erlaubt sofort, durch eine einfache endliche Rechnung die
> sämtlichen umfassendsten Integritätsbereiche über $J_0$ aufzusuchen. **Man
> findet acht solche**, deren jeder die Gestalt hat …"*

Kirmse explicitly claims there are **eight** maximal orders. The correct
number is **seven**. This is "Kirmse's mistake," in his own words — and it
matches the recollection that prompted this investigation, and the Wikipedia
statement that "he thought there were eight maximal orders rather than seven."

## Verification from scratch

Using *only* Kirmse's own data — his table (1) and his definition
$J_1=[J_0,\;a_1,a_2,a_3,a_4]$ with
$a_1=\tfrac12(1{+}i_1{+}i_2{+}i_3)$,
$a_2=\tfrac12(i_4{+}i_5{+}i_6{+}i_7)$,
$a_3=\tfrac12(1{+}i_1{+}i_4{+}i_5)$,
$a_4=\tfrac12(1{+}i_2{+}i_4{+}i_7)$ (p. 70) — the script confirms:

1. Table (1) is a genuine octonion (composition) algebra: $N(uv)=N(u)N(v)$
   holds (300 random pairs, exact arithmetic).
2. **Kirmse's $J_1$ is *not* closed under his own multiplication.** Of the
   144 generator-pair products, 10 leave $J_1$. Explicit witness:
   $$a_1\cdot a_3=\tfrac12(1{+}i_1{+}i_2{+}i_3)\cdot\tfrac12(1{+}i_1{+}i_4{+}i_5)
   =\tfrac12(i_1{+}i_2{+}i_4{+}i_7)\notin J_1.$$
3. Enumerating the 30 $E_8$-type candidate orders (doubly-even self-dual
   $[8,4]$ codes), **exactly 7** are closed under table (1) — seven maximal
   orders, not eight.
4. Kirmse's named $J_1$ is among the 30 candidates but is **one of the
   non-closed ones**. So Kirmse not only miscounted (8 vs 7); the specific
   order he named and carried through §§3–5 is itself not multiplicatively
   closed.

That his §§4–5 conclusions are nonetheless correct is consistent with point
above: those are quadratic-form facts about the $E_8$ lattice; a genuine
maximal order yields the same counts. The flaw is confined to the algebraic
claim "$J_1$ is an order."

## Notes for the digest / wording update

- The manuscript's "Kirmse proposed a labelling under which $D_8^+$ is not
  closed" is **correct in substance** and now verified: Kirmse's $J_1$ is
  exactly such a non-closed lattice.
- The phrase to add, supported by the source: Kirmse's specific error was
  claiming **eight** maximal orders where there are **seven** — his $J_1$
  being a non-closed candidate. (So the answer to "did Kirmse propose seven
  correct orders and err on an eighth?" is: he claimed eight; seven are
  genuine; the verification shows his *named* one was a non-closed one.
  Whether the other seven of "his eight" are exactly the seven genuine
  maximal orders is not determinable from his text — he gives a schematic
  form, not an explicit list — and is left for Coxeter (1946).)
- Biographical: the article's by-line reads "Johannes Kirmse **in
  Hannover**", and it is dated "Hannover, am 15. April 1924" — not Apolda.
  `terminology.md` (following Petersson) says Apolda; minor, worth a check.
- The correction itself (who fixed it, and the precise form of the fix)
  still requires **Coxeter (1946)**; this verification settles Kirmse's
  side only.

## Status

Kirmse (1924): **verified.** Pending: Dickson (1923) — downloaded; Mahler
(1942) — loaded; Coxeter (1946) — to be obtained.
