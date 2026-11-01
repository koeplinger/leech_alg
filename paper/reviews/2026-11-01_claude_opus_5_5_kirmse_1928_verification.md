# Primary-source verification — Kirmse's corrigendum (1928)

Created: 2026-11-01.  Updated: 2026-11-01.

**Reviewer:** Claude Opus 5.5 (Anthropic), at the direction of Jens Köplinger
**Prompts:** 365 (reading), 366 (the paper's footnote).
**Source read:** J. Kirmse, *Berichtigung aus der Abhandlung 1924: Zur
Darstellbarkeit natürlicher ganzer Zahlen als Summen von acht Quadraten*,
Ber. Verh. Sächs. Akad. Wiss. Leipzig, Math.-Phys. Kl. **80** (1928), 33–34
(JFM 54.0164.05).  Read from the author's screenshots of the HathiTrust copy
uc1.b3552809 (University of California, v. 79–80), kept outside the
repository.  Presented by Herglotz; dated "Schmölln (Thüringen), den
10. November 1927"; the session heading reads "Sitzung vom 16. Januar 1918",
evidently for 1928, the year of the volume.
**Question:** does the corrigendum touch the order claim of the 1924 paper
(§2, p. 70: eight maximal orders, among them $J_1$)?  See
`2026-05-20_claude_opus_4_7_kirmse_1924_verification.md` for that claim.

---

## What the corrigendum says

Kirmse describes the 1924 paper as his "Abhandlung über die Idealtheorie
der Cayleyschen Zahlen", in which "noch einige Schreibirrtümer zu
berichtigen und einige begriffstheoretische Ergänzungen anzubringen" are
(some slips of the pen are to be corrected and some conceptual additions
made).  There are two items and a summary.

1. **§1, the multiplication matrices.**  "In § 1 sind die beiden
   Produktzeilen S. 64 unten und S. 65 vertauscht, und so ist die
   rechtsseitige Matrix an die Stelle der linksseitigen gelangt.  Die
   weitere Entwicklung ist leicht zu verfolgen, die Schlußergebnisse
   bleiben."  (In §1 the two product lines at the bottom of p. 64 and on
   p. 65 are interchanged, so that the right-sided matrix has come into the
   place of the left-sided one.  The further development is easy to follow;
   the final results stand.)  The two lines are the lists
   "$a\cdot 1,\ a i_1, \ldots, a i_7$" (p. 64, which defines the left-sided
   matrix (3)) and "$1\cdot a,\ i_1 a, \ldots, i_7 a$" (p. 65, which defines
   the right-sided one).
2. **§3, the definition of ideals.**  "§ 3 ist mit einer neuen
   Idealdefinition auszustatten": ideal moduli are the modules
   $a_1g_1 + \cdots + a_\nu g_\nu$ (or $g_1a_1 + \cdots + g_\nu a_\nu$), with
   the $g_k$ in $J$; ideals are the ideal moduli all of whose numbers have
   relative norms divisible by the relative norm of the module.  The ideals
   of the paper are then exactly the ideals of the integral domain, but
   more of them are distinct: for $(\alpha, \alpha_1) = (\alpha, \alpha_2)$,
   $\bar\alpha_2\alpha_1 \equiv 0 \pmod{(\alpha)}$ is necessary, and
   $\alpha_2 \equiv \rho\,\alpha_1 \pmod{(\alpha)}$ with $\rho$ rational is
   necessary and sufficient (the coprimality condition on $\rho$ prints as
   "$(\rho, 1) = 1$" in the screenshot).
3. **Restated results (p. 34).**  The uniqueness of the primary
   representation, the unique prime-product representation of all numbers
   in $J_1$, the zeta functions of $J_1$ and $J_0$, and the counting
   theorems; in substance, every number of $J_1$ of norm $\alpha$ lies in
   $\omega(\alpha)$ primitive ideals of norm $\alpha$, each of which contains
   it together with $240 - 1$ other numbers of the same norm, not
   necessarily associates.

## Check of item 1

Table (1) of the 1924 paper (p. 64) lists the products of the basis units
explicitly; `python_project/src/verify_kirmse_1924.py` encodes it.  Against
table (1):

- Row $j$ of the printed matrix (3) (p. 65) is the expansion of $i_j a$, not
  of $a i_j$.  Row 1, for example, reads
  $-\alpha_1, \alpha_0, -\alpha_3, \alpha_2, \alpha_5, -\alpha_4, -\alpha_7,
  \alpha_6$, which is $i_1 a$ by table (1) ($i_1i_2 = i_3$,
  $i_1i_3 = -i_2$, $i_1i_4 = -i_5$, $i_1i_5 = i_4$, $i_1i_6 = i_7$,
  $i_1i_7 = -i_6$); $a i_1$ has the opposite signs in the last six places.
  So the printed matrix, called left-sided and defined through the $a i_j$,
  is the right-sided one, as the corrigendum says.
- The product formula (4), printed as $ab$ and obtained from (3) by the rule
  printed with it, gives $ba$ by table (1).  Its $i_1$ component, for
  example, contains $\alpha_3\beta_2 - \alpha_2\beta_3$, where that of $ab$
  contains $\alpha_2\beta_3 - \alpha_3\beta_2$.
- A scratch check at prompt 365 confirmed both for all rows and components,
  and that (4) is a composition product (log 365).

Whether a lattice is closed under a product does not depend on the order of
the factors: the set of all products $xy$ and the set of all products $yx$,
over all pairs from the lattice, coincide.  So the nonclosure of $J_1$ (10
of the 144 generator products leave it) and the count of seven closed
lattices among the thirty candidates hold under either product; the scratch
check recomputed both under the reversed product.  The pair of the
nonclosure witness (Coxeter 1946, §4) leaves $J_1$ in either order:
$a_1 a_3 = \tfrac12(i_1+i_2+i_4+i_7)$ and
$a_3 a_1 = \tfrac12(i_1+i_3+i_5-i_7)$ by table (1).

## Verdict

The corrigendum corrects §1 (an interchange of the left and right
multiplication matrices, with the final results standing) and §3 (the
definition of ideals, and with it the number of distinct ideals).  It does
not touch §2, where p. 70 states the eight maximal orders and defines
$J_1$, and it treats $J_1$ as the integral domain throughout.  The order
claim, and this project's account of it (the paper's Appendix B,
`verify_kirmse_1924.py`, the 1924 and forensic notes), is unaffected.  The
paper cites the corrigendum in a footnote of Appendix B (prompt 366).

## Status

Kirmse (1928): **verified.**
