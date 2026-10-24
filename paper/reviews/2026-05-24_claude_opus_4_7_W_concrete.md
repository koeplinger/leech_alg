# Concrete example: the subalgebra structure of W = σ(Ls)/2L

Created: 2026-05-24.  Updated: 2026-10-24.
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 155
**Script:** `python_project/src/verify_discovery1_W_subalgebra.py`

External feedback from **OpenCode DeepSeek, Hermes** proposed an explicit
basis for W = σ(Ls)/2L and a multiplication table, claiming it gives a
structural reason that replaces the appendix tables. This note records
the verified concrete example (corrected from the proposal in two cells
and in one structural claim). The "F₂-algebra" framing of the proposal
is not used here: the structure on L/2L is the descended octonion
product, which is non-commutative in the L-basis and has no
multiplicative identity, so it does not fit the conventional F₂-algebra
picture cleanly.

The companion finding about V = σ(Ls̄)/2L (totally isotropic, also a
right ideal, hence a Lagrangian) led to the Lagrangian framing that the
main paper's Appendix A, "The closure criterion modulo 2", records. This
note focuses on the W internal structure (idempotent + 3-dim nilpotent
ideal), which is not propagated into the paper.

---

## Setting

From the main paper's Appendix A, "The closure criterion modulo 2":

- The lattice $L = D_8^+$ does **not** contain $e_0$ (the octonion
  identity has odd coordinate sum); so $L$ is a nonunital order and, as a
  lattice, a $\sqrt2$-scaled copy of a maximal order of integral octonions
  (paper Section 2.2).
- $L \cdot L \subseteq L$ and Z-bilinearity descend the octonion product
  to a Z-bilinear map $\bar\mu \colon L/2L \times L/2L \to L/2L$.
- $L/2L \cong \mathbb{F}_2^8$ as a vector space; expressed in the
  L-basis $L_1 = s$, $L_2 = e_0 + e_1$, $L_3 = e_0 - e_1$,
  $L_{k+2} = e_0 + e_k$ ($k = 2, \ldots, 6$) of
  `python_project/src/verify_mod2_quotient.py`, the
  multiplication $\bar\mu$ is **non-commutative** and **non-unital**.
- $V := \sigma(L\bar s)/2L = 2\mathcal O/2L$ is a 4-dim left ideal of
  $L/2L$ (Lemma 4.3 ⟺ this), in fact a two-sided ideal, since
  $2\mathcal O$ is one in $L$ (paper Remark 4.5).
- $W := \sigma(Ls)/2L = p\mathcal O/2L$ is a 4-dim subalgebra of $L/2L$
  (Lemma 4.4 ⟺ this).
- $L/2L = V \oplus W$.

Here $\mathcal O$ is the maximal order of integral octonions of the paper's
Section 2.2 and $p = e_2 - e_4 + e_5 - e_6$; the two identifications are the
paper's display (A.3). The paper proves Lemma 4.4 in Section 4 from
$\sigma(Ls) = p\mathcal O$ (display (4.1)).

## A basis of W in L-coordinates

The following four vectors (in the L-basis, mod 2L) span W:

| | $L_1$ | $L_2$ | $L_3$ | $L_4$ | $L_5$ | $L_6$ | $L_7$ | $L_8$ |
|---|---|---|---|---|---|---|---|---|
| $w_0$ | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| $w_1$ | 0 | 1 | 0 | 1 | 1 | 0 | 1 | 0 |
| $w_2$ | 0 | 0 | 0 | 1 | 0 | 1 | 1 | 1 |
| $w_3$ | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 1 |

Note: $w_0 = L_1 + L_2 + L_3 = s + 2e_0$, since $L_2 + L_3 = 2e_0$. As
$e_0 \notin L$, $2e_0 \notin 2L$ (the paper's Section 5.2: $2L$ "contains
$4e_0$ but not $2e_0$"), so $w_0$ is not congruent to $s$ modulo $2L$:
$w_0$ lies in $\sigma(Ls)$, while $s$ does not. Since $s^2 = -s - 2e_0$,
$w_0 = -s^2$, so $w_0$ represents $s^2$ in $L/2L$, not the (non-existent)
identity of $L/2L$. The script's printed line "e_0 (identity of L/2L) in
W: False" tests the class $(1, 0, \ldots, 0)$ of $L_1 = s$: the class of
$s$ is not in W.

## Multiplication table

Verified directly from $\bar\mu$ on the proposed basis (with $\bar\mu$
inherited from the standard octonion product on $L$):

|   $\cdot$ | $w_0$ | $w_1$ | $w_2$ | $w_3$ |
|----------:|:-----:|:-----:|:-----:|:-----:|
| $w_0$     | $w_0$ |  $0$  | $w_2$ |  $0$  |
| $w_1$     | $w_1$ |  $0$  |  $0$  | $w_2$ |
| $w_2$     |  $0$  |  $0$  |  $0$  |  $0$  |
| $w_3$     | $w_3$ | $w_2$ |  $0$  |  $0$  |

Six non-zero entries. The table is closed in W
(every product is again in W), as required by W being a subalgebra of
$L/2L$ — i.e., the mod-2 image of Lemma 4.4.

## Structural features

1. **$w_0$ is idempotent in W**: $w_0 \cdot w_0 = w_0$. Since
   $w_0 = -s^2$, this records that $s^4 \equiv s^2 \pmod{2L}$; $s$ itself
   is not idempotent modulo $2L$, since $s^2 - s = -(e_0 + \cdots + e_7)$
   lies in $L$ but not in $2L$.
2. **$w_2$ is a universal zero divisor in W**: $w_2 \cdot w = 0$ for all
   $w \in W$, and $w \cdot w_2$ is $0$ or $w_2$ ($w_0 \cdot w_2 = w_2$,
   $w_1 \cdot w_2 = 0$, $w_3 \cdot w_2 = 0$; so $w \cdot w_2 = w_2$
   exactly when $w$ has a $w_0$ component). So $w_2$ is a **left**
   annihilator of W and a **right** annihilator of half of W; precisely,
   $w_2$ is the unique non-zero element in $W \cdot w_2 \cup w_2 \cdot W$.
3. **3-dim nilpotent ideal**: let $N := \mathrm{span}_{\mathbb{F}_2}\{w_1, w_2, w_3\}$.
   Verified: $N$ is closed under both left and right multiplication by
   W; $N^2 = \mathbb{F}_2 w_2$; $N^3 = 0$.
4. **Decomposition**: $W = \mathbb{F}_2 w_0 \oplus N$.
   The idempotent $w_0$ acts on $N$ as: $w_0 \cdot w_1 = 0$,
   $w_0 \cdot w_2 = w_2$, $w_0 \cdot w_3 = 0$;
   and $w_1 \cdot w_0 = w_1$, $w_2 \cdot w_0 = 0$, $w_3 \cdot w_0 = w_3$.
   The Peirce-like decomposition with respect to $w_0$ is asymmetric
   (since $\bar\mu$ is non-commutative).

## What this is *not*

This concrete description is not, by itself, a structural argument for
Lemma 4.4. The mod-2 reading converts "Lemma 4.4 holds" to "W is a
subalgebra of L/2L"; identifying W's internal multiplication table is a
refinement, but verifying closure of W still requires checking the 16
products $\bar\mu(w_i, w_j)$.

The structural argument is the paper's (Section 4): $\sigma(Ls) =
p\mathcal O$ (display (4.1)), with $p$ a purely imaginary element of the
maximal order $\mathcal O$, so $p\mathcal O$ is closed (display (4.3)).
Modulo $2L$ this identifies $W = p\mathcal O/2L$ (display (A.3)) without
reference to σ. Other intrinsic descriptions of W, for instance

- as the kernel (or image) of a natural F₂-linear or F₂-algebra map
  on L/2L,
- as the algebra generated by a canonical set of elements not
  referencing σ,
- as a Peirce component or radical-like piece of L/2L,

are not established. The 16-entry multiplication table above is the
concrete object.

## Cross-references

- Paper Appendix A, "The closure criterion modulo 2": the mod-2 reading of
  Lemmas 4.3 and 4.4, with $V = 2\mathcal O/2L$ and $W = p\mathcal O/2L$
  (display (A.3)); paper Section 4: Lemma 4.4 from
  $\sigma(Ls) = p\mathcal O$ (display (4.1)).
- `python_project/src/verify_mod2_quotient.py`: verifies V is a left
  ideal and W is a subalgebra of L/2L (Lemma 4.3 ⟺ V left ideal,
  Lemma 4.4 ⟺ W subalgebra).
- `python_project/src/verify_discovery1_W_subalgebra.py`: this note's
  verification (proposed basis, table, structural features).
- `prompt_logs/155_W_subalgebra_concrete.txt`: prompt log.

## Three corrections to the originally proposed table/claims

1. The actual multiplication table differs from the proposed one in
   two cells:
   - actual $w_0 \cdot w_0 = w_0$ (proposed: $w_0 + w_1$)
   - actual $w_1 \cdot w_0 = w_1$ (proposed: $w_1 + w_2 + w_3$)
2. The claim "$w_0$ contains the identity of L/2L" is wrong: L/2L
   has no left or right identity (since $e_0 \notin L$; checked over all
   256 classes), and $w_0 = -s^2$ represents $s^2$ modulo $2L$.
   $w_0$ is an idempotent in W, not an inherited identity.
3. The non-commutativity of the multiplication table is consistent
   with $L/2L$'s descended product being non-commutative in the
   L-basis — but the specific asymmetry of the proposed table was
   incorrect.

## Status

The concrete example is verified, with three corrections recorded. No
paper edit follows (per user direction). The script is in
`python_project/src/` and the verification reproduces this note's
table from the L-basis of `python_project/src/verify_mod2_quotient.py`
with exact-arithmetic mod-2L reduction.
