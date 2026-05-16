# Review note — appendix table reduction, and a footnote error (round 9)

**Date:** 2026-05-16 (second of the day)
**Reviewer:** Claude Opus 4.7 (Anthropic), acting as a mathematician at the direction of Jens Köplinger
**Manuscript:** `paper/main.tex` — *An order on the Leech lattice from a $\mathbb{Z}_3$-symmetric triple-octonion product*, v3.
**Prompt:** 096 — (a) can the two integer tables for Lemmas 4.3 and 4.4 be reduced to one; (b) findings from the closer look this required.

This note is investigative, not a full pass. It was triggered by the question of whether Lemmas 4.3/4.4 can rest on fewer explicit tables. Looking closely enough to answer that turned up one genuine error, reported first.

All claims below were checked in an independent implementation (`/tmp/verify_paper.py` plus the closure sweep): octonions from the Fano triples, exact `fractions.Fraction` arithmetic.

---

## A. Substantive item — the §2.3 footnote is wrong (and self-contradictory)

**Location.** The footnote attached to "Both $L\bar s$ and $Ls$ are sublattices of $L$" in §2.3 ([main.tex:260–268](../main.tex#L260-L268)).

**What it says.** It states two properties of $L\bar s$ — (i) it is a $\mathbb{Z}$-subgroup of $L$, and (ii) it "is moreover closed under left-multiplication by $L$" — and then concludes: *"For $Ls$ both properties are immediate from $s \in L$."*

**The error.** Property (ii) is **false for $Ls$**. I tested $L \cdot Ls \subseteq Ls$ at the basis level: **40 of the 64 products $L_i \cdot (L_j s)$ leave $Ls$**. So $Ls$ is *not* closed under left-multiplication by $L$.

This is not a borderline call — it contradicts the paper itself:

- Remark 4.5 (`rem:Ls-not-closed`) states "$Ls$ is *not* closed under the standard octonion product." Since $Ls \subseteq L$, a left ideal $L \cdot Ls \subseteq Ls$ would force $Ls \cdot Ls \subseteq L \cdot Ls \subseteq Ls$ — i.e. exactly the closure Remark 4.5 (and Table A.3's whole reason for existing) denies. The footnote and Remark 4.5 cannot both be true.
- The non-closure of $Ls$ is the **central obstruction the paper exists to resolve**. The footnote inadvertently asserts that obstruction away.

**Why it is not load-bearing.** The footnote only needs to justify the sentence it hangs on — "$Ls$ is a sublattice of $L$" — for which the *only* facts required are $Ls \subseteq L$ and that $Ls$ is a subgroup. Both of those *are* immediate from $s \in L$ ($\ell s \in L\cdot L \subseteq L$). Property (ii) was never needed for $Ls$; it was pulled in because it is genuine and relevant for $L\bar s$ (it is Wilson's condition (2)), and then over-extended to $Ls$ by the phrase "both properties."

**Suggested fix (forward).** Split the two cases instead of saying "both properties":

> *... For $Ls$, the subgroup property is likewise immediate from $s \in L$. Closure under left-multiplication by $L$, however, holds for $L\bar s$ but **fails** for $Ls$ — this is precisely the obstruction resolved by the twist (Remark \ref{rem:Ls-not-closed}).*

This costs one sentence, removes the contradiction, and actually foreshadows the proof's crux.

**Aside (correct, optional).** While checking this I found that $L\bar s$ is in fact a **two-sided** ideal of $L$: $L \cdot L\bar s$, $L\bar s \cdot L$, and $L\bar s \cdot L\bar s$ are all $\subseteq L\bar s$ (0 failures each). The paper only claims, and only needs, the left-multiplication half. No change required; noted in case it is useful.

---

## B. The question asked — can Lemmas 4.3 and 4.4 use fewer tables?

**Yes.** There is a clean reduction. It does better than "one table instead of two": it makes **Table A.1 the single integer table that carries all three lemmas**, with Lemmas 4.3 and 4.4 becoming finite $\mathbb{F}_2$-linear-algebra checks. The reduction is *symbolic* in the sense the friend suggested — it replaces unbounded integer coefficient tables by arithmetic in an 8-dimensional algebra over $\mathbb{F}_2$.

### B.1 The idea: reduce modulo $2L$

The paper already states (§2.3, citing Wilson) that $L\bar s \cap Ls = 2L$. Hence $2L \subseteq Ls$ and $2L \subseteq L\bar s$. Applying the involution $\sigma$ and using $\sigma(2L) = 2L$:
$$2L \subseteq \sigma(Ls), \qquad 2L \subseteq \sigma(L\bar s).$$
(Verified directly.)

Because $L \cdot L \subseteq L$ (Lemma 4.2) and the octonion product is $\mathbb{Z}$-bilinear, $(a + 2c)\cdot b \equiv a\cdot b \pmod{2L}$ for $c \in L$, and likewise in the right factor. So the product **descends** to a well-defined $\mathbb{F}_2$-bilinear product
$$\bar\mu : \bar L \times \bar L \to \bar L, \qquad \bar L := L/2L \cong \mathbb{F}_2^8.$$
The structure constants of $\bar\mu$ are **exactly the entries of Table A.1 reduced mod 2**. No new table.

Since $2L \subseteq \sigma(L\bar s)$ and $2L \subseteq \sigma(Ls)$, those two sublattices are *full preimages* of $\mathbb{F}_2$-subspaces
$$V := \sigma(L\bar s)/2L, \qquad W := \sigma(Ls)/2L.$$
From $\sigma(L\bar s) + \sigma(Ls) = L$ and $\sigma(L\bar s)\cap\sigma(Ls) = 2L$ one gets $\bar L = V \oplus W$ with $\dim V = \dim W = 4$ (index $16 = 2^4$ on each side). Verified: $\dim V = \dim W = 4$, $V \cap W = 0$.

### B.2 The two lemmas become subspace conditions

- **Lemma 4.3** ($L \cdot \sigma(L\bar s) \subseteq \sigma(L\bar s)$) $\iff$ $\bar\mu(\bar L, V) \subseteq V$ — i.e. **$V$ is a left ideal of $(\bar L, \bar\mu)$**.
- **Lemma 4.4** ($\sigma(Ls)\cdot\sigma(Ls) \subseteq \sigma(Ls)$) $\iff$ $\bar\mu(W, W) \subseteq W$ — i.e. **$W$ is a subalgebra of $(\bar L, \bar\mu)$**.

The equivalences are clean both ways. For the substantive ("if") direction of 4.3: given $a \in L$, $w \in \sigma(L\bar s)$, the product $a\cdot w$ lies in $L$ (Lemma 4.2), and its class is $\bar\mu([a],[w]) \in \bar\mu(\bar L, V) \subseteq V$; since $a\cdot w \in L$ projects into $V$ and $\sigma(L\bar s)$ is the *full* preimage of $V$, we get $a\cdot w \in \sigma(L\bar s)$. Lemma 4.4 is identical with $W$. (The mod-2 reductions are exact: I verified $\bar\mu(\bar L,V)\subseteq V$ and $\bar\mu(W,W)\subseteq W$ both hold — PASS — matching the integer lemmas.)

### B.3 What this does to the appendix

| | current | after reduction |
|---|---|---|
| Lemma 4.2 ($L\cdot L\subseteq L$) | Table A.1 (64 integer rows; textbook, Coxeter) | Table A.1 — **kept** as the method demonstration |
| Lemma 4.3 | Table A.2 (64 integer rows) | $\bar\mu(\bar L,V)\subseteq V$: $8\times4 = 32$ products over $\mathbb{F}_2$ |
| Lemma 4.4 | Table A.3 (64 integer rows) | $\bar\mu(W,W)\subseteq W$: $4\times4 = 16$ products over $\mathbb{F}_2$ |

Both integer tables A.2 and A.3 (128 rows over $\mathbb{Z}$, entries unbounded) are replaced by checks against **one object** — the mod-2 reduction of Table A.1 — plus the explicit mod-2 generators of $V$ and $W$ (the $N_k$ and $M_k$ reduced mod 2). The $\mathbb{F}_2$ checks are small enough to display compactly or even inline.

### B.4 Bonus: it answers an open question

The Conclusion lists as open: *"A structural reason for Lemma \ref{lem:sigmaLs-closed}"* — why is $\sigma(Ls)$ closed while $Ls$ is not? The mod-2 reduction gives a partial answer in exactly the requested register: **closure of $\sigma(Ls)$ is the statement that the 4-dimensional subspace $W$ is a subalgebra of the 8-dimensional $\mathbb{F}_2$-algebra $L/2L$**, and the twist's role is to move the condition-3 subspace from $Ls/2L$ (not a subalgebra) to $W$ (a subalgebra). This is more structural than 64 integer rows, though it is still a verification, not yet a representation-theoretic "why".

### B.5 Honest caveats

- This is **a reduction in size and a gain in transparency, not a removal of all computation**. Lemmas 4.3 and 4.4 each still require a finite check; they are genuinely independent of Lemma 4.2 (e.g. $W$ being a subalgebra does not follow formally from $V$ being an ideal). What changes is that the checks are now small $\mathbb{F}_2$ computations against a single shared table, not two separate $\mathbb{Z}$-valued tables.
- The reduction leans on $L\bar s \cap Ls = 2L$, which the paper currently *cites* to Wilson. If the appendix is to be self-contained, that one equality should be stated as the hypothesis the reduction consumes (it is cheap to verify: both lattices have index 16 and $2L \subseteq$ each, so $L\bar s\cap Ls$ has index dividing $16$ and contains $2L$ of index $2^8/2^8\cdot$… — simplest to just verify $2L \subseteq L\bar s\cap Ls$ and the index, which the code already does).
- Whether to *adopt* this is an exposition choice. Keeping Table A.1 as the worked demonstration and recasting 4.3/4.4 in the mod-2 language would shorten the appendix and strengthen §7's structural narrative. Keeping all three integer tables is also defensible (uniform, fully explicit). The mathematics of the paper is unchanged either way.

---

## C. Recommendation

1. **Fix the §2.3 footnote** (item A) — this is a real error and should be corrected before the paper progresses, independently of anything else. The fix is one sentence.
2. **Consider the mod-2 reduction** (item B) — optional but recommended: it cuts the appendix from three integer tables to one, makes Lemmas 4.3/4.4 transparent $\mathbb{F}_2$ statements, and partially answers an open question the paper itself poses.

Item A is the only required change. The central theorem and its proof remain correct and unaffected.
