# Formal referee review of the main paper

**Date:** 2026-04-16
**Reviewer:** Claude Opus 4.7 (Anthropic), acting in the role of a highly trained mathematician at the direction of Jens Köplinger
**Manuscript reviewed:** `paper/main.tex` — *An order on the Leech lattice from transposition-twisted octonions* (Köplinger), version of 14 April 2026, together with the companion note `paper/companion.tex` and the supplemental repository as context.

This document records the referee report verbatim as delivered. Two demonstrable errors in the manuscript were verified computationally against `python_project/src/e8_wilson.py` and the triple-product formula of Definition 3.6 before the report was written; those checks are noted inline.

---

## Recommendation

Major revision. The central theorem is correct, but two claims in the exposition are demonstrably false, the status of the main result relative to existing constructions is not squarely addressed, and several stylistic choices weaken the paper's standing as a self-contained mathematical contribution.

---

## 1. Summary of the claim

The paper asserts that Wilson's Leech lattice $\Lambda\subset L^3$ is closed under a bilinear triple product $\star$ on $\mathbb{R}^{24}$ built from the $\sigma$-twisted octonion multiplication $x\cdot_\sigma y := \sigma(\sigma(x)\cdot\sigma(y))$ with $\mathbb{Z}_3$ cross-block routing. Closure is reduced to four lemmas, three structural and one computational, involving the transposition $\sigma$, the $E_8$ lattice $L$, and the sublattices $Ls$, $L\bar s$, $\sigma(Ls)$, $\sigma(L\bar s)$.

The underlying theorem is plausible and the reduction (Propositions 4.5, 4.7, 4.8) is logically sound *modulo* the four lemmas. The companion note and supplemental code support the computational lemmas. The concerns below are with (i) specific false statements, (ii) the relationship to existing work, and (iii) exposition.

---

## 2. Substantive errors

### 2.1 $\bar s \notin L$ (Section 2.3)

The paper states:

> "Define $s = \tfrac{1}{2}(-e_0+e_1+\cdots+e_7)$ and $\bar s = \tfrac{1}{2}(-e_0-e_1-\cdots-e_7)$. Both lie in $L$ with $N(s)=N(\bar s)=2$."

Under the Wilson/$D_8^+$ convention that the paper evidently intends (the convention stated in the companion, Section 2, and implemented in `e8_wilson.py`: half-integer vectors require *odd* coordinate sum), $\bar s$ has coordinate sum $-4$ and is **not** in $L$. Verified numerically against the supplemental code (`is_in_L(sbar)` returns `False`). The paper's own Section 2.2 description of $L$ — "either all integers or all half-integers, with even coordinate sum" — compounds the confusion: read literally it excludes $s$, read charitably it admits $\bar s$ but contradicts Wilson's convention.

This is not merely cosmetic. The argument uses $L\bar s = \{\ell\cdot\bar s : \ell\in L\}$ as a sublattice of $L$, which is true (verified: all 240 root products $r\cdot\bar s$ lie in $L$), but this has nothing to do with $\bar s$ itself being in $L$. The reader deserves either (a) a corrected definition of $L$ that includes both $s$ and $\bar s$ — which no standard $D_8^+$ convention does — or (b) the honest statement that $L\bar s\subseteq L$ while $\bar s\notin L$, with a one-line justification that $L\bar s$ is nonetheless a sublattice.

### 2.2 No multiplicative identity exists in $(\mathbb{R}^{24},\star)$ (Section 4)

Two assertions must be corrected:

> Property table: "Multiplicative identity — No — $(e_0,e_0,e_0)$ has $N=3$; not in $\Lambda$."
>
> Remark: "The ambient algebra $(\mathbb{R}^{24},+,\star)$ has multiplicative identity $(e_0,e_0,e_0)$, but this element has squared norm 3 and does not lie in $\Lambda$."

Both are false. Using Definition 3.6 with $u=(e_0,e_0,e_0)$:

$$P_1 = e_0\cdot_\sigma x' + e_0\cdot_\sigma y' + e_0\cdot_\sigma z' = x'+y'+z',$$

and similarly $P_2 = P_3 = x'+y'+z'$. So $(e_0,e_0,e_0)\star(x',y',z')=(x'+y'+z',\,x'+y'+z',\,x'+y'+z')$, which is not $v$ in general. Verified on a concrete example.

More importantly, the matrix of left-multiplication by $(a,b,c)$ on $(x',y',z')$ has the form

$$\begin{pmatrix} a & c & b \\ c & b & a \\ b & a & c \end{pmatrix}$$

(entries acting by left-multiplication in $\mathbb{O}$); for this to be the identity map one needs simultaneously $a=e_0, b=0, c=0$ and $b=e_0, a=0, c=0$ and $c=e_0$, which is inconsistent. Hence $(\mathbb{R}^{24},\star)$ admits no multiplicative identity; it is not merely that the identity fails to land in $\Lambda$. The "non-unital order" discussion must be rewritten. Given that the table purports to summarise computer-verified properties, this error also invites closer scrutiny of the methodology's claim to have "independently verified every foundation".

---

## 3. Mathematical substance: is the result genuinely new?

### 3.1 Remark 4.6 and the relationship to $\star_0$

Remark 4.6 correctly observes that $\Sigma := \sigma\oplus\sigma\oplus\sigma$ is an isomorphism $(\mathbb{R}^{24},\star)\to(\mathbb{R}^{24},\star_0)$, so that

$$\Lambda\star\Lambda\subseteq\Lambda \iff \Sigma(\Lambda)\star_0\Sigma(\Lambda)\subseteq\Sigma(\Lambda).$$

This is a candid framing, but it reduces the content of Theorem 1.2 to: *there exists a Leech embedding $\Sigma(\Lambda)\subset\mathbb{O}^3$ on which the standard triple product closes*. The paper then owes the reader an answer to the obvious questions:

- **(Novelty vs. Egan.)** The paper itself cites, in Section 6, that Egan "enumerated at least $17{,}280$ distinct off-diagonal Leech embeddings $L_L\subset\mathbb{O}^3$" on which a (Jordan-doubled) product closes. Is $\Sigma(\Lambda)$ in that list? Related to it? This is not a peripheral question — it is central to whether Theorem 1.2 is new.
- **(Novelty vs. Dixon.)** The cited Dixon (1995, 2010) papers construct $\Lambda$ via an octonionic XY-product. Dixon's XY-product on $\mathbb{O}^3$ (specialised appropriately) furnishes a bilinear product on $\mathbb{R}^{24}$; is the author's $\star$ essentially the same, a conjugate, or genuinely distinct? Footnote 4 of Section 6 compares Dixon–Wilson only historically, not mathematically.
- **(Novelty vs. Lepowsky–Meurman.)** Lepowsky–Meurman (1982) described $\Lambda$ as a sublattice of $E_8^3$ containing $2E_8^3$; do their generators fit into $\Sigma(\Lambda)$ or into a related orbit?

Without a direct comparison, a hostile reviewer could argue that the paper is a reformulation of known material with new bookkeeping. A sympathetic reviewer will still ask for at least one paragraph in Section 6 establishing *what is new*, preferably showing that $\Sigma(\Lambda)$ is (or is not) in the orbit of $\Lambda$ under $\mathrm{Aut}(\Lambda)=\mathrm{Co}_0$, and relating $\star$ to Dixon's XY-product on equal footing.

### 3.2 Status of the "symbolic proof"

Section 5 calls the argument a "symbolic proof", but Lemmas 4.3 and 4.4 are verified by computer (64 basis products each, on specific $\mathbb{Z}$-bases constructed in `symbolic_proof_checks.py`). A proof that invokes a finite computation is acceptable, but the branding is not: a genuine symbolic proof would provide a *structural* reason that $\sigma(Ls)$ closes while $Ls$ does not — for instance, identifying $\sigma(Ls)$ with a left ideal, with an image under a norm-preserving octonion automorphism, or with some module-theoretic structure over $L$. As it stands, the paper demonstrates the phenomenon without explaining it. That is acceptable, but should be labelled honestly ("computer-assisted proof" or "verification on a $\mathbb{Z}$-basis"), and the open question — *why* $\sigma(Ls)$ is closed under $\cdot$ — should be explicitly posed.

### 3.3 Coxeter's transposition vs. $\sigma$

The abstract and Section 6 frame the result as a chain "Kirmse $\to$ Coxeter $\to$ $\sigma$-twist, each step a transposition of two basis indices". This is elegant but conflates two structurally different operations: Coxeter's correction is a permutation of Fano *lines* (equivalently, a relabelling of the multiplication table that need not be induced by a linear map on $\mathbb{R}^8$); $\sigma$ is a linear coordinate involution of $\mathbb{R}^8$. The companion (Section 8, "Could the $\sigma$-twist undo Coxeter's fix?") addresses this clearly and correctly. The main paper must do the same, because a reader of only the main paper will walk away with a misleading picture of what kind of "transposition" is involved at each step.

### 3.4 Routing in Definition 3.6

The routing rule "$\alpha\cdot\beta\to\gamma$ with $\{\alpha,\beta,\gamma\}=\{1,2,3\}$" is under-specified: the three displayed formulas for $P_1,P_2,P_3$ make specific choices of *ordered* pairs ($z\cdot y',\,y\cdot z'$ both go to $P_1$), and the resulting product is not block-symmetric in any obvious sense. The construction is equivalent to a well-known ternary-style routing, but stating this without motivation or uniqueness discussion leaves the reader to wonder whether (i) another routing gives a different closed product, (ii) the closure depends sensitively on this routing, and (iii) the routing has any conceptual content beyond convention. Given Remark 3.5's "essentially unique" up to relabelling, the routing invites the same treatment: is it forced by the $\mathbb{Z}_3$ block symmetry together with $\sigma$-equivariance, or is it a free choice?

---

## 4. Minor mathematical issues

- **"Order" terminology.** Classically, an order is a subring of a finite-dimensional $\mathbb{Q}$-algebra that is a full $\mathbb{Z}$-lattice. Here the ambient algebra is over $\mathbb{R}$, not $\mathbb{Q}$, so calling $\Lambda$ an "order" requires the generalisation used in Remark 1.4 — which is given, but should appear earlier (before Corollary 1.3) and should be pinned to a reference, since "order" in non-associative real algebras is non-standard terminology.

- **Lemma 4.3.** The claim $L\cdot\sigma(L\bar s)\subseteq\sigma(L\bar s)$ under the *standard* (untwisted) product is verified by computation. This is not automatic from any $L$-module structure on $L\bar s$, because the non-associativity of $\mathbb{O}$ breaks the would-be derivation. A sentence making this explicit would help the reader trust the lemma.

- **$L$ as $D_8^+$ (Section 2.2).** The description is ambiguous between the "even sum for both branches" reading (which excludes $s$) and "odd sum for half-integer" (which excludes $\bar s$). Replace with the correct formulation used in the companion, noting explicitly that $\bar s\notin L$ but $L\bar s\subseteq L$ holds nonetheless.

- **Proof of Corollary 1.3.** One line asserts that $\mathrm{Min}(\Lambda)$ generates $\Lambda$ over $\mathbb{Z}$. This is true and standard, but the corollary ultimately requires bilinearity plus closure on *all* of $\Lambda\times\Lambda$, which the theorem already provides; the invocation of $\mathrm{Min}(\Lambda)$ is a red herring unless the preceding computational evidence is being folded in. Either excise the $\mathrm{Min}(\Lambda)$ sentence (the theorem alone suffices) or explain the role it plays.

- **Table of algebraic properties (Section 4).** The percentages are given without sample sizes or test distributions. "$\approx 47\%$" for norm-multiplicativity is a striking figure — is there a structural characterisation of the pairs on which it holds? Without any accompanying analysis, these numbers raise more questions than they answer.

- **Remark 3.5 (uniqueness).** "All transposition twists produce the same multiplication table up to basis relabelling." Correct ($\mathrm{GL}(3,\mathbb{F}_2)$ is 2-transitive on the 7 imaginary points), but a citation or one-line proof would serve the reader. Also, this is uniqueness *among transposition twists*; longer-cycle twists may yield different (possibly closed, possibly not) products.

- **Citations.** Baez (2014) is cited as a sequence of blog posts (the $n$-Category Café), and Petersson (2018) as a workshop lecture. Both are appropriate in context, but the manuscript leans on them for structural claims (the Baez–Egan enumeration; Petersson's modern account of the integral octonions) that probably admit stronger citations.

---

## 5. Exposition and presentation

- **Section 5 ("Research methodology").** The extensive discussion of the human–AI collaboration protocol, including remarks on logging-protocol violations by the AI, is unusual for a mathematics research paper. Substantive in its own right, but not load-bearing for the theorem. Consider moving to an appendix and trimming the autobiographical material.

- **Section 7 ("Outlook").** The speculation linking the construction to Standard Model physics (Furey–Hughes; Köplinger 2023) is prominent but thinly connected to Theorem 1.2. Either strengthen the connection to something concrete (e.g., identifying an order-3 automorphism of $(\mathbb{R}^{24},\star)$ in the spirit of Okubo) or move this material to a separate note.

- **Companion paper authorship.** Listing Claude (an AI system) as author of the companion note is a publication-politics question outside the referee's remit, but the journal will want to be aware of it. The companion is expository and may, depending on venue policy, be better presented as an appendix or supplemental material attributed to the principal author with AI assistance acknowledged.

- **Length.** The paper is currently roughly balanced between the theorem (Sections 1–4, compact and efficient) and meta-material (Sections 5, 7, acknowledgments). For a research note announcing a new closure result, target a tighter manuscript: theorem, proof, properties, a focused related-work section explicitly placing the result against Dixon/Wilson/Egan/Lepowsky–Meurman, and a short conclusion. The methodology and outlook can live in a companion document or postprint.

---

## 6. Summary of requested revisions

**Mandatory (errors):**
1. Correct the claim "$\bar s\in L$" (Section 2.3). State the $L$ membership criterion so it is consistent with the supplemental code's `is_in_L`.
2. Correct the multiplicative-identity statements (Section 4). The ambient algebra $(\mathbb{R}^{24},\star)$ has no multiplicative identity; $(e_0,e_0,e_0)$ is not one.
3. Clarify or rename the "symbolic proof" so that the computational nature of Lemmas 4.3 and 4.4 is honest on the tin.

**Strongly recommended:**
4. Add a paragraph in Section 6 comparing $\Sigma(\Lambda)$ and $\star$ against the Egan enumeration and against Dixon's XY-product, explicitly addressing novelty.
5. Pull the Kirmse–Coxeter vs. $\sigma$ distinction from the companion into Section 6 of the main paper.
6. Move Section 5 to an appendix and tighten Section 7.
7. State or cite the argument that $L\bar s\subseteq L$ is a *sublattice* (closure under addition and $L$-multiplication) so that the reader sees the consistency despite $\bar s\notin L$.

**Nice to have:**
8. A structural (non-computational) reason for Lemma 4.4 — if one exists, it would substantially strengthen the paper.
9. Sample-size and distribution information for the algebraic-property statistics in Section 4.
10. An explicit statement and justification of the routing choice in Definition 3.6.

The core mathematical content — a concrete bilinear product on $\mathbb{R}^{24}$ under which $\Lambda$ is closed — is genuine and, modulo the novelty question against Dixon/Egan, worth publishing. But the manuscript as submitted contains two verifiable false statements in its exposition and leaves unresolved the question of whether the construction is a new closure phenomenon or a presentation of a known one. Major revision and a second referee round recommended.

---

## Disposition by the author (2026-04-16)

Upon receipt of this review, the author directed the following immediate corrections in the main paper (the companion is checked for parallel issues):

- **2.1 (mandatory):** adopt option (b) — state explicitly that $L\bar s\subseteq L$ while $\bar s\notin L$, with a brief justification. Wilson's symbol conventions are preserved.
- **2.2 (mandatory):** state simply that $(\mathbb{R}^{24},\star)$ has no multiplicative identity; remove the incorrect $(e_0,e_0,e_0)$ justification.

The remaining points (3.1–3.4, §§4–5, §6 items 4–10) are deferred to a later revision pass. The author notes specifically:

- The "symbolic proof" re-labelling (3.2) is accepted and will be addressed when a genuinely symbolic proof is supplied.
- The Coxeter-vs-$\sigma$ distinction (3.3) is accepted as important and will be pulled into the main paper; it is not made now because the paper is in early circulation and further edits risk adding confusion.
- The novelty comparison (3.1) will be revisited closer to formal submission; the author's present view is that the simplicity of the $\mathbb{Z}_3$-symmetric triple-octonion product is itself the notable feature, and that making the comparison the centre of the paper invites a "simplicity debate" that is not productive absent a maximality result.

All deferred items remain open for formal review and are not to be considered settled.
