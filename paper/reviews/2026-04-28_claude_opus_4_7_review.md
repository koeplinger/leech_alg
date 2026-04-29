# Formal referee review of the main paper (round 2)

**Date:** 2026-04-28
**Reviewer:** Claude Opus 4.7 (Anthropic), acting in the role of a highly trained mathematician at the direction of Jens Köplinger
**Manuscript reviewed:** `paper/main.tex` — *An order on the Leech lattice from a $\mathbb{Z}_3$-symmetric triple-octonion product* (Köplinger), version of 28 April 2026, together with the supplemental repository as context.

This is a second-round review. The 16 April review flagged two demonstrable errors and a number of recommended changes; this round assesses the manuscript as updated. One residual error survives (the description of $L$ inside the proof of Lemma 4.1), and several items remain partially addressed. The recommendation is minor revision: the mathematical core is sound, the symbolic proof is now genuinely symbolic, and the remaining concerns are smaller in scope than those of the previous round.

---

## Items resolved since the 16 April review

For continuity:

- **2.1** ($\bar s \notin L$): both §2.2 and §2.3 are now consistent with the supplemental code's `is_in_L`.
- **2.2** (multiplicative identity): the block-matrix argument is correct; the table entry and the remark agree.
- **3.2** ("symbolic proof" branding): the proofs of Lemmas 4.2–4.4 in §4 close by reference to explicit basis tables in Appendix A. All 64 integer coefficients per lemma are displayed.
- **3.3** (Coxeter vs $\sigma$): pulled into Remark 3.3 in §3 and into the §6 Kirmse–Coxeter–$\sigma$ paragraph; the structural separation between Coxeter's Fano-triple permutation and $\sigma$'s linear coordinate involution is explicit.
- **§5 → Appendix B**; **§7 (Outlook)** tightened to two short subsections.
- **Sublattice status of $L\bar s$**: explained in a new footnote in §2.3.
- **Comparison paragraph** in §6 distinguishes the present construction from Baez–Egan and Dixon.

The revision is substantial and the headline claims are now better supported.

---

## Recommendation

Minor revision. One residual error (§3.1 below). Six "strongly recommended" items, mostly of clarification or one-sentence justification. No structural rewrite needed.

---

## 1. Summary of the claim

The paper asserts that Wilson's Leech lattice $\Lambda \subset L^3$ is closed under a bilinear $\mathbb{Z}_3$-symmetric triple product $\star$ on $\mathbb{R}^{24} = \mathbb{O}_1 \oplus \mathbb{O}_2 \oplus \mathbb{O}_3$, where each cross-block product uses the same $\sigma$-twisted octonion multiplication $x \cdot_\sigma y := \sigma(\sigma(x)\cdot\sigma(y))$ and the routing is the cyclic $\mathbb{Z}_3$ one of Definition 3.6. Closure reduces to four lemmas, the last two proved by exhibiting the 64 integer coefficients of the basis-by-basis products in Appendix A.

The claim is mathematically clean and the proof, modulo the reservations below, is correct.

---

## 2. The residual error

### 2.1 Lemma 4.1's proof restates the wrong description of $L$ (line 423–424)

The proof of Lemma 4.1 ($\sigma(L) = L$) reads:

> "The lattice $L = D_8^+$ consists of all vectors in $\mathbb{R}^8$ whose coordinates are either all integers or all half-integers, with even coordinate sum."

This is the description that was corrected in §2.2 in the present round (the half-integer branch must have *odd* coordinate sum, not even). Section 2.2 now says it correctly; §4.1's proof still says it wrong.

The conclusion of the proof — that $\sigma(L) \subseteq L$ — is unaffected, because a coordinate transposition preserves both branches' parities and the sum, regardless of which version of the criterion is invoked. So the lemma remains correct. But the paper's own statement of what $L$ *is* now disagrees with itself between §2.2 and §4.1. This is the same class of inconsistency that the 16 April review flagged in §2.3 and that the author intended to fix in whole (Prompt 077, item 1; Prompt 078). The fix needs to be propagated.

Verified by direct comparison of [main.tex:180–185](../main.tex#L180-L185) vs. [main.tex:423–424](../main.tex#L423-L424).

---

## 3. Substantive items

### 3.1 The bilinearity-reduction preamble does not justify "only if" (§4.2, lines 431–445)

The preamble at the head of §4.2 reads:

> "...an inclusion $A \cdot B \subseteq C$ between such lattices holds if and only if it holds on a chosen pair of $\mathbb{Z}$-bases..."

The "if" direction is what the proof uses, and that direction is established by $\mathbb{Z}$-bilinearity of the octonion product (correctly, in the same sentence). The "only if" direction — which would say "if $A \cdot B \subseteq C$ holds on the lattices, it holds on the basis pairs" — is trivial (the basis vectors lie in $A$ and $B$, so the lattice-level inclusion restricts to them), but the reader is asked to take this on inspection.

A parenthetical or a single sentence would close the gap. As stated the equivalence is presented as if it were the load-bearing step in the proof, which it is not (the if-direction is); but a careful reader will pause on the wording.

### 3.2 The $\mathbb{Z}$-basis of $L$ is not stated in closed form (Appendix A)

Appendix A defines:

> "The first vector is Wilson's $s = \tfrac12(-e_0+e_1+\cdots+e_7)$; the remaining seven are type-1 roots of $L$ chosen so that $\{L_1,\ldots,L_8\}$ is linearly independent over $\mathbb{Q}$ and spans $L$ over $\mathbb{Z}$."

Inspection of the displayed coordinates shows the seven type-1 roots are $e_0+e_1$, $e_0-e_1$, $e_0+e_2$, $e_0+e_3$, $e_0+e_4$, $e_0+e_5$, $e_0+e_6$ — six of the form $e_0+e_k$ ($1 \le k \le 6$), plus $e_0-e_1$. This is a specific, canonical basis (and a useful one: $L_2 + L_3 = 2e_0$, $L_2 - L_3 = 2e_1$, $L_{k+2} - L_2 = e_k - e_1$ for $k \ge 2$, etc.).

The appendix should state this explicitly. Reasons:

1. **Reproducibility.** A reader who wants to verify any single entry of any of the three tables — or who wants to translate the tables into another basis — needs the basis explicitly, not implicitly via "chosen so that...".
2. **Independence from the supplemental code.** With the recipe stated, the appendix becomes self-contained: the basis vectors are listed in coordinates, the multiplication is defined by the Fano triples, and the products in the tables can be checked by hand or by any independent computer-algebra system.
3. **Structural visibility.** The basis has structure (six roots through $e_0+e_k$ plus $s$ plus an extra $e_0-e_1$) that the reader will not see otherwise. Whether that structure is incidental (an artifact of the algorithm) or meaningful is itself a question worth raising.

Two lines fix this.

### 3.3 The novelty question against Egan's enumeration is now raised but not addressed (§6 Comparison paragraph)

The new comparison paragraph correctly distinguishes the present construction from Baez–Egan and Dixon, and ends:

> "Whether $\Sigma(\Lambda)$ overlaps with one of the embeddings catalogued by Egan, or sits outside that family, is a question this paper does not address."

In the 16 April review, this was the reviewer's question for the author to consider; in the 28 April version, it has become the author's question for the reader. The escalation is fair, but:

- Egan's catalogue contains $17{,}280$ explicit Leech embeddings — finite, enumerable, and (per Egan) with explicit constructive descriptions.
- $\Sigma(\Lambda)$ has an explicit construction in this paper.
- A direct check (is $\Sigma(\Lambda)$ in the catalogue, or in the same $\mathrm{Aut}(\Lambda)$-orbit as one in the catalogue?) is computationally tractable in principle.

A footnote-level statement on whether (a) the check has been attempted, (b) the check is tractable in the supplemental repository's existing code, or (c) the relevant data is not available in the catalogue's published form, would let the reader assess the gap. As the question now sits, a hostile reviewer could argue the paper is dodging its own central comparison.

### 3.4 The routing of cross-block products in Definition 3.6 is still under-specified

The 16 April review flagged this; the 28 April version is unchanged. The formulas for $P_1, P_2, P_3$ make specific choices of ordered pairs (e.g., $z \cdot_\sigma y'$ and $y \cdot_\sigma z'$ both contribute to $P_1$). The "$\{\alpha,\beta,\gamma\} = \{1,2,3\}$" routing rule is symmetric in two indices, but the displayed formulas are not symmetric in the pair — they fix a particular order.

The reader is left wondering: is the displayed routing the unique $\mathbb{Z}_3$-symmetric one consistent with the chosen labelling? If yes, a one-line remark to that effect. If no, then the closure result depends on the choice and the choice should be motivated, ideally by an equivariance argument or by stating which other routings have been tested and how they behaved.

### 3.5 The fifth open-question bullet (§7 Conclusion) contradicts Remark 3.3

The open questions list ends with:

> "Whether other permutations of the Fano triples (beyond transpositions) also yield closed products on $\Lambda$."

But the new Remark 3.3 carefully distinguishes $\sigma$ — a coordinate involution, not a Fano-triple permutation — from Coxeter's $\kappa$, which is a Fano-triple permutation. Under Remark 3.3's terminology, $\sigma$ is not a "permutation of Fano triples" at all. The bullet is therefore ambiguous between:

- (a) other linear coordinate permutations (beyond simple transpositions $\sigma$) that conjugate the multiplication;
- (b) other Fano-triple permutations (in Coxeter's sense) applied on top of the present convention;
- (c) something more general, e.g., longer permutation cycles in the imaginary index set.

Each of (a), (b), (c) is a different question. Pick one or split into two bullets.

---

## 4. Minor mathematical issues

- **"Order" terminology, first occurrence (Theorem 1.1, Corollary 1.3).** The term is used as if established. Remark 1.4 ("Order-closure, not shell-closure") explains the generalisation, but the terminology itself ("non-unital order in a real algebra") remains non-standard. A footnote at first occurrence pointing to Remark 1.4 would prevent the reader from importing assumptions.

- **Closure-of-shells distinction is duplicated.** Remark 1.4 ("Order-closure, not shell-closure") and Remark 2.4 ("Closure of the lattice, not the shells") make essentially the same point: products of minimal vectors land on higher shells, and lattice closure is independent of shell closure. The duplication is benign but consolidating one of them (or making them complementary rather than overlapping) would tighten the paper.

- **The footnote justifying $L\bar s \subseteq L$ (§2.3) appeals to "implicit use".** It reads: "...closed under left-multiplication by $L$ (Wilson's condition (2) for the untwisted product, used implicitly throughout [Wilson2009])." A precise pointer (page or section in Wilson 2009) or a one-line direct argument — e.g., $L\bar s$ is the conjugate of $Ls$, $L$ is conjugation-stable, $Ls \subseteq L$ since $s \in L$, hence $L\bar s = \overline{Ls} \subseteq \bar L = L$ — would make the footnote self-contained.

- **Lemma 4.5 / rem:Ls-not-closed.** Has an explicit witness for $\sigma(Ls) \ne Ls$: the vector $v = (-1, 0, 1, 0, 0, 1, 0, 1)$. The companion claim $Ls \cdot Ls \not\subseteq Ls$ is asserted ("there exist basis vectors $a, b \in Ls$ with $a \cdot b \notin Ls$") without naming them. The supplemental code names them; the paper could, in a parenthetical.

- **Sample sizes in §5's properties table.** The percentages — $\approx 47\%$, $\approx 8.6\%$, $\approx 29\%$, $<\!0.1\%$ — are still given without sample size, distribution, or test set. Carried over from 16 April. A single sentence ("over $N$ uniformly sampled pairs from $\operatorname{Min}(\Lambda)\times\operatorname{Min}(\Lambda)$, …") would establish the experimental protocol.

- **Validation cross-check ([Appendix B](../main.tex#L1170-L1185), "Symbolic proof and numeric validation").** "Both implementations reproduce the same integer coefficients displayed in the appendix" is a strong assertion. It would be tightened by stating either (a) which test in the supplemental repository performs the equality check between the appendix's tables and the implementations' outputs, or (b) that the appendix's tables were programmatically generated from one implementation and re-verified against the other. As stated, the reader cannot tell whether the cross-check is automated or eyeballed.

- **$\Sigma(\Lambda) \ne \Lambda$ in Remark 4.6.** Asserted, with a sketch ("$\sigma(Ls) \ne Ls$") supplied by the surrounding remarks. A direct one-line: "by Remark 4.5, $\sigma(Ls)$ contains a vector not in $Ls$, hence $\Sigma(\Lambda)$ contains a triple whose third-block-sum lies outside $Ls$" would close the gap.

- **"Essentially unique" in Remark 3.5.** "All transposition twists produce the same multiplication table up to basis relabelling." Correct. A one-line proof or citation (the $\mathrm{GL}(3,\mathbb{F}_2)$ orbit on the 7 imaginary points) was requested in the 16 April review and is still missing.

- **Petersson 2018.** Still cited as a workshop lecture. Carried over.

- **Egan's count of $17{,}280$.** Repeated verbatim from the previous draft. A footnote pointing to the specific n-Category Café part where the count is established would help readers chasing the comparison in §6.

---

## 5. Exposition and presentation

- **Appendix B's logging-protocol-violation sentence.** Reads: "The AI agent violated its own logging protocol on three occasions (prompts 011, 018–019, and 029–031), each time caught by the human researcher on review." For a published mathematics paper, even in an appendix, this level of process self-criticism is unusual. The author may wish to consider whether this sentence belongs in the paper or is better confined to the repository's methodology document.

- **Conclusion's open-questions list overlaps with Outlook.** Six bullets in the Conclusion (§7), two subsections in Outlook (§8). They are different in tone (open questions vs. directions of further research), but both touch on the Conway group, ternary structure, and classification. Some consolidation — perhaps merging the Outlook's "Ternary reformulation" into a Conclusion bullet, leaving Outlook for the physics pointer — would reduce redundancy.

- **§6 Comparison paragraph and the earlier precursors.** The paragraph contrasts the present construction with Baez–Egan and Dixon explicitly. The same section names three earlier precursors (Conway–Sloane 1982, Lepowsky–Meurman 1982, Elkies–Gross 1996) and observes that none defines or tests a bilinear product on $\Lambda$. The Comparison paragraph could either pick these up (the same simplicity argument applies — none gives a bilinear product on $\mathbb{R}^{24}$) or explicitly say so by reference. As it stands the comparison addresses only the two more recent works.

- **Proposition 3.3.** The proof writes $\sigma(x) \cdot_\sigma \sigma(y) = \sigma(\sigma(\sigma(x)) \cdot \sigma(\sigma(y))) = \sigma(x \cdot y)$. A reader may briefly worry about which direction is being shown. One short sentence ("setting $x' = \sigma(x)$, $y' = \sigma(y)$ and using $\sigma^2 = \mathrm{id}$") before the calculation removes the ambiguity; the calculation is immediate after that.

---

## 6. Summary of requested revisions

**Mandatory (errors):**

1. Correct the description of $L$ in the proof of Lemma 4.1 (§4.2) so that it matches §2.2 and the supplemental code: "all integers with even sum, or all half-integers with odd sum". (As in the 16 April–round corrections, the conclusion of the proof is unaffected; only the descriptive sentence needs the fix.)

**Strongly recommended:**

2. Justify the "only if" direction of the bilinearity reduction in §4.2's preamble. One sentence.
3. State the $\mathbb{Z}$-basis of $L$ in closed form in Appendix A: $L_1 = s$, $L_2 = e_0 + e_1$, $L_3 = e_0 - e_1$, $L_4 = e_0 + e_2$, …, $L_8 = e_0 + e_6$.
4. Footnote-level statement in §6 on whether the $\Sigma(\Lambda)$-vs.-Egan question has been examined in the supplemental code, is computationally tractable, or is left as future work; if any work has been done, point to it.
5. State whether the routing of cross-block products in Definition 3.6 is uniquely forced by $\mathbb{Z}_3$-equivariance or one of several possibilities. One line.
6. Clarify the fifth open-question bullet in §7 ("other permutations of the Fano triples") so it does not contradict Remark 3.3.
7. Footnote at first occurrence of "order" pointing to Remark 1.4, ideally with a precedent for non-unital terminology in non-associative real algebras.

**Nice to have:**

8. Strengthen the $L\bar s \subseteq L$ footnote with a precise citation or direct one-line argument in place of "implicit use".
9. State how the cross-check between Appendix A's tables and the two computer-algebra implementations is performed.
10. Add sample sizes and protocol to §5's properties table.
11. Consider whether the logging-protocol-violation sentence in Appendix B belongs in the paper.
12. Consolidate the closure-of-shells remarks (Remark 1.4 and Remark 2.4) into one occurrence.
13. Verify the Egan count of $17{,}280$ against the n-Category Café source and add the precise pointer.
14. One-line orientation in Proposition 3.3's proof identifying which equality is being shown.
15. Extend or annotate the §6 Comparison paragraph to cover the three earlier precursors (Conway–Sloane, Lepowsky–Meurman, Elkies–Gross).

---

## 7. Disposition by the author (2026-04-28)

Upon receipt of this review, the author directed the following immediate corrections in the main paper:

- **§2.1 (mandatory):** Lemma 4.1's proof updated to use the corrected description of $L$ ("all integers with even sum, or all half-integers with odd sum"), with a back-reference to §2.

- **§3.1 (strongly recommended):** the bilinearity-reduction preamble now adds a parenthetical noting that the forward direction uses bilinearity and the reverse direction is immediate.

- **§3.2 (strongly recommended):** Appendix A now states the $\ZZ$-basis of $L$ in closed form: $L_1 = s$, $L_2 = e_0+e_1$, $L_3 = e_0-e_1$, $L_{k+2} = e_0+e_k$ for $k=2,\ldots,6$.

- **§3.5 (strongly recommended):** the fifth open-question bullet rewritten under interpretation (a) — other linear coordinate permutations of $\RR^8$ beyond simple transpositions — to remove the conflict with Remark 3.3.

- **§4 minor:**
  - Sample sizes and the sampling protocol added to §5's properties table.
  - Remark 3.5 ("essentially unique") now contains a one-line proof using $2$-transitivity of $\GL(3, \FF_2)$.
  - Remark 4.6 closes the $\Sigma(\Lambda) \ne \Lambda$ gap by chasing a witness through the block sum.
  - Egan's $17{,}280$ count cross-referenced to part 10 of \cite{Baez2014}.
- **§5 exposition:** the logging-protocol-violations sentence in Appendix B has been replaced by a broader statement that mistakes were made on both sides (human and AI) throughout the investigation, and that each was caught by the systematic cross-checking the methodology requires.

The remaining items (§3.3 Egan/Baez novelty question, §3.4 routing uniqueness, and several minor items) are deferred:

- **§3.3** is now an active research direction: the author intends to (a) independently verify Baez's claim about the doubled Jordan product on $\ZZ^3 \oplus L_L$ by direct symbolic computation, (b) document this verification, and (c) answer a structural question relating Baez–Egan's closure anatomy on $\Lambda$ to the present construction. This will be addressed in follow-up prompts and reflected in the paper as a short, factual mention of the verification result, with details deferred to the supplemental repository.

- **§3.4** (routing uniqueness): the paper makes an existence claim, not a uniqueness claim; the question of whether other $\ZZ_3$-symmetric routings yield closed products is genuine but deferred to mathematical follow-up.

- Minor items 8 (precise $L\bar s$ citation), 9 (validation cross-check protocol description), 11 (logging-violations consideration — partially addressed: the specific prompt numbers were removed; the broader point retained), 12 (consolidating closure-of-shells remarks), 14 (Proposition 3.3 orientation), 15 (extending the §6 Comparison to earlier precursors): not addressed in this round; left open.

All deferred items remain open for formal review and are not to be considered settled.
