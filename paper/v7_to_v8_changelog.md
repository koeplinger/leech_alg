# v7 → v8 revision notes

**v7** — dated 19 July 2026, 24 pages. Frozen at [paper/main_2026-07-19.tex](main_2026-07-19.tex); the last released version, published as [doi:10.13140/RG.2.2.22093.19686](https://doi.org/10.13140/RG.2.2.22093.19686) (ResearchGate).

**v8** — **in progress, not frozen.** Opened 12 October 2026 ([prompt 281](../prompt_logs/281_reopen_project_v8_readthrough_comments.txt)), after the closure of the effort on 27 July 2026 ([CLOSING_NOTE.md](../CLOSING_NOTE.md), with its dated addendum). The working source is [paper/main.tex](main.tex). This changelog is self-contained: it records every v7 → v8 change on its own, and the frozen v7 artifacts (including [v6_to_v7_changelog.md](v6_to_v7_changelog.md) and everything earlier) are not touched.

The planned changes are recorded item by item in [evidence_and_reasoning/2026-10-12_plan.md](../evidence_and_reasoning/2026-10-12_plan.md) and, for item 12, in [evidence_and_reasoning/2026-10-14_plan.md](../evidence_and_reasoning/2026-10-14_plan.md); each enters this changelog when it is applied.

---

## 0. Opening of v8 (12 October 2026)

- **What changed.** Subtitle "19 July 2026 (v7)" → "12 October 2026 (v8)"; set to "13 October 2026 (v8)" at prompt 286 and to "14 October 2026 (v8)" at prompt 289. No change to the text.
- **Where.** Title block.
- **Rationale.** Marks the working source as diverging from the frozen, published v7.

## 1. Abstract (12 October 2026)

*Revised in part at prompt 284 (review response); the wording now in the paper is in item 11.*

- **What changed.** Two sentences. After the four-lemma sentence: "We determine the idempotents, the automorphism group, the square-zero elements, and the span of the image of the algebra." (the author's draft read "We show select algebraic properties of the algebra, including ..."; reworded for idiom, same content). The sentence "Reduced modulo $2L$, the construction yields a Witt decomposition of $L/2L$ into a complementary pair of Lagrangians, and we ask which polarizations lift." is replaced by "We discuss the construction in the context of related work."
- **Where.** Abstract.
- **Rationale.** The abstract now says what Section 5 does. The Witt decomposition is a statement about $L/2L$ (Appendix A.1), an $E_8$-level fact, and the abstract no longer singles it out; Section 6 is where the construction is set against related work. Word count 149 → 155; v7's 150-word cap no longer applies (the envisioned venue recommends 150 to 250 words).

## 2. Section 5.1: lead-in

*Revised in part at prompt 284 (review response); the wording now in the paper is in item 11.*

- **What changed.** "The order $(\Lambda, +, \star)$ has been tested for several standard algebraic identities. Samples are drawn uniformly ..." → "Several standard algebraic identities are easily checked to fail on the order $(\Lambda, +, \star)$. To give an initial idea of the degree of their violation, we have drawn samples uniformly ...". The sampling details that follow are unchanged.
- **Where.** Section 5.1, first paragraph.

## 3. Section 5.1, the no-identity remark: the block entries act in the twisted algebra

- **What changed.** "(entries acting by left-multiplication in $\mathbb{O}$)" → "in $(\mathbb{O}, \cdot_\sigma)$". The entries of the displayed block matrix multiply in the twisted product of Definition 3.2; for the $e_0$ entries of the remark's argument the two readings coincide.
- **Where.** Section 5.1, Remark 5.2.

## 4. Section 5.2: idempotents

*Revised in part at prompt 284 (review response); the wording now in the paper is in item 11.*

- **What changed.** The clause "a symbolic classification carried out and verified in exact arithmetic" is cut; the footnote citing `verify_idempotent_classification.py` stays. After eq. (idempotent-multiples): "On the line spanned by a nonzero idempotent, the points of $\Lambda$ are the integer multiples $kw$ of the corresponding vector $w$ in (idempotent-multiples), of norm $16k^2$, $48k^2$ or $24k^2$; in particular, no real multiple of an idempotent lies on the minimal shell."
- **Where.** Section 5.2.
- **Verification.** $\Lambda \cap \mathbb{R}E = \mathbb{Z}\,t_0 E$ for the least $t_0 > 0$ with $t_0 E \in \Lambda$, and the least multiples with their norms are those of eq. (idempotent-multiples) (`verify_idempotent_lattice_rescaling.py`, ALL PASS, exhaustive over the minimal shell); second pass in `evidence_and_reasoning/2026-10-12_plan.md`, Section 4.

## 5. Section 5.3: the automorphism group

*Revised in part at prompt 284 (review response); the wording now in the paper is in item 11.*

- **What changed.** The sentence "Every automorphism is a blockwise octonion automorphism followed by a permutation of the three blocks; in the explicit demonstration below, blockwise means with respect to the twisted algebra ..." is removed; the second paragraph states that fact in the ambient form. Its load-bearing part, that the displayed generator is an automorphism of $(\mathbb{O}, \cdot_\sigma)$ and not of the standard product (review finding M2, prompts 245, 260, 261), now sits with the generator: "the $C_6$ is generated by the following signed permutation of imaginary units, an automorphism of the twisted algebra $(\mathbb{O}, \cdot_\sigma)$ of Definition 3.2, applied in each block, with $e_0$ fixed:", followed by three numbered equations, $e_1 \mapsto -e_5$, $e_5 \mapsto e_6$, $e_6 \mapsto e_1$; $e_2 \mapsto -e_7$, $e_7 \mapsto -e_3$, $e_3 \mapsto e_2$; $e_4 \mapsto -e_4$. The closing sentence on the companion note is a paragraph of its own, as it applies to both the order and the ambient algebra.
- **Where.** Section 5.3.
- **Verification.** Orbit structure and twisted-versus-standard preservation re-verified exactly (plan, Section 4): the map has order 6 (its cube is $-1$ on the first orbit, $+1$ on the second, $-1$ on $e_4$), and the committed GAP generator is this map.

## 6. Section 6: cross-reference to Appendix A.1

- **What changed.** In the Baez--Egan paragraph, after "equivalently the number of $4$-dimensional complementary isotropic subspaces": "; Section A.1 reads one such pair, $\sigma(L\bar s)/2L$ and $\sigma(Ls)/2L$, as ideal and subalgebra under the product that $L/2L$ inherits from the octonion order $L$."
- **Where.** Section 6, Baez--Egan paragraph.
- **Rationale.** Appendix A.1 pointed to Section 6; the reference back was missing.

## 7. Section 7: the order is not maximal; the span-index question

*Revised in part at prompt 284 (review response); the wording now in the paper is in item 11.*

*The span-index bullet was deleted at prompt 289; see item 12.*

- **What changed.** New paragraph after the opening one: "The order is not maximal. The idempotent $\varepsilon_1 = (e_0, 0, 0)$ of Section 5.2 acts on either side as the swap of blocks 2 and 3, $\varepsilon_1 \star v = v \star \varepsilon_1 = (x', z', y')$ for $v = (x', y', z')$, a swap that preserves $\Lambda$ because Wilson's three conditions are symmetric in the blocks. Moreover $\varepsilon_1 \star \varepsilon_1 = \varepsilon_1$, so $\Lambda + \mathbb{Z}\varepsilon_1$ is $\star$-closed by bilinearity. Since $4\varepsilon_1$ is the least positive multiple of $\varepsilon_1$ in $\Lambda$ (Section 5.2), this overlattice has index 4; it is not contained in $L^3$ because $e_0 \notin L$." The bullet "Whether the order is maximal ..." is replaced by "Whether a closed bilinear product on $\Lambda$ can have a smaller span index than $\star$: the index (Section 5.5) is $2^{25}$ for the Baez--Egan product $\varphi$ (Section 6) and $2^{16}$ for $\star$."
- **Where.** Section 7.
- **Verification.** The witness was checked by hand from Definition 3.3 and in exact arithmetic on two independent code paths (plan, Section 4; key claim 016). Both span indices were re-run: `verify_phi_span_index.py` gives $2^{25}$, `verify_product_span_index.py` gives $2^{16}$.

## 8. Appendix A.1: last clause

- **What changed.** "and whether any such polarization lifts to a closed bilinear product on $\Lambda$ itself" → "and which polarizations, beyond the one constructed here, lift to a closed bilinear product on $\Lambda$ itself". Theorem 1.1 answers the question for the polarization constructed.
- **Where.** Appendix A.1, closing paragraph.

## 9. Appendix C and Acknowledgments: the AI's role; the Hermes attribution; "Python and GAP"

*Revised in part at prompt 284 (review response); the wording now in the paper is in item 11.*

- **What changed.** (a) After "and contributed mathematical analysis.": "It also drafted the manuscript text and made every revision to it, at the author's prompt-by-prompt direction." (b) The Hermes sentence moves from the Acknowledgments to Appendix C, as a paragraph of its own: "Structural proposals from a second AI system (the Hermes agent, Nous Research, running DeepSeek V4 Flash through OpenCode) prompted the Lagrangian framing of Section A.1; the AI agent (Claude Opus 4.7) verified and corrected them." (The author's draft ended "; Claude Opus 4.7 verified and corrected them"; "the AI agent" ties the model name to the agent named at the start of Appendix C.) The Acknowledgments keep the five named people. (c) "independent reimplementation in two computer-algebra systems" → "independent reimplementation in Python and GAP": Python with `fractions.Fraction` is exact arithmetic, not a computer-algebra system. The two other "computer-algebra" wordings (Appendix A, before the basis tables; the "Symbolic proof and numeric validation" paragraph) stay, by the author's decision.
- **Where.** Appendix C; Acknowledgments.
- **Record.** The verifying model and the OpenCode/DeepSeek chain are in `paper/hermes_mod2_structure_supplement.tex` and prompt 156; "Nous Research", "the Hermes agent" and "a second AI system" are stated on the author's authority (prompt 282). The proposals reached the repository as GitHub issue #1 (closed 25 May 2026). The frozen supplement `paper/hermes_mod2_structure_supplement.tex` describes Hermes as "an external contributor identified as Hermes (using OpenCode with DeepSeek)"; it stays as released, and the paper's description is the author's later information.

## 10. Declarations: use of AI

*Revised in part at prompt 284 (review response); the wording now in the paper is in item 11.*

- **What changed.** "the human-run verification protocol" → "the human-designed verification and governance process"; appended: "The author directed the work prompt by prompt and reviewed the AI-assisted output; machine-obtained results were verified under the protocol documented there. The author takes full responsibility for the content, including the accuracy of the references."
- **Where.** Declarations.

## Trackers that follow the paper

`CURRENT_STATE.md` (open questions 1 and 3; "Structure of the algebra": non-maximality), `README.md` (result bullet; "What Remains Open"), `evidence_and_reasoning/research_result.md`, `evidence_and_reasoning/key_claims/README.md` (claim 016, ledger-only row; row 011 wording), `evidence_and_reasoning/README.md` (ledger-row range).

## 11. Journal-quality review response (13 October 2026)

*Revised in part at prompt 289: the Section 7 span-index bullet (MF1) was deleted, and the Section 1 sentence (SF4) was corrected; see item 12.*

A journal-quality review of the assembled v8 ([reviews/2026-10-12_claude_fable_5_1_review_v8.md](reviews/2026-10-12_claude_fable_5_1_review_v8.md); 2 must-fix, 11 should-fix, 5 nits, 1 refuted; every verified finding adversarially checked twice) and the author's rulings (prompt 284) produced these changes:

- **MF1, Section 7 bullet.** The span-index question names its family: "Whether a closed bilinear product on $\Lambda$ built from octonion multiplication, as $\varphi$ and $\star$ both are, can have a smaller span index than $\star$: ...". As posed before, index 1 is attained by a degenerate product.
- **MF2, Section 5.2.** "has exactly eight idempotents, by symbolic classification" re-supplies the antecedent of Section 5.4's "again a symbolic classification".
- **SF1, abstract.** "We determine the idempotents, the automorphism group, and the square-zero elements of the algebra, and the index of the span of the products of the order."
- **SF2, abstract.** "the proof reduces to four lemmas, three of them certified by integer coefficient tables" (Appendix C: three lemmas carry certificates, the fourth is a parity argument).
- **SF3, Section 7.** Non-maximality leads with $L^3$: the proof of Proposition 4.6 uses only Wilson's condition (1), so $L^3 \star L^3 \subseteq L^3$, an even integral overlattice of index $[L^3 : \Lambda] = 2^{12}$; then $\Lambda + \mathbb{Z}\varepsilon_1$ as a second overlattice, of index 4, not contained in $L^3$ and not integral ($\varepsilon_1$ has half-integer inner products with vectors of $\Lambda$). The opening sentence, "The order is not maximal: it is properly contained in the $\star$-closed lattice $L^3$", supplies the sense of "maximal". Verified at prompts 283 and 284: $L$ is closed under $\cdot_\sigma$ on all 64 basis products (prompt 283); the first coordinate over the $\mathbb{Z}$-basis of $\Lambda$ takes the value $3/2$.
- **SF4, Section 1.** "Seven octonion-multiplication conventions on $\mathbb{R}^8$ have this property; the resulting structures are in bijection with the seven maximal orders ..." → "For a fixed octonion multiplication, the integral octonions have seven maximal orders, each an $E_8$ lattice [Coxeter1946, ConwaySmith2003] (Appendix B); the placement of $L$ used here is fixed in Section 2." The old sentence does not hold for the coordinate-symmetric $L$: a signed permutation of the imaginary units maps $L$ to itself or to the mirror lattice $\{-\bar x : x \in L\}$, also closed under the standard product (checked exactly), so $L$ is closed under every one of the 480 multiplication tables. Sections 1 to 4 had been declared read-OK; the author reviews this change in context.
- **SF5, Appendix C.** Footnote after "verified and corrected them": the supplement `paper/hermes_mod2_structure_supplement.pdf` and the two discovery scripts, in the paper's path-citation convention.
- **SF6, Declarations.** "AI assistance is a stated, material part of this work, including the drafting of the manuscript text; the methodology, the role of AI, and the human-designed verification and governance process are documented in Appendix C."
- **SF8, Section 5.1.** The lead-in carries the list ("among them associativity, alternativity, flexibility, power-associativity and commutativity"); the closing paragraph is reduced to "The defining property of the order is closure: the Leech lattice is stable under the product."
- **SF9, Section 5.3.** The generator display uses one alignment point per line with `\quad` separators (standard `align`, nothing journal-specific), so each cycle reads as one sequence.
- **Nits.** Section 7: "$V := \sigma(L\bar s)/2L$ is an ideal of $L/2L$". Section 6: "an order in the sense of Corollary 1.2". "repositories of [repo]" in Appendix C, matching the Declarations.
- **Not applied, by the author's ruling.** SF7 (a script footnote for "over 12,000,000 pairs"); SF10 (the ambient algebra as $(\mathbb{O}, \cdot_\sigma) \otimes \mathbb{R}[\mathbb{Z}_3]$ with the inversion-twisted product; the paper stays on its expository line); SF11 (a page-15 widow; no forced layout in the preprint); the Hermite-normal-form nit (the method is the machine method of the footnote); the Remark 1.3 "distinct" nit (left as settled at prompt 277).

Abstract 165 words. Page count unchanged at 24.

## 12. The span-index findings, and Section 1's seven maximal orders (14 October 2026)

The external review of v8 by Claude Opus 5.5 ([reviews/2026-10-13_claude_opus_5_5_review_v8.md](reviews/2026-10-13_claude_opus_5_5_review_v8.md), with its evidence) and the second pass on it ([reviews/2026-10-14_claude_fable_5_1_index_one_verification.md](reviews/2026-10-14_claude_fable_5_1_index_one_verification.md)) showed that the span index depends on the twist. The author's goal for the wrap-up was "less is more"; the plan is [../evidence_and_reasoning/2026-10-14_plan.md](../evidence_and_reasoning/2026-10-14_plan.md), approved at prompt 289.

- **Sections 5.5 and 5.6 merged** into one subsection, "Other twists: closure and span of the image" (5.5; both labels kept, so every reference resolves). It defines the *span index* in words; states the span index $2^{16}$ of $\star$ and that the span contains $2\Lambda$; adds a span-index column to the census table ($2^{16}$ for the closing odd cycle types, $2^{10}$ for the even ones, n/a for the 4-cycles, none of which closes) and the sentence "The span index depends only on the parity of $\pi$"; and closes with one paragraph: $\sigma$ is the reflection in the root $e_1 - e_2$, the reflection $\rho$ in the root $r = \tfrac12(e_0 + e_1 + \cdots + e_6 - e_7)$ gives a triple product $\star_\rho$ that closes on $\Lambda$ with span index 1, and "The span index is invariant under ring isomorphism (for a unital order it would always be 1), so $\star_\sigma$, $\star_\pi$ for an even closing $\pi$, and $\star_\rho$ make $\Lambda$ into three pairwise non-isomorphic orders." Removed from the old Section 5.5: its opening, which defined the symbol $S$ ("Let $S := \mathbb Z$-span$\{u \star v : u, v \in \Lambda\}$.  By bilinearity, $S$ is generated by the $24 \times 24 = 576$ basis-pair products of a $\mathbb Z$-basis of $\Lambda$, so its index in $\Lambda$ is exactly computable:"; $S$ no longer appears in the paper), the displayed equation (5.7), the coset sentences ("The index counts the cosets ... unreachable"), "span-surjective", "$S$ is the full preimage of an 8-dimensional $\mathbb F_2$-subspace of $\Lambda/2\Lambda$, and $\Lambda/S \cong (\mathbb Z/2)^{16}$", the 2-primary remark tying the index to Section A.1, and "Identifying the subspace $S/2\Lambda$ inside $\Lambda/2\Lambda$ structurally is open" (answered by the review, $S_\sigma = \Delta(Ls) + 2\Lambda$). That opening and the old census lead-in "Testing closure exactly on the 576 basis-pair products ..." became the one bilinearity sentence; the lead-in now names the product $\star_\pi$; the paragraph after the table is kept word for word, with the parity sentence appended. New script citation: `verify_twist_span_index_census.py`.
- **Section 7: the span-index bullet deleted.** It is answered by Section 5.5.
- **Section 7, first bullet: the structural pointer.** "... and the same asymmetry survives reduction modulo $2L$: $V := \sigma(L\bar s)/2L$ is an ideal of $L/2L$, $W$ only a subalgebra.  An integer-level structural reason would replace the appendix tables by a structural argument resting on something other than an ideal property: an identification of $\sigma(Ls)$ with the image under a norm-preserving octonion automorphism, or with some module-theoretic structure over $L$ that lifts the mod-2 picture." → "... and this ideal property has a structural reason: $\sigma(L\bar s) = 2O$ for a maximal order $O \supseteq L$, so Lemma 4.3 follows from $O \cdot O \subseteq O$.  Likewise $\sigma(Ls) = pO$ with $p = e_2 - e_4 + e_5 - e_6$, and the open question becomes why $pO$ is closed." Earlier in the bullet, "$W := \sigma(Ls)/2L$ is a subalgebra" → "$\sigma(Ls)/2L$ is a subalgebra", since $W$ is no longer used in Section 7. New script: `verify_sigma_sublattice_maximal_order.py`.
- **Section 1: the seven maximal orders.** "For a fixed octonion multiplication, the integral octonions have seven maximal orders, each an $E_8$ lattice" → "For a fixed octonion multiplication, exactly seven maximal orders of the rational octonions contain the Gravesian order $\mathbb Z\langle e_0, \ldots, e_7\rangle$, each a scaled copy of the $E_8$ lattice". The old sentence was false as stated: the rational octonion algebra has infinitely many maximal orders, and the classical seven are those containing the Gravesian order (review, Section 6.1(a)). Verified for the multiplication (2.1) by `verify_sigma_sublattice_maximal_order.py`, exhaustively: of all 902 doubly-even codes $C$ of length 8, exactly 16 give lattices $\mathbb Z^8 + \tfrac12 C$ closed under the product, each inside one of seven, and the seven are unimodular for the trace form.
- **Not changed:** "ozone algebra" stays out of the paper (the author's ruling, prompt 289). Appendix B keeps "acting now on a different object" (the author's wording, prompt 143); the plan, Section 4, notes a possible qualification for $L\bar s$ and leaves it to the author. Section 2.2 and the proof of Lemma 4.2 keep "maximal order(s) of the integral octonions": the request covered Section 1 only.
- **Subtitle:** "14 October 2026 (v8)".

Build: 24 pages, no errors, no undefined references, no overfull boxes; two more underfull boxes than before (no action, per the author's generic-layout ruling).

## Record notes (12 October 2026)

- `python_project/src/symbolic_proof_checks.py`: its docstring described doubled-integer arithmetic and a `multiply_doubled()` helper that the file never had; corrected in place at the author's direction, flagged "Updated 12 October 2026" (prompt 283). The script works in exact rational arithmetic (`fractions.Fraction`). v7's abstract wording "four lemmas certified in exact integer arithmetic" (v8, item 11: "three of them certified by integer coefficient tables") rests on the Appendix A certificates, which are integer coefficient tables; the v6 → v7 changelog's justification by way of that docstring is superseded by this note, and a dated addendum in `v6_to_v7_changelog.md` (13 October 2026) points here.
- The frozen supplement `paper/hermes_mod2_structure_supplement.tex` keeps its "external contributor" wording (author's ruling, prompt 283); traceability runs through the paper.

---

*Revision chain: [v3 → v4](v3_to_v4_summary.md), [v4 → v5](v4_to_v5_changelog.md), [v5 → v6](v5_to_v6_changelog.md), [v6 → v7](v6_to_v7_changelog.md), v7 → v8 (this file).*

— Claude (Anthropic), at the direction of Jens Köplinger; opened 2026-10-12, in progress
