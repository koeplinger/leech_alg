# Formal referee review of the main paper (round 3)

**Date:** 2026-04-28 (second of the day)
**Reviewer:** Claude Opus 4.7 (Anthropic), acting in the role of a highly trained mathematician at the direction of Jens Köplinger
**Manuscript reviewed:** `paper/main.tex` — *An order on the Leech lattice from a $\mathbb{Z}_3$-symmetric triple-octonion product* (Köplinger), version of 28 April 2026 (post-Prompt 082 edits), together with the supplemental repository as context.

This is a third-pass review, requested ahead of a full line-by-line read by the author. The objective is to flag every substantive residual error and inconsistency that the author would otherwise discover incrementally. Stylistic preferences are deferred. The recommendation is again minor revision: no demonstrable mathematical errors remain, but several inconsistencies and notational hazards survive the edits and should be addressed in one pass.

---

## Items resolved since the morning review (round 2)

For continuity:

- The residual error in the proof of Lemma 4.1 is fixed; the description of $L = D_8^+$ now matches §2.2 and the supplemental code.
- The bilinearity-reduction preamble in §4.2 has the parenthetical justifying the "only if" direction.
- Appendix A states the $\mathbb{Z}$-basis of $L$ in closed form.
- The fifth open-question bullet in §7 has been rewritten under interpretation (a).
- Sample sizes and the sampling protocol added to §5.
- Remark 3.5 now contains a one-line proof of essential uniqueness.
- Remark 4.6 now closes the $\Sigma(\Lambda) \ne \Lambda$ gap in prose.
- Egan's $17{,}280$ count cross-references part 10 of `\cite{Baez2014}`.
- The logging-protocol-violations sentence in Appendix B has been replaced with a balanced both-sides statement.

The disposition of the morning review records these correctly.

---

## Recommendation

Minor revision. No mathematical errors detected. Three notational/expository inconsistencies, one historical deferral, and a handful of small carry-overs from the previous reviews. A single editing pass should resolve them.

---

## 1. Substantive issues

### 1.1 Variable-name collision: $\sigma = (s\;t)$ vs. Wilson's $s$ (Definition 3.1)

[main.tex:285–297](../main.tex#L285-L297). Definition 3.1 introduces:

> "Let $\sigma = (s\;t)$ be a transposition on $\{1, \ldots, 7\}$."

The symbol $s$ is already defined in §2.3 (line 225) as Wilson's half-integer vector $s = \tfrac12(-e_0 + e_1 + \cdots + e_7) \in L$, and it is used throughout the paper in that meaning ($Ls$, $\sigma(Ls)$, "block-sum sublattice", and so on). Using $s$ as a transposition index in $\sigma = (s\;t)$ creates a notational clash: a reader scanning Definition 3.1 will momentarily parse $(s\;t)$ as "the transposition that exchanges Wilson's vector $s$ with some undefined vector $t$".

The two meanings are formally distinct (one is an integer in $\{1,\ldots,7\}$, the other an octonion in $L$), but a referee will object and the casual reader will lose time. Recommend renaming the transposition indices in Definition 3.1 to $a, b$ or $i, j$ — e.g., "Let $\sigma = (a\;b)$ be a transposition on $\{1,\ldots,7\}$" — and propagating through Remark 3.5 ("conjugation of $\sigma = (a\;b)$ by a Fano automorphism…"). The body of the paper does not use $s, t$ as transposition indices except in these two definitions, so the rename is contained.

### 1.2 Bibliography entry `repo` carries the old subtitle (line 1424–1428)

[main.tex:1424–1428](../main.tex#L1424-L1428). The repository bib entry reads:

> "An order on the Leech lattice from transposition-twisted octonions: research repository"

The paper's title was changed in Prompt 076 to "*An order on the Leech lattice from a $\mathbb{Z}_3$-symmetric triple-octonion product*". The repository's published name should match the paper's published name; the bib entry should be updated, and (separately, outside the manuscript) the repository's README and metadata should follow.

### 1.3 The $\Sigma(\Lambda) \ne \Lambda$ argument in Remark 4.6 elides Wilson surjectivity

[main.tex:522–531](../main.tex#L522-L531). The current argument:

> "an element $(x',y',z') \in \Sigma(\Lambda)$ is by construction $(\sigma(x), \sigma(y), \sigma(z))$ for some $(x,y,z) \in \Lambda$, so its block sum $x'+y'+z' = \sigma(x+y+z)$ lies in $\sigma(Ls)$, not in general in $Ls$. By Remark 4.5, $\sigma(Ls) \ne Ls$, so $\Sigma(\Lambda)$ contains triples whose block sum violates Wilson's condition (3) for membership in $\Lambda$."

The "not in general in $Ls$" hedge depends on the (unstated) fact that the block-sum map $\Lambda \to Ls$ is surjective onto $Ls$ — only then does $\sigma$ applied to the image hit elements of $\sigma(Ls) \setminus Ls$. The required surjectivity is a standard property of Wilson's setup but is not asserted in the paper. As written the argument is morally correct but a careful reader will pause: "why does $\sigma$ applied to the actual block-sum image hit something outside $Ls$, rather than just outside the image-of-$\Lambda$-restricted-to-$\sigma(Ls)$?"

A one-clause fix: "…lies in $\sigma(Ls)$. Since the block-sum map $\Lambda \to Ls$ is surjective (Wilson 2009, Section 3) and $\sigma(Ls) \ne Ls$ (Remark 4.5), there exist $(x,y,z) \in \Lambda$ with $\sigma(x+y+z) \in \sigma(Ls) \setminus Ls$, so the corresponding triple in $\Sigma(\Lambda)$ violates Wilson's condition (3) for membership in $\Lambda$."

### 1.4 Corollary 1.3 proof contains a now-unused step

[main.tex:135–140](../main.tex#L135-L140). The proof:

> "$\operatorname{Min}(\Lambda)$ generates $\Lambda$ over $\mathbb{Z}$ [Conway–Sloane 1999]. The product $\star$ is bilinear (it is a sum of bilinear octonion products). By Theorem 1.1, $\star$ maps $\Lambda \times \Lambda$ into $\Lambda$."

Theorem 1.1 already asserts closure on all of $\Lambda \times \Lambda$, so the $\operatorname{Min}(\Lambda)$-generation citation is dead weight: it would matter only if Theorem 1.1 had been stated on the minimal shell and the corollary needed to extend. As things stand the citation reads as an artifact of an earlier draft where closure was first observed on $\operatorname{Min}(\Lambda)$. Suggest removing the first sentence; the bilinearity sentence and the appeal to Theorem 1.1 are sufficient.

### 1.5 "Conway group" open question stated twice

The relationship between $\star$ and $\Co_0$ is stated as an open question in two places:

- §5 (Algebraic properties), [main.tex:676–678](../main.tex#L676-L678): "The relationship between $\star$ and the Conway group $\Co_0$ (the automorphism group of $\Lambda$ as a lattice) is an open question."
- §7 (Conclusion) open-questions list, [main.tex:803–804](../main.tex#L803-L804): "The automorphism group of $(\Lambda, +, \star)$ and its relationship to the Conway group $\Co_0$."

Two formulations of the same question. The §7 bullet is slightly different (asking about the automorphism group of $(\Lambda, +, \star)$ and its relation to $\Co_0$, rather than the relation of $\star$ itself to $\Co_0$), but a reader will plausibly read both as the same point. Pick one place — §7 is the natural home — and remove the §5 sentence (or rewrite it to give a properly different observation, e.g., that $\star$ is not invariant under all of $\Co_0$).

---

## 2. Cosmetic/structural inconsistencies in the proof of Lemmas 4.2–4.4

[main.tex:454–484](../main.tex#L454-L484). Lemma 4.2's proof concludes with "By $\mathbb{Z}$-bilinearity, $L \cdot L \subseteq L$. This recovers the classical fact…". Lemmas 4.3 and 4.4 proofs end with the table reference and stop:

> Lemma 4.3 proof: "The 64 products $L_i \cdot N_j$ are integer combinations of $\{N_1,\ldots,N_8\}$; the explicit coefficients are Table 2."
>
> Lemma 4.4 proof: "The 64 products $M_i \cdot M_j$ are integer combinations of $\{M_1,\ldots,M_8\}$; the explicit coefficients are Table 3."

The conclusion "By $\mathbb{Z}$-bilinearity, [the lemma's inclusion]" is implicit but not stated. For symmetry of presentation across the three lemmas — and for a reader who does not chase back to the §4.2 preamble — adding the corresponding closing sentence to Lemmas 4.3 and 4.4 would help. (One sentence each.)

---

## 3. Carry-overs from earlier reviews

These were noted in the morning review and remain "deferred" in the disposition. Listing them here so a single line-by-line pass can see the full picture:

- **Petersson 2018 cited as a workshop lecture only.** [main.tex:1414–1417](../main.tex#L1414-L1417). If a published version exists (Petersson did publish on integral octonions in this period), it would be preferred over the lecture.
- **Closure-of-shells distinction duplicated.** Remark 1.4 (`rem:order-closure`, [main.tex:142–154](../main.tex#L142-L154)) and Remark 2.4 (`rem:shells`, [main.tex:210–221](../main.tex#L210-L221)) make essentially the same point. They are placed in different sections but a reader sees both.
- **Footnote on $L\bar s \subseteq L$ appeals to "implicit use".** [main.tex:235–243](../main.tex#L235-L243). A precise pointer or a short direct argument would be cleaner than "used implicitly throughout [Wilson 2009]".
- **Remark 4.5 names a witness for $\sigma(Ls) \ne Ls$ but not for $Ls \cdot Ls \not\subseteq Ls$.** [main.tex:486–501](../main.tex#L486-L501). The supplemental code names both; a parenthetical naming the second pair would close the asymmetry.
- **"Comparison" paragraph in §6 omits the earlier precursors.** [main.tex:766–776](../main.tex#L766-L776). The paragraph discusses Baez–Egan and Dixon explicitly but does not extend the comparison to Conway–Sloane 1982, Lepowsky–Meurman 1982, or Elkies–Gross 1996, which are named in the same section as additive/geometric constructions that do not define a bilinear product.
- **Validation cross-check between appendix tables and the two implementations not specified.** [main.tex:1230–1246](../main.tex#L1230-L1246). The reader cannot tell whether the cross-check is automated or eyeballed.
- **First-occurrence footnote on the "order" terminology not yet added.** Theorem 1.1 / Corollary 1.3 use the term as if standard; Remark 1.4 explains the generalisation but does not get a forward-pointer.

---

## 4. Things checked and confirmed correct

For the record (and because the author asked for thoroughness):

- The seven Fano triples [eq:fano] form a Fano plane (every pair of imaginary indices appears in exactly one triple).
- The block-matrix for left-multiplication in §5's multiplicative-identity remark is correct.
- The "30 of 64 entries change" count in Remark 3.4 is consistent with $\sigma$ fixing $e_0$ (15 entries unchanged), squares $e_i \cdot_\sigma e_i = -e_0$ (7 unchanged), leaving 42 off-diagonal entries with both indices $> 0$ to apportion.
- The argument that $L \cdot_\sigma L = \sigma(L \cdot L) \subseteq L$ in §6 (Kirmse–Coxeter–$\sigma$ paragraph) is correct.
- The Lemma 4.1 proof, after the morning fix, is consistent with §2.2's description of $D_8^+$.
- Equation labels (`eq:fano`, `eq:twist-formula`, `eq:iso`, `eq:triple`, `eq:P1`–`eq:P3`, `eq:P1-std`–`eq:P3-std`) all resolve.
- All theorem/lemma/proposition cross-references (`thm:main`, `lem:sigma-L`, `lem:L-order`, `lem:L-sigmaLsbar`, `lem:sigmaLs-closed`, `prop:cond1`, `prop:cond2`, `prop:cond3`, `cor:order`, `def:twist`, `def:triple`, `def:wilson-leech`, `rem:*`) resolve.
- Bibliography entries are well-formed with no missing fields; no undefined citations.

---

## 5. Summary of requested revisions

**Mandatory (errors):** none.

**Strongly recommended (in order of mechanical cost):**

1. Rename the transposition indices in Definition 3.1 (and Remark 3.5, where they recur) from $(s\;t)$ to $(a\;b)$ or $(i\;j)$ to remove the collision with Wilson's $s$.
2. Update the `\bibitem{repo}` title to match the paper's current title ("from a $\mathbb{Z}_3$-symmetric triple-octonion product"), and propagate the rename to the repository's README/metadata.
3. Tighten the $\Sigma(\Lambda) \ne \Lambda$ argument in Remark 4.6 by invoking Wilson surjectivity onto $Ls$ explicitly (one clause).
4. Remove the unused $\operatorname{Min}(\Lambda)$ sentence from Corollary 1.3's proof.
5. De-duplicate the Conway-group open-question statement: keep §7's bullet, remove §5's sentence.
6. Add the "By $\mathbb{Z}$-bilinearity, $A \cdot B \subseteq C$" closing sentence to Lemmas 4.3 and 4.4 proofs for parallel form with Lemma 4.2.

**Carry-overs (from earlier reviews) — list above, reproduced for completeness:**

7. Petersson 2018 reference quality.
8. Consolidate the closure-of-shells remarks (1.4 / 2.4).
9. Strengthen the $L\bar s$ footnote with a precise citation or short direct argument.
10. Name a witness for $Ls \cdot Ls \not\subseteq Ls$ in Remark 4.5.
11. Extend the §6 "Comparison" paragraph to cover the earlier precursors.
12. Specify the validation cross-check protocol in Appendix B.
13. Footnote at first occurrence of "order" pointing to Remark 1.4.

---

## 6. Disposition by the author (to be added)

Author to mark each item as accepted, deferred, or declined, before the next investigative step (Egan/Baez deep dive).
