# Journal-quality review of the assembled v8 (12 October 2026)

**Date:** 2026-10-12
**Reviewer:** Claude Fable 5.1, at the direction of Jens Köplinger
**Manuscript:** `paper/main.tex`, v8 (12 October 2026), 24 pp., in progress
**Prompt:** 283.
**Method:** six read-only referee dimensions (mathematical correctness of the
assembled text; internal consistency and cross-references; reproducibility and
the experiment-versus-proof relationship; presentation against the editorial
standards; attribution, references and the research-record statements; a
referee's verdict for a journal in algebra or number theory with a
computational component), then deduplication and two independent adversarial
refuters per finding. 54 agents. 50 raw findings, 43 after deduplication; the
24 of highest severity were verified (**9 confirmed, 14 split, 1 refuted**),
the remaining 19 are listed below as unverified. No agent wrote to the
repository. The orchestrator re-verified the consequential confirmed and split
findings first-hand, in exact arithmetic where they are mathematical (the
items marked ✓ below).

Line numbers refer to the v8 working tree at review time. The author's
standing rulings (prompts 282 and 283) were not re-litigated: the two
remaining "computer-algebra" wordings, the Hermes sentence and the frozen
supplement's wording, the abstract length, Remark 1.3's sampled/exhaustive
wording, and the σ-free statement of Theorem 1.1 (prompt 175).

---

## Headline

The v7 → v8 edits are mathematically sound: every changed statement was
re-derived by hand and in exact arithmetic (the idempotent rays, the ε₁
witness, the C₆ generator against the committed GAP generator, the Section 6
clause on the pair σ(Ls̄)/2L, σ(Ls)/2L, Remark 5.2's twisted-product reading,
both span indices re-run). Cross-references, labels, citations and the
bibliography are clean at 24 pages; the registry matches every cited work in
the changed passages; sampled and exhaustive results are labeled as such.

Two **must-fix** items, both consequences of the v8 edits: the new Section 7
open question is settled trivially as posed, and the Section 5.2 cut left a
dangling "again" in Section 5.4. The most substantive **should-fix** is one a
referee will raise: L³ itself is a ⋆-closed overlattice of Λ (the proof of
Proposition 4.6 uses only Wilson's condition (1)), so non-maximality is
immediate from Section 4, and the ε₁ paragraph should lead with that and
present ε₁ as the second example, not contained in L³. A second substantive
item is pre-existing: the Introduction's "seven octonion-multiplication
conventions" sentence does not hold for the paper's coordinate-symmetric L.

---

## Must-fix

### MF1. The new span-index question is settled trivially as posed (§7, lines 1088–1091) ✓

> \item Whether a closed bilinear product on~$\LL$ can have a smaller
> span index than~$\star$: the index (Section~\ref{sec:span}) is
> $2^{25}$ for the Baez--Egan product~$\varphi$
> (Section~\ref{sec:related}) and $2^{16}$ for~$\star$.

Read with the paper's own definition of order (Corollary 1.2 and Remark 1.3:
bilinearity and Λ ⋆ Λ ⊆ Λ, nothing more), the quantifier ranges over all
bilinear products under which Λ closes, and index 1 is attained: with a the
dual vector of a basis vector b₁ (2a ∈ Λ, so ⟨u, a⟩ ∈ ℤ on Λ), the product
u ∘ v := ⟨u, a⟩ v is bilinear, Λ-valued, and b₁ ∘ v = v, so its products span
Λ. Both refuters reproduced this on the repository's ℤ-basis. The comparison of
2²⁵ and 2¹⁶ carries content only for products of octonionic origin, which is
what the author's comment (plan item 7) meant.

**Fix:** name the family. For example: "Whether a closed bilinear product
on~$\LL$ built from octonion multiplication, as $\varphi$ and $\star$ both
are, can have a smaller span index than~$\star$: ...". Any equivalent
restriction (the permutation-twisted triple products of Section 5.6 together
with the projected Jordan products) serves.

### MF2. "again a symbolic classification" has lost its antecedent (§5.4, line 766) ✓

> to zero, again a symbolic classification.

The "again" pointed to the Section 5.2 clause "a symbolic classification
carried out and verified in exact arithmetic", cut in v8 (changelog item 4);
the two clauses were added together at prompt 248 to resolve a presubmission
should-fix asking that the two classifications be named as symbolic (they
cannot be established by enumeration over ℝ²⁴). "symbolic classification" now
occurs only here.

**Fix options:** (a) delete ", again a symbolic classification" (the footnote
cites the script); or (b) restore a two-word status label in 5.2, "has exactly
eight idempotents, by symbolic classification", which re-supplies the
antecedent and keeps the in-text statement that the count is a theorem, not a
search result. Option (b) revisits the prompt-281 cut on new grounds (the
status label was the clause's purpose); the author's call.

---

## Should-fix

### SF1. Abstract: "the span of the image of the algebra" (line 106) ✓

The rewording at prompt 282 attached "of the algebra" to all four
determinations. For the ambient algebra the span of the image is trivially
ℝ²⁴ (already ε₁ ⋆ v = (x′, z′, y′) runs over ℝ²⁴); what Section 5.5 determines
is the index [Λ : S] = 2¹⁶ of the ℤ-span of the products, a property of the
order. **Fix:** "We determine the idempotents, the automorphism group, and
the square-zero elements of the algebra, and the index of the span of the
products of the order." (One refuter adds that "determine the span" slightly
overclaims, since S/2Λ is left open; "the index of the span" is exact.)

### SF2. Abstract: "four lemmas certified in exact integer arithmetic" (line 104) ✓

Appendix C (lines 1728–1731): "three of the four lemmas are established by
exhibiting a finite list of integer coefficients (Appendix A); the fourth,
σ(L) = L, is a parity argument." Section 4 scopes the coefficient tables to
"the next three lemmas". The v6 → v7 justification of the phrase rested on the
script docstring corrected today (the script runs in exact rationals), and the
v8 record note rests it on the Appendix A tables, which cover three lemmas.
**Fix options:** "the proof reduces to four lemmas, three of them certified by
integer coefficient tables"; or "four lemmas certified in exact arithmetic"
(Section 1's own wording). v7 text; the author's call.

### SF3. Non-maximality: lead with L³; placement; the word "maximal" (§7, lines 1045–1054) ✓

Content. The proof of Proposition 4.6 uses only Wilson's condition (1) with
Lemmas 4.1 and 4.2, so it proves L³ ⋆ L³ ⊆ L³: L³ is a ⋆-closed overlattice of
Λ, even and integral, of index [L³ : Λ] = 2¹² (re-checked: L is closed under
·σ on all 64 basis products; the index is the 4096 of the ℤ-basis
construction). Non-maximality is therefore immediate from Section 4, and the
paragraph's closing "it is not contained in L³" reads as if L³ were a boundary
the witness must leave, when L³ is itself a closed overlattice. The record
(key claim 016, CURRENT_STATE, the plan's Section 4) nowhere notes this.
**Fix:** state L³ first, then ε₁ as a second overlattice not contained in L³,
of index 4 and not integral; adjust key claim 016 and the trackers.

Placement. A referee reads the Conclusion as recalling results, not proving
them; the result could sit in Section 5.2 as a proposition after the ray
sentence, with the Conclusion recalling it in one sentence. The author placed
it in Section 7 deliberately (prompt 281); recorded as the referee's
preference, not a defect.

The word. With the v7 bullet gone, "maximal" is defined nowhere for the
paper's nonunital order (its other uses are the classical maximal orders of
the integral octonions). A clause such as "it is properly contained in the
⋆-closed lattice L³" supplies the sense; both refuters rate this a nit.

### SF4. "Seven octonion-multiplication conventions" (§1, lines 136–139) ✓

> Seven octonion-multiplication conventions on~$\RR^8$ have this property;
> the resulting structures are in bijection with the seven maximal orders of
> the integral octonions \cite{Coxeter1946, ConwaySmith2003}.

For the paper's L the count is not seven. The 480 octonion multiplication
tables on e₀, …, e₇ are the signed-permutation transports of the Fano triples
(2.1); a signed permutation of the imaginary units maps L to itself or to the
mirror lattice L′ = {−x̄ : x ∈ L} (a single sign flip does the latter, checked),
and L′ is closed under the standard product because x ↦ −x̄ is an
anti-automorphism (checked on the 64 basis products; L′ ≠ L, e₀ ∉ L′). So L is
closed under every one of the 480 tables, and "seven" has no reading under
which it holds for L; "in bijection" has none either. The classical fact
(Appendix B) is the other way round: for a fixed multiplication there are
seven maximal orders, related by transpositions of imaginary units. Sections
1 to 4 were declared read-OK at prompt 281; this is raised because it is a
factual statement a referee will stop at. **Fix:** replace by the classical
statement, e.g. "for a fixed octonion multiplication the integral octonions
have seven maximal orders, each an E₈ lattice, related by transpositions of
imaginary basis units \cite{Coxeter1946, ConwaySmith2003} (Appendix B)."
(The v4 review had accepted the sentence under a reading the computation now
excludes.)

### SF5. The Hermes verification claim has no locator (Appendix C, lines 1714–1717)

"verified and corrected them" is the paper's one verification claim without a
pointer to its record; the paper never names
`paper/hermes_mod2_structure_supplement.tex`. **Fix:** a footnote in the
paper's path-citation convention (as for the companion note in Section 5.3):
"\cite{repo} \nolinkurl{paper/hermes_mod2_structure_supplement.pdf}". The
author excluded new references from the v8 edits; this one is a locator for a
claim the edits introduced, so it is put for ruling.

### SF6. The Declarations omit the drafting disclosure (lines 1156–1163)

Appendix C now states that the AI drafted the manuscript text and made every
revision; the journal-facing "Use of AI" paragraph does not carry it, while
every other Appendix C statement about the AI's role has its counterpart
there. **Fix:** add the clause, e.g. "AI assistance is a stated, material part
of this work, including the drafting of the manuscript text; ...". Author's
wording.

### SF7. "over 12,000,000 random pairs" has no script footnote (Appendix C, lines 1722–1724; abstract line 102)

The only computational claim without a `\scriptref`. The record composes the
figure as 593,412 + 3 × 4,000,000 = 12,593,412 (key claim 008); the three
4,000,000-pair runs are reproducible (`trial_007_scaled_test.py`,
`trial_007_fast.py`, `trial_007_exhaust.py`), the 593,412-pair harness is not
retained. **Fix:** footnote the three scripts after "zero failures", with a
clause that the earlier 593,412-pair run is recorded in key claim 008.

### SF8. Section 5.1 now states its verdict twice (lines 603–604 and 685–688)

The v8 lead-in "are easily checked to fail" delivers the verdict the closing
paragraph restates ("The order does not satisfy any of the classical algebraic
identities ..."). **Fix:** fold the parenthetical list into the lead-in and
reduce the closing paragraph to "The defining property is closure: ...". One
refuter rates this a nit.

### SF9. The generator display spreads each cycle across the text width (§5.3, lines 727–733)

With three alignment pairs per line, `align` places the three maps of each
cycle at the left margin, the center and the right margin, so a line reads as
three column entries rather than one cyclic sequence, and "e₄ ↦ −e₄" sits
under the first column. Content verified correct. **Fix:** one alignment point
per line with `\quad` separators:
`&e_1 \mapsto -e_5, \quad e_5 \mapsto e_6, \quad e_6 \mapsto e_1,` etc.

### SF10. The ambient algebra has a one-line description the paper does not state (§5.3, §7 classification bullet)

By Definition 3.3 the routing sends blocks (α, β) to −(α + β) mod 3 with the
left factor from u, so (ℝ²⁴, ⋆) = (𝕆, ·σ) ⊗ A with A = ℝ[ℤ₃] under the product
a ⋆ b = ι(ab), ι the group inversion; A = ℝf ⊕ A₀ gives D = 𝕆 ⊗ ℝf and
T = 𝕆 ⊗ A₀ ≅ 𝕆_ℂ with product conj(ab), consistent with Aut = G₂ × G₂ × S₃.
Both refuters verified this exactly on all 576 basis pairs. The referee
dimension suggests one sentence at the start of the ambient-algebra paragraph
and re-scoping the Conclusion's classification bullet to the order. New
content, outside the surgical scope of v8; recorded for the author's decision.

### SF11. A one-word widow opens page 15 ("ℤ." before display (A.1)–(A.8))

New in v8, a pagination shift from the insertions. Cosmetic; re-check after
the fixes above, which move text again.

---

## Nits (confirmed)

- §7, line 1074: V is used before any definition in the main text (W is
  defined inline at 1064–1065). "$V := \sigma(L\bar s)/2L$ is an ideal ...".
- §6, line 957: "an order in the sense of Theorem 1.1" points one step
  before the definition; "in the sense of Corollary 1.2".
- Declarations line 1154 "repositories" versus Appendix C line 1777
  "repository" for the same `\cite{repo}`; one form.
- Remark 1.3, line 193, "sampled distinct pairs" versus Remark 5.1 "sampled
  pairs": the script draws u, v independently (three coincident draws in the
  recorded seeds); "distinct" is the odd one out. Touches the prompt-277
  wording only by that word.
- §6, line 969: "the same Hermite-normal-form method as for ⋆ (Section 5.5)",
  but Section 5.5 names no method; either name it there or reword here.

## Refuted

- Theorem 1.1 omits σ: the σ-free statement is a recorded author ruling
  (prompt 175), and Definition 3.3 pins the twist; both refuters refuted.

## Reported, not verified (beyond the verification cap; severity nit)

Presentation: serial comma in the Section 5.2 norms (line 714); the Section 6
semicolon splice "reads ... under our reading" (916); "this overlattice has
index 4" without both terms (1053); "prompt-by-prompt" twice in one Appendix C
paragraph (1707, 1712) and "prompt by prompt" in the Declarations;
"process" / "protocol" / "methodology" for one thing in the Declarations
(1157–1161); the abstract's related-work sentence asserts nothing checkable
(106; requested by the author); two page-break positions (census table
lead-in, Conclusion heading); "left-multiplication" hyphenation (678 versus
238). Reproducibility: the rescaling script finds the least multiple by a
grid scan (705; the ray argument of the record makes it exhaustive); the Dixon
footnote cites library modules without a main entry point (884); the
square-zero count's citation placement (778). Attribution and registry: no
record pointer for the Hermes paragraph (same as SF5); Baez2014 year range
in the registry; the Elduque 2023 venue stale in the registry; the `[repo]`
bibitem unregistered; GAP and LOOPS registered but uncited in the paper. Tone:
§6 line 989 "none defines or tests a bilinear product on Λ" characterizes what
prior work did not do, against the author's own tone standard. The registry
items belong to the documentation sweep planned after the paper is final.

---

## Verdict

**Minor revision.** The assembled v8 holds up: the edits requested at prompt
281 are correct and consistent with the body, and the paper's reproducibility
discipline is intact. Before freezing: MF1 and MF2 are required; SF3's content
point (L³) is the one item that changes a mathematical statement's framing and
should be settled; SF4 is a pre-existing factual error in the Introduction;
SF1, SF2 and SF5 to SF9 are cheap. SF10 is a referee's gift, not a defect.
Nothing has been applied to `main.tex`; dispositions are the author's.
