# v3 → v4 changes (top-level summary)

**v3** — frozen 29 April 2026, 12 pages. Source preserved at [paper/main_2026-04-29.tex](main_2026-04-29.tex).

**v4** — frozen 25 May 2026, 19 pages. Source at [paper/main.tex](main.tex) (working copy) and [paper/main_2026-05-25.tex](main_2026-05-25.tex) (frozen snapshot).

The page-count growth from 12 → 19 is concentrated in three additions: the historical synthesis (new Appendix B), the mod-2-quotient structural reading (new Appendix A.1), and the expanded algebraic-properties test at N = 10⁶ (§5). The rest of the change is editorial refinement.

---

## A. Mathematical content additions

### A.1 Section A.1 — *The mod-2-quotient reading of Tables A.2 and A.3* (new)

The integer-coefficient tables that certify Lemmas 4.3 and 4.4 now have an explicit structural reading. From Wilson's identity $L\bar s \cap Ls = 2L$, the octonion product on $L$ descends to an $\mathbb{F}_2$-bilinear product $\bar\mu$ on $L/2L \cong \mathbb{F}_2^8$, and the two sublattices in play project to a 4 + 4 direct-sum decomposition:
- $V := \sigma(L\bar s)/2L$ is a **two-sided ideal** of $(L/2L, \bar\mu)$ (Lemma 4.3 mod 2L).
- $W := \sigma(Ls)/2L$ is a **subalgebra** of $(L/2L, \bar\mu)$ (Lemma 4.4 mod 2L).

Under the natural plus-type quadratic form $q(v + 2L) := N(v)/2 \pmod 2$ on $L/2L$, both V and W are **maximal totally isotropic subspaces (Lagrangians)**, and $L/2L = V \oplus W$ is a **Witt decomposition into a complementary pair of Lagrangians** — a polarization of the orthogonal $\mathbb{F}_2$-space. The lattice-level picture of L/2L with this plus-type form and the 270 × 64 = 17,280 enumeration of complementary Lagrangian pairs is already in Baez/Egan (2014); what v4 adds is the descended product $\bar\mu$ and the identification of V and W as ideal and subalgebra of $(L/2L, \bar\mu)$.

§7 open question 1 is reframed as a **research program**: which polarizations of $L/2L$ lift to closed bilinear products on $\Lambda$ itself.

### A.2 Appendix B — *Historical note: integral Cayley numbers, 1923–1946* (new)

A 1-page synthesis of the four primary papers, all re-derived in exact arithmetic in this project:
- **Dickson (1923)** — first verified maximal order ($O_D$); his Theorem XV undercounts to three because of an implicit Hurwitz-quaternion containment assumption.
- **Kirmse (1924)** — independent treatment; recovered the full count of seven, but exhibited a non-closed candidate $J_1$ among them.
- **Mahler (1942)** — Euclidean property $N(X - G) \le 1/2$ and principal-ideal theorem in a maximal order.
- **Coxeter (1946)** — identified Kirmse's error and the Bruck transposition that repairs $J_1$; first correct account of all seven maximal orders.

This places the paper's σ as **exactly the Coxeter–Bruck transposition** — the same linear involution of $\mathbb{R}^8$, but acting on Wilson's coordinate-symmetric sublattice $Ls$ rather than on Kirmse's axis-asymmetric $J_1$. The abstract and Remark 2.2 (since removed) both make this identification.

### A.3 §5 — N = 10⁶ algebraic-properties test

§5 in v3 reported short-sample percentages with hand-waved precision. v4 reruns the eight properties on N = 10⁶ samples (per property) with proper standard-error columns. Two new rows:
- **Cube power-associativity** split out from a general "power-associativity" row.
- **Quartic power-associativity** (all five degree-4 parenthesisations equal) — a strictly stronger identity than cube.
- **Symmetric composition (Elduque)** — added as a new row.

Empirically, every quartic-passing sample also satisfies the cube identity (verified in a joint exact-arithmetic test on 3×10⁵ samples), while the cube-passing set is ~2.4× larger. None of the eight properties holds globally; the paper does not overclaim composition or alternativity.

### A.4 Section A.1 closing — *Polarization* and *Structural reframing* paragraphs

The two closing paragraphs of Section A.1 record:
- **V is also a right ideal**, not only a left ideal (the right-ideal direction is not formally implied by the left-ideal direction since $\bar\mu$ is non-commutative in the L-basis; footnote explains why).
- The Lagrangian/Witt picture (V and W both totally isotropic).
- The research-program reframing of §7's open question 1.

(Originating feedback credited in the acknowledgments to "OpenCode DeepSeek, Hermes".)

### A.5 Remark 4.9 — *Consecutive twists* (moved from Theorem 1.1 footnote)

The empirical enumeration of products of two transpositions:
- $\pi = $ identity (the untwisted product) — does not close.
- 70 three-cycles — 35 close; one of two orientations per 3-subset of imaginary indices.
- 105 (2,2)-doubles — 42 close; the 21 whose fixed three indices form a Fano line fail uniformly, the other 84 split half-half.

Verified by exact-arithmetic enumeration on a Z-basis of $\Lambda$. In v3 this was a footnote attached to Theorem 1.1; in v4 it is Remark 4.9 at the end of §4.

---

## B. Structural refinements (without new mathematics)

### B.1 §3 — *The construction* tightened

v3 had Definition 3.1 introduce "the transposition-twisted octonion algebra" with the twist formula folded in, followed by Proposition 3.2 stating and proving the algebra-isomorphism property.

v4 cleanly separates σ from its action on the product:
- **Definition 3.1**: σ is just the transposition (no algebra).
- **Definition 3.2**: $\cdot_\sigma$ is the σ-twisted variant of the standard octonion product.
- A single sentence then states σ is the algebra isomorphism (no proposition — it follows by $\sigma^2 = \mathrm{id}$).
- **Definition 3.3**: the triple product with $\mathbb{Z}_3$ routing (unchanged from v3).

### B.2 §4 — *Proof of closure* reordered and decluttered

In v3, the blockwise expansion of $u \star v$ came first; in v4, the four lemmas come first (R⁸-internal facts about σ, $L$, $Ls$, $L\bar s$), then the blockwise form is introduced where it is needed (to set up the three condition propositions). The intro paragraph was rewritten to identify the load-bearing effect of σ (it *moves* $L\bar s$ and $Ls$) in one declarative sentence. Two remarks were merged into one (`rem:Ls-not-closed`), and one was removed entirely (`rem:two-framings`).

### B.3 §6 — *Related work* historical material collapsed

The full historical synthesis (Dickson 1923 / Kirmse 1924 / Mahler 1942 / Coxeter 1946) was moved into Appendix B. §6 keeps a one-line pointer. The Corradetti 2026 reference is added as a thematically adjacent recent work on E₈ via the Okubo algebra.

### B.4 §1 / §2 / §5 / §7 / §8 — wording

- **Abstract** rewritten to open with the closure + polarization story; the historical paragraph now lists Dickson (first), Kirmse, Bruck (transposition), Coxeter (full account) in correct sequence.
- **§1 intro**: σ is a transposition; preview matches §3.
- **§2.1**: "division algebra" scoped to *unital* (Hurwitz uniqueness); the relaxation (non-unital, e.g. Okubo) is explicitly noted.
- **§2.2**: maximal-order line sharpened — "closed under multiplication *and admits no enlargement*".
- **Remarks 2.1 and 2.2** removed (redundant with §1 / abstract).
- **§5 table**: "Multiplicative identity" row reads "n/a (structural)" (no longer "0%" implying a sampling rate).
- **§7 bullet 1**: parenthetical pointer to Section A.1's polarization-reframing.
- **§8 Outlook**: two subsection headers removed; the two paragraphs flow as a single forward-looking statement.

### B.5 Theorem 1.1, Corollary 1.2

- **Theorem 1.1**: simplified to reference only the triple-product definition; the "in Wilson's representation" qualifier added explicitly.
- **Corollary 1.2**: trimmed of conflated "Z-subring of full rank" language; the non-unital-order disclosure is now a footnote on "order".

---

## C. References, attribution, and editorial

### C.1 Bibliography

Added in v4:
- **Dickson1923** — *A new simple theory of hypercomplex integers*, J. Math. Pures Appl. (9) 2 (1923), 281–326.
- **Mahler1942** — *On ideals in the Cayley–Dickson algebra*, Proc. Roy. Irish Acad. Sect. A 48 (1942), 123–133.
- **Corradetti2026** — *Integral elements of Okubo algebra and the $E_8$-lattice*, arXiv:2605.09333.
- **MarraniCorradettiZucconi2025** — *Physics with non-unital algebras? An invitation to the Okubo algebra*, J. Phys. A: Math. Theor. 58 (2025), 075202.

Petersson 2018 URL added; repo reference now cites both the Bitbucket primary and GitHub mirror.

### C.2 First-name additions

On first body occurrence, 15 authors gained full first names: Robert A. Wilson, Richard H. Bruck, Johannes Kirmse, H. S. M. ("Donald") Coxeter, Geoffrey M. Dixon, John C. Baez, Greg Egan, Susumu Okubo, Alberto Elduque, John H. Conway, Neil J. A. Sloane, James Lepowsky, Arne Meurman, Noam D. Elkies, Benedict H. Gross, Leonard E. Dickson, Kurt Mahler, Daniele Corradetti, Jonathan I. Hall.

### C.3 Wilson citations consolidated

§2 went from 5 Wilson 2009 citations to 1 (at the start of §2.2). Subsequent §2.3 sublattice identities and Λ definition derive Wilson's attribution from the section title "Wilson's Leech lattice".

### C.4 Acknowledgments

Added in v4: Matthew Barley, Petr Vojtěchovský. A second sentence credits **OpenCode DeepSeek, Hermes** for the feedback that prompted the Lagrangian-decomposition framing in Section A.1.

### C.5 Editorial standards

Adopted during v4 development (`evidence_and_reasoning/editorial_standards.md`): five standards governing prose changes to durable artifacts. Most visibly: en-dash (`--`) only, never em-dash (`---`). Other standards (lead programmatically, introduce one thing at a time, cut every redundancy, one consistent vocabulary, precision over good-sounding vagueness) were applied during v4 assembly.

### C.6 Research methodology (Appendix C)

Kept and lightly updated. Records the human–AI collaboration model, the prompt-log audit trail, the cross-checking discipline, and is candid about mistakes made on both sides.

---

## D. Verification artefacts

All claims in v4 have been independently verified in exact arithmetic. Key scripts:
- `python_project/src/symbolic_proof_checks.py` — Lemmas 4.1–4.4.
- `gap_project/tests/test_lemmas.g` — same lemmas, independent CAS.
- `python_project/src/verify_section5_properties.py` — §5 N = 10⁶ run.
- `python_project/src/verify_mod2_quotient.py` — Section A.1 V/W ideal/subalgebra.
- `python_project/src/verify_discovery1_W_subalgebra.py` — W multiplication table.
- `python_project/src/verify_discovery2_V_isotropic.py` — V totally isotropic, two-sided.
- `python_project/src/verify_consecutive_twists_exact.py` — Remark 4.9 enumeration.
- `python_project/src/verify_all_cycles_exact.py` — extends to pure n-cycles n ∈ {3,4,5,6,7}.
- `python_project/src/verify_{dickson_1923, kirmse_1924, mahler_1942, coxeter_1946}.py` — historical re-derivations.

Total: 27 bibliography entries, all cross-checked against the project's central registry in `evidence_and_reasoning/references/` (paper/reviews/2026-05-24_claude_opus_4_7_reference_check.md).

---

## E. Companion artefacts updated for v4

- [paper/companion.tex](companion.tex) — re-versioned as v4 companion (Claude Opus 4.7, 23 May 2026); the σ ≡ Kirmse-twist framing aligned.
- [paper/historical_appendix.tex](historical_appendix.tex), [paper/historical_appendix_expanded.tex](historical_appendix_expanded.tex) — standalone short and expanded forms of Appendix B.
- [paper/hermes_mod2_structure_supplement.tex](hermes_mod2_structure_supplement.tex) — supplement on the W subalgebra structure (Hermes Discovery 1 and 2 with verification and corrections).
- [paper/Baez-Egan_Leech_Verification_addendum_2026-05-24.tex](Baez-Egan_Leech_Verification_addendum_2026-05-24.tex) — addendum to the original 16 April 2026 Baez–Egan verification PDF, explicitly recording the L/2L parallel.

---

## F. Page-count attribution

| Region | v3 | v4 | Δ |
|---|---|---|---|
| Front matter + §1–§3 + Theorem 1.1 | ~3 pp | ~3 pp | +0 |
| §4 proof | ~3 pp | ~4 pp | +1 |
| §5 algebraic properties | ~1 pp | ~1.5 pp | +0.5 |
| §6 related work | ~2 pp | ~2.5 pp | +0.5 (Corradetti + collapsed history) |
| §7 conclusion + §8 outlook | ~1 pp | ~1 pp | +0 |
| Appendix A (basis tables + mod-2 reading) | ~2 pp | ~4 pp | +2 (new Section A.1) |
| Appendix B (historical) | — | ~1 pp | +1 (entirely new) |
| Appendix C (methodology) | ~1 pp | ~1 pp | +0 |
| Bibliography | ~1 pp | ~1 pp | +0 |
| **Total** | **~12 pp** | **~19 pp** | **+7 pp** |

---

## G. v4 → v3 cross-references

- Frozen v3 source: [paper/main_2026-04-29.tex](main_2026-04-29.tex).
- Frozen v4 source: [paper/main_2026-05-25.tex](main_2026-05-25.tex) (this freeze).
- Working source (continues in this file going forward): [paper/main.tex](main.tex) — currently identical to the v4 freeze.

— Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger, 2026-05-25
