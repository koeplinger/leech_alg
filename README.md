# An Order on the Leech Lattice

Created: 2026-04-06.  Updated: 2026-10-31.

The Leech lattice Λ admits an **order**: it is closed under a bilinear product
on R²⁴, a Z₃-symmetric triple-octonion product assembled from nine
octonion-multiplication blocks under cyclic Z₃ permutation of three octonion
factors.  In Robert A. Wilson's sublattice description of Λ, closure turns on a
transposition σ of two imaginary basis units: σ leaves the E₈ lattice invariant
and is an octonion-algebra isomorphism, but it moves Wilson's sublattices Ls̄
and Ls.

**The result is published:**

> J. Köplinger, *An order on the Leech lattice from a Z₃-symmetric
> triple-octonion product*,
> [doi:10.13140/RG.2.2.22093.19686](https://doi.org/10.13140/RG.2.2.22093.19686)
> (ResearchGate).

**A further revision of the paper (v8) is in progress** in
[paper/main.tex](paper/main.tex).

The published version is
[paper/main_2026-07-19.tex](paper/main_2026-07-19.tex) (v7, 24 pages; compiled:
[paper/main_2026-07-19.pdf](paper/main_2026-07-19.pdf)); the working source is
[paper/main.tex](paper/main.tex) (v8, in progress; compiled:
[paper/main.pdf](paper/main.pdf)).  [paper/](paper/) also holds the earlier
released versions with their changelogs, the companion material, and the review
record.

This repository is the complete record of how the result was obtained.  The
research was conducted with AI assistance under a strict integrity protocol
([MANIFESTO.md](MANIFESTO.md)): every prompt is logged, every claim is backed
by evidence, and every conclusion is traceable.  The construction was proposed
in [prompt 025](prompt_logs/025_triple_kirmse_twisted_octonions.txt)
(11 April 2026); the prompts that follow trace how the result was checked and
verified.

## Repository Structure

| Path | Purpose |
|---|---|
| [CURRENT_STATE.md](CURRENT_STATE.md) | Technical summary: what was established, what was ruled out, what remains open |
| [MANIFESTO.md](MANIFESTO.md) | Operating rules for AI-assisted research |
| [DOCUMENT_GENRES.md](DOCUMENT_GENRES.md) | What each artifact is for, and how it may be changed: immutable, public record, or current state |
| [tools/](tools/) | Project scaffolding; `lint_docs.py` enforces the document genres mechanically |
| [TRIAL_METHODOLOGY.md](TRIAL_METHODOLOGY.md) | Structure and philosophy for trial files |
| [paper/](paper/) | Formal write-up (`main.tex`, `main.pdf`), released versions and changelogs, companion material, and the review record |
| [evidence_and_reasoning/](evidence_and_reasoning/) | Key claims, trial results, and references |
| [python_project/](python_project/) | Python code: shared tools, trial experiments, symbolic-proof verification |
| [gap_project/](gap_project/) | GAP / LOOPS independent re-derivation of the paper's verification tests |
| [prompt_logs/](prompt_logs/) | Chronological log of all AI interaction prompts |
| [source_documents/](source_documents/) | Primary source PDFs (freely redistributable only; others listed by DOI) |
| [LICENSE.md](LICENSE.md), [LICENSE-CODE](LICENSE-CODE) | CC BY 4.0 for text / evidence / paper; MIT for source code |

## The Result and Its Evidence

**Result (trial 007):** The Leech lattice Λ admits an order under a
Z₃-symmetric triple-octonion product on R²⁴: the same octonion product in
all three blocks, with Z₃-symmetric cross-block routing.  The transposition σ
exhibits the fit of this octonion product with Wilson's representation.

- **Symbolic proof** of closure on Wilson's three sublattice conditions, via
  four lemmas on the interplay between the transposition σ and Wilson's
  sublattices Ls̄ and Ls of L (Section 4 of the paper).  Lemma 4.1, σ(L) = L,
  is a parity argument; Lemmas 4.2 to 4.4 follow from the closure of one
  maximal order of integral octonions, 𝒪 = ½(e₁ + e₂)L, through the lattice
  identities σ(Ls̄) = 2𝒪 and σ(Ls) = p𝒪 (p = e₂ − e₄ + e₅ − e₆) and the
  identity (xa)(xb) = x((ax)b) for purely imaginary x
  ([python_project/src/verify_structural_proof.py](python_project/src/verify_structural_proof.py)
  checks each of these facts in exact arithmetic).  The lemma inclusions are
  also checked on Z-bases: Lemmas 4.1, 4.3 and 4.4 in
  [python_project/src/symbolic_proof_checks.py](python_project/src/symbolic_proof_checks.py),
  Lemmas 4.2 to 4.4 in `verify_structural_proof.py`, and all four in
  [gap_project/tests/test_lemmas.g](gap_project/tests/test_lemmas.g).
- **Computational verification** on 12M+ random pairs of minimal vectors
  with zero failures, plus tests from first principles on every foundation
  (octonion properties, Wilson's construction).  198 tests pass.
- **Independent GAP / LOOPS re-derivation** of the same verification tests
  ([gap_project/](gap_project/), 110 checks; uses the LOOPS package to
  cross-verify the octonion multiplication as a nonassociative Moufang
  loop of order 16, isomorphic to the standard octonion loop
  MoufangLoop(16, 3)).  All arithmetic exact (rational); no floating point.
- **Structure of the algebra**, exact and exhaustive (Section 5 and Appendix E
  of the paper): the ambient algebra splits as R²⁴ = D ⊕ T into mutually
  annihilating two-sided ideals (D the diagonal octonions, T the sum-zero
  triples); it has exactly eight idempotents, of which only 0 lies in Λ;
  the order is not maximal (Section 5.3: L³ and Λ + Zε₁ are ⋆-closed
  overlattices, of index 2¹² and 4); the span index of the products depends
  on the twist (Appendix E: 2¹⁶ for σ, proved by hand in Proposition E.1,
  2¹⁰ for the even closing permutations, 1 for a reflection twist); Λ
  contains square-zero vectors of norms in 12Z, exactly 4,032 of them of
  norm 12; and Aut(Λ, +, ⋆) ≅ C₆ × S₃ is
  finite of order 36, with Aut(Λ, +, ⋆) ⊊ Co₀
  ([python_project/src/verify_aut_lambda_star.py](python_project/src/verify_aut_lambda_star.py)
  and [gap_project/aut_lambda_star.g](gap_project/aut_lambda_star.g)).
- **Which twists close** (Section 5.6, Proposition 5.3, named in the
  abstract and the conclusion): for every invertible linear map g of R⁸ with
  g(L) = L, the triple product built from the g-twisted octonion product
  closes on Λ exactly when Wilson's sublattices Ls̄ and Ls are closed under
  that product.  Appendix A reads this modulo 2: every polarization of L/2L
  is (g(Ls̄)/2L, g(Ls)/2L) for an isometry g of L, and the product closes
  exactly when both halves are subalgebras.  Section 6 applies it to Egan's
  construction of Λ from E₈ alone, of which Wilson's Λ and σ(Λ) are two
  instances.
- **Key claims:** [007](evidence_and_reasoning/key_claims/007_triple_octonion_ruled_out.txt)
  (the triple product over Wilson's own multiplication convention does
  not close on Λ) and
  [008](evidence_and_reasoning/key_claims/008_transposition_twist_order.txt)
  (over the σ-related convention it does: (Λ, +, ⋆) is an order).

## What Remains Open

These are the questions the paper leaves open (Sections 7 and 8), with further
detail in [CURRENT_STATE.md](CURRENT_STATE.md).

- Where Aut(Λ, +, ⋆) ≅ C₆ × S₃ sits inside Co₀ = Aut(Λ) up to conjugacy,
  and whether the order-6 octonion automorphism generating the C₆ has an
  intrinsic description.
- A reason for the closure pattern of the census of Appendix E.  All
  5,040 permutations of S₇ have been enumerated exactly (1,764 of the
  5,039 nonidentity permutations close), and Proposition 5.3 decides
  closure for each; the 2-cycles close by Theorem 1.1, and the pattern is
  described for the 3-cycles.  Open: why every 4-cycle fails and why
  several cycle types close at exactly one half.
- Whether a **ternary** reformulation (via composition algebras in the
  sense of Elduque, or via Okubo's ternary structure) gives a more
  natural classification: the rigid Z₃ cross-block routing and the
  order-3 automorphism content of the construction both point that
  way.  This is the principal direction suggested for future work.

## Running the Code

### Python suite

```bash
cd python_project
python3 -m pytest tests/ -v                     # Run all foundation tests
cd src
python3 symbolic_proof_checks.py                # Lemmas 4.1, 4.3, 4.4
python3 verify_structural_proof.py              # Lemmas 4.2 to 4.4: the maximal order, the lattice identities, the Moufang identity
python3 verify_appendix_tables_v7.py            # the Appendix A tables of v7, recomputed
python3 trial_007_kirmse_twist.py               # The winning trial
python3 trial_001_triple_octonion.py            # Any earlier trial
python3 verify_aut_lambda_star.py               # Aut(Λ, +, ⋆) = C₆ × S₃, order 36
python3 verify_idempotent_classification.py     # The eight idempotents
python3 verify_square_zero_classification.py    # Square-zero classification (norm-12 stratum: 4,032)
python3 verify_all_permutations_exact.py        # Closure census over all of S₇
```

Requires Python 3.x with NumPy and SymPy (see [python_project/requirements.txt](python_project/requirements.txt)).

### GAP / LOOPS suite

```bash
gap -q gap_project/run_all.g                    # 110 checks across 8 files
(cd gap_project && gap -q -b aut_lambda_star.g < /dev/null)   # Automorphism group, second pass
(cd gap_project && gap -q -b octonion_stabilisers.g < /dev/null)   # Stabilizers of L, Ls, Ls̄
```

Requires GAP 4.x (`apt install gap` on Debian/Ubuntu) with the LOOPS package,
which the Debian/Ubuntu archive does not carry: unpack its release archive into
`~/.gap/pkg/`, or run `LoadPackage("PackageManager"); InstallPackage("loops");`
in a GAP that has the PackageManager package.
The GAP suite re-derives the same verification claims as the Python suite,
with all arithmetic exact (rational), and uses LOOPS to cross-verify the
octonion multiplication as a Moufang loop of order 16.  See
[gap_project/README.md](gap_project/README.md) for the full layout.

### Compiling the paper

```bash
cd paper
pdflatex main.tex && pdflatex main.tex && pdflatex main.tex   # three times: refs, citations, PDF outline
pdflatex companion.tex && pdflatex companion.tex
pdflatex automorphism_group_2026-07-12.tex && pdflatex automorphism_group_2026-07-12.tex
```

The bibliography is inline via `\bibitem`, so no BibTeX pass is needed; the
second `pdflatex` run resolves cross-references and citations, and `main.tex`
needs a third for its PDF outline (the two notes settle in two runs).
