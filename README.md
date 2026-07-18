# Exploring an Order on the Leech Lattice

This repository documents the exploration of a specific mathematical hypothesis:
whether the Leech lattice Λ admits an **order** structure — a bilinear
multiplication making (Λ, +, ·) an order in a 24-dimensional real algebra.

The research is conducted with AI assistance under a strict integrity protocol
([MANIFESTO.md](MANIFESTO.md)).  Every step is logged, every claim is backed by
evidence, and every conclusion is traceable.

**Start here: [CURRENT_STATE.md](CURRENT_STATE.md)** — the entry point for
anyone continuing this research, with or without prior context.

The write-up is in [paper/main.tex](paper/main.tex) (compiled:
[paper/main.pdf](paper/main.pdf)), at v6, frozen 18 July 2026 as
[paper/main_2026-07-18.tex](paper/main_2026-07-18.tex).  [paper/](paper/)
also holds the frozen versions with their changelogs, the
companion material, and the review record.

## Repository Structure

| Path | Purpose |
|---|---|
| [CURRENT_STATE.md](CURRENT_STATE.md) | **Entry point** — what's established, what's ruled out, what's next |
| [MANIFESTO.md](MANIFESTO.md) | Operating rules for AI-assisted research |
| [DOCUMENT_GENRES.md](DOCUMENT_GENRES.md) | What each artifact is for, and how it may be changed: immutable, frozen, or current state |
| [tools/](tools/) | Project scaffolding; `lint_docs.py` enforces the document genres mechanically |
| [TRIAL_METHODOLOGY.md](TRIAL_METHODOLOGY.md) | Structure and philosophy for trial files |
| [paper/](paper/) | Formal write-up (`main.tex`, `main.pdf`), earlier frozen versions and changelogs, companion material, and the review record |
| [evidence_and_reasoning/](evidence_and_reasoning/) | Key claims, trial results, and references |
| [python_project/](python_project/) | Python code: shared tools, trial experiments, symbolic-proof verification |
| [gap_project/](gap_project/) | GAP / LOOPS independent re-derivation of the paper's verification tests |
| [prompt_logs/](prompt_logs/) | Chronological log of all AI interaction prompts |
| [source_documents/](source_documents/) | Primary source PDFs (freely redistributable only; others listed by DOI) |
| [LICENSE.md](LICENSE.md), [LICENSE-CODE](LICENSE-CODE) | CC BY 4.0 for text / evidence / paper; MIT for source code |

## Current Status

**Research goal:** Find a highly symmetric order on the Leech lattice.

**Result (trial 007):** The Leech lattice Λ admits an order under a
Z₃-symmetric triple-octonion product on R²⁴ — the same octonion product in
all three blocks, with Z₃-symmetric cross-block routing.  In Wilson's
representation of Λ, the fit of this octonion product with the
representation is exhibited by a transposition σ of two imaginary basis
units: σ leaves the E₈ lattice invariant and is an octonion-algebra
isomorphism, but it moves Wilson's sublattices Ls and Ls̄.

**Evidence:**
- **Symbolic proof** of closure on Wilson's three sublattice conditions, via
  four lemmas on the interplay between the transposition σ and the Leech
  sublattices Ls, Ls̄ ([python_project/src/symbolic_proof_checks.py](python_project/src/symbolic_proof_checks.py)
  executes the finite case verifications with exact rational arithmetic).
- **Computational verification** on 12M+ random pairs of minimal vectors
  with zero failures, plus tests from first principles on every foundation
  (octonion properties, Wilson's construction).  197 tests pass.
- **Independent GAP / LOOPS re-derivation** of the same verification tests
  ([gap_project/](gap_project/), 110 checks; uses the LOOPS package to
  cross-verify the octonion multiplication as a nonassociative Moufang
  loop of order 16, isomorphic to the standard octonion loop
  MoufangLoop(16, 3)).  All arithmetic exact (rational); no floating point.
- **Structure of the algebra**, exact and exhaustive (Section 5 of the
  paper): the ambient algebra splits as R²⁴ = D ⊕ T into mutually
  annihilating two-sided ideals (D the diagonal octonions, T the sum-zero
  triples); it has exactly eight idempotents, of which only 0 lies in Λ; Λ
  contains square-zero vectors of norms in 12Z, exactly 4,032 of them of
  norm 12; and Aut(Λ, +, ⋆) ≅ C₆ × S₃ is
  finite of order 36, with Aut(Λ, +, ⋆) ⊊ Co₀
  ([python_project/src/verify_aut_lambda_star.py](python_project/src/verify_aut_lambda_star.py)
  and [gap_project/aut_lambda_star.g](gap_project/aut_lambda_star.g)).
- **Formal write-up** with full proof, related work, and methodology in
  [paper/main.tex](paper/main.tex).
- **Key claims:** [007](evidence_and_reasoning/key_claims/007_triple_octonion_ruled_out.txt)
  (the triple product over Wilson's own multiplication convention does
  not close on Λ) and
  [008](evidence_and_reasoning/key_claims/008_transposition_twist_order.txt)
  (over the σ-related convention it does: (Λ, +, ⋆) is an order).

**Open questions:**
- A structural reason for Lemma 4.4 (σ(Ls) is closed under the standard
  octonion product, while Ls is not).  One candidate shape is excluded:
  σ(Ls) is not an ideal of L on either side.
- Maximality of the order: is there a lattice Λ' ⊋ Λ of finite index with
  Λ' ⋆ Λ' ⊆ Λ'?
- Where Aut(Λ, +, ⋆) ≅ C₆ × S₃ sits inside Co₀ = Aut(Λ) up to conjugacy,
  and whether the order-6 octonion automorphism generating the C₆ has an
  intrinsic description.
- A structural characterization of which permutation cycles close.  All
  5,040 permutations of S₇ have been enumerated exactly (1,764 of the
  5,039 non-identity permutations close), but the pattern is described
  only for the 3-cycles.
- Whether a **ternary** reformulation (via composition algebras in the
  sense of Elduque, or via Okubo's ternary structure) gives a more
  natural classification — the rigid Z₃ cross-block routing and the
  order-3 automorphism content of the construction both point that
  way.  This is the principal direction suggested for future work.

## Running the Code

### Python suite

```bash
cd python_project
python3 -m pytest tests/ -v                     # Run all foundation tests
cd src
python3 symbolic_proof_checks.py                # Verify the four lemmas
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
gap -q gap_project/aut_lambda_star.g            # Automorphism group, second pass
gap -q gap_project/octonion_stabilisers.g       # Stabilizers of L, Ls, Ls̄
```

Requires GAP 4.x with the LOOPS package
(`apt install gap gap-pkg-loops` on Debian/Ubuntu, or via GAP's PackageManager).
The GAP suite re-derives the same verification claims as the Python suite,
with all arithmetic exact (rational), and uses LOOPS to cross-verify the
octonion multiplication as a Moufang loop of order 16.  See
[gap_project/README.md](gap_project/README.md) for the full layout.

### Compiling the paper

```bash
cd paper
pdflatex main.tex
pdflatex companion.tex
pdflatex automorphism_group_2026-07-12.tex
```
