# Exploring an Order on the Leech Lattice

This repository documents the exploration of a specific mathematical hypothesis:
whether the Leech lattice Λ admits an **order** structure — a bilinear
multiplication making (Λ, +, ·) an order in a 24-dimensional real algebra.

The research is conducted with AI assistance under a strict integrity protocol
([MANIFESTO.md](MANIFESTO.md)).  Every step is logged, every claim is backed by
evidence, and every conclusion is traceable.

**Start here: [CURRENT_STATE.md](CURRENT_STATE.md)** — the entry point for
anyone continuing this research, with or without prior context.

The finished write-up is in [paper/main.tex](paper/main.tex)
(compiled: [paper/main.pdf](paper/main.pdf)), now at v6
(10 July 2026; source snapshot: [paper/main_2026-07-10.tex](paper/main_2026-07-10.tex)).
Revision chain: [v3 → v4](paper/v3_to_v4_summary.md),
[v4 → v5](paper/v4_to_v5_changelog.md),
[v5 → v6](paper/v5_to_v6_changelog.md).  Prior freezes:
v4 at [paper/main_2026-05-25.tex](paper/main_2026-05-25.tex),
v5 at [paper/main_2026-06-07.tex](paper/main_2026-06-07.tex).

## Repository Structure

| Path | Purpose |
|---|---|
| [CURRENT_STATE.md](CURRENT_STATE.md) | **Entry point** — what's established, what's ruled out, what's next |
| [MANIFESTO.md](MANIFESTO.md) | Operating rules for AI-assisted research |
| [TRIAL_METHODOLOGY.md](TRIAL_METHODOLOGY.md) | Structure and philosophy for trial files |
| [paper/main.tex](paper/main.tex), [paper/main.pdf](paper/main.pdf) | Formal write-up: LaTeX source (bibliography inline via \bibitem) and compiled PDF |
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
  executes the finite case verifications with exact integer arithmetic).
- **Computational verification** on 12M+ random pairs of minimal vectors
  with zero failures, plus tests from first principles on every foundation
  (octonion properties, Wilson's construction).  197 tests pass.
- **Independent GAP / LOOPS re-derivation** of the same verification tests
  ([gap_project/](gap_project/), 110 checks; uses the LOOPS package to
  cross-verify the octonion multiplication as the unique non-associative
  Moufang loop of order 16).  All arithmetic exact (rational); no floating point.
- **Formal write-up** with full proof, related work, and methodology in
  [paper/main.tex](paper/main.tex).
- **Key claims:** [007](evidence_and_reasoning/key_claims/007_triple_octonion_ruled_out.txt)
  (the triple product over Wilson's own multiplication convention does
  not close on Λ) and
  [008](evidence_and_reasoning/key_claims/008_transposition_twist_order.txt)
  (over the σ-related convention it does: (Λ, +, ⋆) is an order).

**Open questions:**
- Algebraic properties as a binary product have been characterised on
  Min(Λ): the order is non-unital, non-commutative, not
  norm-multiplicative, and fails alternativity, flexibility, and
  power-associativity (see Section 5 of the paper and
  [CURRENT_STATE.md](CURRENT_STATE.md)).  Open: whether any of these
  negative findings admits a tighter structural statement on Λ
  itself rather than on random samples of Min(Λ).
- Maximality of the order.
- Relationship to the Conway group Co₀ = Aut(Λ).
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
```

Requires Python 3.x with NumPy and SciPy.

### GAP / LOOPS suite

```bash
gap -q gap_project/run_all.g                    # 110 checks across 8 files
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
```
