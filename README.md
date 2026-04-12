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
(compiled: [paper/main.pdf](paper/main.pdf)).

## Repository Structure

| Path | Purpose |
|---|---|
| [CURRENT_STATE.md](CURRENT_STATE.md) | **Entry point** — what's established, what's ruled out, what's next |
| [MANIFESTO.md](MANIFESTO.md) | Operating rules for AI-assisted research |
| [TRIAL_METHODOLOGY.md](TRIAL_METHODOLOGY.md) | Structure and philosophy for trial files |
| [paper/main.tex](paper/main.tex), [paper/main.pdf](paper/main.pdf), [paper/references.bib](paper/references.bib) | Formal write-up: LaTeX source, compiled PDF, BibTeX database |
| [evidence_and_reasoning/](evidence_and_reasoning/) | Key claims, trial results, and references |
| [python_project/](python_project/) | Python code: shared tools, trial experiments, symbolic-proof verification |
| [prompt_logs/](prompt_logs/) | Chronological log of all AI interaction prompts |
| [source_documents/](source_documents/) | Primary source PDFs (freely redistributable only; others listed by DOI) |
| [LICENSE.md](LICENSE.md), [LICENSE-CODE](LICENSE-CODE) | CC BY 4.0 for text / evidence / paper; MIT for source code |

## Current Status

**Research goal:** Find a highly symmetric order on the Leech lattice.

**Result (trial 007):** The Leech lattice Λ admits an order under a bilinear
product built from three identical copies of a *transposition-twisted* octonion
algebra with Z₃-symmetric cross-block routing.  The twist is a single
transposition of two imaginary basis elements of the standard Fano-plane
multiplication.

**Evidence:**
- **Symbolic proof** of closure on Wilson's three sublattice conditions, via
  five lemmas on the interplay between the twist σ and the Leech sublattices
  Ls, Ls̄ ([python_project/src/symbolic_proof_checks.py](python_project/src/symbolic_proof_checks.py)
  executes the finite case verifications with exact integer arithmetic).
- **Computational verification** on 4M+ random pairs of minimal vectors with
  zero failures, plus tests from first principles on every foundation
  (octonion properties, Wilson's construction).  157+ tests pass.
- **Formal write-up** with full proof, related work, and methodology in
  [paper/main.tex](paper/main.tex).
- **Key claims:** [007](evidence_and_reasoning/key_claims/007_triple_octonion_ruled_out.txt)
  (untwisted triple-octonion ruled out) and
  [008](evidence_and_reasoning/key_claims/008_transposition_twist_order.txt)
  (the transposition-twisted product is an order on Λ).

**Open questions:**
- Algebraic properties (identity, alternativity, norm multiplicativity).
- Maximality of the order.
- Relationship to the Conway group Co₀ = Aut(Λ).
- Whether a ternary reformulation (via composition algebras in the sense of
  Elduque) gives a more natural classification.

## Running the Code

```bash
cd python_project
python3 -m pytest tests/ -v                     # Run all foundation tests
cd src
python3 symbolic_proof_checks.py                # Verify the five lemmas
python3 trial_007_kirmse_twist.py               # The winning trial
python3 trial_001_triple_octonion.py            # Any earlier trial
```

Requires Python 3.x with NumPy and SciPy.

Compiling the paper:

```bash
cd paper
pdflatex main.tex
```
