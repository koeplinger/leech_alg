# Mathematical Hypothesis Exploration

This repository documents the exploration of a mathematical hypothesis: whether
the Leech lattice Λ admits an **order** structure — a bilinear multiplication
making (Λ, +, ·) an order in a 24-dimensional real algebra.

The research is conducted with AI assistance under a strict integrity protocol
([MANIFESTO.md](MANIFESTO.md)).  Every step is logged, every claim is backed by
evidence, and every conclusion is traceable.

**Start here: [CURRENT_STATE.md](CURRENT_STATE.md)** — the entry point for
anyone continuing this research, with or without prior context.

## Repository Structure

| Path | Purpose |
|---|---|
| [CURRENT_STATE.md](CURRENT_STATE.md) | **Entry point** — what's established, what's ruled out, what's next |
| [MANIFESTO.md](MANIFESTO.md) | Operating rules for AI-assisted research |
| [TRIAL_METHODOLOGY.md](TRIAL_METHODOLOGY.md) | Structure and philosophy for trial files |
| [evidence_and_reasoning/](evidence_and_reasoning/) | Key claims, trial results, and references |
| [python_project/](python_project/) | Python code: shared tools and trial experiments |
| [prompt_logs/](prompt_logs/) | Chronological log of all AI interaction prompts |
| [source_documents/](source_documents/) | Primary source PDFs |
| [research_output/](research_output/) | Reserved for future experimental output |

## Current Status

**Research goal:** Find a highly symmetric order on the Leech lattice.

**Progress:**
- Foundations established: octonion algebra, E8 lattice (Wilson + Dixon),
  Leech lattice (Wilson + Dixon), Okubo algebra.  157 tests pass.
- **Triple-octonion algebra (O₁ ⊕ O₂ ⊕ O₃): RULED OUT** (4 trials,
  1,536+ variants tested, structural obstruction identified).
  See [key claim 007](evidence_and_reasoning/key_claims/007_triple_octonion_ruled_out.txt).
- Next candidates: Okubo-based products, non-block-decomposed products,
  Conway-group-informed constructions.

## Running the Code

```bash
cd python_project
python3 -m pytest tests/ -v          # Run all 157 foundation tests
cd src
python3 trial_001_triple_octonion.py  # Run any trial
```

Requires Python 3.x with NumPy and SciPy.
