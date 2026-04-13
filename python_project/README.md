# Python Project

Self-contained Python project for the computational work of this repository: lattice constructions, octonion / Okubo algebras, trial experiments testing candidate orders on the Leech lattice Λ, exact-arithmetic verification of the symbolic-proof lemmas, and the foundation test suite.

## Structure

```
python_project/
├── src/                # Source code (see src/README.md)
├── tests/              # Test suite (see tests/README.md)
├── conftest.py         # pytest configuration (puts src/ on sys.path)
├── requirements.txt    # numpy, pytest
└── README.md           # This file
```

## Setup

Requires Python 3.x.  Install dependencies:

```bash
pip install -r requirements.txt
```

The only external dependencies are `numpy` (used for fast array arithmetic in the trial harnesses) and `pytest` (test runner).  All other imports are standard-library (`fractions`, `itertools`, `collections`, `typing`, `multiprocessing`, `argparse`).

## Running the tests

From this directory:

```bash
python3 -m pytest tests/ -v
```

The 197-test suite verifies every foundation (octonion algebra, Wilson's E8 and Leech constructions, Dixon's E8 and Leech constructions, Okubo algebra) before any trial relies on it.

## Running individual scripts

Each trial and verification script is self-contained and runnable directly:

```bash
cd src
python3 symbolic_proof_checks.py        # Verify the five lemmas behind the paper
python3 trial_007_kirmse_twist.py       # Initial test of the transposition-twisted product
python3 trial_007_fast.py               # 4M-pair random closure check
python3 trial_007_exhaust.py            # Multiprocessor harness for full 38.6B-pair sweep
python3 trial_001_triple_octonion.py    # Any earlier trial; see src/README.md for the catalogue
python3 consistency_checks.py           # Pre-paper consistency sweep (10 checks)
```

## Design Principles

- Self-contained and runnable without internet access after `pip install`.
- All computations are deterministic and reproducible (fixed seeds wherever randomness is used).
- The symbolic-proof lemmas (`symbolic_proof_checks.py`) use exact rational arithmetic via `fractions.Fraction` — no floating-point.
- Test coverage is sufficient to validate every mathematical claim used in the paper.
- No personal or sensitive information is stored in this project.
