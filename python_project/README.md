# Python Project

Self-contained Python project for the computational work of this repository: lattice constructions, octonion / Okubo algebras, trial experiments testing candidate orders on the Leech lattice Λ, exact-arithmetic verification of the symbolic-proof lemmas and of the algebraic structure of the resulting order (ideals, idempotents, square-zero elements, automorphism group), and the foundation test suite.

## Structure

```
python_project/
├── src/                # Source code (see src/README.md)
├── tests/              # Test suite (see tests/README.md)
├── conftest.py         # pytest configuration (puts src/ on sys.path)
├── requirements.txt    # numpy, sympy, pytest
└── README.md           # This file
```

## Setup

Requires Python 3.x.  Install dependencies:

```bash
pip install -r requirements.txt
```

The external dependencies are `numpy` (fast array arithmetic in the trial harnesses), `sympy` (exact linear algebra over Q and Q(ζ₃) in the structure and automorphism verification scripts), and `pytest` (test runner).  All other imports are standard-library (`fractions`, `itertools`, `collections`, `typing`, `multiprocessing`, `argparse`).

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
python3 symbolic_proof_checks.py        # Verify the four lemmas behind the paper
python3 trial_007_kirmse_twist.py       # The winning trial: triple-octonion closure on Λ
python3 trial_007_fast.py               # 4M-pair random closure check
python3 trial_007_exhaust.py            # Multiprocessor harness for full 38.6B-pair sweep
python3 trial_001_triple_octonion.py    # Any earlier trial; see src/README.md for the catalogue
python3 consistency_checks.py           # Pre-paper consistency sweep (checks 1-3, 5-10)
```

The algebraic properties reported in Section 5 of the paper are certified by their own scripts, each stating what it proves:

```bash
cd src
python3 verify_ideal_decomposition.py          # R^24 = D + T, two annihilating two-sided ideals
python3 verify_idempotent_classification.py    # exactly eight idempotents; no identity element
python3 verify_idempotent_lattice_rescaling.py # no nonzero idempotent lies in Λ; least lattice multiples
python3 verify_square_zero_classification.py   # the square-zero cone; norm-12 stratum of Λ is 4,032
python3 verify_star_algebra_structure.py       # Aut(R^24, ⋆) = G_2 x G_2 x S_3, dimension 28
python3 verify_aut_lambda_star.py              # Aut(Λ, +, ⋆): order 36, C_6 x S_3; inside Co_0
python3 verify_aut_octonion_crosscheck.py      # independent re-derivation of the stabilizers
python3 verify_sigma_Ls_ideal_exclusion.py     # σ(Ls) is not an ideal of L (excludes a candidate)
python3 verify_all_permutations_exact.py       # closure census over all of S₇ (1,764 of 5,039)
```

See [`src/README.md`](src/README.md) for what each one certifies, and [`../gap_project/README.md`](../gap_project/README.md) for the independent GAP determination of the automorphism group's order and structure.

## Design Principles

- Self-contained and runnable without internet access after `pip install`.
- All computations are deterministic and reproducible (fixed seeds wherever randomness is used).
- The symbolic-proof lemmas (`symbolic_proof_checks.py`) use exact rational arithmetic via `fractions.Fraction`, with no floating-point.
- The pytest suite covers the foundations (octonions, E8, Leech, Okubo) that everything else is built on.  Each claim the paper makes on top of them names the script that reproduces it; those scripts are run directly, not through pytest.
- No personal or sensitive information is stored in this project.
