# Source Code

This folder contains the Python source code for mathematical computations,
lattice constructions, algebra implementations, and trial experiments.

## Shared tools (foundations)

These modules are verified by 157 tests in `../tests/` and used by all trials.

| Module | Purpose |
|---|---|
| `octonions.py` | Octonion algebra: standard product and Dixon X-product |
| `e8_wilson.py` | Wilson's E8 lattice L (240 roots, membership test, structural element s) |
| `e8_dixon.py` | Dixon's E8 constructions (Ξ^even, A^odd) |
| `leech_wilson.py` | Wilson's Leech lattice Λ (196,560 minimal vectors, membership test) |
| `leech_dixon.py` | Dixon's Leech lattice construction |
| `okubo.py` | Okubo algebra (Petersson construction via order-3 automorphism τ) |
| `okubo_samples.py` | Instructive examples of Petersson τ constructions |

## Trial files

Each trial is a self-contained Python script testing a specific algebra against
the Leech lattice.  Run any trial directly:

```
cd python_project/src && python3 trial_NNN_*.py
```

| Trial | File | What it tests |
|---|---|---|
| 001 | `trial_001_triple_octonion.py` | Base triple-octonion algebra |
| 002 | `trial_002_scaled_triple_octonion.py` | Per-block scaling search |
| 003 | `trial_003_discrete_variants.py` | Conjugation × sign × routing (1,536 variants) |
| 004 | `trial_004_basis_automorphisms.py` | E8 automorphism basis changes |
| 005 | `trial_005_triple_okubo.py` | Triple Okubo/para-octonion (base + discrete variants) |
| 006 | `trial_006_triple_okubo_automorphisms.py` | Triple Okubo/para-octonion + E8 automorphisms |

All trials use fixed random seeds for reproducibility.  Results are recorded
in `../../evidence_and_reasoning/trial_NNN_results.md`.

## Guidelines

- Each module should have a clear, single responsibility.
- All non-trivial algorithms must reference the mathematical source they
  implement (e.g., a paper or textbook section).
- Trial files should be readable top-to-bottom by a human without requiring
  context from prior trials.
- No hardcoded personal information or environment-specific paths.
