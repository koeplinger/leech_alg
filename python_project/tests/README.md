# Tests

Test suite for the Python project.  197 tests verify every shared tool used by the trials and by the symbolic proof.

Run the full suite from `python_project/`:

```bash
python3 -m pytest tests/ -v
```

Scope: this suite covers the **foundations** (octonion algebra, the E8 and Leech constructions of Wilson and Dixon, the Okubo algebra), so that no trial or proof rests on an unverified tool.  It does not cover the paper's own results.  Those are certified one by one by the standalone verification scripts in [`../src/`](../src/README.md), which are run directly and print their own pass/fail report, and independently in [`../../gap_project/`](../../gap_project/README.md).

## Contents

| File | Subject |
|---|---|
| `test_octonions.py` | Octonion algebra: structure constants, Fano-plane rule, composition (norm-multiplicativity), Dixon X-product |
| `test_e8_wilson.py` | Wilson's E8 lattice L: 240-root enumeration, membership test, structural element s |
| `test_e8_dixon.py` | Dixon's E8 constructions Ξ^even and A^odd: enumeration, lattice properties, integrality |
| `test_leech_wilson.py` | Wilson's Leech lattice Λ: 196,560 minimal vectors, three-condition membership test, type-1/2/3 decomposition |
| `test_leech_dixon.py` | Dixon's Leech lattice: 196,560 minimal vectors, comparison with Wilson's embedding (17,232 shared vectors) |
| `test_okubo.py` | Okubo algebra: Petersson construction, order-3 automorphism τ validation, structure constants |

## Categories

- **Unit tests**: individual functions and modules in isolation.
- **Mathematical validation tests**: computed results match known properties (root counts, lattice determinants, composition law).
- **Cross-construction tests**: Wilson and Dixon implementations agree where they should and differ where they should.

## Guidelines

- Tests are self-contained and runnable without external services.
- Each test file corresponds to a source module.
- Test names describe the property or behavior being verified.
- Random tests use fixed seeds for reproducibility.
