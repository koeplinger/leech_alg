"""
trial_007_explanation.py — Demonstrate the transposition-twisted triple
octonion product on the Leech lattice.

This file is a self-contained, minimal demonstration of the finding from
trial 007: applying a transposition of two imaginary basis elements to the
standard Fano triples produces a different multiplication on R^8 that,
when used in the triple product with Z_3 routing, closes on the Leech
lattice.

The construction
================
1. Start with the standard Fano triples:
       (1,2,4), (2,3,5), (3,4,6), (4,5,7), (5,6,1), (6,7,2), (7,1,3)

2. Apply the transposition (e_1 <-> e_2) to every triple:
       (2,1,4), (1,3,5), (3,4,6), (4,5,7), (5,6,2), (6,7,1), (7,2,3)

   This changes the multiplication table.  For example:
       Standard:   e_1 * e_2 = +e_4
       Swapped:    e_1 * e_2 = -e_4    (sign flip!)

   In total, 30 of the 64 table entries change (sign or target index).

3. Build three copies of this swapped algebra on R^24:
       Block 0: indices  0– 7
       Block 1: indices  8–15
       Block 2: indices 16–23

4. Define the triple product with Z_3 cross-block routing:
       Same block:  O_a × O_a → O_a
       Cross block: O_a × O_b → O_c   where {a,b,c} = {0,1,2}

   All block products (same and cross) use the SAME swapped multiplication.

5. Test: for any two vectors u, v in Min(Λ), is u · v in Λ?

Why "different but isomorphic" matters
======================================
The standard and swapped algebras are both octonion algebras over R.
Since the real octonion algebra is unique up to isomorphism, there exists
a linear map φ: R^8 → R^8 such that φ(x ·_std y) = φ(x) ·_swap φ(y).

But φ moves the lattice: φ(L) ≠ L in general.  The Leech lattice is a
fixed geometric object in R^24.  The swapped multiplication happens to
map Λ × Λ → Λ; the standard multiplication does not.

Analogy: consider the integers Z ⊂ R.  Define two multiplications on R:
  (a) standard: x · y = xy
  (b) scaled:   x * y = 2xy
Both are isomorphic as R-algebras (via φ(x) = 2x).  But Z is closed
under (a) and not under (b).  A different lattice (½Z) is closed under (b).

References
==========
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              J. Algebra 322 (2009) 2186–2190.
"""

import numpy as np
from collections import Counter
from typing import Dict, Tuple

from octonions import OctonionAlgebra, STANDARD_FANO_TRIPLES, STANDARD_ALGEBRA
from leech_wilson import (
    leech_type1_vectors,
    leech_type2_vectors,
    leech_type3_vectors,
    is_in_leech,
)
from trial_001_triple_octonion import (
    multiply_24,
    check_leech_membership,
    flatten_leech_vector,
)


# ---------------------------------------------------------------------------
# Step 1: Build the transposition-twisted algebra
# ---------------------------------------------------------------------------

def build_swap_algebra(s: int, t: int) -> OctonionAlgebra:
    """
    Build the octonion algebra obtained by swapping imaginary basis
    elements e_s and e_t in the standard Fano triples.

    Parameters
    ----------
    s, t : int in {1,...,7}, s ≠ t

    Returns
    -------
    OctonionAlgebra with the transposed Fano triples.
    """
    perm = {i: i for i in range(1, 8)}
    perm[s] = t
    perm[t] = s
    new_triples = tuple(
        (perm[a], perm[b], perm[c]) for (a, b, c) in STANDARD_FANO_TRIPLES
    )
    return OctonionAlgebra(new_triples, name=f"swap({s},{t})")


# The specific algebra used throughout:
SWAP_ALGEBRA = build_swap_algebra(1, 2)


# ---------------------------------------------------------------------------
# Step 2: Show the multiplication table differences
# ---------------------------------------------------------------------------

def show_table_differences():
    """Print every entry where the standard and swapped tables differ."""
    std = STANDARD_ALGEBRA
    swp = SWAP_ALGEBRA

    print("Multiplication table differences (standard vs swap(1,2)):")
    print(f"{'':>12s}  {'standard':>12s}  {'swap(1,2)':>12s}")
    print("-" * 42)

    count = 0
    for i in range(8):
        for j in range(8):
            s_s, k_s = std._table[i][j]
            s_w, k_w = swp._table[i][j]
            if s_s != s_w or k_s != k_w:
                count += 1
                print(f"  e_{i}*e_{j}:  {s_s:+d}·e_{k_s}        {s_w:+d}·e_{k_w}")

    print(f"\nTotal differing entries: {count}/64")


# ---------------------------------------------------------------------------
# Step 3: Closure test
# ---------------------------------------------------------------------------

def run_closure_test(verbose: bool = True):
    """
    Test Leech lattice closure for the transposition-twisted triple product.

    Uses swap(1,2) throughout.  Tests all type combinations with sample
    sizes chosen for confidence without excessive runtime.
    """
    alg = SWAP_ALGEBRA

    if verbose:
        print("=" * 60)
        print("Transposition-twisted triple octonion product")
        print(f"Algebra: {alg.name}")
        print(f"Fano triples: {alg.fano_triples}")
        print("=" * 60)
        print()
        show_table_differences()
        print()

    # Generate minimal-shell vectors
    type1 = [flatten_leech_vector(t) for t in leech_type1_vectors()]
    type2 = [flatten_leech_vector(t) for t in leech_type2_vectors()]
    type3_raw = leech_type3_vectors()
    rng = np.random.RandomState(42)
    idx = rng.choice(len(type3_raw), 2000, replace=False)
    type3 = [flatten_leech_vector(type3_raw[i]) for i in idx]

    if verbose:
        print(f"Type 1: {len(type1)} vectors")
        print(f"Type 2: {len(type2)} vectors")
        print(f"Type 3: {len(type3)} vectors (sampled from {len(type3_raw)})")
        print()

    # Test sets
    tests = [
        ("t1×t1", type1, type1, 5000),
        ("t1×t2", type1[:100], type2[:100], 5000),
        ("t2×t2", type2[:200], type2[:200], 5000),
        ("t1×t3", type1[:100], type3[:100], 5000),
        ("t2×t3", type2[:100], type3[:100], 5000),
        ("t3×t3", type3[:200], type3[:200], 5000),
    ]

    total_tested = 0
    total_fail = 0

    for label, sa, sb, limit in tests:
        tested = 0
        fails = 0
        for i, a in enumerate(sa):
            if tested >= limit:
                break
            for j, b in enumerate(sb):
                if tested >= limit:
                    break
                prod = multiply_24(a, b, alg)
                pn = float(np.dot(prod, prod))
                tested += 1
                if pn < 1e-12:
                    continue
                d = check_leech_membership(prod)
                if not d['in_leech']:
                    fails += 1

        total_tested += tested
        total_fail += fails
        rate = fails / max(tested, 1) * 100
        if verbose:
            print(f"  {label}: {tested:>6d} tested, {fails} failures ({rate:.1f}%)")

    if verbose:
        print(f"\n  TOTAL: {total_tested} tested, {total_fail} failures")
        if total_fail == 0:
            print("  *** ALL products lie in Λ ***")


if __name__ == "__main__":
    run_closure_test(verbose=True)
