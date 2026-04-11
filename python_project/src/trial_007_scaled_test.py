"""
trial_007_scaled_test.py — Scaled closure test for the transposition-twisted
triple octonion product.

Purpose
=======
Run 4,000,000 random pairs from Min(Λ) through the swap(1,2) multiplication
and check Leech lattice closure.  Measure wall-clock time to extrapolate how
long the fully exhaustive test (196,560² ≈ 3.86 × 10¹⁰ pairs) would take.

This addresses caveat 1 from trial_007_results.md: the previous test covered
593,400 pairs, which is suggestive but not exhaustive.

References
==========
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              J. Algebra 322 (2009) 2186–2190.
"""

import numpy as np
import time
import sys

from octonions import OctonionAlgebra, STANDARD_FANO_TRIPLES
from leech_wilson import (
    leech_type1_vectors,
    leech_type2_vectors,
    leech_type3_vectors,
)
from trial_001_triple_octonion import (
    multiply_24,
    check_leech_membership,
    flatten_leech_vector,
)


def build_swap_algebra(s: int, t: int) -> OctonionAlgebra:
    """Build the octonion algebra with imaginary basis elements s,t swapped."""
    perm = {i: i for i in range(1, 8)}
    perm[s] = t
    perm[t] = s
    new_triples = tuple(
        (perm[a], perm[b], perm[c]) for (a, b, c) in STANDARD_FANO_TRIPLES
    )
    return OctonionAlgebra(new_triples, name=f"swap({s},{t})")


def run_scaled_test(n_pairs: int = 4_000_000, seed: int = 2026):
    """
    Test n_pairs random pairs from Min(Λ) for Leech closure under swap(1,2).

    Draws pairs uniformly at random from the full set of 196,560 minimal
    vectors (all three types combined).  Reports timing and extrapolation.
    """
    alg = build_swap_algebra(1, 2)

    print("=" * 70)
    print(f"SCALED CLOSURE TEST — swap(1,2), {n_pairs:,} pairs")
    print("=" * 70)
    print()

    # ------------------------------------------------------------------
    # Load all minimal vectors
    # ------------------------------------------------------------------
    print("Loading minimal vectors...")
    t0 = time.time()

    type1_raw = leech_type1_vectors()
    type2_raw = leech_type2_vectors()
    type3_raw = leech_type3_vectors()

    n1 = len(type1_raw)
    n2 = len(type2_raw)
    n3 = len(type3_raw)
    n_total = n1 + n2 + n3

    print(f"  Type 1: {n1:>8,} vectors")
    print(f"  Type 2: {n2:>8,} vectors")
    print(f"  Type 3: {n3:>8,} vectors")
    print(f"  Total:  {n_total:>8,} vectors")
    print()

    # Flatten all into a single array for fast random access
    print("Flattening to 24-vectors...")
    all_vecs = np.empty((n_total, 24), dtype=np.float64)
    idx = 0
    for v in type1_raw:
        all_vecs[idx] = flatten_leech_vector(v)
        idx += 1
    for v in type2_raw:
        all_vecs[idx] = flatten_leech_vector(v)
        idx += 1
    for v in type3_raw:
        all_vecs[idx] = flatten_leech_vector(v)
        idx += 1
    assert idx == n_total

    t_load = time.time() - t0
    print(f"  Loading + flattening took {t_load:.1f}s")
    print()

    exhaustive_pairs = n_total * n_total
    print(f"Exhaustive pool: {n_total:,} × {n_total:,} = {exhaustive_pairs:,.0f} pairs")
    print(f"Testing: {n_pairs:,} random pairs ({n_pairs/exhaustive_pairs*100:.4f}%)")
    print()

    # ------------------------------------------------------------------
    # Draw random pairs and test
    # ------------------------------------------------------------------
    rng = np.random.RandomState(seed)

    failures = 0
    zero_products = 0
    tested = 0

    # Report progress every this many pairs
    report_interval = 200_000

    print("Running closure test...")
    print(f"  {'Pairs':>12s}  {'Failures':>10s}  {'Zero':>8s}  {'Rate':>10s}  {'Elapsed':>10s}  {'Per pair':>10s}")
    print("  " + "-" * 66)

    t_start = time.time()

    for batch_start in range(0, n_pairs, report_interval):
        batch_end = min(batch_start + report_interval, n_pairs)
        batch_size = batch_end - batch_start

        # Draw random index pairs
        idx_a = rng.randint(0, n_total, size=batch_size)
        idx_b = rng.randint(0, n_total, size=batch_size)

        for k in range(batch_size):
            a = all_vecs[idx_a[k]]
            b = all_vecs[idx_b[k]]

            prod = multiply_24(a, b, alg)
            prod_norm_sq = float(np.dot(prod, prod))

            tested += 1

            if prod_norm_sq < 1e-12:
                zero_products += 1
                continue

            diag = check_leech_membership(prod)
            if not diag['in_leech']:
                failures += 1
                # On first failure, print details
                if failures <= 3:
                    print(f"\n  *** FAILURE #{failures} at pair {tested}: ***")
                    print(f"      idx_a={idx_a[k]}, idx_b={idx_b[k]}")
                    print(f"      ||prod||² = {prod_norm_sq}")
                    for key in ['x_in_L', 'y_in_L', 'z_in_L',
                                'xy_in_Lsbar', 'xz_in_Lsbar', 'yz_in_Lsbar',
                                'xyz_in_Ls']:
                        if not diag[key]:
                            print(f"      VIOLATED: {key}")
                    print()

        elapsed = time.time() - t_start
        per_pair_us = elapsed / tested * 1e6
        fail_rate = failures / tested * 100
        print(f"  {tested:>12,}  {failures:>10,}  {zero_products:>8,}  "
              f"{fail_rate:>9.4f}%  {elapsed:>9.1f}s  {per_pair_us:>8.1f}µs")

    t_test = time.time() - t_start

    # ------------------------------------------------------------------
    # Results and extrapolation
    # ------------------------------------------------------------------
    print()
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    print()
    print(f"  Pairs tested:    {tested:>12,}")
    print(f"  Failures:        {failures:>12,}")
    print(f"  Zero products:   {zero_products:>12,}")
    print(f"  Failure rate:    {failures/tested*100:>12.6f}%")
    print()
    print(f"  Wall-clock time: {t_test:>12.1f}s  ({t_test/60:.1f} min)")
    print(f"  Per pair:        {t_test/tested*1e6:>12.1f}µs")
    print()

    # Extrapolation
    per_pair_s = t_test / tested
    exhaustive_s = per_pair_s * exhaustive_pairs
    exhaustive_h = exhaustive_s / 3600
    exhaustive_d = exhaustive_h / 24

    print("EXTRAPOLATION TO EXHAUSTIVE TEST")
    print(f"  Exhaustive pairs:  {exhaustive_pairs:>18,.0f}")
    print(f"  Estimated time:    {exhaustive_s:>18,.0f}s")
    print(f"                     {exhaustive_h:>18,.1f}h")
    print(f"                     {exhaustive_d:>18,.1f}d")
    print()

    if failures == 0:
        print(f"  *** ZERO FAILURES in {tested:,} pairs ***")
        print(f"  This is {tested/exhaustive_pairs*100:.4f}% of the exhaustive pool.")
    else:
        print(f"  *** {failures} FAILURES — closure does NOT hold ***")

    return {
        'tested': tested,
        'failures': failures,
        'zero_products': zero_products,
        'wall_time_s': t_test,
        'per_pair_us': t_test / tested * 1e6,
        'exhaustive_pairs': exhaustive_pairs,
        'exhaustive_time_s': exhaustive_s,
        'exhaustive_time_h': exhaustive_h,
        'exhaustive_time_d': exhaustive_d,
    }


if __name__ == "__main__":
    n = 4_000_000
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
    run_scaled_test(n_pairs=n)
