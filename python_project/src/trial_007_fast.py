"""
trial_007_fast.py — Performance-optimized Leech closure test for swap(1,2).

NOT intended as readable reference code.  The clean, verifiable implementations
are in trial_001_triple_octonion.py, trial_007_kirmse_twist.py, and
trial_007_explanation.py.  This file trades clarity for speed: it precomputes
a dense structure-constant tensor, batches all multiplications via np.einsum,
and vectorises all Wilson membership checks with bulk NumPy operations.

The purpose is to push toward exhaustive verification of the 38.6 billion
pairs in Min(Λ) × Min(Λ), or at least to test a much larger sample than
the 4 million pairs that took 14 minutes in pure Python.

Usage
=====
    python3 trial_007_fast.py [N_PAIRS]       # default 4,000,000
    python3 trial_007_fast.py exhaustive      # run all 38.6B pairs
"""

import numpy as np
import time
import sys

# ---------------------------------------------------------------------------
# We import ONLY what we need for vector generation.  All hot-path code is
# reimplemented here in vectorised NumPy.
# ---------------------------------------------------------------------------
from octonions import OctonionAlgebra, STANDARD_FANO_TRIPLES
from leech_wilson import (
    leech_type1_vectors,
    leech_type2_vectors,
    leech_type3_vectors,
    _L_zbasis,
    _Ls_bar_basis_inv,
    _Ls_basis_inv,
)
from trial_001_triple_octonion import flatten_leech_vector


# ---------------------------------------------------------------------------
# Swap(1,2) algebra — structure constant tensor
# ---------------------------------------------------------------------------

def _build_swap_algebra():
    perm = {i: i for i in range(1, 8)}
    perm[1], perm[2] = 2, 1
    triples = tuple((perm[a], perm[b], perm[c])
                    for (a, b, c) in STANDARD_FANO_TRIPLES)
    return OctonionAlgebra(triples, name="swap(1,2)")

_SWAP = _build_swap_algebra()

def _build_structure_tensor(alg):
    """8×8×8 tensor M where M[i,j,k] = sign for e_i*e_j = sign*e_k."""
    M = np.zeros((8, 8, 8), dtype=np.float64)
    for i in range(8):
        for j in range(8):
            sign, k = alg._table[i][j]
            M[i, j, k] = float(sign)
    return M

# Precompute once at import time
_M = _build_structure_tensor(_SWAP)

# Block routing: target_block[bi, bj] = bt
_TARGET = np.array([
    [0, 2, 1],   # block 0 × {0,1,2} → {0,2,1}
    [2, 1, 0],   # block 1 × {0,1,2} → {2,1,0}
    [1, 0, 2],   # block 2 × {0,1,2} → {1,0,2}
], dtype=np.int32)


# ---------------------------------------------------------------------------
# Batch 24-dim multiplication
# ---------------------------------------------------------------------------

def batch_multiply_24(A, B):
    """
    Multiply pairs of 24-vectors using the swap(1,2) triple-octonion product.

    A, B : (N, 24) float64
    Returns: (N, 24) float64
    """
    N = A.shape[0]
    result = np.zeros((N, 24), dtype=np.float64)

    for bi in range(3):
        s0i = bi * 8
        Ai = A[:, s0i:s0i+8]
        for bj in range(3):
            s0j = bj * 8
            Bj = B[:, s0j:s0j+8]
            bt = _TARGET[bi, bj]
            s0t = bt * 8
            # Core: einsum over the 8×8×8 structure tensor
            # prod[n, k] = sum_{i,j} M[i,j,k] * Ai[n,i] * Bj[n,j]
            result[:, s0t:s0t+8] += np.einsum('ijk,ni,nj->nk', _M, Ai, Bj,
                                               optimize=True)
    return result


# ---------------------------------------------------------------------------
# Batch Wilson membership checks
# ---------------------------------------------------------------------------

# Precompute L Z-basis inverse for is_in_L check via lattice coordinates.
_L_basis_inv = np.linalg.inv(_L_zbasis)  # 8×8


def batch_is_in_L(X, tol=1e-9):
    """
    Check whether each row of X (shape N×8) is in Wilson's E8 lattice L = D8+.

    Returns bool array of shape (N,).
    """
    c2 = X * 2.0
    c2r = np.round(c2)
    # All coords must be half-integer multiples
    ok_half = np.max(np.abs(c2 - c2r), axis=1) < tol
    c2i = c2r.astype(np.int64)
    parities = np.abs(c2i) % 2  # 0 or 1; use abs to handle negatives safely
    # Actually in NumPy, (-3) % 2 == 1 (Python semantics), so parities is fine
    # without abs.  But let's be safe for clarity.
    all_even = np.all(parities == 0, axis=1)
    all_odd = np.all(parities == 1, axis=1)
    s = np.sum(c2i, axis=1)
    # Integer coords: sum(x_i) even  ⟺  sum(c2i) % 4 == 0
    case_int = all_even & (s % 4 == 0)
    # Half-integer coords: sum(x_i) odd integer  ⟺  sum(c2i) % 4 == 2
    # (sum of 8 odd numbers is always even, so s%4 ∈ {0,2})
    case_half = all_odd & (np.abs(s) % 4 == 2)
    return ok_half & (case_int | case_half)


def batch_is_in_sublattice(V, basis_inv, tol=1e-9):
    """
    Check whether each row of V (shape N×8) lies in the Z-span of a lattice
    whose basis inverse is given.

    V @ basis_inv should yield integer coordinates.
    Returns bool array of shape (N,).
    """
    coords = V @ basis_inv   # (N, 8)
    rounded = np.round(coords)
    return np.max(np.abs(coords - rounded), axis=1) < tol


def batch_check_leech(prods, tol=1e-9):
    """
    Check Wilson's 3 conditions on each product vector.

    prods : (N, 24) float64
    Returns : bool array (N,) — True if in Λ.
    """
    X = prods[:, 0:8]
    Y = prods[:, 8:16]
    Z = prods[:, 16:24]

    # Condition 1: x, y, z ∈ L
    ok = batch_is_in_L(X, tol)
    ok &= batch_is_in_L(Y, tol)
    ok &= batch_is_in_L(Z, tol)

    # Early exit if all fail (unlikely but free)
    if not np.any(ok):
        return ok

    # Condition 2: x+y, x+z, y+z ∈ Ls̄
    ok &= batch_is_in_sublattice(X + Y, _Ls_bar_basis_inv, tol)
    ok &= batch_is_in_sublattice(X + Z, _Ls_bar_basis_inv, tol)
    ok &= batch_is_in_sublattice(Y + Z, _Ls_bar_basis_inv, tol)

    # Condition 3: x+y+z ∈ Ls
    ok &= batch_is_in_sublattice(X + Y + Z, _Ls_basis_inv, tol)

    return ok


# ---------------------------------------------------------------------------
# Verification: compare fast path against slow path on a small sample
# ---------------------------------------------------------------------------

def verify_against_reference(all_vecs, n_check=1000):
    """Spot-check that fast multiply + membership agrees with reference code."""
    from trial_001_triple_octonion import multiply_24, check_leech_membership

    rng = np.random.RandomState(9999)
    N = all_vecs.shape[0]
    idxA = rng.randint(0, N, size=n_check)
    idxB = rng.randint(0, N, size=n_check)

    A = all_vecs[idxA]
    B = all_vecs[idxB]

    # Fast path
    prods_fast = batch_multiply_24(A, B)
    in_leech_fast = batch_check_leech(prods_fast)

    # Slow path
    mismatches_mul = 0
    mismatches_mem = 0
    for i in range(n_check):
        prod_slow = multiply_24(A[i], B[i], _SWAP)
        if not np.allclose(prods_fast[i], prod_slow, atol=1e-12):
            mismatches_mul += 1

        pn = float(np.dot(prod_slow, prod_slow))
        if pn < 1e-12:
            # Zero product — fast path should also show True (vacuously in Λ)
            # Actually we need to check: batch_check_leech would check
            # conditions on the zero vector.  Zero ∈ L, zero ∈ Ls̄, zero ∈ Ls.
            continue
        diag = check_leech_membership(prod_slow)
        if diag['in_leech'] != in_leech_fast[i]:
            mismatches_mem += 1

    return mismatches_mul, mismatches_mem


# ---------------------------------------------------------------------------
# Main test
# ---------------------------------------------------------------------------

def run_test(n_pairs=4_000_000, seed=2026, batch_size=50_000, exhaustive=False):
    print("=" * 70)
    print("FAST CLOSURE TEST — swap(1,2), vectorised NumPy")
    print("=" * 70)
    print()

    # ------------------------------------------------------------------
    # Load vectors
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

    print(f"  Type 1: {n1:>8,}")
    print(f"  Type 2: {n2:>8,}")
    print(f"  Type 3: {n3:>8,}")
    print(f"  Total:  {n_total:>8,}")

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

    t_load = time.time() - t0
    print(f"  Load time: {t_load:.1f}s")
    print()

    # ------------------------------------------------------------------
    # Verify fast path matches slow path
    # ------------------------------------------------------------------
    print("Verifying fast path against reference (1000 pairs)...")
    t_ver = time.time()
    mm_mul, mm_mem = verify_against_reference(all_vecs, n_check=1000)
    t_ver = time.time() - t_ver
    print(f"  Multiplication mismatches: {mm_mul}")
    print(f"  Membership mismatches:     {mm_mem}")
    print(f"  Verification time:         {t_ver:.1f}s")
    if mm_mul > 0 or mm_mem > 0:
        print("  *** FAST PATH DOES NOT MATCH REFERENCE — ABORTING ***")
        return
    print("  Fast path verified OK.")
    print()

    # ------------------------------------------------------------------
    # Determine pair count
    # ------------------------------------------------------------------
    exhaustive_total = n_total * n_total

    if exhaustive:
        n_pairs = exhaustive_total
        print(f"EXHAUSTIVE MODE: testing all {n_pairs:,} pairs")
    else:
        print(f"Exhaustive pool: {exhaustive_total:,} pairs")
        print(f"Testing: {n_pairs:,} random pairs ({n_pairs/exhaustive_total*100:.6f}%)")
    print()

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------
    rng = np.random.RandomState(seed)

    total_tested = 0
    total_failures = 0
    total_zero = 0

    print(f"Batch size: {batch_size:,}")
    print()
    print(f"  {'Pairs':>14s}  {'Fail':>8s}  {'Zero':>8s}  "
          f"{'Elapsed':>10s}  {'Per pair':>10s}  {'ETA total':>12s}")
    print("  " + "-" * 70)

    t_start = time.time()
    pairs_remaining = n_pairs

    while pairs_remaining > 0:
        bs = min(batch_size, pairs_remaining)

        if exhaustive:
            # Sequential enumeration: row = total_tested // n_total,
            # col = total_tested % n_total
            start = total_tested
            row_indices = np.arange(start, start + bs) // n_total
            col_indices = np.arange(start, start + bs) % n_total
            idxA = row_indices.astype(np.int64)
            idxB = col_indices.astype(np.int64)
        else:
            idxA = rng.randint(0, n_total, size=bs)
            idxB = rng.randint(0, n_total, size=bs)

        A = all_vecs[idxA]
        B = all_vecs[idxB]

        prods = batch_multiply_24(A, B)

        # Count zero products
        norms_sq = np.sum(prods * prods, axis=1)
        zero_mask = norms_sq < 1e-12
        n_zero = int(np.sum(zero_mask))

        # Check non-zero products for Leech membership
        nonzero_mask = ~zero_mask
        n_nonzero = int(np.sum(nonzero_mask))

        if n_nonzero > 0:
            in_leech = batch_check_leech(prods[nonzero_mask])
            n_fail = int(np.sum(~in_leech))
        else:
            n_fail = 0

        total_tested += bs
        total_failures += n_fail
        total_zero += n_zero
        pairs_remaining -= bs

        # Progress report every ~500K pairs or at end
        if total_tested % (batch_size * 10) < batch_size or pairs_remaining == 0:
            elapsed = time.time() - t_start
            per_pair = elapsed / total_tested * 1e6
            if n_pairs > total_tested:
                eta = per_pair * (n_pairs - total_tested) / 1e6
                eta_str = f"{eta:.0f}s" if eta < 3600 else f"{eta/3600:.1f}h"
            else:
                eta_str = "done"
            print(f"  {total_tested:>14,}  {total_failures:>8,}  {total_zero:>8,}  "
                  f"{elapsed:>9.1f}s  {per_pair:>8.1f}µs  {eta_str:>12s}")

        # Early report on first failure
        if n_fail > 0 and total_failures <= 3:
            fail_indices = np.where(~in_leech)[0]
            for fi in fail_indices[:3]:
                global_idx = np.where(nonzero_mask)[0][fi]
                print(f"\n  *** FAILURE at pair ({idxA[global_idx]}, {idxB[global_idx]}) ***")
                print(f"      ||prod||² = {norms_sq[global_idx]:.1f}")

    t_total = time.time() - t_start

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print()
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    print()
    print(f"  Pairs tested:    {total_tested:>14,}")
    print(f"  Failures:        {total_failures:>14,}")
    print(f"  Zero products:   {total_zero:>14,}")
    print(f"  Failure rate:    {total_failures/total_tested*100:>14.8f}%")
    print()
    print(f"  Wall-clock time: {t_total:>14.1f}s  ({t_total/60:.1f} min)")
    print(f"  Per pair:        {t_total/total_tested*1e6:>14.1f}µs")
    print()

    # Speedup vs slow path
    slow_per_pair = 213.3  # µs, from trial_007_scaled_test.py
    fast_per_pair = t_total / total_tested * 1e6
    speedup = slow_per_pair / fast_per_pair
    print(f"  Speedup vs pure Python: {speedup:.1f}×")
    print()

    # Extrapolation
    if not exhaustive:
        est_s = fast_per_pair * exhaustive_total / 1e6
        est_h = est_s / 3600
        est_d = est_h / 24
        print("EXTRAPOLATION TO EXHAUSTIVE TEST")
        print(f"  Exhaustive pairs:  {exhaustive_total:>18,}")
        print(f"  Estimated time:    {est_s:>18,.0f}s")
        print(f"                     {est_h:>18,.1f}h")
        print(f"                     {est_d:>18,.1f}d")
        print()

    if total_failures == 0:
        print(f"  *** ZERO FAILURES in {total_tested:,} pairs ***")
    else:
        print(f"  *** {total_failures} FAILURES ***")

    return {
        'tested': total_tested,
        'failures': total_failures,
        'zero_products': total_zero,
        'wall_time_s': t_total,
        'per_pair_us': fast_per_pair,
    }


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "exhaustive":
        run_test(exhaustive=True, batch_size=100_000)
    else:
        n = 4_000_000
        if len(sys.argv) > 1:
            n = int(sys.argv[1])
        run_test(n_pairs=n)
