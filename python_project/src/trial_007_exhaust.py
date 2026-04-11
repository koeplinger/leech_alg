"""
trial_007_exhaust.py — Exhaustive / large-scale Leech closure verification.

Performance-optimized for the exhaustive test: all 196,560² = 38,635,833,600
pairs of Min(Λ) vectors through the swap(1,2) triple octonion product.

NOT reference code.  The clean implementations are in trial_001_triple_octonion.py
and trial_007_explanation.py.  This file trades readability for speed.

Key optimisations
=================
1. Fixed-first-vector trick.  For exhaustive row-by-row enumeration, the first
   vector a is constant across all 196,560 second vectors.  The bilinear product
   result[k] = Σ_{ij} M[i,j,k] a[i] b[j]  factors as  result = B @ C  where
   C[j,k] = Σ_i M[i,j,k] a[i]  is precomputed once per row.  This replaces
   N×8×8×8 = 512N multiply-adds with N×8×8 = 64N — an 8× reduction.

2. Vectorised membership checks.  All Wilson conditions (is_in_L, is_in_Ls̄,
   is_in_Ls) are batch-evaluated over NumPy/CuPy arrays.

3. CPU multiprocessing.  Rows are distributed across cores via mp.Pool.

4. GPU via CuPy (optional, auto-detected).  All computation stays on-device;
   only scalar failure counts transfer back.

Usage
=====
    python3 trial_007_exhaust.py                    # 4M random, all CPU cores
    python3 trial_007_exhaust.py --random 10000000  # 10M random pairs
    python3 trial_007_exhaust.py --exhaustive       # all 38.6B pairs
    python3 trial_007_exhaust.py --gpu              # use GPU (requires CuPy)
    python3 trial_007_exhaust.py --workers 8        # limit CPU workers

Install CuPy for GPU support:
    pip install cupy-cuda12x    # for CUDA 12.x
    pip install cupy-cuda11x    # for CUDA 11.x
"""

import numpy as np
import multiprocessing as mp
import time
import sys
import argparse
import os

# ---------------------------------------------------------------------------
# Imports from the clean codebase (used only for setup, not in hot path)
# ---------------------------------------------------------------------------
from octonions import OctonionAlgebra, STANDARD_FANO_TRIPLES
from leech_wilson import (
    leech_type1_vectors,
    leech_type2_vectors,
    leech_type3_vectors,
    _Ls_bar_basis_inv,
    _Ls_basis_inv,
)
from trial_001_triple_octonion import flatten_leech_vector

# ---------------------------------------------------------------------------
# GPU auto-detection
# ---------------------------------------------------------------------------
try:
    import cupy as cp
    # Test that a GPU is actually usable
    cp.array([1.0])
    HAS_GPU = True
except Exception:
    HAS_GPU = False

# ---------------------------------------------------------------------------
# Structure tensor for swap(1,2)
# ---------------------------------------------------------------------------

def _build_swap_algebra():
    perm = {i: i for i in range(1, 8)}
    perm[1], perm[2] = 2, 1
    triples = tuple((perm[a], perm[b], perm[c])
                    for (a, b, c) in STANDARD_FANO_TRIPLES)
    return OctonionAlgebra(triples, name="swap(1,2)")

_SWAP = _build_swap_algebra()

def _build_structure_tensor():
    M = np.zeros((8, 8, 8), dtype=np.float64)
    for i in range(8):
        for j in range(8):
            sign, k = _SWAP._table[i][j]
            M[i, j, k] = float(sign)
    return M

_M = _build_structure_tensor()

# Block routing: target_block[bi, bj] = bt
_TARGET = np.array([[0, 2, 1], [2, 1, 0], [1, 0, 2]], dtype=np.int32)

# Sublattice basis inverses (precomputed in leech_wilson.py)
_LSBAR_INV = np.ascontiguousarray(_Ls_bar_basis_inv)
_LS_INV = np.ascontiguousarray(_Ls_basis_inv)


# ===================================================================
# CPU path — NumPy
# ===================================================================

def multiply_fixed_first(fixed_vec, all_vecs):
    """
    Multiply a fixed 24-vector against all rows of all_vecs.

    Uses the factorisation: for fixed a and varying b,
      prod[k] = Σ_{ij} M[i,j,k] a[i] b[j] = Σ_j b[j] C[j,k]
    where C[j,k] = Σ_i M[i,j,k] a[i] is precomputed once.

    Parameters: fixed_vec (24,), all_vecs (N, 24)
    Returns: (N, 24) products
    """
    N = all_vecs.shape[0]
    result = np.zeros((N, 24), dtype=np.float64)

    # Precompute C matrices: one per source block of the fixed vector
    C = [np.tensordot(fixed_vec[bi*8:(bi+1)*8], _M, axes=([0], [0]))
         for bi in range(3)]  # each (8, 8)

    for bi in range(3):
        for bj in range(3):
            bt = _TARGET[bi, bj]
            result[:, bt*8:(bt+1)*8] += all_vecs[:, bj*8:(bj+1)*8] @ C[bi]

    return result


def batch_multiply_24(A, B):
    """General batch multiply for random pair sampling. A, B: (N, 24)."""
    N = A.shape[0]
    result = np.zeros((N, 24), dtype=np.float64)
    for bi in range(3):
        Ai = A[:, bi*8:(bi+1)*8]
        for bj in range(3):
            Bj = B[:, bj*8:(bj+1)*8]
            bt = _TARGET[bi, bj]
            result[:, bt*8:(bt+1)*8] += np.einsum(
                'ijk,ni,nj->nk', _M, Ai, Bj, optimize=True)
    return result


def batch_is_in_L(X, tol=1e-9):
    """Check each row of X (N×8) for E8 lattice membership. Returns bool (N,)."""
    c2 = X * 2.0
    c2r = np.round(c2)
    ok = np.max(np.abs(c2 - c2r), axis=1) < tol
    c2i = c2r.astype(np.int64)
    parities = c2i % 2
    all_even = np.all(parities == 0, axis=1)
    all_odd = np.all(parities == 1, axis=1)
    s = np.sum(c2i, axis=1)
    return ok & ((all_even & (s % 4 == 0)) | (all_odd & (s % 4 == 2)))


def batch_is_in_sublattice(V, basis_inv, tol=1e-9):
    """Check each row of V (N×8) for sublattice membership. Returns bool (N,)."""
    coords = V @ basis_inv
    return np.max(np.abs(coords - np.round(coords)), axis=1) < tol


def batch_check_leech(prods, tol=1e-9):
    """Check Wilson's 3 conditions on product vectors. prods: (N, 24). Returns bool (N,)."""
    X = prods[:, 0:8]
    Y = prods[:, 8:16]
    Z = prods[:, 16:24]
    ok = batch_is_in_L(X, tol)
    ok &= batch_is_in_L(Y, tol)
    ok &= batch_is_in_L(Z, tol)
    ok &= batch_is_in_sublattice(X + Y, _LSBAR_INV, tol)
    ok &= batch_is_in_sublattice(X + Z, _LSBAR_INV, tol)
    ok &= batch_is_in_sublattice(Y + Z, _LSBAR_INV, tol)
    ok &= batch_is_in_sublattice(X + Y + Z, _LS_INV, tol)
    return ok


# ===================================================================
# GPU path — CuPy
# ===================================================================

def gpu_multiply_fixed_first(fixed_vec_gpu, all_vecs_gpu, M_gpu):
    """GPU version of multiply_fixed_first."""
    N = all_vecs_gpu.shape[0]
    result = cp.zeros((N, 24), dtype=cp.float64)
    C = [cp.tensordot(fixed_vec_gpu[bi*8:(bi+1)*8], M_gpu, axes=([0], [0]))
         for bi in range(3)]
    for bi in range(3):
        for bj in range(3):
            bt = _TARGET[bi, bj]
            result[:, bt*8:(bt+1)*8] += all_vecs_gpu[:, bj*8:(bj+1)*8] @ C[bi]
    return result


def gpu_batch_multiply(A, B, M_gpu):
    """GPU general batch multiply."""
    N = A.shape[0]
    result = cp.zeros((N, 24), dtype=cp.float64)
    for bi in range(3):
        Ai = A[:, bi*8:(bi+1)*8]
        for bj in range(3):
            Bj = B[:, bj*8:(bj+1)*8]
            bt = _TARGET[bi, bj]
            result[:, bt*8:(bt+1)*8] += cp.einsum(
                'ijk,ni,nj->nk', M_gpu, Ai, Bj, optimize=True)
    return result


def gpu_batch_is_in_L(X, tol=1e-9):
    """GPU version of batch_is_in_L."""
    c2 = X * 2.0
    c2r = cp.round(c2)
    ok = cp.max(cp.abs(c2 - c2r), axis=1) < tol
    c2i = c2r.astype(cp.int64)
    parities = c2i % 2
    all_even = cp.all(parities == 0, axis=1)
    all_odd = cp.all(parities == 1, axis=1)
    s = cp.sum(c2i, axis=1)
    return ok & ((all_even & (s % 4 == 0)) | (all_odd & (s % 4 == 2)))


def gpu_batch_is_in_sublattice(V, basis_inv, tol=1e-9):
    """GPU version of batch_is_in_sublattice."""
    coords = V @ basis_inv
    return cp.max(cp.abs(coords - cp.round(coords)), axis=1) < tol


def gpu_batch_check_leech(prods, lsbar_inv, ls_inv, tol=1e-9):
    """GPU version of batch_check_leech."""
    X = prods[:, 0:8]
    Y = prods[:, 8:16]
    Z = prods[:, 16:24]
    ok = gpu_batch_is_in_L(X, tol)
    ok &= gpu_batch_is_in_L(Y, tol)
    ok &= gpu_batch_is_in_L(Z, tol)
    ok &= gpu_batch_is_in_sublattice(X + Y, lsbar_inv, tol)
    ok &= gpu_batch_is_in_sublattice(X + Z, lsbar_inv, tol)
    ok &= gpu_batch_is_in_sublattice(Y + Z, lsbar_inv, tol)
    ok &= gpu_batch_is_in_sublattice(X + Y + Z, ls_inv, tol)
    return ok


# ===================================================================
# Verification against reference implementation
# ===================================================================

def verify_fast_path(all_vecs, use_gpu=False, n_check=2000):
    """Spot-check fast multiplication + membership against reference code."""
    from trial_001_triple_octonion import multiply_24, check_leech_membership

    rng = np.random.RandomState(9999)
    N = all_vecs.shape[0]
    idxA = rng.randint(0, N, size=n_check)
    idxB = rng.randint(0, N, size=n_check)

    # --- Test 1: general batch multiply ---
    A_np = all_vecs[idxA]
    B_np = all_vecs[idxB]

    if use_gpu:
        A_g = cp.asarray(A_np)
        B_g = cp.asarray(B_np)
        M_g = cp.asarray(_M)
        lsbar_g = cp.asarray(_LSBAR_INV)
        ls_g = cp.asarray(_LS_INV)
        prods_fast = cp.asnumpy(gpu_batch_multiply(A_g, B_g, M_g))
        in_leech_fast = cp.asnumpy(gpu_batch_check_leech(
            cp.asarray(prods_fast), lsbar_g, ls_g))
    else:
        prods_fast = batch_multiply_24(A_np, B_np)
        in_leech_fast = batch_check_leech(prods_fast)

    mm_mul = 0
    mm_mem = 0
    for i in range(n_check):
        prod_slow = multiply_24(A_np[i], B_np[i], _SWAP)
        if not np.allclose(prods_fast[i], prod_slow, atol=1e-12):
            mm_mul += 1
        pn = float(np.dot(prod_slow, prod_slow))
        if pn < 1e-12:
            continue
        diag = check_leech_membership(prod_slow)
        if diag['in_leech'] != in_leech_fast[i]:
            mm_mem += 1

    # --- Test 2: fixed-first-vector multiply ---
    n_row_check = min(50, N)
    row_indices = rng.choice(N, n_row_check, replace=False)
    mm_row = 0
    # Use a small subset as second vectors
    subset_size = min(500, N)
    subset_idx = rng.choice(N, subset_size, replace=False)
    subset = all_vecs[subset_idx]

    for ri in row_indices:
        if use_gpu:
            prods_row = cp.asnumpy(gpu_multiply_fixed_first(
                cp.asarray(all_vecs[ri]), cp.asarray(subset), M_g))
        else:
            prods_row = multiply_fixed_first(all_vecs[ri], subset)
        prods_ref = batch_multiply_24(
            np.broadcast_to(all_vecs[ri:ri+1], subset.shape).copy(), subset)
        if not np.allclose(prods_row, prods_ref, atol=1e-12):
            mm_row += 1

    return mm_mul, mm_mem, mm_row


# ===================================================================
# CPU multiprocessing — worker
# ===================================================================

# Global state set by pool initializer (shared via fork COW)
_g_vecs = None
_g_M = None
_g_lsbar_inv = None
_g_ls_inv = None

def _init_worker(vecs, M, lsbar_inv, ls_inv):
    global _g_vecs, _g_M, _g_lsbar_inv, _g_ls_inv
    _g_vecs = vecs
    _g_M = M
    _g_lsbar_inv = lsbar_inv
    _g_ls_inv = ls_inv


def _worker_exhaustive(row_range):
    """Process a range of rows for the exhaustive test."""
    start, end = row_range
    N = _g_vecs.shape[0]
    failures = 0
    rows_done = 0

    for i in range(start, end):
        # Precompute C matrices for this fixed first vector
        fixed = _g_vecs[i]
        C = [np.tensordot(fixed[bi*8:(bi+1)*8], _g_M, axes=([0], [0]))
             for bi in range(3)]

        result = np.zeros((N, 24), dtype=np.float64)
        for bi in range(3):
            for bj in range(3):
                bt = _TARGET[bi, bj]
                result[:, bt*8:(bt+1)*8] += _g_vecs[:, bj*8:(bj+1)*8] @ C[bi]

        # Check membership
        norms_sq = np.sum(result * result, axis=1)
        nonzero = norms_sq >= 1e-12
        if np.any(nonzero):
            prods_nz = result[nonzero]
            X = prods_nz[:, 0:8]
            Y = prods_nz[:, 8:16]
            Z = prods_nz[:, 16:24]
            ok = batch_is_in_L(X)
            ok &= batch_is_in_L(Y)
            ok &= batch_is_in_L(Z)
            ok &= batch_is_in_sublattice(X + Y, _g_lsbar_inv)
            ok &= batch_is_in_sublattice(X + Z, _g_lsbar_inv)
            ok &= batch_is_in_sublattice(Y + Z, _g_lsbar_inv)
            ok &= batch_is_in_sublattice(X + Y + Z, _g_ls_inv)
            n_fail = int(np.sum(~ok))
            failures += n_fail

        rows_done += 1

    return (start, end, failures, rows_done)


def _worker_random(args):
    """Process a batch of random pairs."""
    seed, n_pairs = args
    rng = np.random.RandomState(seed)
    N = _g_vecs.shape[0]
    failures = 0
    batch_size = 50_000

    for offset in range(0, n_pairs, batch_size):
        bs = min(batch_size, n_pairs - offset)
        idxA = rng.randint(0, N, size=bs)
        idxB = rng.randint(0, N, size=bs)
        A = _g_vecs[idxA]
        B = _g_vecs[idxB]

        result = np.zeros((bs, 24), dtype=np.float64)
        for bi in range(3):
            Ai = A[:, bi*8:(bi+1)*8]
            for bj in range(3):
                Bj = B[:, bj*8:(bj+1)*8]
                bt = _TARGET[bi, bj]
                result[:, bt*8:(bt+1)*8] += np.einsum(
                    'ijk,ni,nj->nk', _g_M, Ai, Bj, optimize=True)

        norms_sq = np.sum(result * result, axis=1)
        nonzero = norms_sq >= 1e-12
        if np.any(nonzero):
            ok = batch_check_leech(result[nonzero])
            failures += int(np.sum(~ok))

    return failures


# ===================================================================
# GPU exhaustive run
# ===================================================================

def run_gpu_exhaustive(all_vecs_np):
    """Run the exhaustive test on GPU."""
    N = all_vecs_np.shape[0]
    total_pairs = N * N

    print(f"GPU exhaustive: {N:,} × {N:,} = {total_pairs:,} pairs")
    print()

    all_vecs_g = cp.asarray(all_vecs_np)
    M_g = cp.asarray(_M)
    lsbar_g = cp.asarray(_LSBAR_INV)
    ls_g = cp.asarray(_LS_INV)

    total_failures = 0
    report_interval = max(1, N // 200)  # ~200 progress reports

    t_start = time.time()

    for i in range(N):
        fixed_g = all_vecs_g[i]

        # Precompute C matrices on GPU
        C = [cp.tensordot(fixed_g[bi*8:(bi+1)*8], M_g, axes=([0], [0]))
             for bi in range(3)]

        result = cp.zeros((N, 24), dtype=cp.float64)
        for bi in range(3):
            for bj in range(3):
                bt = _TARGET[bi, bj]
                result[:, bt*8:(bt+1)*8] += all_vecs_g[:, bj*8:(bj+1)*8] @ C[bi]

        # Membership check on GPU
        norms_sq = cp.sum(result * result, axis=1)
        nonzero = norms_sq >= 1e-12
        if cp.any(nonzero):
            prods_nz = result[nonzero]
            ok = gpu_batch_check_leech(prods_nz, lsbar_g, ls_g)
            n_fail = int(cp.sum(~ok))
            total_failures += n_fail

            if n_fail > 0:
                elapsed = time.time() - t_start
                print(f"\n  *** FAILURE at row {i}, {n_fail} failures ***")
                print(f"      Total failures so far: {total_failures}")
                print(f"      Elapsed: {elapsed:.1f}s")

        if i % report_interval == 0 or i == N - 1:
            elapsed = time.time() - t_start
            pairs_done = (i + 1) * N
            per_pair = elapsed / pairs_done * 1e6
            pct = pairs_done / total_pairs * 100
            if i > 0:
                eta = (elapsed / (i + 1)) * (N - i - 1)
                eta_str = f"{eta:.0f}s" if eta < 3600 else f"{eta/3600:.1f}h"
            else:
                eta_str = "..."
            print(f"  Row {i:>7,}/{N:,}  ({pct:>6.2f}%)  "
                  f"fail={total_failures}  "
                  f"{elapsed:>8.1f}s  {per_pair:.2f}µs/pair  ETA {eta_str}")

    t_total = time.time() - t_start
    return total_failures, t_total


def run_gpu_random(all_vecs_np, n_pairs, seed=2026):
    """Run random sampling on GPU."""
    N = all_vecs_np.shape[0]
    batch_size = 200_000

    print(f"GPU random: {n_pairs:,} pairs from pool of {N*N:,}")
    print()

    all_vecs_g = cp.asarray(all_vecs_np)
    M_g = cp.asarray(_M)
    lsbar_g = cp.asarray(_LSBAR_INV)
    ls_g = cp.asarray(_LS_INV)

    rng = np.random.RandomState(seed)
    total_failures = 0
    tested = 0
    report_interval = max(1, n_pairs // (batch_size * 20)) * batch_size

    t_start = time.time()

    for offset in range(0, n_pairs, batch_size):
        bs = min(batch_size, n_pairs - offset)
        idxA = rng.randint(0, N, size=bs)
        idxB = rng.randint(0, N, size=bs)

        A = all_vecs_g[idxA]
        B = all_vecs_g[idxB]

        result = cp.zeros((bs, 24), dtype=cp.float64)
        for bi in range(3):
            Ai = A[:, bi*8:(bi+1)*8]
            for bj in range(3):
                Bj = B[:, bj*8:(bj+1)*8]
                bt = _TARGET[bi, bj]
                result[:, bt*8:(bt+1)*8] += cp.einsum(
                    'ijk,ni,nj->nk', M_g, Ai, Bj, optimize=True)

        norms_sq = cp.sum(result * result, axis=1)
        nonzero = norms_sq >= 1e-12
        if cp.any(nonzero):
            ok = gpu_batch_check_leech(result[nonzero], lsbar_g, ls_g)
            total_failures += int(cp.sum(~ok))

        tested += bs

        if tested % report_interval < batch_size or offset + bs >= n_pairs:
            elapsed = time.time() - t_start
            per_pair = elapsed / tested * 1e6
            eta = per_pair * (n_pairs - tested) / 1e6
            eta_str = f"{eta:.0f}s" if eta < 3600 else f"{eta/3600:.1f}h"
            print(f"  {tested:>14,}  fail={total_failures}  "
                  f"{elapsed:>8.1f}s  {per_pair:.1f}µs/pair  ETA {eta_str}")

    t_total = time.time() - t_start
    return total_failures, t_total


# ===================================================================
# Load vectors
# ===================================================================

def load_all_vectors():
    """Load and flatten all 196,560 minimal vectors."""
    print("Loading minimal vectors...")
    t0 = time.time()

    type1_raw = leech_type1_vectors()
    type2_raw = leech_type2_vectors()
    type3_raw = leech_type3_vectors()

    n1, n2, n3 = len(type1_raw), len(type2_raw), len(type3_raw)
    n_total = n1 + n2 + n3

    print(f"  Type 1: {n1:>8,}")
    print(f"  Type 2: {n2:>8,}")
    print(f"  Type 3: {n3:>8,}")
    print(f"  Total:  {n_total:>8,}")

    all_vecs = np.empty((n_total, 24), dtype=np.float64)
    idx = 0
    for raw in (type1_raw, type2_raw, type3_raw):
        for v in raw:
            all_vecs[idx] = flatten_leech_vector(v)
            idx += 1

    print(f"  Load time: {time.time() - t0:.1f}s")
    print(f"  Memory: {all_vecs.nbytes / 1e6:.1f} MB")
    print()
    return all_vecs


# ===================================================================
# Main
# ===================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Exhaustive Leech closure test for swap(1,2)")
    parser.add_argument('--exhaustive', action='store_true',
                        help='Run all 38.6B pairs')
    parser.add_argument('--random', type=int, default=4_000_000,
                        help='Number of random pairs (default: 4M)')
    parser.add_argument('--gpu', action='store_true',
                        help='Use GPU via CuPy')
    parser.add_argument('--workers', type=int, default=0,
                        help='CPU workers (0=auto, default: all cores)')
    parser.add_argument('--skip-verify', action='store_true',
                        help='Skip reference verification')
    parser.add_argument('--seed', type=int, default=2026)
    args = parser.parse_args()

    if args.gpu and not HAS_GPU:
        print("ERROR: --gpu requested but CuPy is not available.")
        print("Install with: pip install cupy-cuda12x")
        sys.exit(1)

    use_gpu = args.gpu
    n_workers = args.workers or os.cpu_count()

    print("=" * 70)
    if args.exhaustive:
        mode_str = "EXHAUSTIVE"
    else:
        mode_str = f"RANDOM {args.random:,} pairs"
    backend_str = "GPU (CuPy)" if use_gpu else f"CPU ({n_workers} workers)"
    print(f"CLOSURE TEST — swap(1,2) — {mode_str} — {backend_str}")
    print("=" * 70)
    print()

    if use_gpu:
        dev = cp.cuda.Device()
        print(f"GPU: {dev.id} — {cp.cuda.runtime.getDeviceProperties(dev.id)['name'].decode()}")
        mem = dev.mem_info
        print(f"     VRAM: {mem[1]/1e9:.1f} GB total, {mem[0]/1e9:.1f} GB free")
        print()

    # Load vectors
    all_vecs = load_all_vectors()
    N = all_vecs.shape[0]
    exhaustive_total = N * N

    # Verify
    if not args.skip_verify:
        print("Verifying fast path against reference (2000 pairs)...")
        t_ver = time.time()
        mm_mul, mm_mem, mm_row = verify_fast_path(all_vecs, use_gpu=use_gpu)
        t_ver = time.time() - t_ver
        print(f"  Multiply mismatches:     {mm_mul}")
        print(f"  Membership mismatches:   {mm_mem}")
        print(f"  Fixed-row mismatches:    {mm_row}")
        if mm_mul > 0 or mm_mem > 0 or mm_row > 0:
            print("  *** MISMATCH — ABORTING ***")
            sys.exit(1)
        print(f"  Verified OK ({t_ver:.1f}s)")
        print()

    # ---------------------------------------------------------------
    # GPU path
    # ---------------------------------------------------------------
    if use_gpu:
        t_start = time.time()
        if args.exhaustive:
            failures, t_total = run_gpu_exhaustive(all_vecs)
            tested = exhaustive_total
        else:
            failures, t_total = run_gpu_random(all_vecs, args.random, args.seed)
            tested = args.random

    # ---------------------------------------------------------------
    # CPU multiprocessing path
    # ---------------------------------------------------------------
    else:
        t_start = time.time()

        if args.exhaustive:
            # Distribute rows across workers
            chunk_size = max(1, N // (n_workers * 20))  # ~20 chunks per worker
            row_ranges = []
            for s in range(0, N, chunk_size):
                row_ranges.append((s, min(s + chunk_size, N)))

            print(f"Distributing {N:,} rows across {n_workers} workers "
                  f"({len(row_ranges)} chunks of ~{chunk_size:,})")
            print()

            total_failures = 0
            rows_completed = 0

            with mp.Pool(n_workers,
                         initializer=_init_worker,
                         initargs=(all_vecs, _M, _LSBAR_INV, _LS_INV)) as pool:

                for result in pool.imap_unordered(_worker_exhaustive, row_ranges):
                    start, end, fails, n_rows = result
                    total_failures += fails
                    rows_completed += n_rows
                    elapsed = time.time() - t_start
                    pairs_done = rows_completed * N
                    pct = rows_completed / N * 100
                    per_pair = elapsed / pairs_done * 1e6
                    if rows_completed < N:
                        eta = (elapsed / rows_completed) * (N - rows_completed)
                        eta_str = (f"{eta:.0f}s" if eta < 3600
                                   else f"{eta/3600:.1f}h")
                    else:
                        eta_str = "done"
                    print(f"  Rows {rows_completed:>7,}/{N:,} ({pct:>6.2f}%)  "
                          f"fail={total_failures}  "
                          f"{elapsed:>8.1f}s  {per_pair:.2f}µs/pair  "
                          f"ETA {eta_str}")

            failures = total_failures
            tested = exhaustive_total

        else:
            # Random sampling: split across workers
            per_worker = args.random // n_workers
            remainder = args.random % n_workers
            worker_args = []
            for w in range(n_workers):
                n = per_worker + (1 if w < remainder else 0)
                worker_args.append((args.seed + w, n))

            print(f"Distributing {args.random:,} random pairs across "
                  f"{n_workers} workers")
            print()

            with mp.Pool(n_workers,
                         initializer=_init_worker,
                         initargs=(all_vecs, _M, _LSBAR_INV, _LS_INV)) as pool:
                results = pool.map(_worker_random, worker_args)

            failures = sum(results)
            tested = args.random

        t_total = time.time() - t_start

    # ---------------------------------------------------------------
    # Results
    # ---------------------------------------------------------------
    print()
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    print()
    print(f"  Mode:            {'exhaustive' if args.exhaustive else 'random'}")
    print(f"  Backend:         {'GPU' if use_gpu else f'CPU ({n_workers} workers)'}")
    print(f"  Pairs tested:    {tested:>18,}")
    print(f"  Failures:        {failures:>18,}")
    print(f"  Failure rate:    {failures/tested*100:>18.10f}%")
    print()
    print(f"  Wall-clock time: {t_total:>18.1f}s  ({t_total/60:.1f} min)")
    print(f"  Per pair:        {t_total/tested*1e6:>18.2f}µs")
    print()

    if not args.exhaustive:
        per_pair_us = t_total / tested * 1e6
        est_s = per_pair_us * exhaustive_total / 1e6
        est_h = est_s / 3600
        print("EXTRAPOLATION TO EXHAUSTIVE")
        print(f"  Exhaustive pairs:  {exhaustive_total:>18,}")
        print(f"  Estimated time:    {est_h:>18.1f}h ({est_h/24:.1f}d)")
        print()

    if failures == 0:
        print(f"  *** ZERO FAILURES in {tested:,} pairs ***")
        if args.exhaustive:
            print()
            print("  EXHAUSTIVE VERIFICATION COMPLETE.")
            print("  The swap(1,2) triple octonion product closes on Min(Λ).")
            print("  By bilinearity and the fact that Min(Λ) generates Λ over Z,")
            print("  ALL products of Λ vectors lie in Λ.")
            print("  Therefore (Λ, +, ·) is an order under this product.")
    else:
        print(f"  *** {failures} FAILURES — closure does NOT hold ***")

    print()


if __name__ == "__main__":
    main()
