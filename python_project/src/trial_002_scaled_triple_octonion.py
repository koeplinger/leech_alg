"""
trial_002_scaled_triple_octonion.py — Per-block scaling of the triple-octonion algebra.

Algebra
=======
The same triple-octonion multiplication as trial 001 (three copies O₁, O₂, O₃
at indices 0–7, 8–15, 16–23 with Z₃-symmetric cross-block mixing), but now with
per-block scaling.

Each block α ∈ {1,2,3} has a real scaling factor sα.  The basis vectors for
block α are bᵢ = sα · eᵢ (the standard unit vectors scaled by sα).  The
multiplication of two vectors expressed in this basis uses the same ±1 octonion
structure constants, but when converted back to standard R²⁴ coordinates the
product acquires scale factors:

    Block k of product = (1/sₖ) · same_k + sₖ/(sᵢsⱼ) · cross_k

where same_k = aₖ * bₖ  (octonion product of same-block components)
      cross_k = aᵢ * bⱼ + aⱼ * bᵢ  (sum of cross-block octonion products
                                      that land in block k, with {i,j,k}={1,2,3})

Key simplification
------------------
Scaling all sα by a common factor c multiplies the entire product by 1/c.
So the product depends on only TWO independent ratios:

    u = s₂/s₁,   v = s₃/s₁

with the overall scale s₁ determined by the norm constraint:

    ||product||² = 8   ⟹   s₁ = ||q(u,v)|| / √8

where q(u,v) is the "shape" of the product (before the 1/s₁ overall factor):

    q_block_0 = same₀ + cross₀/(uv)
    q_block_1 = same₁/u + u·cross₁/v
    q_block_2 = same₂/v + v·cross₂/u

The normalized product is then (√8/||q||) · q.  We check whether this lies in Λ.

Search strategy
===============
1. Precompute same_k and cross_k building blocks for a sample of Min(Λ) pairs.
2. Sweep a 2D grid over (u, v), plus the Z₃-symmetric point u=v=1.
3. For each (u, v), compute the normalized product for each pair and batch-check
   Leech lattice membership.
4. Report the (u, v) values with the highest closure rate on Min(Λ).

References
==========
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186–2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
"""

import numpy as np
from collections import Counter
from typing import List, Tuple, Dict

from octonions import STANDARD_ALGEBRA, OctonionAlgebra
from leech_wilson import (
    leech_type1_vectors,
    leech_type2_vectors,
    leech_type3_vectors,
    is_in_leech,
)
import leech_wilson
from e8_wilson import is_in_L

# Access precomputed sublattice inverse bases from leech_wilson
_Ls_bar_basis_inv = leech_wilson._Ls_bar_basis_inv
_Ls_basis_inv = leech_wilson._Ls_basis_inv

BLOCK_SLICES = [slice(0, 8), slice(8, 16), slice(16, 24)]
ALG = STANDARD_ALGEBRA


# ---------------------------------------------------------------------------
# Vectorized Leech lattice membership
# ---------------------------------------------------------------------------

def batch_is_in_L(vecs: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    """
    Check if each row of vecs (N×8) is in Wilson's E8 lattice L = D8+.

    Returns a boolean array of shape (N,).
    """
    c2 = np.round(vecs * 2.0)
    close = np.all(np.abs(vecs * 2.0 - c2) < tol, axis=1)
    c2i = c2.astype(np.int64)
    parities = c2i % 2  # 0 or 1 (Python % is non-negative for int)
    all_even = np.all(parities == 0, axis=1)
    all_odd = np.all(parities == 1, axis=1)
    s = np.sum(c2i, axis=1)
    even_check = (s % 4 == 0)
    odd_check = (s % 4 == 2)
    return close & ((all_even & even_check) | (all_odd & odd_check))


def batch_in_sublattice(vecs: np.ndarray, inv_basis: np.ndarray,
                        tol: float = 1e-9) -> np.ndarray:
    """
    Check if each row of vecs (N×8) is in the Z-span of the given basis.

    Parameters
    ----------
    vecs : (N, 8) array
    inv_basis : (8, 8) array — inverse of the basis matrix.

    Returns boolean array of shape (N,).
    """
    coords = vecs @ inv_basis  # (N, 8)
    return np.all(np.abs(coords - np.round(coords)) < tol, axis=1)


def batch_is_in_leech(prods: np.ndarray) -> np.ndarray:
    """
    Check if each row of prods (N×24) is in the Leech lattice Λ.

    Implements Wilson's three conditions in batch.
    Returns boolean array of shape (N,).
    """
    x = prods[:, 0:8]
    y = prods[:, 8:16]
    z = prods[:, 16:24]

    # Condition 1: x, y, z ∈ L
    ok = batch_is_in_L(x) & batch_is_in_L(y) & batch_is_in_L(z)

    # Condition 2: x+y, x+z, y+z ∈ Ls̄
    ok &= batch_in_sublattice(x + y, _Ls_bar_basis_inv)
    ok &= batch_in_sublattice(x + z, _Ls_bar_basis_inv)
    ok &= batch_in_sublattice(y + z, _Ls_bar_basis_inv)

    # Condition 3: x+y+z ∈ Ls
    ok &= batch_in_sublattice(x + y + z, _Ls_basis_inv)

    return ok


# ---------------------------------------------------------------------------
# Precompute building blocks
# ---------------------------------------------------------------------------

def flatten(trip: Tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    return np.concatenate(trip)


def precompute_pair_blocks(a: np.ndarray, b: np.ndarray) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """
    Precompute the same-block and cross-block octonion products for a pair.

    Returns (same, cross) where:
        same[k] = aₖ * bₖ                          (8-vector)
        cross[k] = aᵢ * bⱼ + aⱼ * bᵢ  with {i,j,k}={0,1,2}  (8-vector)
    """
    a_blocks = [a[s] for s in BLOCK_SLICES]
    b_blocks = [b[s] for s in BLOCK_SLICES]

    same = [ALG._mul_coords(a_blocks[k], b_blocks[k]) for k in range(3)]

    cross = [None, None, None]
    perm = [(1, 2), (0, 2), (0, 1)]  # for target block k, the other two
    for k in range(3):
        i, j = perm[k]
        cross[k] = (ALG._mul_coords(a_blocks[i], b_blocks[j])
                     + ALG._mul_coords(a_blocks[j], b_blocks[i]))

    return same, cross


def compute_q(same: List[np.ndarray], cross: List[np.ndarray],
              u: float, v: float) -> np.ndarray:
    """
    Compute the unnormalized product shape q(u, v).

    q_block_0 = same₀ + cross₀/(u*v)
    q_block_1 = same₁/u + u*cross₁/v
    q_block_2 = same₂/v + v*cross₂/u
    """
    q = np.zeros(24)
    q[0:8] = same[0] + cross[0] / (u * v)
    q[8:16] = same[1] / u + u * cross[1] / v
    q[16:24] = same[2] / v + v * cross[2] / u
    return q


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def run_trial(verbose: bool = True):
    """
    Execute trial 002: per-block scaling search.
    """
    if verbose:
        print("=" * 72)
        print("TRIAL 002: Per-block scaled triple-octonion algebra")
        print("=" * 72)

    # ------------------------------------------------------------------
    # Step 1: Select and precompute pairs
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 1: Generating Min(Λ) vectors and precomputing blocks...")

    type1 = [flatten(t) for t in leech_type1_vectors()]
    type2_raw = leech_type2_vectors()
    type3_raw = leech_type3_vectors()

    rng = np.random.RandomState(42)

    # Sample pairs from each type combination
    pair_configs = [
        ("t1×t1", type1, type1, 50),
        ("t1×t2", type1, [flatten(t) for t in type2_raw[:500]], 50),
        ("t2×t2", [flatten(t) for t in type2_raw[:500]],
                   [flatten(t) for t in type2_raw[:500]], 50),
        ("t1×t3", type1, [flatten(type3_raw[i]) for i in rng.choice(len(type3_raw), 500, replace=False)], 30),
        ("t2×t3", [flatten(t) for t in type2_raw[:200]],
                   [flatten(type3_raw[i]) for i in rng.choice(len(type3_raw), 500, replace=False)], 30),
        ("t3×t3", [flatten(type3_raw[i]) for i in rng.choice(len(type3_raw), 1000, replace=False)],
                   [flatten(type3_raw[i]) for i in rng.choice(len(type3_raw), 1000, replace=False)], 90),
    ]

    all_pairs = []  # list of (label, same, cross)
    for label, set_a, set_b, n_pairs in pair_configs:
        count = 0
        for _ in range(n_pairs * 5):  # oversample then take first n_pairs
            if count >= n_pairs:
                break
            i = rng.randint(len(set_a))
            j = rng.randint(len(set_b))
            same, cross = precompute_pair_blocks(set_a[i], set_b[j])
            all_pairs.append((label, same, cross))
            count += 1

    if verbose:
        labels = [p[0] for p in all_pairs]
        for lbl in dict.fromkeys(labels):
            print(f"  {lbl}: {labels.count(lbl)} pairs")
        print(f"  Total: {len(all_pairs)} pairs")

    # ------------------------------------------------------------------
    # Step 2: Z₃-symmetric sweep (u = v = 1, vary overall scale)
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 2: Z₃-symmetric case (u=v=1)...")
        print("  Product = (√8/||p||) · p where p is the trial-001 product.")

    # For u=v=1, q = same+cross = original trial-001 product
    z3_hits = 0
    z3_total = 0
    z3_by_type = Counter()
    z3_hits_by_type = Counter()
    for label, same, cross in all_pairs:
        q = np.zeros(24)
        for k in range(3):
            q[BLOCK_SLICES[k]] = same[k] + cross[k]
        norm_sq = np.dot(q, q)
        if norm_sq < 1e-12:
            continue
        z3_total += 1
        z3_by_type[label] += 1
        p_norm = q * np.sqrt(8.0 / norm_sq)
        # Check Leech membership
        p_norm_batch = p_norm.reshape(1, 24)
        if batch_is_in_leech(p_norm_batch)[0]:
            z3_hits += 1
            z3_hits_by_type[label] += 1

    if verbose:
        print(f"  Total nonzero products: {z3_total}")
        print(f"  Products on Min(Λ) after normalization: {z3_hits}")
        print(f"  Closure rate: {z3_hits}/{z3_total} = "
              f"{z3_hits/max(z3_total,1)*100:.1f}%")
        print(f"  By type:")
        for lbl in dict.fromkeys(p[0] for p in all_pairs):
            total = z3_by_type[lbl]
            hits = z3_hits_by_type[lbl]
            print(f"    {lbl}: {hits}/{total} "
                  f"({hits/max(total,1)*100:.1f}%)")

    # ------------------------------------------------------------------
    # Step 3: 2D grid search over (u, v)
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 3: 2D grid search over (u, v)...")

    # Grid in log space
    log_vals = np.linspace(-1.0, 1.0, 60)  # 60 values from 0.1 to 10
    uv_vals = 10.0 ** log_vals

    n_grid = len(uv_vals)
    n_pairs = len(all_pairs)

    if verbose:
        print(f"  Grid: {n_grid}×{n_grid} = {n_grid**2} (u,v) points")
        print(f"  Pairs: {n_pairs}")
        print(f"  Total products to check: {n_grid**2 * n_pairs:,}")

    # Precompute same/cross as arrays for vectorized access
    same_arr = np.zeros((n_pairs, 3, 8))
    cross_arr = np.zeros((n_pairs, 3, 8))
    for idx, (label, same, cross) in enumerate(all_pairs):
        for k in range(3):
            same_arr[idx, k] = same[k]
            cross_arr[idx, k] = cross[k]

    # Search
    best_hits = -1
    best_uv = (1.0, 1.0)
    results_grid = np.zeros((n_grid, n_grid), dtype=int)

    for iu, u in enumerate(uv_vals):
        for iv, v in enumerate(uv_vals):
            # Compute q for all pairs: shape (n_pairs, 24)
            q = np.zeros((n_pairs, 24))
            q[:, 0:8] = same_arr[:, 0] + cross_arr[:, 0] / (u * v)
            q[:, 8:16] = same_arr[:, 1] / u + u * cross_arr[:, 1] / v
            q[:, 16:24] = same_arr[:, 2] / v + v * cross_arr[:, 2] / u

            # Norm of each q
            norm_sq = np.sum(q * q, axis=1)

            # Skip zero products
            nonzero = norm_sq > 1e-12
            if not np.any(nonzero):
                continue

            # Normalize to ||p|| = √8
            scale = np.zeros(n_pairs)
            scale[nonzero] = np.sqrt(8.0 / norm_sq[nonzero])
            p_norm = q * scale[:, np.newaxis]  # (n_pairs, 24)

            # Only check nonzero products
            mask = nonzero
            if not np.any(mask):
                continue

            in_leech = batch_is_in_leech(p_norm[mask])
            hits = int(np.sum(in_leech))
            results_grid[iu, iv] = hits

            if hits > best_hits:
                best_hits = hits
                best_uv = (u, v)

    if verbose:
        print(f"\n  Best (u, v) = ({best_uv[0]:.4f}, {best_uv[1]:.4f}) "
              f"with {best_hits}/{n_pairs} products on Min(Λ)")

        # Report top-20 (u, v) values
        flat_idx = np.argsort(results_grid.ravel())[::-1]
        print(f"\n  Top-20 scaling parameters:")
        print(f"  {'Rank':<5} {'u':>8} {'v':>8} {'s₂/s₁':>8} {'s₃/s₁':>8} "
              f"{'Hits':>6} {'Rate':>8}")
        print(f"  {'-'*51}")
        seen_hits = set()
        reported = 0
        for flat in flat_idx:
            if reported >= 20:
                break
            iu, iv = divmod(flat, n_grid)
            h = results_grid[iu, iv]
            if h == 0:
                break
            u_val = uv_vals[iu]
            v_val = uv_vals[iv]
            reported += 1
            print(f"  {reported:<5} {u_val:8.4f} {v_val:8.4f} "
                  f"{u_val:8.4f} {v_val:8.4f} "
                  f"{h:6d} {h/n_pairs*100:7.1f}%")

    # ------------------------------------------------------------------
    # Step 4: Detailed analysis of best scaling
    # ------------------------------------------------------------------
    if verbose and best_hits > 0:
        u_best, v_best = best_uv
        print(f"\n{'='*72}")
        print(f"DETAILED ANALYSIS at (u, v) = ({u_best:.4f}, {v_best:.4f})")
        print(f"{'='*72}")

        # Per-type analysis
        type_hits = Counter()
        type_total = Counter()
        type_norms_orig = {}  # original (unscaled) product norms

        for idx, (label, same, cross) in enumerate(all_pairs):
            q = compute_q(same, cross, u_best, v_best)
            norm_sq = np.dot(q, q)
            if norm_sq < 1e-12:
                continue
            type_total[label] += 1

            if label not in type_norms_orig:
                type_norms_orig[label] = []
            type_norms_orig[label].append(norm_sq)

            p_norm = q * np.sqrt(8.0 / norm_sq)
            p_batch = p_norm.reshape(1, 24)
            if batch_is_in_leech(p_batch)[0]:
                type_hits[label] += 1

        print(f"\n  Per-type closure rates:")
        for lbl in dict.fromkeys(p[0] for p in all_pairs):
            total = type_total[lbl]
            hits = type_hits[lbl]
            norms = type_norms_orig.get(lbl, [])
            norm_vals = sorted(set(round(n, 1) for n in norms))
            print(f"    {lbl}: {hits}/{total} ({hits/max(total,1)*100:.1f}%)")
            if len(norm_vals) <= 10:
                print(f"      Unnormalized ||q||² values: {norm_vals}")

        # Failure analysis at best scaling
        print(f"\n  Failure analysis (products NOT on Min(Λ)):")
        cond_fails = Counter()
        n_fail = 0
        for idx, (label, same, cross) in enumerate(all_pairs):
            q = compute_q(same, cross, u_best, v_best)
            norm_sq = np.dot(q, q)
            if norm_sq < 1e-12:
                continue
            p_norm = q * np.sqrt(8.0 / norm_sq)
            x = ALG.element(p_norm[0:8])
            y = ALG.element(p_norm[8:16])
            z = ALG.element(p_norm[16:24])

            if not is_in_leech(x, y, z):
                n_fail += 1
                if n_fail <= 200:  # track conditions for first 200 failures
                    if not is_in_L(x): cond_fails['x∉L'] += 1
                    if not is_in_L(y): cond_fails['y∉L'] += 1
                    if not is_in_L(z): cond_fails['z∉L'] += 1
                    from leech_wilson import is_in_Ls_bar, is_in_Ls
                    if not is_in_Ls_bar(ALG.element(x.coords + y.coords)):
                        cond_fails['x+y∉Ls̄'] += 1
                    if not is_in_Ls_bar(ALG.element(x.coords + z.coords)):
                        cond_fails['x+z∉Ls̄'] += 1
                    if not is_in_Ls_bar(ALG.element(y.coords + z.coords)):
                        cond_fails['y+z∉Ls̄'] += 1
                    if not is_in_Ls(ALG.element(x.coords + y.coords + z.coords)):
                        cond_fails['x+y+z∉Ls'] += 1

        print(f"    Total failures: {n_fail}")
        if cond_fails:
            analyzed = min(n_fail, 200)
            print(f"    Condition violations (of {analyzed} analyzed):")
            for cond, count in cond_fails.most_common():
                print(f"      {cond}: {count} ({count/analyzed*100:.1f}%)")

    # ------------------------------------------------------------------
    # Step 5: Random search (finer exploration)
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n{'='*72}")
        print("RANDOM SEARCH: 50,000 random (u, v) samples")
        print(f"{'='*72}")

    n_random = 50000
    rng2 = np.random.RandomState(123)
    log_u = rng2.uniform(-1.0, 1.0, n_random)
    log_v = rng2.uniform(-1.0, 1.0, n_random)
    u_samples = 10.0 ** log_u
    v_samples = 10.0 ** log_v

    random_hits = np.zeros(n_random, dtype=int)

    for trial_idx in range(n_random):
        u = u_samples[trial_idx]
        v = v_samples[trial_idx]

        q = np.zeros((n_pairs, 24))
        q[:, 0:8] = same_arr[:, 0] + cross_arr[:, 0] / (u * v)
        q[:, 8:16] = same_arr[:, 1] / u + u * cross_arr[:, 1] / v
        q[:, 16:24] = same_arr[:, 2] / v + v * cross_arr[:, 2] / u

        norm_sq = np.sum(q * q, axis=1)
        nonzero = norm_sq > 1e-12
        if not np.any(nonzero):
            continue

        scale = np.zeros(n_pairs)
        scale[nonzero] = np.sqrt(8.0 / norm_sq[nonzero])
        p_norm = q * scale[:, np.newaxis]

        in_leech = batch_is_in_leech(p_norm[nonzero])
        random_hits[trial_idx] = int(np.sum(in_leech))

    if verbose:
        hit_dist = Counter(random_hits.tolist())
        print(f"\n  Distribution of Min(Λ) hit counts:")
        for h in sorted(hit_dist.keys()):
            print(f"    {h} hits: {hit_dist[h]} samples "
                  f"({hit_dist[h]/n_random*100:.2f}%)")

        # Top random scalings
        top_idx = np.argsort(random_hits)[::-1][:20]
        print(f"\n  Top-20 random scalings:")
        print(f"  {'Rank':<5} {'u':>8} {'v':>8} {'Hits':>6} {'Rate':>8}")
        print(f"  {'-'*37}")
        for rank, idx in enumerate(top_idx):
            h = random_hits[idx]
            if h == 0:
                break
            print(f"  {rank+1:<5} {u_samples[idx]:8.4f} {v_samples[idx]:8.4f} "
                  f"{h:6d} {h/n_pairs*100:7.1f}%")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    if verbose:
        print(f"\n{'='*72}")
        print("SUMMARY")
        print(f"{'='*72}")
        max_random = int(np.max(random_hits))
        max_grid = int(np.max(results_grid))
        overall_max = max(max_random, max_grid)
        print(f"\n  Best closure rate found: {overall_max}/{n_pairs} "
              f"= {overall_max/n_pairs*100:.1f}%")
        print(f"  Best from grid search: {max_grid}/{n_pairs}")
        print(f"  Best from random search: {max_random}/{n_pairs}")
        print(f"  Z₃-symmetric (u=v=1): {z3_hits}/{z3_total}")

        if overall_max == 0:
            print(f"\n  No scaling found that places ANY product on Min(Λ).")
            print(f"  The normalization step maps products to vectors of the")
            print(f"  correct norm (√8), but these vectors do not satisfy")
            print(f"  Wilson's lattice conditions for any tested scaling.")
        elif overall_max < n_pairs:
            print(f"\n  Partial closure: some products land on Min(Λ) but not all.")
        else:
            print(f"\n  FULL closure on tested sample!")

    return {
        'z3_hits': z3_hits,
        'z3_total': z3_total,
        'best_grid_uv': best_uv,
        'best_grid_hits': best_hits,
        'results_grid': results_grid,
        'random_hits': random_hits,
        'u_samples': u_samples,
        'v_samples': v_samples,
        'n_pairs': n_pairs,
    }


if __name__ == "__main__":
    run_trial(verbose=True)
