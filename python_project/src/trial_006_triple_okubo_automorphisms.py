"""
trial_006_triple_okubo_automorphisms.py — E8 automorphism basis changes
for the triple symmetric-composition algebra.

This is the automorphism companion to trial 005, analogous to how trial 004
was the automorphism companion to trial 003.

Trial 005 showed that the triple Okubo/para-octonion algebra fails
catastrophically: 91.5% of products leave the E8 lattice entirely (Wilson
condition 1 fails).  This is because the Petersson construction introduces
irrational structure constants (√3/2 from the 2π/3 rotation in τ).

A priori, an E8 automorphism basis change COULD absorb some of these
irrational factors if the automorphism happens to align with the rotation
planes of τ.  This trial checks whether any such alignment exists.

Expected outcome: FAIL.  The √3 factors are intrinsic to the Petersson
construction and cannot be eliminated by a change of basis that preserves
the E8 lattice.

References
==========
Same as trial 005.
"""

import numpy as np
from collections import Counter

from octonions import STANDARD_ALGEBRA
from leech_wilson import (
    leech_type1_vectors,
    leech_type2_vectors,
    leech_type3_vectors,
    is_in_leech,
)
from e8_wilson import is_in_L

from trial_005_triple_okubo import (
    build_three_algebras,
    flatten_triple,
    BLOCK_SLICES,
    target_block,
    OkuboAlgebra,
)

# ---------------------------------------------------------------------------
# Multiplication with per-block basis automorphisms
# ---------------------------------------------------------------------------

def multiply_24_with_auts(a, b, algebras, T_list, T_inv_list):
    """
    Like trial_005.multiply_24 but with per-block basis automorphisms.

    For each block pair (i,j) → t:
      1. Transform to "algebra frame": a'_i = T_inv[i] @ a_i, b'_j = T_inv[j] @ b_j
      2. Compute the 8-dim product: p' = algebras[t]._mul_coords(a'_i, b'_j)
      3. Transform back: contribution to block t += T[t] @ p'
    """
    result = np.zeros(24, dtype=np.float64)
    a_blocks = [a[s] for s in BLOCK_SLICES]
    b_blocks = [b[s] for s in BLOCK_SLICES]

    for bi in range(3):
        for bj in range(3):
            bt = target_block(bi, bj)
            ai_frame = T_inv_list[bi] @ a_blocks[bi]
            bj_frame = T_inv_list[bj] @ b_blocks[bj]
            prod_frame = algebras[bt]._mul_coords(ai_frame, bj_frame)
            result[BLOCK_SLICES[bt]] += T_list[bt] @ prod_frame

    return result


def check_in_leech(v):
    alg = STANDARD_ALGEBRA
    x = alg.element(v[0:8])
    y = alg.element(v[8:16])
    z = alg.element(v[16:24])
    return is_in_leech(x, y, z)


# ---------------------------------------------------------------------------
# Random E8 automorphisms (same as trial 004)
# ---------------------------------------------------------------------------

def random_sign_change(rng):
    signs = np.ones(8)
    flips = rng.randint(0, 2, size=8)
    signs[flips == 1] = -1.0
    if np.sum(flips) % 2 != 0:
        signs[rng.randint(0, 8)] *= -1.0
    return np.diag(signs)


def random_permutation_matrix(rng):
    perm = rng.permutation(8)
    P = np.zeros((8, 8))
    for i, j in enumerate(perm):
        P[i, j] = 1.0
    return P


def random_e8_automorphism(rng):
    return random_permutation_matrix(rng) @ random_sign_change(rng)


# ---------------------------------------------------------------------------
# Main trial
# ---------------------------------------------------------------------------

def run_trial(n_automorphisms=300, verbose=True):
    if verbose:
        print("=" * 72)
        print("TRIAL 006: E8 automorphism basis changes for triple Okubo")
        print("=" * 72)

    rng = np.random.RandomState(456)
    algebras = list(build_three_algebras())
    I8 = np.eye(8)

    # Generate test vectors
    if verbose:
        print("\nStep 1: Generating Min(Λ) vectors...")

    type1 = [flatten_triple(t) for t in leech_type1_vectors()]
    type2 = [flatten_triple(t) for t in leech_type2_vectors()]
    type3_raw = leech_type3_vectors()
    type3_idx = rng.choice(len(type3_raw), 500, replace=False)
    type3 = [flatten_triple(type3_raw[i]) for i in type3_idx]

    # Build test pairs
    test_pairs = []
    pair_labels = []
    for i in range(0, 100, 2):
        test_pairs.append((type3[i], type3[i+1]))
        pair_labels.append("t3×t3")
    idx2 = rng.choice(len(type2), 40, replace=False)
    for i in range(0, 40, 2):
        test_pairs.append((type2[idx2[i]], type2[idx2[i+1]]))
        pair_labels.append("t2×t2")
    for i in range(10):
        test_pairs.append((type1[i], type2[i]))
        pair_labels.append("t1×t2")

    n_pairs = len(test_pairs)
    if verbose:
        tc = Counter(pair_labels)
        print(f"  Test pairs: {n_pairs} ({dict(tc)})")

    # Baseline
    if verbose:
        print("\nStep 2: Baseline (identity)...")
    n_ok_base = 0
    n_nz_base = 0
    for (a, b), pl in zip(test_pairs, pair_labels):
        from trial_005_triple_okubo import multiply_24
        prod = multiply_24(a, b, algebras)
        ns = float(np.dot(prod, prod))
        if ns < 1e-12:
            continue
        n_nz_base += 1
        if check_in_leech(prod):
            n_ok_base += 1
    base_rate = n_ok_base / max(1, n_nz_base)
    if verbose:
        print(f"  Baseline: {n_ok_base}/{n_nz_base} ({base_rate*100:.1f}%)")

    # Single-block automorphisms
    if verbose:
        print(f"\nStep 3: Single-block automorphisms ({n_automorphisms} random)...")

    best_rate = 0.0
    rate_hist = Counter()

    for trial_i in range(n_automorphisms):
        T = random_e8_automorphism(rng)
        T_inv = np.linalg.inv(T)
        T_list = [T, I8, I8]
        T_inv_list = [T_inv, I8, I8]

        n_ok = 0
        n_nz = 0
        for (a, b), pl in zip(test_pairs, pair_labels):
            prod = multiply_24_with_auts(a, b, algebras, T_list, T_inv_list)
            ns = float(np.dot(prod, prod))
            if ns < 1e-12:
                continue
            n_nz += 1
            if check_in_leech(prod):
                n_ok += 1

        rate = n_ok / max(1, n_nz)
        rate_hist[round(rate, 2)] += 1
        if rate > best_rate:
            best_rate = rate

        if verbose and (trial_i + 1) % 100 == 0:
            print(f"  ... {trial_i+1}/{n_automorphisms} (best: {best_rate*100:.1f}%)")

    if verbose:
        print(f"\n  Best single-block rate: {best_rate*100:.1f}%")

    # Three-block automorphisms
    if verbose:
        print(f"\nStep 4: Three-block automorphisms ({n_automorphisms} random)...")

    best_rate_3 = 0.0

    for trial_i in range(n_automorphisms):
        Ts = [random_e8_automorphism(rng) for _ in range(3)]
        T_invs = [np.linalg.inv(T) for T in Ts]

        n_ok = 0
        n_nz = 0
        for (a, b), pl in zip(test_pairs, pair_labels):
            prod = multiply_24_with_auts(a, b, algebras, Ts, T_invs)
            ns = float(np.dot(prod, prod))
            if ns < 1e-12:
                continue
            n_nz += 1
            if check_in_leech(prod):
                n_ok += 1

        rate = n_ok / max(1, n_nz)
        if rate > best_rate_3:
            best_rate_3 = rate

        if verbose and (trial_i + 1) % 100 == 0:
            print(f"  ... {trial_i+1}/{n_automorphisms} (best: {best_rate_3*100:.1f}%)")

    if verbose:
        print(f"\n  Best three-block rate: {best_rate_3*100:.1f}%")

    # Summary
    if verbose:
        print("\n" + "=" * 72)
        print("SUMMARY")
        print("=" * 72)
        overall = max(base_rate, best_rate, best_rate_3)
        print(f"\n  Baseline:         {base_rate*100:.1f}%")
        print(f"  Best single-block: {best_rate*100:.1f}%")
        print(f"  Best three-block:  {best_rate_3*100:.1f}%")

        if overall < 1.0:
            print(f"\n  ✗ No E8 automorphism rescues the triple Okubo algebra.")
            print(f"    The irrational structure constants (√3 from Petersson's")
            print(f"    2π/3 rotation) cannot be absorbed by E8 automorphisms.")

    return {
        'baseline': base_rate,
        'best_single': best_rate,
        'best_three': best_rate_3,
    }


if __name__ == "__main__":
    run_trial(n_automorphisms=300, verbose=True)
