"""
trial_004_basis_automorphisms.py — Basis changes via E8 lattice automorphisms.

Algebra family tested
=====================
Same base multiplication rule as trial 001: three copies of the octonion
algebra on R²⁴ with Z₃-symmetric cross-block mixing.  But now we apply an
independent change of basis within each block using an automorphism of the E8
lattice L.

If T₁, T₂, T₃ are invertible 8×8 real matrices, we can define a new
multiplication on R²⁴ by:

    (a₁, a₂, a₃) ★ (b₁, b₂, b₃) = (c₁, c₂, c₃)

where each cₖ is computed by:
  1. Transform factors into the "octonion frame":  aᵢ' = Tᵢ⁻¹ aᵢ,  bⱼ' = Tⱼ⁻¹ bⱼ
  2. Compute the octonion product:  pₖ' = aᵢ' · bⱼ'  (summed over pairs)
  3. Transform back:  cₖ = Tₖ pₖ'

For a Min(Λ) vector v = (v₁, v₂, v₃), each vᵢ lives in L (the E8 lattice).
If Tᵢ is an automorphism of L (Tᵢ ∈ Aut(L)), then Tᵢ⁻¹vᵢ ∈ L, so the
factors are still valid E8 vectors.  The octonion product of two E8 vectors
is well-defined.  The question is whether Tₖ maps the result back into L
in a way that the full 24-vector satisfies Wilson's conditions.

Key insight: if T₁ = T₂ = T₃ = T for some T ∈ Aut(O) ∩ Aut(L) (an
octonion automorphism that also preserves L), then ★ is identical to the
original product (conjugation by an algebra automorphism doesn't change the
algebra).  So only automorphisms OUTSIDE the octonion automorphism group G₂
give genuinely different algebras.

Strategy
========
The full Aut(L) = W(E₈) has order 696,729,600 — far too large to sweep.
Instead, we sample a manageable set of E8 automorphisms:

  1. SIGN CHANGES: flipping signs of individual coordinates.  These are
     diagonal matrices with ±1 entries.  Must preserve L, so the number of
     sign flips must be even (to preserve sum parity).  This gives 2⁷ = 128
     sign patterns (we fix the sign of e₀ to reduce redundancy).

  2. COORDINATE PERMUTATIONS: permuting the 8 coordinates.  Not all
     permutations preserve L; we need the sum-parity to be maintained.  All
     permutations preserve the "all-integer, even sum" condition.  For the
     half-integer vectors, permutations also preserve "all half-integer, odd
     sum."  So ALL permutations of {0,...,7} are L-automorphisms.  S₈ has
     order 40,320.

  3. COMBINED: sign changes × permutations.  This subgroup has order up to
     128 × 40,320 = 5,160,960.  Still large.

We sample 500 random automorphisms from this family and apply each one
to ONE block (T₁ = T, T₂ = T₃ = I), testing whether the modified product
closes on the Leech lattice.  We also test a few cases with different
automorphisms on different blocks.

References
==========
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186–2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
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

# ---------------------------------------------------------------------------
# Block layout
# ---------------------------------------------------------------------------

BLOCK_SLICES = [slice(0, 8), slice(8, 16), slice(16, 24)]


# ---------------------------------------------------------------------------
# E8 lattice automorphisms: sign changes + permutations
# ---------------------------------------------------------------------------

def random_sign_change(rng):
    """
    Random sign-flip matrix preserving L.

    Flip an EVEN number of coordinates (to preserve the parity of the sum).
    Returns an 8×8 diagonal matrix with ±1 entries.
    """
    signs = np.ones(8)
    # Flip each coordinate with probability 0.5, then fix parity
    flips = rng.randint(0, 2, size=8)
    signs[flips == 1] = -1.0
    # Ensure even number of flips
    if np.sum(flips) % 2 != 0:
        # Flip one more random coordinate
        idx = rng.randint(0, 8)
        signs[idx] *= -1.0
    return np.diag(signs)


def random_permutation_matrix(rng):
    """Random permutation of the 8 coordinates.  Always preserves L."""
    perm = rng.permutation(8)
    P = np.zeros((8, 8))
    for i, j in enumerate(perm):
        P[i, j] = 1.0
    return P


def random_e8_automorphism(rng):
    """Random E8 automorphism: permutation × sign change."""
    P = random_permutation_matrix(rng)
    S = random_sign_change(rng)
    return P @ S


def verify_preserves_L(T, sample_vectors):
    """Verify that T maps L-vectors to L-vectors on a sample."""
    alg = STANDARD_ALGEBRA
    for v in sample_vectors[:20]:
        Tv = T @ v
        if not is_in_L(alg.element(Tv)):
            return False
    return True


# ---------------------------------------------------------------------------
# Multiplication with basis automorphisms
# ---------------------------------------------------------------------------

def multiply_24_with_auts(a, b, T_list, T_inv_list):
    """
    Multiply two R²⁴ vectors using the triple-octonion rule with per-block
    basis automorphisms.

    T_list[k] is the 8×8 matrix for block k.
    T_inv_list[k] = T_list[k]⁻¹.

    The product is:
      For blocks i, j with target block t:
        a_i' = T_inv[i] @ a_block_i
        b_j' = T_inv[j] @ b_block_j
        p'   = octonion_product(a_i', b_j')
        contribution to block t += T[t] @ p'
    """
    alg = STANDARD_ALGEBRA
    result = np.zeros(24, dtype=np.float64)

    a_blocks = [a[s] for s in BLOCK_SLICES]
    b_blocks = [b[s] for s in BLOCK_SLICES]

    for bi in range(3):
        for bj in range(3):
            # Target block
            if bi == bj:
                bt = bi
            else:
                bt = 3 - bi - bj

            # Transform to octonion frame
            ai_oct = T_inv_list[bi] @ a_blocks[bi]
            bj_oct = T_inv_list[bj] @ b_blocks[bj]

            # Octonion product
            prod_oct = alg._mul_coords(ai_oct, bj_oct)

            # Transform back
            result[BLOCK_SLICES[bt]] += T_list[bt] @ prod_oct

    return result


# ---------------------------------------------------------------------------
# Leech membership check
# ---------------------------------------------------------------------------

def check_in_leech(v):
    """Check if a 24-vector is in the Leech lattice."""
    alg = STANDARD_ALGEBRA
    x = alg.element(v[0:8])
    y = alg.element(v[8:16])
    z = alg.element(v[16:24])
    return is_in_leech(x, y, z)


def flatten_triple(trip):
    return np.concatenate(trip)


# ---------------------------------------------------------------------------
# Main trial
# ---------------------------------------------------------------------------

def run_trial(n_automorphisms=500, verbose=True):
    """
    Execute trial 004: basis automorphism search.
    """
    if verbose:
        print("=" * 72)
        print("TRIAL 004: E8 automorphism basis changes")
        print("=" * 72)

    rng = np.random.RandomState(123)
    alg = STANDARD_ALGEBRA

    # ------------------------------------------------------------------
    # Step 1: Generate test vectors
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 1: Generating Min(Λ) vectors...")

    type1 = [flatten_triple(t) for t in leech_type1_vectors()]
    type2 = [flatten_triple(t) for t in leech_type2_vectors()]
    type3_raw = leech_type3_vectors()
    type3_idx = rng.choice(len(type3_raw), 500, replace=False)
    type3 = [flatten_triple(type3_raw[i]) for i in type3_idx]

    if verbose:
        print(f"  Type 1: {len(type1)}, Type 2: {len(type2)}, "
              f"Type 3: {len(type3)} (sampled)")

    # ------------------------------------------------------------------
    # Step 2: Build test pairs (focused on the critical type3×type3)
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 2: Building test pairs...")

    test_pairs = []
    pair_labels = []

    # type3 × type3: 80 pairs (the critical case)
    for i in range(0, 160, 2):
        test_pairs.append((type3[i], type3[i+1]))
        pair_labels.append("t3×t3")

    # type2 × type2: 30 pairs
    idx2 = rng.choice(len(type2), 60, replace=False)
    for i in range(0, 60, 2):
        test_pairs.append((type2[idx2[i]], type2[idx2[i+1]]))
        pair_labels.append("t2×t2")

    # type1 × type1: 20 pairs
    idx1 = rng.choice(len(type1), 40, replace=False)
    for i in range(0, 40, 2):
        test_pairs.append((type1[idx1[i]], type1[idx1[i+1]]))
        pair_labels.append("t1×t1")

    # type2 × type3: 20 pairs
    for i in range(20):
        test_pairs.append((type2[i], type3[i]))
        pair_labels.append("t2×t3")

    n_pairs = len(test_pairs)
    if verbose:
        tc = Counter(pair_labels)
        print(f"  Total: {n_pairs} pairs ({dict(tc)})")

    # Collect sample E8 vectors for verification
    sample_e8 = [type1[0][0:8], type1[1][0:8], type2[0][0:8]]

    # ------------------------------------------------------------------
    # Step 3: Identity baseline
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 3: Baseline (identity automorphism = trial 001)...")

    I8 = np.eye(8)
    T_id = [I8, I8, I8]

    n_ok_baseline = 0
    n_nz_baseline = 0
    for (a, b), pl in zip(test_pairs, pair_labels):
        prod = multiply_24_with_auts(a, b, T_id, T_id)
        ns = float(np.dot(prod, prod))
        if ns < 1e-12:
            continue
        n_nz_baseline += 1
        if check_in_leech(prod):
            n_ok_baseline += 1

    if verbose:
        print(f"  Baseline: {n_ok_baseline}/{n_nz_baseline} in Λ "
              f"({n_ok_baseline/max(1,n_nz_baseline)*100:.1f}%)")

    # ------------------------------------------------------------------
    # Step 4a: Single-block automorphisms (T on block 0, I on blocks 1,2)
    # ------------------------------------------------------------------
    if verbose:
        print(f"\nStep 4a: Single-block automorphisms "
              f"(T on block 0, I on blocks 1,2)...")
        print(f"  Testing {n_automorphisms} random E8 automorphisms...")

    best_rate = 0.0
    best_results = []
    rate_histogram = Counter()

    for trial_i in range(n_automorphisms):
        T = random_e8_automorphism(rng)
        T_inv = np.linalg.inv(T)

        # Quick check: does T preserve L?
        if not verify_preserves_L(T, sample_e8):
            continue

        T_list = [T, I8, I8]
        T_inv_list = [T_inv, I8, I8]

        n_ok = 0
        n_nz = 0
        fail_by_type = Counter()

        for (a, b), pl in zip(test_pairs, pair_labels):
            prod = multiply_24_with_auts(a, b, T_list, T_inv_list)
            ns = float(np.dot(prod, prod))
            if ns < 1e-12:
                continue
            n_nz += 1
            if check_in_leech(prod):
                n_ok += 1
            else:
                fail_by_type[pl] += 1

        rate = n_ok / max(1, n_nz)
        rate_histogram[round(rate, 2)] += 1

        if rate > best_rate:
            best_rate = rate
            best_results = [(T.copy(), n_ok, n_nz, dict(fail_by_type))]
        elif rate == best_rate:
            best_results.append((T.copy(), n_ok, n_nz, dict(fail_by_type)))

        if verbose and (trial_i + 1) % 100 == 0:
            print(f"  ... {trial_i+1}/{n_automorphisms} done "
                  f"(best rate: {best_rate*100:.1f}%)")

    if verbose:
        print(f"\n  Best single-block rate: {best_rate*100:.1f}%")
        print(f"  Rate distribution:")
        for r, c in sorted(rate_histogram.items(), reverse=True)[:10]:
            print(f"    {r*100:5.1f}%: {c} automorphisms")

        if best_results:
            T_best, n_ok, n_nz, fbt = best_results[0]
            print(f"  Best result: {n_ok}/{n_nz} in Λ, failures: {fbt}")

    # ------------------------------------------------------------------
    # Step 4b: Independent automorphisms on all three blocks
    # ------------------------------------------------------------------
    if verbose:
        print(f"\nStep 4b: Independent automorphisms on all three blocks...")
        print(f"  Testing {n_automorphisms} random (T₁, T₂, T₃) triples...")

    best_rate_3 = 0.0
    best_results_3 = []
    rate_histogram_3 = Counter()

    for trial_i in range(n_automorphisms):
        Ts = [random_e8_automorphism(rng) for _ in range(3)]
        T_invs = [np.linalg.inv(T) for T in Ts]

        # Verify all preserve L
        ok = True
        for T in Ts:
            if not verify_preserves_L(T, sample_e8):
                ok = False
                break
        if not ok:
            continue

        n_ok = 0
        n_nz = 0
        fail_by_type = Counter()

        for (a, b), pl in zip(test_pairs, pair_labels):
            prod = multiply_24_with_auts(a, b, Ts, T_invs)
            ns = float(np.dot(prod, prod))
            if ns < 1e-12:
                continue
            n_nz += 1
            if check_in_leech(prod):
                n_ok += 1
            else:
                fail_by_type[pl] += 1

        rate = n_ok / max(1, n_nz)
        rate_histogram_3[round(rate, 2)] += 1

        if rate > best_rate_3:
            best_rate_3 = rate
            best_results_3 = [(n_ok, n_nz, dict(fail_by_type))]
        elif rate == best_rate_3:
            best_results_3.append((n_ok, n_nz, dict(fail_by_type)))

        if verbose and (trial_i + 1) % 100 == 0:
            print(f"  ... {trial_i+1}/{n_automorphisms} done "
                  f"(best rate: {best_rate_3*100:.1f}%)")

    if verbose:
        print(f"\n  Best three-block rate: {best_rate_3*100:.1f}%")
        print(f"  Rate distribution:")
        for r, c in sorted(rate_histogram_3.items(), reverse=True)[:10]:
            print(f"    {r*100:5.1f}%: {c} automorphisms")

        if best_results_3:
            n_ok, n_nz, fbt = best_results_3[0]
            print(f"  Best result: {n_ok}/{n_nz} in Λ, failures: {fbt}")

    # ------------------------------------------------------------------
    # Step 4c: Sign-only automorphisms (exhaustive for single block)
    # ------------------------------------------------------------------
    if verbose:
        print(f"\nStep 4c: Sign-only automorphisms on block 0 (exhaustive)...")

    # All even-parity sign patterns: 2⁷ = 128
    best_rate_sign = 0.0
    n_sign_tested = 0

    for bits in range(128):
        signs = np.ones(8)
        flip_count = 0
        for k in range(7):
            if bits & (1 << k):
                signs[k + 1] = -1.0  # flip coordinates 1–7
                flip_count += 1
        # Fix parity: if odd number of flips, also flip coordinate 0
        if flip_count % 2 != 0:
            signs[0] = -1.0

        T = np.diag(signs)
        T_inv = np.diag(signs)  # self-inverse for sign matrices

        T_list = [T, I8, I8]
        T_inv_list = [T_inv, I8, I8]

        n_ok = 0
        n_nz = 0

        for (a, b), pl in zip(test_pairs, pair_labels):
            prod = multiply_24_with_auts(a, b, T_list, T_inv_list)
            ns = float(np.dot(prod, prod))
            if ns < 1e-12:
                continue
            n_nz += 1
            if check_in_leech(prod):
                n_ok += 1

        rate = n_ok / max(1, n_nz)
        if rate > best_rate_sign:
            best_rate_sign = rate
        n_sign_tested += 1

    if verbose:
        print(f"  Tested all {n_sign_tested} even-parity sign patterns")
        print(f"  Best sign-only rate: {best_rate_sign*100:.1f}%")

    # ------------------------------------------------------------------
    # Step 5: Summary
    # ------------------------------------------------------------------
    if verbose:
        print("\n" + "=" * 72)
        print("SUMMARY")
        print("=" * 72)

        overall_best = max(best_rate, best_rate_3, best_rate_sign,
                           n_ok_baseline / max(1, n_nz_baseline))

        print(f"\n  Baseline (identity):          "
              f"{n_ok_baseline/max(1,n_nz_baseline)*100:.1f}%")
        print(f"  Best single-block random:     {best_rate*100:.1f}%")
        print(f"  Best three-block random:      {best_rate_3*100:.1f}%")
        print(f"  Best sign-only (exhaustive):  {best_rate_sign*100:.1f}%")

        if overall_best >= 1.0:
            print("\n✓ Found basis change(s) with 100% closure! "
                  "Further testing needed.")
        else:
            print(f"\n✗ No E8 automorphism basis change achieves closure.")
            print(f"  Best overall rate: {overall_best*100:.1f}%")
            print(f"  This provides strong evidence that the triple-octonion")
            print(f"  algebra cannot be rescued by changing the basis within")
            print(f"  each E8 block.")

    return {
        'baseline_rate': n_ok_baseline / max(1, n_nz_baseline),
        'best_single_block': best_rate,
        'best_three_block': best_rate_3,
        'best_sign_only': best_rate_sign,
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_trial(n_automorphisms=500, verbose=True)
