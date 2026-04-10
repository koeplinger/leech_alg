"""
trial_003_discrete_variants.py — Exhaustive sweep over discrete modifications
of the triple-octonion algebra.

Algebra family tested
=====================
Same base structure as trials 001–002: three copies of the octonion algebra on
R²⁴ (O₁ at indices 0–7, O₂ at 8–15, O₃ at 16–23).  But now we vary THREE
discrete degrees of freedom in the cross-block terms:

  (A) CONJUGATION — what operation is applied to the left/right factor before
      computing the octonion product in a cross-block term.

      For each cross-block pair (Oα, Oβ) → Oγ, the product could be any of:
        a · b       (plain)
        ā · b       (conjugate left)
        a · b̄       (conjugate right)
        ā · b̄       (conjugate both)

      Octonion conjugation: ā = (a₀, -a₁, -a₂, ..., -a₇).  It swaps the
      relationship between Ls and Ls̄ sublattices — directly relevant since
      trial 001 failed exclusively on Wilson condition 3 (x+y+z ∈ Ls).

  (B) SIGN — multiply each cross-block contribution by +1 or -1.

      For each of the 3 cross-block pairs, we can flip the sign of the
      entire contribution.

  (C) ROUTING — which block receives the cross-block product.

      Trial 001 uses "third-block" routing: Oα × Oβ → Oγ where {α,β,γ} = {0,1,2}.
      Alternatives:
        "left"   routing: Oα × Oβ → Oα  (product stays in left factor's block)
        "right"  routing: Oα × Oβ → Oβ  (product stays in right factor's block)

      These break Z₃ symmetry but are legitimate algebraic choices.

Parameter space
===============
  Conjugation: 4 choices per cross-block pair, 3 pairs → 4³ = 64
  Sign:        2 choices per pair → 2³ = 8
  Routing:     3 choices (third, left, right)

  Total: 64 × 8 × 3 = 1,536 variants.

  We can reduce by Z₃ symmetry for the "third-block" routing (but not for
  left/right routing).  However, 1,536 is small enough to sweep exhaustively,
  so we test all of them without symmetry reduction.

For each variant, we test a SAMPLE of Min(Λ) pairs and count how many products
land in Λ.  Variants with 100% closure on the sample are tested more
thoroughly.

Motivation
==========
Trial 001 showed that the triple-octonion algebra fails exclusively on
condition 3 (x+y+z ∈ Ls) for type3 × type3 products.  Conjugation directly
modifies the Ls/Ls̄ structure.  Sign flips and routing changes alter the
interference pattern in the three-block sum.  Together, these cover all "simple
discrete tweaks" to the cross-block rule.

If ALL 1,536 variants fail, this rules out the entire family of
"three-octonion algebras with simple cross-block modifications" as candidates
for an order on Λ.

References
==========
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186–2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
"""

import numpy as np
from collections import Counter
from itertools import product as iterproduct

from octonions import STANDARD_ALGEBRA
from leech_wilson import (
    leech_type1_vectors,
    leech_type2_vectors,
    leech_type3_vectors,
    is_in_leech,
    is_in_Ls_bar,
    is_in_Ls,
)
from e8_wilson import is_in_L

# ---------------------------------------------------------------------------
# Block layout (same as trial 001)
# ---------------------------------------------------------------------------

BLOCK_SLICES = [slice(0, 8), slice(8, 16), slice(16, 24)]


def conjugate_octonion(v):
    """Octonion conjugation: (a₀, a₁, ..., a₇) → (a₀, -a₁, ..., -a₇)."""
    c = v.copy()
    c[1:] = -c[1:]
    return c


# ---------------------------------------------------------------------------
# Multiplication with discrete parameters
# ---------------------------------------------------------------------------

# Conjugation codes: 0 = plain, 1 = conj left, 2 = conj right, 3 = conj both
CONJ_LABELS = ["a·b", "ā·b", "a·b̄", "ā·b̄"]

# Routing codes: 0 = third-block, 1 = left-block, 2 = right-block
ROUTING_LABELS = ["third", "left", "right"]


def _apply_conj(v, code):
    """Apply conjugation according to code (0=none, 1=conjugate)."""
    if code:
        return conjugate_octonion(v)
    return v


def multiply_24_variant(a, b, conj_codes, sign_codes, routing):
    """
    Multiply two R²⁴ vectors using a variant triple-octonion rule.

    Parameters
    ----------
    a, b : ndarray of shape (24,)
    conj_codes : list of 3 tuples (conj_left, conj_right)
        For each cross-block pair index k (0,1,2), whether to conjugate
        the left and/or right factor.  Each element is 0 or 1.
    sign_codes : list of 3 ints (+1 or -1)
        Sign multiplier for each cross-block pair.
    routing : int
        0 = third-block, 1 = left-block, 2 = right-block.

    The 3 cross-block pairs are indexed as:
        pair 0: (O₀, O₁)  — blocks 0 and 1
        pair 1: (O₀, O₂)  — blocks 0 and 2
        pair 2: (O₁, O₂)  — blocks 1 and 2
    """
    alg = STANDARD_ALGEBRA
    result = np.zeros(24, dtype=np.float64)

    a_blocks = [a[s].copy() for s in BLOCK_SLICES]
    b_blocks = [b[s].copy() for s in BLOCK_SLICES]

    # Same-block products: always Oα × Oα → Oα, no modifications
    for bi in range(3):
        product = alg._mul_coords(a_blocks[bi], b_blocks[bi])
        result[BLOCK_SLICES[bi]] += product

    # Cross-block products
    # The 6 ordered cross-block pairs, grouped by their unordered pair index:
    cross_pairs = [
        # pair_index, left_block, right_block
        (0, 0, 1),  # O₀ × O₁
        (0, 1, 0),  # O₁ × O₀
        (1, 0, 2),  # O₀ × O₂
        (1, 2, 0),  # O₂ × O₀
        (2, 1, 2),  # O₁ × O₂
        (2, 2, 1),  # O₂ × O₁
    ]

    for pair_idx, left_blk, right_blk in cross_pairs:
        conj_left, conj_right = conj_codes[pair_idx]

        left_vec = _apply_conj(a_blocks[left_blk], conj_left)
        right_vec = _apply_conj(b_blocks[right_blk], conj_right)

        product = alg._mul_coords(left_vec, right_vec)
        product = sign_codes[pair_idx] * product

        # Determine target block
        if routing == 0:  # third-block
            target = 3 - left_blk - right_blk
        elif routing == 1:  # left-block
            target = left_blk
        else:  # right-block
            target = right_blk

        result[BLOCK_SLICES[target]] += product

    return result


# ---------------------------------------------------------------------------
# Leech membership check for 24-vectors
# ---------------------------------------------------------------------------

def check_in_leech(v):
    """Check if a 24-vector is in the Leech lattice."""
    alg = STANDARD_ALGEBRA
    x = alg.element(v[0:8])
    y = alg.element(v[8:16])
    z = alg.element(v[16:24])
    return is_in_leech(x, y, z)


def flatten_triple(trip):
    """Convert (x, y, z) triple to 24-vector."""
    return np.concatenate(trip)


# ---------------------------------------------------------------------------
# Enumerate all 1,536 variants
# ---------------------------------------------------------------------------

def enumerate_variants():
    """
    Yield all (conj_codes, sign_codes, routing, label) tuples.

    conj_codes: list of 3 tuples, each (conj_left, conj_right) ∈ {0,1}²
    sign_codes: list of 3 ints ∈ {+1, -1}
    routing: int ∈ {0, 1, 2}
    label: human-readable string
    """
    # Conjugation: for each of 3 unordered cross-block pairs, 4 choices
    conj_options = [(0, 0), (1, 0), (0, 1), (1, 1)]

    for routing in range(3):
        for c0 in conj_options:
            for c1 in conj_options:
                for c2 in conj_options:
                    for s0 in (+1, -1):
                        for s1 in (+1, -1):
                            for s2 in (+1, -1):
                                conj_codes = [c0, c1, c2]
                                sign_codes = [s0, s1, s2]

                                conj_strs = []
                                for ci, (cl, cr) in enumerate(conj_codes):
                                    idx = cl * 2 + cr
                                    conj_strs.append(CONJ_LABELS[idx])

                                sign_str = "".join("+" if s > 0 else "-"
                                                   for s in sign_codes)
                                label = (f"route={ROUTING_LABELS[routing]}, "
                                         f"conj=[{','.join(conj_strs)}], "
                                         f"sign={sign_str}")

                                yield conj_codes, sign_codes, routing, label


# ---------------------------------------------------------------------------
# Main trial
# ---------------------------------------------------------------------------

def run_trial(verbose=True):
    """
    Execute trial 003: exhaustive sweep over all 1,536 discrete variants.
    """
    if verbose:
        print("=" * 72)
        print("TRIAL 003: Discrete variants of triple-octonion cross-block rule")
        print("=" * 72)

    # ------------------------------------------------------------------
    # Step 1: Generate test vectors
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 1: Generating Min(Λ) vectors...")

    type1 = [flatten_triple(t) for t in leech_type1_vectors()]
    type2 = [flatten_triple(t) for t in leech_type2_vectors()]
    type3_raw = leech_type3_vectors()

    rng = np.random.RandomState(42)
    type3_idx = rng.choice(len(type3_raw), 500, replace=False)
    type3 = [flatten_triple(type3_raw[i]) for i in type3_idx]

    if verbose:
        print(f"  Type 1: {len(type1)}")
        print(f"  Type 2: {len(type2)}")
        print(f"  Type 3: {len(type3)} (sampled from {len(type3_raw)})")

    # ------------------------------------------------------------------
    # Step 2: Build test pairs
    # ------------------------------------------------------------------
    # We use a modest sample that covers all type combinations.
    # The critical combination is type3×type3 (the only one that failed
    # in trial 001), but we test all to catch new failure modes.
    if verbose:
        print("\nStep 2: Building test pairs...")

    test_pairs = []
    pair_labels = []

    # type1 × type1: 50 pairs
    idx1 = rng.choice(len(type1), min(50, len(type1)), replace=False)
    for i in range(0, len(idx1) - 1, 2):
        test_pairs.append((type1[idx1[i]], type1[idx1[i+1]]))
        pair_labels.append("t1×t1")

    # type2 × type2: 50 pairs
    idx2 = rng.choice(len(type2), min(100, len(type2)), replace=False)
    for i in range(0, min(100, len(idx2)) - 1, 2):
        test_pairs.append((type2[idx2[i]], type2[idx2[i+1]]))
        pair_labels.append("t2×t2")

    # type3 × type3: 100 pairs (the critical case)
    for i in range(0, min(200, len(type3)) - 1, 2):
        test_pairs.append((type3[i], type3[i+1]))
        pair_labels.append("t3×t3")

    # type1 × type2: 25 pairs
    for i in range(25):
        test_pairs.append((type1[i % len(type1)], type2[i % len(type2)]))
        pair_labels.append("t1×t2")

    # type1 × type3: 25 pairs
    for i in range(25):
        test_pairs.append((type1[i % len(type1)], type3[i % len(type3)]))
        pair_labels.append("t1×t3")

    # type2 × type3: 25 pairs
    for i in range(25):
        test_pairs.append((type2[i % len(type2)], type3[i % len(type3)]))
        pair_labels.append("t2×t3")

    n_pairs = len(test_pairs)
    if verbose:
        type_counts = Counter(pair_labels)
        print(f"  Total test pairs: {n_pairs}")
        for k, v in sorted(type_counts.items()):
            print(f"    {k}: {v}")

    # ------------------------------------------------------------------
    # Step 3: Sweep all variants
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 3: Sweeping all 1,536 variants...")
        print("  (testing each variant against all pairs)\n")

    variant_results = []  # (label, n_in_leech, n_nonzero, n_fail, fail_by_type)
    best_closure_rate = 0.0
    best_variants = []
    n_variants = 0

    for conj_codes, sign_codes, routing, label in enumerate_variants():
        n_in_leech = 0
        n_nonzero = 0
        n_fail = 0
        fail_by_type = Counter()

        for (a, b), plabel in zip(test_pairs, pair_labels):
            prod = multiply_24_variant(a, b, conj_codes, sign_codes, routing)
            norm_sq = float(np.dot(prod, prod))

            if norm_sq < 1e-12:
                continue

            n_nonzero += 1
            if check_in_leech(prod):
                n_in_leech += 1
            else:
                n_fail += 1
                fail_by_type[plabel] += 1

        if n_nonzero > 0:
            closure_rate = n_in_leech / n_nonzero
        else:
            closure_rate = 0.0

        variant_results.append((label, n_in_leech, n_nonzero, n_fail,
                                dict(fail_by_type), closure_rate))

        if closure_rate > best_closure_rate:
            best_closure_rate = closure_rate
            best_variants = [(label, n_in_leech, n_nonzero, n_fail,
                              dict(fail_by_type))]
        elif closure_rate == best_closure_rate and closure_rate > 0:
            best_variants.append((label, n_in_leech, n_nonzero, n_fail,
                                  dict(fail_by_type)))

        n_variants += 1
        if verbose and n_variants % 200 == 0:
            print(f"  ... tested {n_variants}/1536 variants "
                  f"(best closure so far: {best_closure_rate*100:.1f}%)")

    # ------------------------------------------------------------------
    # Step 4: Results
    # ------------------------------------------------------------------
    if verbose:
        print("\n" + "=" * 72)
        print("RESULTS")
        print("=" * 72)

        print(f"\nTotal variants tested: {n_variants}")
        print(f"Test pairs per variant: {n_pairs}")

        # Distribution of closure rates
        rates = [r[5] for r in variant_results]
        rate_bins = Counter(round(r, 3) for r in rates)

        print(f"\nClosure rate distribution:")
        for rate_val, count in sorted(rate_bins.items(), reverse=True)[:20]:
            print(f"  {rate_val*100:6.1f}%: {count} variants")

        # Best variants
        print(f"\nBest closure rate: {best_closure_rate*100:.1f}%")
        print(f"Number of variants achieving best rate: {len(best_variants)}")

        if best_closure_rate < 1.0:
            print("\n⚠ NO variant achieves 100% closure on the sample.")

        print(f"\nTop variants (up to 20):")
        # Sort all by closure rate descending
        variant_results.sort(key=lambda x: x[5], reverse=True)
        for i, (label, n_ok, n_nz, n_f, fbt, rate) in enumerate(
                variant_results[:20]):
            fail_str = ", ".join(f"{k}:{v}" for k, v in sorted(fbt.items()))
            print(f"  [{rate*100:5.1f}%] {label}")
            if n_f > 0:
                print(f"         {n_ok}/{n_nz} in Λ, failures: {fail_str}")

        # Check: do ANY variants have zero type3×type3 failures?
        zero_t3_count = sum(1 for _, _, _, _, fbt, _ in variant_results
                           if fbt.get("t3×t3", 0) == 0)
        print(f"\nVariants with zero type3×type3 failures: {zero_t3_count}")

        # Analyze whether any variant fixes condition 3
        # For the best variants, do a detailed Wilson-condition breakdown
        if best_variants:
            print(f"\n--- Detailed analysis of best variant ---")
            best_label = best_variants[0][0]
            # Re-parse the label to get parameters
            # Run a detailed check on the best variant
            # Find the best variant's parameters
            for conj_codes, sign_codes, routing, label in enumerate_variants():
                if label == best_label:
                    _detailed_check(test_pairs, pair_labels, conj_codes,
                                    sign_codes, routing, label)
                    break

    # ------------------------------------------------------------------
    # Step 5: Summary
    # ------------------------------------------------------------------
    if verbose:
        print("\n" + "=" * 72)
        print("SUMMARY")
        print("=" * 72)

        if best_closure_rate >= 1.0:
            print("\n✓ Found variant(s) with 100% closure! Further testing needed.")
        else:
            print(f"\n✗ ALL {n_variants} variants fail.")
            print(f"  Best closure rate: {best_closure_rate*100:.1f}%")
            print(f"  This rules out the entire family of triple-octonion")
            print(f"  algebras with conjugation/sign/routing modifications")
            print(f"  on cross-block terms.")

    return variant_results


def _detailed_check(test_pairs, pair_labels, conj_codes, sign_codes,
                    routing, label):
    """Run a detailed Wilson-condition breakdown for one variant."""
    alg = STANDARD_ALGEBRA

    print(f"  Variant: {label}")

    condition_fails = Counter()
    fail_norms = []
    fail_types = Counter()

    for (a, b), plabel in zip(test_pairs, pair_labels):
        prod = multiply_24_variant(a, b, conj_codes, sign_codes, routing)
        norm_sq = float(np.dot(prod, prod))

        if norm_sq < 1e-12:
            continue

        x = alg.element(prod[0:8])
        y = alg.element(prod[8:16])
        z = alg.element(prod[16:24])

        in_L_x = is_in_L(x)
        in_L_y = is_in_L(y)
        in_L_z = is_in_L(z)
        xy_Lsb = is_in_Ls_bar(alg.element(x.coords + y.coords))
        xz_Lsb = is_in_Ls_bar(alg.element(x.coords + z.coords))
        yz_Lsb = is_in_Ls_bar(alg.element(y.coords + z.coords))
        xyz_Ls = is_in_Ls(alg.element(x.coords + y.coords + z.coords))

        ok = (in_L_x and in_L_y and in_L_z and
              xy_Lsb and xz_Lsb and yz_Lsb and xyz_Ls)

        if not ok:
            fail_types[plabel] += 1
            fail_norms.append(norm_sq)
            if not in_L_x: condition_fails['x∈L'] += 1
            if not in_L_y: condition_fails['y∈L'] += 1
            if not in_L_z: condition_fails['z∈L'] += 1
            if not xy_Lsb: condition_fails['x+y∈Ls̄'] += 1
            if not xz_Lsb: condition_fails['x+z∈Ls̄'] += 1
            if not yz_Lsb: condition_fails['y+z∈Ls̄'] += 1
            if not xyz_Ls:  condition_fails['x+y+z∈Ls'] += 1

    total_fails = sum(fail_types.values())
    if total_fails > 0:
        print(f"  Failures by type: {dict(fail_types)}")
        print(f"  Wilson condition violations:")
        for cond, count in condition_fails.most_common():
            print(f"    {cond}: {count}/{total_fails} "
                  f"({count/total_fails*100:.1f}%)")
        if fail_norms:
            norm_dist = Counter(round(n, 0) for n in fail_norms)
            print(f"  Norm² distribution of failures: "
                  f"{dict(sorted(norm_dist.items()))}")
    else:
        print(f"  No failures detected in sample!")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_trial(verbose=True)
