"""
trial_001_triple_octonion.py — Triple-octonion Z-module on R^24.

Algebra definition
==================
Take the standard octonion algebra O (Fano-plane rule: e_a · e_{a+1} = e_{a+3},
indices mod 7 in {1,…,7}) and make three copies:

    O₁ occupies coordinate indices  0– 7
    O₂ occupies coordinate indices  8–15
    O₃ occupies coordinate indices 16–23

A vector v ∈ R²⁴ decomposes as v = (v₁, v₂, v₃) with vᵢ ∈ R⁸.

The multiplication rule distributes across these blocks. Each basis vector
f_i (i = 0,…,23) lives in exactly one block. When multiplying f_i · f_j, the
two blocks determine where the octonion product is written:

    Same block:      O_α × O_α → O_α    (product stays in the same block)
    Different blocks: O_α × O_β → O_γ    where {α, β, γ} = {1, 2, 3}

Concretely, if a ∈ {0,…,7} is the octonion index within its block:

    f_{a} · f_{b}       = (e_a · e_b)  written into O₁   [both in O₁]
    f_{8+a} · f_{8+b}   = (e_a · e_b)  written into O₂   [both in O₂]
    f_{16+a} · f_{16+b} = (e_a · e_b)  written into O₃   [both in O₃]

    f_{a} · f_{8+b}     = (e_a · e_b)  written into O₃   [O₁ × O₂ → O₃]
    f_{8+a} · f_{b}     = (e_a · e_b)  written into O₃   [O₂ × O₁ → O₃]

    f_{a} · f_{16+b}    = (e_a · e_b)  written into O₂   [O₁ × O₃ → O₂]
    f_{16+a} · f_{b}    = (e_a · e_b)  written into O₂   [O₃ × O₁ → O₂]

    f_{8+a} · f_{16+b}  = (e_a · e_b)  written into O₁   [O₂ × O₃ → O₁]
    f_{16+a} · f_{8+b}  = (e_a · e_b)  written into O₁   [O₃ × O₂ → O₁]

This is a Z₃-symmetric mixing: cross-block products always go to the third
block. The algebra is manifestly symmetric under cyclic permutations of the
three blocks.

Note: the octonion product e_a · e_b in the cross-block terms uses the SAME
structure constants as the within-block terms. No sign changes or conjugations
are introduced for cross terms.

Motivation
==========
This is the simplest conceivable way to build a 24-dimensional algebra from
three copies of the octonions. It serves as a baseline: we expect it to fail,
and the failure pattern will inform more sophisticated constructions.

References
==========
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186–2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
[Wikipedia]   https://en.wikipedia.org/wiki/Octonion
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
    is_in_Ls_bar,
    is_in_Ls,
)
from e8_wilson import is_in_L

# ---------------------------------------------------------------------------
# Block layout
# ---------------------------------------------------------------------------

BLOCK_SLICES = [slice(0, 8), slice(8, 16), slice(16, 24)]
BLOCK_OFFSETS = [0, 8, 16]

# Target block for (block_i, block_j) product.
# {α, β, γ} = {0, 1, 2}: same block → same block; different → third.
def _target_block(bi: int, bj: int) -> int:
    """Return the block index where the product of blocks bi and bj lands."""
    if bi == bj:
        return bi
    return 3 - bi - bj  # the remaining element of {0, 1, 2}


# ---------------------------------------------------------------------------
# 24-dimensional multiplication
# ---------------------------------------------------------------------------

def multiply_24(a: np.ndarray, b: np.ndarray, alg: OctonionAlgebra = STANDARD_ALGEBRA) -> np.ndarray:
    """
    Multiply two R^24 vectors using the triple-octonion rule.

    Parameters
    ----------
    a, b : np.ndarray of shape (24,)
        The two factors.
    alg : OctonionAlgebra
        The octonion algebra to use for each 8-dimensional block product.

    Returns
    -------
    np.ndarray of shape (24,)
        The product vector.
    """
    result = np.zeros(24, dtype=np.float64)

    # Extract the three 8-component blocks from each factor
    a_blocks = [a[s] for s in BLOCK_SLICES]
    b_blocks = [b[s] for s in BLOCK_SLICES]

    # Sum over all 9 block pairs
    for bi in range(3):
        for bj in range(3):
            bt = _target_block(bi, bj)
            product = alg._mul_coords(a_blocks[bi], b_blocks[bj])
            result[BLOCK_SLICES[bt]] += product

    return result


# ---------------------------------------------------------------------------
# Leech lattice membership for 24-vectors
# ---------------------------------------------------------------------------

def check_leech_membership(v: np.ndarray) -> Dict:
    """
    Check whether a 24-vector lies in the Leech lattice Λ.

    Returns a dict with detailed diagnostics:
      - 'in_leech': bool
      - 'norm_sq': squared norm
      - 'x_in_L', 'y_in_L', 'z_in_L': whether each octonion component ∈ L
      - 'xy_in_Lsbar', 'xz_in_Lsbar', 'yz_in_Lsbar': Wilson condition 2
      - 'xyz_in_Ls': Wilson condition 3
    """
    alg = STANDARD_ALGEBRA
    x = alg.element(v[0:8])
    y = alg.element(v[8:16])
    z = alg.element(v[16:24])

    norm_sq = float(np.dot(v, v))

    x_in_L = is_in_L(x)
    y_in_L = is_in_L(y)
    z_in_L = is_in_L(z)

    xy_in_Lsbar = is_in_Ls_bar(alg.element(x.coords + y.coords))
    xz_in_Lsbar = is_in_Ls_bar(alg.element(x.coords + z.coords))
    yz_in_Lsbar = is_in_Ls_bar(alg.element(y.coords + z.coords))

    xyz_in_Ls = is_in_Ls(alg.element(x.coords + y.coords + z.coords))

    in_leech = (x_in_L and y_in_L and z_in_L
                and xy_in_Lsbar and xz_in_Lsbar and yz_in_Lsbar
                and xyz_in_Ls)

    return {
        'in_leech': in_leech,
        'norm_sq': norm_sq,
        'x_in_L': x_in_L,
        'y_in_L': y_in_L,
        'z_in_L': z_in_L,
        'xy_in_Lsbar': xy_in_Lsbar,
        'xz_in_Lsbar': xz_in_Lsbar,
        'yz_in_Lsbar': yz_in_Lsbar,
        'xyz_in_Ls': xyz_in_Ls,
    }


def flatten_leech_vector(trip: Tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    """Convert a (x, y, z) triple of 8-vectors into a single 24-vector."""
    return np.concatenate(trip)


# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------

def sanity_check_norms(vectors: List[np.ndarray], alg: OctonionAlgebra = STANDARD_ALGEBRA) -> Dict:
    """
    Compute product norms for a sample of vector pairs.

    Picks a small sample and computes ||a · b||² for each pair, to check
    whether the products are in the right norm range for the minimal shell.
    """
    n = min(len(vectors), 50)
    sample = vectors[:n]
    norms_sq = []

    for i in range(n):
        for j in range(n):
            prod = multiply_24(sample[i], sample[j], alg)
            ns = float(np.dot(prod, prod))
            if ns > 1e-12:  # skip zero products
                norms_sq.append(ns)

    if not norms_sq:
        return {'min': 0, 'max': 0, 'mean': 0, 'unique_norms': set(), 'count': 0}

    unique = sorted(set(round(x, 6) for x in norms_sq))
    return {
        'min': min(norms_sq),
        'max': max(norms_sq),
        'mean': np.mean(norms_sq),
        'unique_norms': unique,
        'count': len(norms_sq),
        'nonzero_count': len([x for x in norms_sq if x > 1e-12]),
    }


# ---------------------------------------------------------------------------
# Main trial execution
# ---------------------------------------------------------------------------

def run_trial(max_pairs: int = 0, verbose: bool = True):
    """
    Execute trial 001.

    Parameters
    ----------
    max_pairs : int
        Maximum number of (a, b) pairs to test. 0 means test all pairs within
        each vector type and a sample of cross-type pairs.
    verbose : bool
        Print progress and results.

    Returns
    -------
    dict with trial results and failure analysis.
    """
    alg = STANDARD_ALGEBRA

    if verbose:
        print("=" * 72)
        print("TRIAL 001: Triple-octonion Z-module (O₁ ⊕ O₂ ⊕ O₃)")
        print("=" * 72)

    # ------------------------------------------------------------------
    # Step 1: Generate minimal shell vectors
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 1: Generating minimal-shell vectors of Λ...")

    type1 = [flatten_leech_vector(t) for t in leech_type1_vectors()]
    type2 = [flatten_leech_vector(t) for t in leech_type2_vectors()]
    type3_raw = leech_type3_vectors()

    if verbose:
        print(f"  Type 1: {len(type1):,} vectors")
        print(f"  Type 2: {len(type2):,} vectors")
        print(f"  Type 3: {len(type3_raw):,} vectors (will sample)")

    # For type 3, use a sample for efficiency (184,320 is large)
    type3_sample_size = min(len(type3_raw), 2000)
    rng = np.random.RandomState(42)
    type3_indices = rng.choice(len(type3_raw), type3_sample_size, replace=False)
    type3 = [flatten_leech_vector(type3_raw[i]) for i in type3_indices]

    if verbose:
        print(f"  Type 3 sample: {len(type3):,} vectors")

    # ------------------------------------------------------------------
    # Step 2: Sanity check — product norms
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 2: Sanity check — product norms (type 1 × type 1)...")

    norm_info = sanity_check_norms(type1, alg)
    if verbose:
        print(f"  Nonzero products: {norm_info['nonzero_count']} / {norm_info['count']}")
        print(f"  Norm² range: [{norm_info['min']:.1f}, {norm_info['max']:.1f}]")
        print(f"  Distinct norm² values: {norm_info['unique_norms'][:20]}")
        if norm_info['min'] > 8.0 + 1e-6:
            print("  ⚠ WARNING: All product norms exceed minimal shell norm (8).")
            print("    Products never land on the minimal shell.")

    # Also check type 2 × type 2
    if verbose:
        print("\n  Sanity check — product norms (type 2 × type 2, sample)...")
    norm_info_t2 = sanity_check_norms(type2[:50], alg)
    if verbose:
        print(f"  Nonzero products: {norm_info_t2['nonzero_count']} / {norm_info_t2['count']}")
        print(f"  Norm² range: [{norm_info_t2['min']:.1f}, {norm_info_t2['max']:.1f}]")
        print(f"  Distinct norm² values: {norm_info_t2['unique_norms'][:20]}")

    # ------------------------------------------------------------------
    # Step 3: Closure test
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 3: Closure test — checking if products lie in Λ...")

    # Strategy: test systematically within and across types.
    # For efficiency, limit pair counts.
    test_sets = [
        ("type1×type1", type1, type1, 200),
        ("type1×type2", type1[:50], type2[:50], 200),
        ("type2×type2", type2[:100], type2[:100], 500),
        ("type1×type3", type1[:30], type3[:30], 100),
        ("type2×type3", type2[:30], type3[:30], 100),
        ("type3×type3", type3[:50], type3[:50], 200),
    ]

    results = {}
    total_tested = 0
    total_in_leech = 0
    total_zero = 0
    total_fail = 0

    failure_conditions = Counter()  # which Wilson conditions fail
    failure_norms = []
    failure_examples = []

    for label, set_a, set_b, pair_limit in test_sets:
        if verbose:
            print(f"\n  Testing {label}...")

        tested = 0
        in_leech = 0
        zero_products = 0
        failures = 0

        for i, a in enumerate(set_a):
            if tested >= pair_limit:
                break
            for j, b in enumerate(set_b):
                if tested >= pair_limit:
                    break
                prod = multiply_24(a, b, alg)
                prod_norm_sq = float(np.dot(prod, prod))

                if prod_norm_sq < 1e-12:
                    zero_products += 1
                    tested += 1
                    continue

                diag = check_leech_membership(prod)
                tested += 1

                if diag['in_leech']:
                    in_leech += 1
                else:
                    failures += 1

                    # Record which conditions fail
                    for cond in ['x_in_L', 'y_in_L', 'z_in_L',
                                 'xy_in_Lsbar', 'xz_in_Lsbar', 'yz_in_Lsbar',
                                 'xyz_in_Ls']:
                        if not diag[cond]:
                            failure_conditions[cond] += 1

                    failure_norms.append(prod_norm_sq)

                    if len(failure_examples) < 10:
                        failure_examples.append({
                            'label': label,
                            'a_index': i,
                            'b_index': j,
                            'prod_norm_sq': prod_norm_sq,
                            'diagnostics': diag,
                        })

        total_tested += tested
        total_in_leech += in_leech
        total_zero += zero_products
        total_fail += failures

        results[label] = {
            'tested': tested,
            'in_leech': in_leech,
            'zero': zero_products,
            'failures': failures,
        }

        if verbose:
            print(f"    Tested: {tested}, In Λ: {in_leech}, "
                  f"Zero: {zero_products}, Failures: {failures}")

    # ------------------------------------------------------------------
    # Step 4: Failure analysis
    # ------------------------------------------------------------------
    if verbose:
        print("\n" + "=" * 72)
        print("FAILURE ANALYSIS")
        print("=" * 72)

        print(f"\nOverall: {total_tested} pairs tested")
        print(f"  Products in Λ:  {total_in_leech}")
        print(f"  Zero products:  {total_zero}")
        print(f"  Failures:       {total_fail}")

        if total_fail > 0:
            print(f"\nFailure rate: {total_fail}/{total_tested - total_zero} "
                  f"= {total_fail/(total_tested - total_zero)*100:.1f}% "
                  f"(of nonzero products)")

            print("\nWilson condition violations (count of failures where each is violated):")
            for cond, count in failure_conditions.most_common():
                print(f"  {cond}: {count}/{total_fail} "
                      f"({count/total_fail*100:.1f}%)")

            print(f"\nNorm² distribution of failing products:")
            if failure_norms:
                norm_counter = Counter(round(n, 2) for n in failure_norms)
                for norm_val, cnt in sorted(norm_counter.items()):
                    print(f"  ||prod||² = {norm_val}: {cnt} products")

            print(f"\nFirst few failure examples:")
            for ex in failure_examples[:5]:
                d = ex['diagnostics']
                conds_failed = [c for c in ['x_in_L', 'y_in_L', 'z_in_L',
                                            'xy_in_Lsbar', 'xz_in_Lsbar',
                                            'yz_in_Lsbar', 'xyz_in_Ls']
                                if not d[c]]
                print(f"  {ex['label']} pair ({ex['a_index']},{ex['b_index']}): "
                      f"||prod||²={ex['prod_norm_sq']:.2f}, "
                      f"failed: {conds_failed}")

    # ------------------------------------------------------------------
    # Step 5: Summary
    # ------------------------------------------------------------------
    if verbose:
        print("\n" + "=" * 72)
        print("SUMMARY")
        print("=" * 72)

        if total_fail == 0 and total_in_leech > 0:
            print("\n✓ ALL nonzero products lie in Λ. Closure holds on tested pairs.")
        elif total_fail > 0:
            print(f"\n✗ Closure FAILS. {total_fail} of {total_tested - total_zero} "
                  f"nonzero products lie outside Λ.")
            print("\nKey observations:")
            # Identify the dominant failure mode
            if failure_conditions:
                top_cond, top_count = failure_conditions.most_common(1)[0]
                print(f"  - Dominant failure: {top_cond} "
                      f"(violated in {top_count/total_fail*100:.1f}% of failures)")
            if failure_norms:
                unique_norms = sorted(set(round(n, 2) for n in failure_norms))
                if len(unique_norms) <= 5:
                    print(f"  - Failing products have norm² ∈ {unique_norms}")
                else:
                    print(f"  - Failing product norms range from "
                          f"{min(failure_norms):.2f} to {max(failure_norms):.2f}")
        else:
            print("\n? All products are zero. The algebra is trivial on the tested pairs.")

    return {
        'total_tested': total_tested,
        'total_in_leech': total_in_leech,
        'total_zero': total_zero,
        'total_fail': total_fail,
        'results_by_type': results,
        'failure_conditions': dict(failure_conditions),
        'failure_norms': failure_norms,
        'failure_examples': failure_examples,
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_trial(verbose=True)
