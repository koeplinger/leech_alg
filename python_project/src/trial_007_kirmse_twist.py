"""
trial_007_kirmse_twist.py — Triple Kirmse-twisted octonion algebra on R^24.

Preliminary investigation
=========================
The "Kirmse twist" (index-doubling permutation) corrects Kirmse's non-closed
integral octonions to a genuine maximal order (Coxeter–Dickson integers).
There are 7 maximal orders in the octonions, one for each choice of
distinguished imaginary basis element.

Key question: does Wilson's E8 lattice L, as defined in e8_wilson.py, form
a maximal order under the standard octonion multiplication?  If yes, then
trial 001 already used a maximal-order multiplication and this trial's
proposed construction is structurally identical to trial 001.

This trial proceeds in phases:
  Phase 0: Verify whether L is closed under the standard multiplication.
  Phase 1: Test 7 permuted multiplication tables (one per imaginary basis
           element) to determine which make L a maximal order.
  Phase 2: For each multiplication that closes on L, run the triple-product
           Leech lattice closure test (same as trial 001 but with the
           twisted multiplication).

Algebra definition
==================
Same as trial 001 (triple octonion product with Z₃ cross-block routing),
but the octonion multiplication may differ from the standard Fano-plane rule.
Each "twist" is obtained by applying a specific permutation of {1,...,7} to
the Fano triples, potentially changing the multiplication table.

References
==========
[Wilson2009]   R.A. Wilson, "Octonions and the Leech lattice",
               J. Algebra 322 (2009) 2186–2190.
[Conway-Smith] J.H. Conway, D.A. Smith, "On Quaternions and Octonions" (2003).
[Coxeter1946]  H.S.M. Coxeter, "Integral Cayley numbers",
               Duke Math. J. 13 (1946), 561–578.
[Kirmse1924]   J. Kirmse, "Über die Darstellbarkeit natürlicher ganzer Zahlen
               als Summen von acht Quadraten und über ein mit diesem Problem
               zusammenhängendes nichtkommutatives und nichtassoziatives
               Zahlensystem", Ber. Verh. Sächs. Akad. Wiss. Leipzig, Math.-
               Phys. Kl. 76 (1924), 63–82.
"""

import numpy as np
from collections import Counter
from itertools import permutations
from typing import List, Tuple, Dict, Optional

from octonions import OctonionAlgebra, STANDARD_FANO_TRIPLES, STANDARD_ALGEBRA
from e8_wilson import wilson_e8_roots, is_in_L
from leech_wilson import (
    leech_type1_vectors,
    leech_type2_vectors,
    leech_type3_vectors,
    is_in_leech,
    is_in_Ls_bar,
    is_in_Ls,
)
from trial_001_triple_octonion import (
    multiply_24,
    check_leech_membership,
    flatten_leech_vector,
)


# ---------------------------------------------------------------------------
# Phase 0: Verify L closure under standard multiplication
# ---------------------------------------------------------------------------

def verify_L_closure(alg: OctonionAlgebra, verbose: bool = True) -> Dict:
    """
    Test whether Wilson's E8 lattice L is closed under the given octonion
    multiplication.

    Tests all 240 × 240 = 57,600 pairs of E8 roots.  For each pair (λ, μ),
    computes the product λ · μ and checks whether the result is in L.

    Returns
    -------
    dict with:
      'closed': bool — whether ALL products lie in L
      'total_pairs': int
      'in_L': int — count of products in L
      'not_in_L': int — count of products not in L
      'failures': list of (i, j, product_coords) for first few failures
    """
    roots = wilson_e8_roots()
    n = len(roots)

    total = 0
    in_L_count = 0
    not_in_L_count = 0
    failures = []

    for i, r1 in enumerate(roots):
        for j, r2 in enumerate(roots):
            # Compute product using the given algebra
            prod_coords = alg._mul_coords(r1.coords, r2.coords)
            prod = alg.element(prod_coords)
            total += 1

            if is_in_L(prod):
                in_L_count += 1
            else:
                not_in_L_count += 1
                if len(failures) < 20:
                    failures.append((i, j, prod_coords.copy()))

    closed = (not_in_L_count == 0)

    if verbose:
        print(f"\n  L closure test for algebra '{alg.name}':")
        print(f"    Tested: {total} root pairs")
        print(f"    Products in L: {in_L_count}")
        print(f"    Products NOT in L: {not_in_L_count}")
        print(f"    L is {'CLOSED' if closed else 'NOT CLOSED'} under this multiplication.")

        if failures:
            print(f"\n    First few failures:")
            for idx, (i, j, pc) in enumerate(failures[:5]):
                print(f"      roots[{i}] * roots[{j}] = {pc}")
                # Check which coordinates are non-half-integer
                c2 = pc * 2.0
                is_int = np.allclose(c2, np.round(c2), atol=1e-9)
                print(f"        All coords half-integer: {is_int}")
                if is_int:
                    c2i = np.round(c2).astype(int)
                    parities = c2i % 2
                    print(f"        Coord parities (×2): {parities}")
                    print(f"        Sum (×2): {np.sum(c2i)}")

    return {
        'closed': closed,
        'total_pairs': total,
        'in_L': in_L_count,
        'not_in_L': not_in_L_count,
        'failures': failures,
    }


# ---------------------------------------------------------------------------
# Phase 1: Permuted multiplication tables
# ---------------------------------------------------------------------------

def apply_permutation_to_fano(
    perm: Dict[int, int],
    fano_triples: Tuple[Tuple[int, int, int], ...] = STANDARD_FANO_TRIPLES,
) -> Tuple[Tuple[int, int, int], ...]:
    """
    Apply a permutation of {1,...,7} to the Fano triples.

    Parameters
    ----------
    perm : dict mapping i -> perm[i] for i in {1,...,7}
    fano_triples : the base Fano triples

    Returns
    -------
    New Fano triples with permuted indices.
    """
    return tuple(
        (perm[a], perm[b], perm[c]) for (a, b, c) in fano_triples
    )


def build_transposition_algebras() -> List[Tuple[str, OctonionAlgebra]]:
    """
    Build 21 octonion algebras by applying each transposition (s, t) with
    s, t in {1,...,7} to the standard Fano triples.

    A transposition that is NOT a Fano-plane automorphism will produce a
    multiplication table with different structure constants (different signs
    on some products), potentially changing which lattices are closed.

    Returns list of (description, algebra) pairs.
    """
    algebras = []

    for s in range(1, 8):
        for t in range(s + 1, 8):
            perm = {i: i for i in range(1, 8)}
            perm[s] = t
            perm[t] = s

            new_triples = apply_permutation_to_fano(perm)
            name = f"swap({s},{t})"

            try:
                alg = OctonionAlgebra(new_triples, name=name)
                algebras.append((name, alg))
            except ValueError as e:
                print(f"  WARNING: {name} failed: {e}")

    return algebras


def build_index_doubling_algebras() -> List[Tuple[str, OctonionAlgebra]]:
    """
    Build 7 algebras using the index-doubling permutation σ_s(k) = 2k - s (mod 7)
    for each distinguished element s ∈ {1,...,7}.

    In F_7 arithmetic (with 7 ≡ 0), σ_s fixes s and permutes the other 6
    elements in two 3-cycles.

    These are the permutations directly associated with the Kirmse twist.
    """
    algebras = []

    for s_our in range(1, 8):
        # Convert to F_7: our index i corresponds to F_7 element i % 7
        s_f7 = s_our % 7  # 7 -> 0, 1->1, ..., 6->6

        perm = {}
        for k_our in range(1, 8):
            k_f7 = k_our % 7
            img_f7 = (2 * k_f7 - s_f7) % 7
            img_our = img_f7 if img_f7 != 0 else 7
            perm[k_our] = img_our

        # Verify it's a valid permutation
        assert sorted(perm.values()) == list(range(1, 8)), \
            f"σ_{s_our} is not a valid permutation: {perm}"

        new_triples = apply_permutation_to_fano(perm)
        name = f"σ_{s_our}(k)=2k-{s_our}"

        # Check if this is the identity permutation
        is_identity = all(perm[k] == k for k in range(1, 8))
        if is_identity:
            name += " [identity]"

        try:
            alg = OctonionAlgebra(new_triples, name=name)
            algebras.append((name, alg))
        except ValueError as e:
            print(f"  WARNING: {name} failed: {e}")

    return algebras


# ---------------------------------------------------------------------------
# Phase 2: Leech lattice closure test (reused from trial 001)
# ---------------------------------------------------------------------------

def run_leech_closure_test(
    alg: OctonionAlgebra,
    verbose: bool = True,
    pair_limit: int = 200,
) -> Dict:
    """
    Run the Leech lattice closure test using the given algebra in the
    triple-product construction with Z₃ cross-block routing.

    Same methodology as trial 001 but with a configurable algebra.
    """
    if verbose:
        print(f"\n  Leech closure test for algebra '{alg.name}':")

    type1 = [flatten_leech_vector(t) for t in leech_type1_vectors()]
    type2 = [flatten_leech_vector(t) for t in leech_type2_vectors()]
    type3_raw = leech_type3_vectors()
    rng = np.random.RandomState(42)
    type3_indices = rng.choice(len(type3_raw), min(len(type3_raw), 500), replace=False)
    type3 = [flatten_leech_vector(type3_raw[i]) for i in type3_indices]

    test_sets = [
        ("t1×t1", type1[:50], type1[:50], pair_limit),
        ("t1×t2", type1[:30], type2[:30], 100),
        ("t2×t2", type2[:50], type2[:50], pair_limit),
        ("t1×t3", type1[:20], type3[:20], 100),
        ("t2×t3", type2[:20], type3[:20], 100),
        ("t3×t3", type3[:30], type3[:30], pair_limit),
    ]

    results = {}
    total_tested = 0
    total_in_leech = 0
    total_fail = 0
    failure_conditions = Counter()

    for label, set_a, set_b, plimit in test_sets:
        tested = 0
        in_leech = 0
        failures = 0

        for i, a in enumerate(set_a):
            if tested >= plimit:
                break
            for j, b in enumerate(set_b):
                if tested >= plimit:
                    break
                prod = multiply_24(a, b, alg)
                prod_norm_sq = float(np.dot(prod, prod))
                if prod_norm_sq < 1e-12:
                    tested += 1
                    continue

                diag = check_leech_membership(prod)
                tested += 1

                if diag['in_leech']:
                    in_leech += 1
                else:
                    failures += 1
                    for cond in ['x_in_L', 'y_in_L', 'z_in_L',
                                 'xy_in_Lsbar', 'xz_in_Lsbar', 'yz_in_Lsbar',
                                 'xyz_in_Ls']:
                        if not diag[cond]:
                            failure_conditions[cond] += 1

        total_tested += tested
        total_in_leech += in_leech
        total_fail += failures
        results[label] = {'tested': tested, 'in_leech': in_leech, 'failures': failures}

        if verbose:
            rate = failures / max(tested, 1) * 100
            print(f"    {label}: {tested} tested, {in_leech} in Λ, "
                  f"{failures} failures ({rate:.1f}%)")

    closure_rate = total_in_leech / max(total_tested, 1) * 100
    if verbose:
        print(f"    TOTAL: {total_tested} tested, {total_in_leech} in Λ "
              f"({closure_rate:.1f}%), {total_fail} failures")
        if failure_conditions:
            print(f"    Failure conditions: {dict(failure_conditions.most_common())}")

    return {
        'total_tested': total_tested,
        'total_in_leech': total_in_leech,
        'total_fail': total_fail,
        'closure_rate_pct': closure_rate,
        'results_by_type': results,
        'failure_conditions': dict(failure_conditions),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_trial(verbose: bool = True):
    """Execute trial 007."""

    if verbose:
        print("=" * 72)
        print("TRIAL 007: Triple Kirmse-twisted octonion algebra")
        print("=" * 72)

    # ==================================================================
    # Phase 0: Is Wilson's L closed under standard octonion multiplication?
    # ==================================================================
    if verbose:
        print("\n" + "-" * 72)
        print("PHASE 0: Verify L closure under standard multiplication")
        print("-" * 72)

    std_result = verify_L_closure(STANDARD_ALGEBRA, verbose=verbose)

    if std_result['closed']:
        if verbose:
            print("\n  *** FINDING: Wilson's L IS a maximal order under the")
            print("  standard octonion multiplication. This means trial 001")
            print("  already used the Coxeter–Dickson integral octonions.")
            print("  The 'Kirmse twist' is already implicitly applied.")
    else:
        if verbose:
            print("\n  *** FINDING: Wilson's L is NOT closed under the standard")
            print("  multiplication. A Kirmse twist is needed.")

    # ==================================================================
    # Phase 1: Test index-doubling permutations
    # ==================================================================
    if verbose:
        print("\n" + "-" * 72)
        print("PHASE 1: Index-doubling permutations σ_s(k) = 2k − s (mod 7)")
        print("-" * 72)
        print("\n  Testing 7 permutations, one per imaginary basis element...")

    idx_dbl_algebras = build_index_doubling_algebras()

    idx_dbl_results = {}
    for name, alg in idx_dbl_algebras:
        result = verify_L_closure(alg, verbose=verbose)
        idx_dbl_results[name] = result

    # Count how many close
    n_closed = sum(1 for r in idx_dbl_results.values() if r['closed'])
    if verbose:
        print(f"\n  Summary: {n_closed}/7 index-doubling permutations "
              f"make L closed.")

    # ==================================================================
    # Phase 1b: Test transposition permutations
    # ==================================================================
    if verbose:
        print("\n" + "-" * 72)
        print("PHASE 1b: Transposition permutations (s ↔ t)")
        print("-" * 72)
        print("\n  Testing 21 transpositions of imaginary basis elements...")

    trans_algebras = build_transposition_algebras()

    trans_results = {}
    for name, alg in trans_algebras:
        result = verify_L_closure(alg, verbose=False)
        trans_results[name] = result
        if verbose:
            status = "CLOSED" if result['closed'] else f"NOT CLOSED ({result['not_in_L']} failures)"
            print(f"    {name}: {status}")

    n_closed_trans = sum(1 for r in trans_results.values() if r['closed'])
    if verbose:
        print(f"\n  Summary: {n_closed_trans}/21 transpositions make L closed.")

    # ==================================================================
    # Phase 2: Leech lattice closure for all algebras that close on L
    # ==================================================================
    if verbose:
        print("\n" + "-" * 72)
        print("PHASE 2: Leech lattice closure tests")
        print("-" * 72)

    # Collect all algebras where L is closed
    closed_algebras = []

    if std_result['closed']:
        closed_algebras.append(("standard", STANDARD_ALGEBRA))

    for name, alg in idx_dbl_algebras:
        if idx_dbl_results[name]['closed']:
            # Skip if this is the identity permutation (= standard)
            if "[identity]" not in name:
                closed_algebras.append((name, alg))

    for name, alg in trans_algebras:
        if trans_results[name]['closed']:
            closed_algebras.append((name, alg))

    if not closed_algebras:
        if verbose:
            print("\n  No algebra closes on L. Cannot proceed to Leech test.")
            print("  This would mean Wilson's L is NOT a maximal order under")
            print("  ANY of the tested multiplications — a surprising finding.")
        return {
            'phase0': std_result,
            'phase1_idx_dbl': idx_dbl_results,
            'phase1_trans': trans_results,
            'phase2': {},
        }

    if verbose:
        print(f"\n  Testing {len(closed_algebras)} algebra(s) for Leech closure:")
        for name, _ in closed_algebras:
            print(f"    - {name}")

    leech_results = {}
    for name, alg in closed_algebras:
        result = run_leech_closure_test(alg, verbose=verbose)
        leech_results[name] = result

    # ==================================================================
    # Summary
    # ==================================================================
    if verbose:
        print("\n" + "=" * 72)
        print("SUMMARY")
        print("=" * 72)

        print(f"\n  Phase 0: L {'IS' if std_result['closed'] else 'is NOT'} "
              f"closed under standard multiplication.")
        print(f"  Phase 1: {n_closed}/7 index-doubling permutations close on L; "
              f"{n_closed_trans}/21 transpositions close on L.")

        if leech_results:
            print(f"\n  Phase 2: Leech closure rates:")
            for name, r in leech_results.items():
                print(f"    {name}: {r['closure_rate_pct']:.1f}% "
                      f"({r['total_in_leech']}/{r['total_tested']})")

            # Find best
            best_name = max(leech_results, key=lambda n: leech_results[n]['closure_rate_pct'])
            best_rate = leech_results[best_name]['closure_rate_pct']
            print(f"\n  Best: {best_name} at {best_rate:.1f}%")

            if best_rate < 100:
                print(f"\n  VERDICT: FAIL — no tested multiplication achieves "
                      f"full Leech closure.")
            else:
                print(f"\n  VERDICT: SUCCESS — {best_name} achieves full Leech closure!")
        else:
            print(f"\n  Phase 2: No algebras to test.")
            print(f"\n  VERDICT: Cannot proceed — L is not closed under any "
                  f"tested multiplication.")

    return {
        'phase0': std_result,
        'phase1_idx_dbl': idx_dbl_results,
        'phase1_trans': trans_results,
        'phase2': leech_results,
    }


if __name__ == "__main__":
    run_trial(verbose=True)
