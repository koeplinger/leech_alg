"""
trial_005_triple_okubo.py — Triple symmetric-composition algebra on R²⁴
with Z₃-shifted Petersson isotopes.

Algebra definition
==================

Background: the Petersson isotope construction
----------------------------------------------
Given an octonion algebra (O, ·) and an order-3 automorphism τ, the
Petersson isotope with parameter τ^k (k = 0, 1, 2) is:

    x *_k y = τ^k(x̄) · τ^{2k}(ȳ)

where x̄ denotes octonion conjugation.  This gives three different
8-dimensional symmetric composition algebras:

    k=0:  x *₀ y = x̄ · ȳ           (the para-octonion algebra)
    k=1:  x *₁ y = τ(x̄) · τ²(ȳ)   (an Okubo algebra)
    k=2:  x *₂ y = τ²(x̄) · τ(ȳ)   (a different Okubo algebra, isomorphic
                                       to the first as abstract algebra but
                                       with different structure constants)

All three are symmetric composition algebras:
    n(x * y) = n(x) n(y)       (composition norm)
    (x * y) * x = x * (y * x) = n(x) y

The three algebras form a natural Z₃ orbit under the powers of τ.

Why not three "shifted mediator" Okubo algebras?
-------------------------------------------------
The user's original idea was to pick one Fano line (associative triple)
and build three Okubo algebras by cyclically rotating which element of
the triple serves as mediator.  Computational verification shows that for
EVERY Fano line in the standard octonion algebra, only 1 or 2 of the 3
cyclic mediators yield valid automorphisms (see prompt_logs/023).  No
Fano line admits all three.

The Z₃-orbit of Petersson isotopes {τ⁰, τ¹, τ²} is the closest valid
construction to the user's intent: same Fano line, same τ, three
"shifts" that are genuinely different algebras.

The 24-dimensional algebra
--------------------------
Three blocks: B₀ (indices 0–7), B₁ (indices 8–15), B₂ (indices 16–23).

    Block k uses the Petersson isotope *_k:
        B₀: para-octonion  (τ⁰)
        B₁: Okubo from τ   (τ¹)
        B₂: Okubo from τ²  (τ²)

The multiplication rule mirrors the triple-octonion structure:
    Same-block:  Bα × Bα → Bα using *_α
    Cross-block: Bα × Bβ → Bγ using *_γ  (product of the target block)

The choice to use the TARGET block's product for cross-block terms is the
most natural one: the result lives in block γ, so it should be expressed
in block γ's algebra.

The Fano line and automorphism τ
---------------------------------
Fano line: {e₁, e₃, e₇}  (the standard choice from [MarraniCorradettiZucconi2025]).
Valid mediator: e₃.
Rotation pairs: (e₂, e₅) and (e₄, e₆), angle 2π/3.

This trial includes the base test AND all discrete variants (conjugation,
sign, routing) as in trials 001 + 003, for exhaustiveness.

References
==========
[MarraniCorradettiZucconi2025]  eq. (1.5)–(1.6) for the Petersson isotope.
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              J. Algebra 322 (2009) 2186–2190.
"""

import numpy as np
from collections import Counter

from octonions import STANDARD_ALGEBRA, OctonionAlgebra
from okubo import OkuboAlgebra, from_octonions, _validate_order3, _validate_automorphism
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


# ---------------------------------------------------------------------------
# Build the three 8-dimensional algebras
# ---------------------------------------------------------------------------

def build_tau_matrix():
    """
    The standard order-3 automorphism τ of the octonion algebra.

    Fixes: {e₀, e₁, e₃, e₇}  (quaternion subalgebra from Fano line {1,3,7}).
    Rotates: (e₂, e₅) and (e₄, e₆) by 2π/3.
    Mediator: e₃.

    Reference: [MarraniCorradettiZucconi2025] eq. (1.5).
    """
    T = np.eye(8)
    c, s = -0.5, np.sqrt(3.0) / 2.0   # cos 2π/3, sin 2π/3
    # (e₂, e₅) rotation
    T[2, 2] = c;  T[2, 5] = -s
    T[5, 2] = s;  T[5, 5] = c
    # (e₄, e₆) rotation
    T[4, 4] = c;  T[4, 6] = -s
    T[6, 4] = s;  T[6, 6] = c
    return T


def build_para_octonion_struct():
    """
    Structure constants for the para-octonion algebra: x * y = x̄ · ȳ.
    This is the Petersson isotope with τ⁰ = identity.
    """
    alg = STANDARD_ALGEBRA
    conj = np.diag([1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0])
    struct = np.zeros((8, 8, 8))
    for i in range(8):
        ei = np.zeros(8); ei[i] = 1.0
        ei_bar = conj @ ei
        for j in range(8):
            ej = np.zeros(8); ej[j] = 1.0
            ej_bar = conj @ ej
            struct[i, j, :] = alg._mul_coords(ei_bar, ej_bar)
    return struct


def build_three_algebras():
    """
    Build the three 8-dimensional algebras for the Z₃ orbit.

    Returns (alg0, alg1, alg2):
        alg0: para-octonion  (τ⁰)
        alg1: Okubo from τ   (τ¹)
        alg2: Okubo from τ²  (τ²)
    """
    tau = build_tau_matrix()
    _validate_order3(tau)
    _validate_automorphism(STANDARD_ALGEBRA, tau)

    tau2 = tau @ tau
    _validate_order3(tau2)
    _validate_automorphism(STANDARD_ALGEBRA, tau2)

    alg0 = OkuboAlgebra(build_para_octonion_struct(), name="para-octonion (τ⁰)")
    alg1 = from_octonions(STANDARD_ALGEBRA, tau)
    alg1.name = "Okubo from τ¹"
    alg2 = from_octonions(STANDARD_ALGEBRA, tau2)
    alg2.name = "Okubo from τ²"

    return alg0, alg1, alg2


# ---------------------------------------------------------------------------
# Verify the three algebras
# ---------------------------------------------------------------------------

def verify_algebra(alg, name):
    """Check composition norm and symmetric composition on all basis pairs."""
    comp_ok = True
    sym_ok = True
    for i in range(8):
        ei = alg.basis_element(i)
        nx = ei.norm_sq()
        for j in range(8):
            ej = alg.basis_element(j)
            # Composition: n(ei*ej) = n(ei)*n(ej)
            if abs((ei * ej).norm_sq() - nx * ej.norm_sq()) > 1e-10:
                comp_ok = False
            # Symmetric: (ei*ej)*ei = n(ei)*ej
            if not ((ei * ej) * ei).is_close(nx * ej, atol=1e-10):
                sym_ok = False
            if not (ei * (ej * ei)).is_close(nx * ej, atol=1e-10):
                sym_ok = False
    return comp_ok, sym_ok


# ---------------------------------------------------------------------------
# 24-dimensional multiplication
# ---------------------------------------------------------------------------

def target_block(bi, bj):
    """Return target block for product of blocks bi and bj."""
    if bi == bj:
        return bi
    return 3 - bi - bj


def multiply_24(a, b, algebras):
    """
    Multiply two R²⁴ vectors using the triple symmetric-composition rule.

    algebras: list of 3 OkuboAlgebra objects, one per block.
    Cross-block products use the TARGET block's algebra.
    """
    result = np.zeros(24, dtype=np.float64)
    a_blocks = [a[s] for s in BLOCK_SLICES]
    b_blocks = [b[s] for s in BLOCK_SLICES]

    for bi in range(3):
        for bj in range(3):
            bt = target_block(bi, bj)
            product = algebras[bt]._mul_coords(a_blocks[bi], b_blocks[bj])
            result[BLOCK_SLICES[bt]] += product

    return result


# ---------------------------------------------------------------------------
# Conjugation variant multiplication
# ---------------------------------------------------------------------------

CONJ_LABELS = ["a*b", "ā*b", "a*b̄", "ā*b̄"]
ROUTING_LABELS = ["third", "left", "right"]


def conjugate_octonion(v):
    """Octonion conjugation: (a₀, a₁, ..., a₇) → (a₀, -a₁, ..., -a₇)."""
    c = v.copy()
    c[1:] = -c[1:]
    return c


def multiply_24_variant(a, b, algebras, conj_codes, sign_codes, routing):
    """
    Multiply with discrete variant parameters.

    conj_codes: list of 3 tuples (conj_left, conj_right) ∈ {0,1}²
    sign_codes: list of 3 ints (+1 or -1)
    routing: 0=third, 1=left, 2=right
    """
    result = np.zeros(24, dtype=np.float64)
    a_blocks = [a[s].copy() for s in BLOCK_SLICES]
    b_blocks = [b[s].copy() for s in BLOCK_SLICES]

    # Same-block products (no modifications)
    for bi in range(3):
        product = algebras[bi]._mul_coords(a_blocks[bi], b_blocks[bi])
        result[BLOCK_SLICES[bi]] += product

    # Cross-block products with modifications
    cross_pairs = [
        (0, 0, 1), (0, 1, 0),  # pair 0: blocks 0,1
        (1, 0, 2), (1, 2, 0),  # pair 1: blocks 0,2
        (2, 1, 2), (2, 2, 1),  # pair 2: blocks 1,2
    ]

    for pair_idx, left_blk, right_blk in cross_pairs:
        conj_left, conj_right = conj_codes[pair_idx]
        left_vec = conjugate_octonion(a_blocks[left_blk]) if conj_left else a_blocks[left_blk]
        right_vec = conjugate_octonion(b_blocks[right_blk]) if conj_right else b_blocks[right_blk]

        # Target block
        if routing == 0:
            tgt = 3 - left_blk - right_blk
        elif routing == 1:
            tgt = left_blk
        else:
            tgt = right_blk

        product = algebras[tgt]._mul_coords(left_vec, right_vec)
        product = sign_codes[pair_idx] * product
        result[BLOCK_SLICES[tgt]] += product

    return result


# ---------------------------------------------------------------------------
# Leech membership
# ---------------------------------------------------------------------------

def check_in_leech(v):
    alg = STANDARD_ALGEBRA
    x = alg.element(v[0:8])
    y = alg.element(v[8:16])
    z = alg.element(v[16:24])
    return is_in_leech(x, y, z)


def check_wilson_conditions(v):
    """Return dict of which Wilson conditions pass."""
    alg = STANDARD_ALGEBRA
    x = alg.element(v[0:8])
    y = alg.element(v[8:16])
    z = alg.element(v[16:24])
    return {
        'x∈L': is_in_L(x),
        'y∈L': is_in_L(y),
        'z∈L': is_in_L(z),
        'x+y∈Ls̄': is_in_Ls_bar(alg.element(x.coords + y.coords)),
        'x+z∈Ls̄': is_in_Ls_bar(alg.element(x.coords + z.coords)),
        'y+z∈Ls̄': is_in_Ls_bar(alg.element(y.coords + z.coords)),
        'x+y+z∈Ls': is_in_Ls(alg.element(x.coords + y.coords + z.coords)),
    }


def flatten_triple(trip):
    return np.concatenate(trip)


# ---------------------------------------------------------------------------
# Main trial
# ---------------------------------------------------------------------------

def run_trial(verbose=True):
    if verbose:
        print("=" * 72)
        print("TRIAL 005: Triple symmetric-composition algebra")
        print("             (para-octonion + Okubo_τ + Okubo_τ²)")
        print("=" * 72)

    # ------------------------------------------------------------------
    # Step 1: Build and verify the three algebras
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 1: Building three 8-dim symmetric composition algebras...")

    alg0, alg1, alg2 = build_three_algebras()
    algebras = [alg0, alg1, alg2]

    for i, alg in enumerate(algebras):
        comp_ok, sym_ok = verify_algebra(alg, alg.name)
        if verbose:
            print(f"  Block {i} ({alg.name}):")
            print(f"    Composition norm: {'PASS' if comp_ok else 'FAIL'}")
            print(f"    Symmetric composition: {'PASS' if sym_ok else 'FAIL'}")

    # Verify they are genuinely different
    diff_01 = not np.allclose(alg0._struct, alg1._struct)
    diff_02 = not np.allclose(alg0._struct, alg2._struct)
    diff_12 = not np.allclose(alg1._struct, alg2._struct)
    if verbose:
        print(f"\n  All three algebras different: "
              f"{diff_01 and diff_02 and diff_12}")

    # ------------------------------------------------------------------
    # Step 2: Generate test vectors
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 2: Generating Min(Λ) vectors...")

    type1 = [flatten_triple(t) for t in leech_type1_vectors()]
    type2 = [flatten_triple(t) for t in leech_type2_vectors()]
    type3_raw = leech_type3_vectors()
    rng = np.random.RandomState(42)
    type3_idx = rng.choice(len(type3_raw), 2000, replace=False)
    type3 = [flatten_triple(type3_raw[i]) for i in type3_idx]

    if verbose:
        print(f"  Type 1: {len(type1)}")
        print(f"  Type 2: {len(type2)}")
        print(f"  Type 3: {len(type3)} (sampled from {len(type3_raw)})")

    # ------------------------------------------------------------------
    # Step 3: Base algebra test (like trial 001)
    # ------------------------------------------------------------------
    if verbose:
        print("\nStep 3: Base algebra closure test...")

    test_configs = [
        ("t1×t1", type1, type1, 200),
        ("t1×t2", type1[:50], type2[:50], 200),
        ("t2×t2", type2[:100], type2[:100], 500),
        ("t1×t3", type1[:30], type3[:30], 100),
        ("t2×t3", type2[:30], type3[:30], 100),
        ("t3×t3", type3[:100], type3[:100], 500),
    ]

    total_tested = 0
    total_ok = 0
    total_fail = 0
    total_zero = 0
    fail_by_type = Counter()
    condition_fails = Counter()
    fail_norms = []

    for label, set_a, set_b, limit in test_configs:
        tested = 0
        ok_count = 0
        zero_count = 0
        fail_count = 0

        for i, a in enumerate(set_a):
            if tested >= limit:
                break
            for j, b in enumerate(set_b):
                if tested >= limit:
                    break
                prod = multiply_24(a, b, algebras)
                ns = float(np.dot(prod, prod))
                tested += 1

                if ns < 1e-12:
                    zero_count += 1
                    continue

                if check_in_leech(prod):
                    ok_count += 1
                else:
                    fail_count += 1
                    fail_by_type[label] += 1
                    fail_norms.append(ns)

                    conds = check_wilson_conditions(prod)
                    for cname, cval in conds.items():
                        if not cval:
                            condition_fails[cname] += 1

        total_tested += tested
        total_ok += ok_count
        total_fail += fail_count
        total_zero += zero_count

        if verbose:
            rate = ok_count / max(1, tested - zero_count)
            print(f"  {label}: {tested} tested, {ok_count} in Λ, "
                  f"{zero_count} zero, {fail_count} fail "
                  f"({rate*100:.1f}% closure)")

    if verbose:
        print(f"\n  TOTAL: {total_tested} tested, {total_ok} in Λ, "
              f"{total_zero} zero, {total_fail} fail")
        nonzero = total_tested - total_zero
        if nonzero > 0:
            print(f"  Closure rate: {total_ok/nonzero*100:.1f}%")

        if total_fail > 0:
            print(f"\n  Failures by type: {dict(fail_by_type)}")
            print(f"  Wilson condition violations:")
            for cname, count in condition_fails.most_common():
                print(f"    {cname}: {count}/{total_fail} "
                      f"({count/total_fail*100:.1f}%)")
            if fail_norms:
                norm_dist = Counter(round(n, 0) for n in fail_norms)
                top_norms = sorted(norm_dist.items())[:15]
                print(f"  Norm² distribution (top 15): "
                      f"{dict(top_norms)}")

    base_closure = total_ok / max(1, total_tested - total_zero)

    # ------------------------------------------------------------------
    # Step 4: Discrete variant sweep (like trial 003)
    # ------------------------------------------------------------------
    if verbose:
        print("\n" + "=" * 72)
        print("Step 4: Discrete variant sweep (conjugation × sign × routing)")
        print("=" * 72)

    # Build a smaller test set for the sweep
    sweep_pairs = []
    sweep_labels = []

    # 100 type3×type3 (critical), 50 type2×type2, 25 each of the rest
    for i in range(0, 200, 2):
        sweep_pairs.append((type3[i], type3[i+1]))
        sweep_labels.append("t3×t3")
    idx2 = rng.choice(len(type2), 100, replace=False)
    for i in range(0, 100, 2):
        sweep_pairs.append((type2[idx2[i]], type2[idx2[i+1]]))
        sweep_labels.append("t2×t2")
    for i in range(25):
        sweep_pairs.append((type1[i], type2[i]))
        sweep_labels.append("t1×t2")
    for i in range(25):
        sweep_pairs.append((type1[i], type3[i]))
        sweep_labels.append("t1×t3")
    for i in range(25):
        sweep_pairs.append((type2[i], type3[i]))
        sweep_labels.append("t2×t3")

    n_sweep_pairs = len(sweep_pairs)
    if verbose:
        tc = Counter(sweep_labels)
        print(f"  Sweep pairs: {n_sweep_pairs} ({dict(tc)})")

    conj_options = [(0, 0), (1, 0), (0, 1), (1, 1)]
    best_rate = 0.0
    best_label = ""
    n_variants = 0
    rate_hist = Counter()

    for routing in range(3):
        for c0 in conj_options:
            for c1 in conj_options:
                for c2 in conj_options:
                    for s0 in (+1, -1):
                        for s1 in (+1, -1):
                            for s2 in (+1, -1):
                                conj_codes = [c0, c1, c2]
                                sign_codes = [s0, s1, s2]

                                n_ok = 0
                                n_nz = 0

                                for (a, b), pl in zip(sweep_pairs,
                                                      sweep_labels):
                                    prod = multiply_24_variant(
                                        a, b, algebras,
                                        conj_codes, sign_codes, routing)
                                    ns = float(np.dot(prod, prod))
                                    if ns < 1e-12:
                                        continue
                                    n_nz += 1
                                    if check_in_leech(prod):
                                        n_ok += 1

                                rate = n_ok / max(1, n_nz)
                                rate_hist[round(rate, 3)] += 1

                                if rate > best_rate:
                                    best_rate = rate
                                    conj_strs = [CONJ_LABELS[cl*2+cr]
                                                 for cl, cr in conj_codes]
                                    sign_str = "".join("+" if s > 0 else "-"
                                                       for s in sign_codes)
                                    best_label = (
                                        f"route={ROUTING_LABELS[routing]}, "
                                        f"conj=[{','.join(conj_strs)}], "
                                        f"sign={sign_str}")

                                n_variants += 1

        if verbose:
            print(f"  ... routing={ROUTING_LABELS[routing]} done "
                  f"({n_variants} total, best={best_rate*100:.1f}%)")

    if verbose:
        print(f"\n  Total variants: {n_variants}")
        print(f"  Best closure rate: {best_rate*100:.1f}%")
        print(f"  Best variant: {best_label}")
        print(f"\n  Closure rate distribution (top 10):")
        for r, c in sorted(rate_hist.items(), reverse=True)[:10]:
            print(f"    {r*100:6.1f}%: {c} variants")

    # ------------------------------------------------------------------
    # Step 5: Detailed analysis of best variant
    # ------------------------------------------------------------------
    if verbose and best_rate < 1.0 and total_fail > 0:
        print(f"\n  --- Detailed analysis of best variant ---")
        print(f"  {best_label}")
        # Re-run the best variant with Wilson condition breakdown
        # Parse best_label to extract parameters (or just use the base algebra
        # since it's typically the best)
        print(f"  (See base algebra analysis above for condition breakdown)")

    # ------------------------------------------------------------------
    # Step 6: Summary
    # ------------------------------------------------------------------
    if verbose:
        print("\n" + "=" * 72)
        print("SUMMARY")
        print("=" * 72)

        print(f"\n  Base algebra closure: {base_closure*100:.1f}%")
        print(f"  Best discrete variant: {best_rate*100:.1f}%")
        print(f"  Variants tested: {n_variants}")

        if best_rate >= 1.0:
            print("\n  ✓ Found variant(s) with 100% closure!")
        elif total_fail > 0:
            print(f"\n  ✗ ALL variants fail.")
            print(f"    Failures by type: {dict(fail_by_type)}")
            if condition_fails:
                top_cond = condition_fails.most_common(1)[0]
                print(f"    Dominant condition: {top_cond[0]} "
                      f"({top_cond[1]}/{total_fail})")
        else:
            print("\n  ? All products are zero.")

    return {
        'base_closure': base_closure,
        'best_variant_rate': best_rate,
        'best_variant_label': best_label,
        'n_variants': n_variants,
        'fail_by_type': dict(fail_by_type),
        'condition_fails': dict(condition_fails),
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_trial(verbose=True)
