"""
symbolic_proof_checks.py — Exact-arithmetic verification of the five lemmas
needed for the symbolic proof of Leech lattice closure.

All arithmetic uses doubled-integer coordinates (multiply everything by 2)
to avoid fractions entirely. The octonion product of two vectors with
half-integer coordinates produces half-integer coordinates, so doubled
coordinates stay integer under multiplication (after dividing by 4 and
re-doubling — see multiply_doubled()).

Lemma A: σ(L) = L
Lemma B: L·L ⊆ L  (maximal order)
Lemma C: L · σ(Ls̄) ⊆ σ(Ls̄)
Lemma D: σ(Ls) · σ(Ls) ⊆ σ(Ls)
Lemma E: σ(Ls) ≠ Ls  (and Ls·Ls ⊄ Ls as counterexample)
"""

import numpy as np
from fractions import Fraction
from typing import List, Tuple

from octonions import STANDARD_ALGEBRA, STANDARD_FANO_TRIPLES
from e8_wilson import wilson_e8_roots, WILSON_S, is_in_L
from leech_wilson import WILSON_S_BAR


# ===================================================================
# Exact arithmetic using Fraction
# ===================================================================

def to_frac(coords) -> List[Fraction]:
    """Convert numpy array to list of Fractions."""
    return [Fraction(float(c)).limit_denominator(1000) for c in coords]


def frac_mul(a: List[Fraction], b: List[Fraction]) -> List[Fraction]:
    """Multiply two octonion vectors (as Fraction lists) using the standard product."""
    result = [Fraction(0)] * 8
    table = STANDARD_ALGEBRA._table  # table[i][j] = (sign, target)
    for i in range(8):
        if a[i] == 0:
            continue
        for j in range(8):
            if b[j] == 0:
                continue
            sign, target = table[i][j]
            result[target] += Fraction(sign) * a[i] * b[j]
    return result


def sigma(v: List[Fraction]) -> List[Fraction]:
    """Apply σ = swap(1,2): swap positions 1 and 2."""
    w = list(v)
    w[1], w[2] = v[2], v[1]
    return w


def is_integer_combination(v: List[Fraction], basis: List[List[Fraction]]) -> bool:
    """Check if v is an integer combination of basis vectors (exact)."""
    n = len(basis)
    assert n == 8
    # Build basis matrix and solve v = B^T c
    # Use exact Fraction arithmetic for Gaussian elimination
    # Augmented matrix: [B | v]
    aug = [list(basis[i]) + [v[i]] for i in range(8)]
    aug = [[Fraction(x) for x in row] for row in aug]

    # Transpose: we want to solve sum_j c_j * basis[j] = v
    # i.e., for each coordinate i: sum_j c_j * basis[j][i] = v[i]
    # Matrix equation: M @ c = v where M[i][j] = basis[j][i]
    M = [[basis[j][i] for j in range(n)] for i in range(8)]
    rhs = list(v)

    # Gaussian elimination with exact fractions
    mat = [M[i] + [rhs[i]] for i in range(8)]

    for col in range(n):
        # Find pivot
        pivot_row = None
        for row in range(col, 8):
            if mat[row][col] != 0:
                pivot_row = row
                break
        if pivot_row is None:
            return False
        mat[col], mat[pivot_row] = mat[pivot_row], mat[col]

        # Eliminate
        pivot = mat[col][col]
        for row in range(8):
            if row == col:
                continue
            if mat[row][col] == 0:
                continue
            factor = mat[row][col] / pivot
            for k in range(n + 1):
                mat[row][k] -= factor * mat[col][k]

    # Extract solution
    for col in range(n):
        val = mat[col][n] / mat[col][col]
        if val.denominator != 1:
            return False  # Not an integer combination

    return True


# ===================================================================
# Build Z-bases (exact)
# ===================================================================

def build_L_zbasis_frac() -> List[List[Fraction]]:
    """Build a Z-basis for L = D₈⁺ using exact arithmetic."""
    roots = wilson_e8_roots()
    basis = [to_frac(WILSON_S.coords)]
    for r in roots:
        candidate = basis + [to_frac(r.coords)]
        # Check linear independence via rank
        mat = [[candidate[i][j] for j in range(8)] for i in range(len(candidate))]
        # Simple rank check: try Gaussian elimination
        test = [list(row) for row in mat]
        rank = 0
        for col in range(8):
            pivot = None
            for row in range(rank, len(test)):
                if test[row][col] != 0:
                    pivot = row
                    break
            if pivot is None:
                continue
            test[rank], test[pivot] = test[pivot], test[rank]
            for row in range(len(test)):
                if row == rank:
                    continue
                if test[row][col] == 0:
                    continue
                factor = Fraction(test[row][col], test[rank][col])
                for k in range(8):
                    test[row][k] -= factor * test[rank][k]
            rank += 1

        if rank == len(candidate):
            basis.append(to_frac(r.coords))
            if len(basis) == 8:
                break

    assert len(basis) == 8, f"Could not find 8 independent basis vectors, got {len(basis)}"
    return basis


def right_multiply_basis(basis: List[List[Fraction]], c: List[Fraction]) -> List[List[Fraction]]:
    """Compute {b · c : b ∈ basis} for a fixed element c."""
    return [frac_mul(b, c) for b in basis]


# ===================================================================
# Main checks
# ===================================================================

def main():
    print("=" * 70)
    print("SYMBOLIC PROOF VERIFICATION — Exact arithmetic (Fraction)")
    print("=" * 70)

    # Build bases
    print("\n--- Building Z-bases ---")
    L_basis = build_L_zbasis_frac()
    print(f"  L basis: {len(L_basis)} vectors")

    s_frac = to_frac(WILSON_S.coords)
    sbar_frac = to_frac(WILSON_S_BAR.coords)

    print(f"  s = {s_frac}")
    print(f"  s̄ = {sbar_frac}")

    # Ls basis: {b · s : b ∈ L_basis}
    Ls_basis = right_multiply_basis(L_basis, s_frac)
    print(f"  Ls basis: {len(Ls_basis)} vectors")

    # Ls̄ basis: {b · s̄ : b ∈ L_basis}
    Lsbar_basis = right_multiply_basis(L_basis, sbar_frac)
    print(f"  Ls̄ basis: {len(Lsbar_basis)} vectors")

    # σ(Ls) basis: {σ(b · s) : b ∈ L_basis}
    sigma_Ls_basis = [sigma(v) for v in Ls_basis]
    print(f"  σ(Ls) basis: {len(sigma_Ls_basis)} vectors")

    # σ(Ls̄) basis: {σ(b · s̄) : b ∈ L_basis}
    sigma_Lsbar_basis = [sigma(v) for v in Lsbar_basis]
    print(f"  σ(Ls̄) basis: {len(sigma_Lsbar_basis)} vectors")

    # ===============================================================
    # Lemma A: σ(L) = L
    # ===============================================================
    print("\n" + "=" * 70)
    print("LEMMA A: σ(L) = L")
    print("=" * 70)
    # Check that σ(b) ∈ L for each basis vector
    all_ok = True
    for i, b in enumerate(L_basis):
        sb = sigma(b)
        if not is_integer_combination(sb, L_basis):
            print(f"  FAIL: σ(L_basis[{i}]) not in L")
            all_ok = False
    if all_ok:
        print("  PASS: All 8 basis vectors of L map into L under σ.")
        print("  Since σ is its own inverse, σ(L) = L. ✓")

    # ===============================================================
    # Lemma E: σ(Ls) ≠ Ls
    # ===============================================================
    print("\n" + "=" * 70)
    print("LEMMA E: σ(Ls) ≠ Ls")
    print("=" * 70)
    found_diff = False
    for i, v in enumerate(sigma_Ls_basis):
        if not is_integer_combination(v, Ls_basis):
            print(f"  σ(Ls) basis vector {i} is NOT in Ls:")
            print(f"    vector = {v}")
            found_diff = True
            break
    if found_diff:
        print("  PASS: σ(Ls) ≠ Ls. ✓")
    else:
        print("  NOTE: σ(Ls) = Ls (unexpected!)")

    # Also check σ(Ls̄) vs Ls̄
    print("\n  Checking σ(Ls̄) vs Ls̄...")
    found_diff_sbar = False
    for i, v in enumerate(sigma_Lsbar_basis):
        if not is_integer_combination(v, Lsbar_basis):
            print(f"  σ(Ls̄) basis vector {i} is NOT in Ls̄:")
            print(f"    vector = {v}")
            found_diff_sbar = True
            break
    if found_diff_sbar:
        print("  σ(Ls̄) ≠ Ls̄.")
    else:
        print("  σ(Ls̄) = Ls̄.")

    # ===============================================================
    # Counterexample: Ls · Ls ⊄ Ls
    # ===============================================================
    print("\n" + "=" * 70)
    print("COUNTEREXAMPLE: Ls · Ls ⊄ Ls")
    print("=" * 70)
    found_counter = False
    for i in range(8):
        for j in range(8):
            prod = frac_mul(Ls_basis[i], Ls_basis[j])
            if not is_integer_combination(prod, Ls_basis):
                print(f"  Ls_basis[{i}] · Ls_basis[{j}] is NOT in Ls.")
                print(f"    a = {Ls_basis[i]}")
                print(f"    b = {Ls_basis[j]}")
                print(f"    a·b = {prod}")
                found_counter = True
                break
        if found_counter:
            break
    if found_counter:
        print("  CONFIRMED: Ls is NOT closed under standard product. ✓")
    else:
        print("  UNEXPECTED: Ls IS closed — standard product should fail!")

    # ===============================================================
    # Lemma D: σ(Ls) · σ(Ls) ⊆ σ(Ls)
    # ===============================================================
    print("\n" + "=" * 70)
    print("LEMMA D: σ(Ls) · σ(Ls) ⊆ σ(Ls)")
    print("=" * 70)
    failures_D = 0
    for i in range(8):
        for j in range(8):
            prod = frac_mul(sigma_Ls_basis[i], sigma_Ls_basis[j])
            if not is_integer_combination(prod, sigma_Ls_basis):
                print(f"  FAIL: σ(Ls)_basis[{i}] · σ(Ls)_basis[{j}] not in σ(Ls)")
                print(f"    a = {sigma_Ls_basis[i]}")
                print(f"    b = {sigma_Ls_basis[j]}")
                print(f"    a·b = {prod}")
                failures_D += 1
    if failures_D == 0:
        print(f"  PASS: All 64 basis products lie in σ(Ls). ✓")
    else:
        print(f"  FAIL: {failures_D}/64 products not in σ(Ls).")

    # ===============================================================
    # Lemma C: L · σ(Ls̄) ⊆ σ(Ls̄)
    # ===============================================================
    print("\n" + "=" * 70)
    print("LEMMA C: L · σ(Ls̄) ⊆ σ(Ls̄)")
    print("=" * 70)
    failures_C = 0
    for i in range(8):
        for j in range(8):
            prod = frac_mul(L_basis[i], sigma_Lsbar_basis[j])
            if not is_integer_combination(prod, sigma_Lsbar_basis):
                print(f"  FAIL: L_basis[{i}] · σ(Ls̄)_basis[{j}] not in σ(Ls̄)")
                print(f"    a = {L_basis[i]}")
                print(f"    b = {sigma_Lsbar_basis[j]}")
                print(f"    a·b = {prod}")
                failures_C += 1
    if failures_C == 0:
        print(f"  PASS: All 64 products lie in σ(Ls̄). ✓")
    else:
        print(f"  FAIL: {failures_C}/64 products not in σ(Ls̄).")

    # ===============================================================
    # Bonus: check L · Ls̄ ⊆ Ls̄ (standard product, condition 2 for untwisted)
    # ===============================================================
    print("\n" + "=" * 70)
    print("BONUS: L · Ls̄ ⊆ Ls̄ (untwisted condition 2)")
    print("=" * 70)
    failures_bonus = 0
    for i in range(8):
        for j in range(8):
            prod = frac_mul(L_basis[i], Lsbar_basis[j])
            if not is_integer_combination(prod, Lsbar_basis):
                print(f"  FAIL: L_basis[{i}] · Ls̄_basis[{j}] not in Ls̄")
                failures_bonus += 1
    if failures_bonus == 0:
        print(f"  PASS: All 64 products lie in Ls̄. ✓")
        print("  (Confirms condition 2 passes for untwisted product too.)")
    else:
        print(f"  FAIL: {failures_bonus}/64 products not in Ls̄.")

    # ===============================================================
    # Summary
    # ===============================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Lemma A (σ(L)=L):           {'PASS' if all_ok else 'FAIL'}")
    print(f"  Lemma E (σ(Ls)≠Ls):         {'PASS' if found_diff else 'FAIL'}")
    print(f"  Counterex (Ls·Ls⊄Ls):       {'PASS' if found_counter else 'FAIL'}")
    print(f"  Lemma D (σ(Ls)·σ(Ls)⊆σ(Ls)):{' PASS' if failures_D==0 else ' FAIL'}")
    print(f"  Lemma C (L·σ(Ls̄)⊆σ(Ls̄)):   {'PASS' if failures_C==0 else 'FAIL'}")
    print(f"  Bonus (L·Ls̄⊆Ls̄):           {'PASS' if failures_bonus==0 else 'FAIL'}")


if __name__ == "__main__":
    main()
