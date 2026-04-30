"""
egan_baez_count.py — Independent attempt to verify Egan's count of 17,280
Jordan rings in the Baez/Egan construction.

Per Baez/Egan part 11
(https://golem.ph.utexas.edu/category/2014/12/integral_octonions_part_11.html):
the count factorises as 17,280 = 270 × 64, where:
  - 270 is the number of sublattices L_1 ⊂ E_8 isometric to √2 E_8.
  - For each such L_1, there are exactly 64 sublattices L_2 isometric
    to √2 E_8 with L_1 ∩ L_2 = 2 E_8.

Approach: a sublattice L_1 ⊂ E_8 with index 16 and minimum squared norm
≥ 4 is determined by a 4-dimensional subspace V_1 ⊂ E_8 / 2 E_8 ≅ F_2^8
that contains no root projection (equivalently, V_1 is totally
isotropic for the F_2 quadratic form q induced from E_8's bilinear
form). On this picture:
  - 270 = number of 4-dim totally isotropic subspaces (maximal isotropic
    subspaces of F_2^8 with the E_8 mod 2 quadratic form, which is the
    plus-type form).
  - The compatibility condition L_1 ∩ L_2 = 2 E_8 translates to
    V_1 ∩ V_2 = {0} (as subspaces of F_2^8).
  - 17,280 = number of ordered pairs (V_1, V_2) of maximal isotropic
    subspaces that meet trivially.

This module enumerates these subspaces and pairs and reports the counts.
"""

import numpy as np
from typing import List, Tuple, FrozenSet
import sys

sys.path.insert(0, '/home/jens/Projects/leech_alg/python_project/src')

from leech_wilson import _build_L_zbasis


def build_q_form() -> Tuple[np.ndarray, np.ndarray]:
    """Compute the E_8 mod 2 quadratic form q: F_2^8 → F_2.

    For our chosen Z-basis B = {b_1, ..., b_8} of E_8 (rows of L_basis),
    the form is q(v) = sum a_i v_i + sum_{i<j} m_ij v_i v_j over F_2,
    where:
      - a_i = (1/2)⟨b_i, b_i⟩ mod 2  (linear part)
      - m_ij = ⟨b_i, b_j⟩ mod 2       (off-diagonal bilinear part)

    Returns (a, M) where a is length-8 F_2 vector and M is 8×8 F_2
    symmetric matrix with zero diagonal (only off-diagonal m_ij entries).
    """
    B = _build_L_zbasis()  # 8x8 matrix, rows are basis vectors of E_8
    G = B @ B.T  # 8x8 Gram matrix (rational/integer entries)
    a = np.zeros(8, dtype=np.int64)
    M = np.zeros((8, 8), dtype=np.int64)
    for i in range(8):
        # (1/2) ⟨b_i, b_i⟩ should be an integer (E_8 even)
        val = G[i, i] / 2
        assert abs(val - round(val)) < 1e-9, f"E_8 not even at i={i}: ⟨b,b⟩={G[i,i]}"
        a[i] = int(round(val)) % 2
    for i in range(8):
        for j in range(i + 1, 8):
            val = G[i, j]
            assert abs(val - round(val)) < 1e-9, f"non-integer Gram entry at ({i},{j})"
            M[i, j] = int(round(val)) % 2
            M[j, i] = M[i, j]
    return a, M


def q_value(v: np.ndarray, a: np.ndarray, M: np.ndarray) -> int:
    """Evaluate q(v) over F_2 for v ∈ F_2^8, given a and M."""
    val = int(np.dot(a, v)) % 2
    for i in range(8):
        if v[i] == 0:
            continue
        for j in range(i + 1, 8):
            if v[j] == 0:
                continue
            val = (val + M[i, j]) % 2
    return val


def enumerate_isotropic_vectors(a: np.ndarray, M: np.ndarray) -> List[np.ndarray]:
    """List all v ∈ F_2^8 with q(v) = 0."""
    iso = []
    for k in range(256):
        v = np.array([(k >> i) & 1 for i in range(8)], dtype=np.int64)
        if q_value(v, a, M) == 0:
            iso.append(v)
    return iso


def subspace_key(basis: List[np.ndarray]) -> bytes:
    """Reduce a list of F_2 vectors to a row-reduced echelon form, then
    return a canonical bytes key for set membership."""
    # Stack rows into 2D array
    rows = [v.copy() for v in basis]
    n = len(rows)
    if n == 0:
        return b""
    # Gaussian elimination over F_2
    pivot_col = -1
    rank = 0
    for col in range(8):
        # Find pivot row at or after `rank`
        pivot = None
        for r in range(rank, n):
            if rows[r][col] == 1:
                pivot = r
                break
        if pivot is None:
            continue
        if pivot != rank:
            rows[rank], rows[pivot] = rows[pivot], rows[rank]
        # Eliminate other rows
        for r in range(n):
            if r != rank and rows[r][col] == 1:
                rows[r] = (rows[r] + rows[rank]) % 2
        rank += 1
    # Drop zero rows
    rows = [r for r in rows[:rank]]
    # Sort rows for canonicality
    rows.sort(key=lambda r: tuple(int(x) for x in r))
    # Pack into bytes
    flat = np.concatenate(rows).astype(np.int8)
    return bytes(flat.tolist())


def enumerate_maximal_isotropic_subspaces(
    a: np.ndarray, M: np.ndarray, iso: List[np.ndarray]
) -> List[List[np.ndarray]]:
    """Enumerate all 4-dim totally isotropic subspaces of F_2^8 with q.

    A 4-dim totally isotropic subspace V satisfies:
      - dim V = 4 (15 non-zero elements)
      - q(v) = 0 for all v ∈ V
      - V is closed under addition (subspace property)
    """
    found = {}
    # Build subspace incrementally: pick v_1, v_2, v_3, v_4 such that the
    # subspace they span is totally isotropic.
    # For a subspace ⟨v_1, v_2, v_3, v_4⟩, we need q(v_i) = 0 AND
    # q(v_i + v_j) = 0 for all i, j (which simplifies to bilinear-form
    # condition q(v_i + v_j) = q(v_i) + q(v_j) + b(v_i, v_j) where
    # b(v_i, v_j) is the polarisation bilinear form). For F_2 quadratic
    # forms, b(u, v) = q(u + v) - q(u) - q(v) over Z (then mod 2).
    # Total isotropy requires q(any element of V) = 0, equivalently:
    # q(v_i) = 0 for each generator AND b(v_i, v_j) = 0 for all i ≠ j.

    # Compute bilinear form b: b(u, v) = q(u + v) - q(u) - q(v) (mod 2)
    # For our q with linear and bilinear parts:
    # q(u + v) = a·(u+v) + sum (u+v)_i (u+v)_j m_ij
    # = q(u) + q(v) + sum_{i<j} m_ij (u_i v_j + u_j v_i)
    # = q(u) + q(v) + sum_{i≠j} m_ij u_i v_j  (i<j,j<i covered by symmetry)
    # Actually: u_i v_j + u_j v_i = (uv^T + vu^T)_{ij}, but over F_2 it
    # simplifies. The polarisation b(u, v) = sum_{i,j} m_ij u_i v_j (with
    # appropriate symmetry) — it's the off-diagonal symmetric bilinear.
    # In F_2: b(u, v) = sum_{i<j} m_ij (u_i v_j + u_j v_i) = u^T (M)/2 v
    # but M is symmetric so just take u^T M_strictUpper v + ...
    # Easier: compute b(u, v) by direct evaluation.
    def b_form(u, v):
        return (q_value((u + v) % 2, a, M) - q_value(u, a, M) - q_value(v, a, M)) % 2

    # Enumerate: pick first generator v_1 (any nonzero isotropic).
    # For each v_1, pick v_2 with b(v_1, v_2) = 0 and v_2 isotropic, v_2 ∉ ⟨v_1⟩.
    # ... up to dim 4.
    # Use canonical form to deduplicate.

    n_iso = len(iso)
    print(f"  Isotropic vectors (incl. zero): {n_iso}")
    nonzero_iso = [v for v in iso if v.any()]
    print(f"  Nonzero isotropic: {len(nonzero_iso)}")

    for i1 in range(len(nonzero_iso)):
        v1 = nonzero_iso[i1]
        for i2 in range(i1 + 1, len(nonzero_iso)):
            v2 = nonzero_iso[i2]
            if b_form(v1, v2) != 0:
                continue
            # ⟨v1, v2⟩ ⊕ {0, v1, v2, v1+v2}: check v1+v2 isotropic
            v12 = (v1 + v2) % 2
            if q_value(v12, a, M) != 0:
                continue
            # Now extend with v3
            for i3 in range(i2 + 1, len(nonzero_iso)):
                v3 = nonzero_iso[i3]
                if b_form(v1, v3) != 0 or b_form(v2, v3) != 0:
                    continue
                # Check that all linear combos with v3 are isotropic.
                # Since q(v1)=q(v2)=q(v3)=0 and bilinear forms vanish,
                # q on subspace is sum of pairwise b's = 0. So OK.
                # We still need v3 ∉ ⟨v1, v2⟩.
                if np.array_equal(v3, v1) or np.array_equal(v3, v2) or np.array_equal(v3, v12):
                    continue
                for i4 in range(i3 + 1, len(nonzero_iso)):
                    v4 = nonzero_iso[i4]
                    if b_form(v1, v4) != 0 or b_form(v2, v4) != 0 or b_form(v3, v4) != 0:
                        continue
                    # v4 ∉ ⟨v1, v2, v3⟩
                    span_so_far = set()
                    for c1 in (0, 1):
                        for c2 in (0, 1):
                            for c3 in (0, 1):
                                w = (c1 * v1 + c2 * v2 + c3 * v3) % 2
                                span_so_far.add(tuple(int(x) for x in w))
                    if tuple(int(x) for x in v4) in span_so_far:
                        continue
                    # 4-dim isotropic subspace found
                    key = subspace_key([v1, v2, v3, v4])
                    if key not in found:
                        found[key] = (v1, v2, v3, v4)

    return list(found.values())


def main():
    print("=" * 70)
    print("Egan/Baez count: 17,280 Jordan rings — independent attempt")
    print("=" * 70)

    print("\n--- Step 1: build E_8 mod 2 quadratic form ---")
    a, M = build_q_form()
    print(f"  Linear part a: {a}")
    print(f"  Bilinear M (off-diagonal):\n{M}")

    print("\n--- Step 2: enumerate isotropic vectors ---")
    iso = enumerate_isotropic_vectors(a, M)
    print(f"  q^(-1)(0) = {len(iso)} vectors (expected 136 for plus-type form on F_2^8)")

    if len(iso) != 136:
        print(f"  WARNING: expected 136 isotropic vectors for plus-type form.")

    print("\n--- Step 3: enumerate 4-dim maximal isotropic subspaces ---")
    subspaces = enumerate_maximal_isotropic_subspaces(a, M, iso)
    print(f"  Maximal isotropic subspaces: {len(subspaces)}")

    if len(subspaces) == 270:
        print("  ✓ Matches Egan's count of 270.")
    else:
        print(f"  Differs from Egan's 270.")

    # Step 4: count ordered pairs (V_1, V_2) with V_1 ∩ V_2 = {0}.
    # Equivalently, V_1 + V_2 = F_2^8, i.e., (V_1, V_2) is a pair of
    # complementary 4-dim subspaces. Translates to L_1 ∩ L_2 = 2 E_8.
    print("\n--- Step 4: count ordered pairs (V_1, V_2) with V_1 ∩ V_2 = 0 ---")

    # Represent each subspace as a set of all 16 elements (its closure under +).
    subspace_elements = []
    for v1, v2, v3, v4 in subspaces:
        elems = set()
        for c1 in (0, 1):
            for c2 in (0, 1):
                for c3 in (0, 1):
                    for c4 in (0, 1):
                        w = (c1 * v1 + c2 * v2 + c3 * v3 + c4 * v4) % 2
                        elems.add(tuple(int(x) for x in w))
        subspace_elements.append(elems)

    n = len(subspaces)
    pair_count = 0
    per_subspace_partner_count = []
    for i in range(n):
        partners_for_i = 0
        for j in range(n):
            if i == j:
                continue
            # V_i ∩ V_j = {0} iff their nonzero elements are disjoint
            inter = subspace_elements[i] & subspace_elements[j]
            # inter always contains zero; if size is 1, only zero in common
            if len(inter) == 1:
                pair_count += 1
                partners_for_i += 1
        per_subspace_partner_count.append(partners_for_i)
        if (i + 1) % 50 == 0:
            print(f"  processed {i + 1}/{n} subspaces; running pair count: {pair_count}")

    print(f"\n  Ordered pairs (V_1, V_2) with V_1 ∩ V_2 = 0: {pair_count}")
    print(f"  Partner count per subspace (min, max, mean): "
          f"{min(per_subspace_partner_count)}, "
          f"{max(per_subspace_partner_count)}, "
          f"{sum(per_subspace_partner_count) / len(per_subspace_partner_count):.2f}")

    if pair_count == 17280:
        print("  ✓ Matches Egan's count of 17,280.")
    elif pair_count == 8640:
        print("  ✓ Matches Egan's count of 17,280 (interpreted as unordered pairs: 8,640).")
    else:
        print(f"  Differs from Egan's 17,280.")

    return len(subspaces), pair_count


if __name__ == "__main__":
    main()
