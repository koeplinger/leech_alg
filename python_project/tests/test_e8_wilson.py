"""
test_e8_wilson.py — Verification of Wilson's E8 lattice construction.

Wilson [Wilson2009] defines an E8 lattice L as the Z-span of 240 octonions
(Section 2).  This module verifies that the construction has all expected
properties of the E8 root lattice, and additionally checks the structural
identities involving Wilson's element s that Wilson himself proves.

Properties verified
-------------------
 1.  Type-1 root count: 112
 2.  Type-2 root count: 128
 3.  Total root count: 240
 4.  No duplicate roots
 5.  Type-1 coordinate structure: exactly two non-zero entries, each ±1
 6.  Type-2 coordinate structure: all entries ±1/2
 7.  All type-1 roots have squared norm 2
 8.  All type-2 roots have squared norm 2
 9.  All 240 roots have squared norm 2
10.  All roots satisfy the L membership criterion (is_in_L)
11.  Single basis elements e_0,...,e_7 are NOT in L
12.  Scaled basis elements 2*e_i ARE in L
13.  The element (1/2)(e_0+...+e_7) (0 minus signs, even) is NOT in L
14.  Full Gram matrix (240×240) has all-integer entries
     [Required for L to be an integer lattice]
15.  Rank 8: 8 linearly independent roots exist
16.  Gram matrix of a basis has determinant 1 (L is unimodular / self-dual)
     [Wilson2009 Section 2: "L is self-dual"]
17.  Wilson's element s = (1/2)(-e_0 + e_1 + ... + e_7):
       a. Coordinates are (-1/2, 1/2, ..., 1/2)
       b. N(s) = 2
       c. s is a type-2 root (1 minus sign = odd)
       d. s-bar has N(s-bar) = 2
18.  Even lattice: sample of sums and differences of roots have even norm
19.  Ls ⊆ L: r*s ∈ L for all 240 roots r
     [Wilson2009 Section 2: "Ls ⊆ LR = 2B ⊂ L"]
20.  2L ⊆ Ls: r*s-bar ∈ L for all 240 roots r
     [Wilson2009 Section 2: "2L ⊂ Ls"; by right-alternativity, 2r ∈ Ls iff
      (2r)*s^{-1} ∈ L; since s^{-1} = s-bar/N(s) = s-bar/2, this becomes
      2r*(s-bar/2) = r*s-bar ∈ L]

References
----------
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186-2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
"""

import numpy as np
import pytest

from octonions import STANDARD_ALGEBRA
from e8_wilson import (
    wilson_e8_roots,
    wilson_type1_roots,
    wilson_type2_roots,
    WILSON_S,
    inner_product,
    is_in_L,
)


# ---------------------------------------------------------------------------
# Fixtures / cached root lists
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def all_roots():
    return wilson_e8_roots()


@pytest.fixture(scope="module")
def type1_roots():
    return wilson_type1_roots()


@pytest.fixture(scope="module")
def type2_roots():
    return wilson_type2_roots()


@pytest.fixture(scope="module")
def gram_matrix(all_roots):
    """Full 240×240 Gram matrix G[i][j] = <root_i, root_j>."""
    R = np.array([r.coords for r in all_roots])  # 240 x 8
    return R @ R.T                                # 240 x 240


# ---------------------------------------------------------------------------
# 1–4: Root counts and distinctness
# ---------------------------------------------------------------------------

class TestRootCounts:
    """Wilson describes exactly 112 + 128 = 240 roots. [Wilson2009 Section 2]"""

    def test_type1_count(self, type1_roots):
        """112 = C(8,2) * 4 type-1 roots."""
        assert len(type1_roots) == 112

    def test_type2_count(self, type2_roots):
        """128 = 2^8 / 2 type-2 roots (odd number of minus signs)."""
        assert len(type2_roots) == 128

    def test_total_count(self, all_roots):
        """Total: 112 + 128 = 240 roots."""
        assert len(all_roots) == 240


class TestRootsDistinct:
    def test_no_duplicate_roots(self, all_roots):
        """All 240 coordinate tuples are distinct."""
        coord_set = {tuple(r.coords) for r in all_roots}
        assert len(coord_set) == 240


# ---------------------------------------------------------------------------
# 5–6: Coordinate structure
# ---------------------------------------------------------------------------

class TestRootCoordinateStructure:

    def test_type1_have_exactly_two_nonzero_entries(self, type1_roots):
        """
        Type-1 roots ±e_a ± e_b have exactly two nonzero coordinates,
        each equal to ±1.
        """
        for r in type1_roots:
            nonzero = np.where(np.abs(r.coords) > 1e-12)[0]
            assert len(nonzero) == 2, f"Expected 2 non-zero coords, got {nonzero} in {r}"
            assert np.allclose(np.abs(r.coords[nonzero]), 1.0)

    def test_type2_all_coords_are_pm_half(self, type2_roots):
        """
        Type-2 roots (1/2)(±e_0 ... ±e_7) have all eight coordinates equal
        to ±1/2.
        """
        for r in type2_roots:
            assert np.allclose(np.abs(r.coords), 0.5), \
                f"Expected all |coords| = 1/2, got {r.coords}"


# ---------------------------------------------------------------------------
# 7–9: Root norms
# ---------------------------------------------------------------------------

class TestRootNorms:
    """
    All roots have squared norm 2.

    Type-1: N(±e_a ± e_b) = 1^2 + 1^2 = 2.
    Type-2: N((1/2)(±e_0...±e_7)) = 8 * (1/2)^2 = 2.

    This is the defining property of minimal vectors (roots) of E8.
    Reference: [Wilson2009] Section 2 — "vectors of minimal norm".
    """

    def test_type1_roots_have_norm_2(self, type1_roots):
        for r in type1_roots:
            assert abs(r.norm_sq() - 2.0) < 1e-12

    def test_type2_roots_have_norm_2(self, type2_roots):
        for r in type2_roots:
            assert abs(r.norm_sq() - 2.0) < 1e-12

    def test_all_roots_have_norm_2(self, all_roots):
        for r in all_roots:
            assert abs(r.norm_sq() - 2.0) < 1e-12


# ---------------------------------------------------------------------------
# 10–13: Lattice membership criterion
# ---------------------------------------------------------------------------

class TestLMembership:
    """
    Verifies the is_in_L predicate on known elements.

    L = D8+ = { all-integer coords, even sum }
            u { all half-integer coords, odd sum }.
    """

    def test_all_roots_are_in_L(self, all_roots):
        for r in all_roots:
            assert is_in_L(r), f"Root {r} not recognised as being in L"

    def test_single_basis_elements_not_in_L(self):
        """
        Each individual basis element e_i has integer coords summing to 1
        (odd), so it is NOT in L.
        """
        for i in range(8):
            e = STANDARD_ALGEBRA.basis_element(i)
            assert not is_in_L(e), f"e_{i} should not be in L"

    def test_doubled_basis_elements_are_in_L(self):
        """
        2*e_i has integer coords summing to 2 (even), so it IS in L.
        (E.g. 2*e_1 = (e_1+e_2) + (e_1-e_2) ∈ L since both summands are roots.)
        """
        for i in range(8):
            e2 = STANDARD_ALGEBRA.element(
                2.0 * STANDARD_ALGEBRA.basis_element(i).coords
            )
            assert is_in_L(e2), f"2*e_{i} should be in L"

    def test_all_plus_halfint_not_in_L(self):
        """
        (1/2)(e_0 + e_1 + ... + e_7): all-half-integer coords summing to 4
        (even integer), so NOT in L.
        """
        v = STANDARD_ALGEBRA.element(np.full(8, 0.5))
        assert not is_in_L(v)


# ---------------------------------------------------------------------------
# 14: Integer Gram matrix
# ---------------------------------------------------------------------------

class TestIntegerInnerProducts:
    """
    The full 240×240 Gram matrix has all-integer entries.

    This is necessary (and, together with all norms being even, sufficient) for
    L to be an even integer lattice.  By bilinearity, integer inner products on
    the root generators imply integer inner products on all of L.
    Reference: standard requirement for a lattice to be integral.
    """

    def test_full_gram_matrix_is_integer(self, gram_matrix):
        """All 240*240 inner products are integers (checked via rounding error)."""
        assert np.allclose(gram_matrix, np.round(gram_matrix), atol=1e-10), \
            "Gram matrix contains non-integer entries"


# ---------------------------------------------------------------------------
# 15–16: Rank and self-duality
# ---------------------------------------------------------------------------

class TestRankAndSelfDuality:

    def test_rank_8(self, all_roots):
        """
        The 240 roots span an 8-dimensional space.
        Reference: [Wilson2009] Section 2 — roots "span the lattice".
        """
        R = np.array([r.coords for r in all_roots])
        assert np.linalg.matrix_rank(R) == 8

    def test_gram_matrix_det_is_1(self, all_roots, type2_roots):
        """
        A Z-basis of L has Gram matrix with determinant 1 (unimodular).

        L = D8+ = E8 is self-dual (unimodular), so the Gram matrix of any
        Z-basis has det = 1.  Reference: [Wilson2009] Section 2 — "L is self-dual".

        Construction of a valid Z-basis:
        L = D8+ consists of two cosets of D8 (the integer even-sum sublattice).
        A real basis found from type-1 roots alone spans only D8 (a sublattice
        of index 2, Gram det = 4).  Including a type-2 root (here WILSON_S) as
        the first element ensures the basis covers the full D8+ lattice.

        Verification that det(B)=1:
        With basis {WILSON_S, e0+e1, e0-e1, e0+e2, ..., e0+e6}, column-expansion
        along the e7 column (only non-zero entry from WILSON_S) gives det(B) = 1.
        """
        # Build the basis: WILSON_S first (type-2, covers the D8+ coset),
        # then 7 independent type-1 roots.
        basis_coords = [WILSON_S.coords.copy()]
        for r in all_roots:
            candidate = basis_coords + [r.coords]
            if np.linalg.matrix_rank(np.array(candidate)) == len(candidate):
                basis_coords.append(r.coords)
                if len(basis_coords) == 8:
                    break
        assert len(basis_coords) == 8, "Could not find 8 independent roots"
        B = np.array(basis_coords)   # 8 × 8
        G = B @ B.T
        det = np.linalg.det(G)
        assert abs(abs(det) - 1.0) < 1e-6, \
            f"Gram matrix det = {det:.6f}, expected 1 (self-dual / unimodular lattice)"


# ---------------------------------------------------------------------------
# 17: Wilson's element s
# ---------------------------------------------------------------------------

class TestWilsonElementS:
    """
    Tests on Wilson's special element s = (1/2)(-e_0 + e_1 + ... + e_7).

    Wilson writes s = (1/2)(-1 + i_0 + ... + i_6) and proves s ∈ L, s-bar ∈ R.
    Under our index map this becomes s = (1/2)(-e_0 + e_1 + ... + e_7).
    Reference: [Wilson2009] Section 2.
    """

    def test_s_coords(self):
        """s has coordinates (-1/2, 1/2, 1/2, ..., 1/2)."""
        expected = np.array([-0.5] + [0.5] * 7)
        assert np.allclose(WILSON_S.coords, expected)

    def test_s_has_norm_2(self):
        """N(s) = 2."""
        assert abs(WILSON_S.norm_sq() - 2.0) < 1e-12

    def test_s_is_type2_root(self, type2_roots):
        """s has exactly 1 minus sign (odd), so it is a type-2 root."""
        type2_set = {tuple(r.coords) for r in type2_roots}
        assert tuple(WILSON_S.coords) in type2_set

    def test_s_bar_has_norm_2(self):
        """N(s-bar) = N(s) = 2 (conjugation preserves the norm)."""
        assert abs(WILSON_S.conjugate().norm_sq() - 2.0) < 1e-12


# ---------------------------------------------------------------------------
# 18: Even lattice
# ---------------------------------------------------------------------------

class TestEvenLattice:
    """
    L is an even lattice: N(v) is even for all v in L.

    By the polarisation identity  N(r1+r2) = N(r1) + N(r2) + 2<r1,r2>
    and the fact that N(root) = 2 and <r1,r2> ∈ Z (verified separately),
    every integer combination of roots has even norm.  We verify this for
    a sample of lattice sums and differences.
    """

    def test_sums_of_roots_have_even_norm(self, all_roots):
        sample = all_roots[:20]
        for r1 in sample:
            for r2 in sample:
                n = (r1 + r2).norm_sq()
                n_int = round(n)
                assert abs(n - n_int) < 1e-12 and n_int % 2 == 0, \
                    f"Odd norm {n} for sum of roots"

    def test_differences_of_roots_have_even_norm(self, all_roots):
        sample = all_roots[:20]
        for r1 in sample:
            for r2 in sample:
                n = (r1 - r2).norm_sq()
                n_int = round(n)
                assert abs(n - n_int) < 1e-12 and n_int % 2 == 0, \
                    f"Odd norm {n} for difference of roots"


# ---------------------------------------------------------------------------
# 19: Ls ⊆ L
# ---------------------------------------------------------------------------

class TestLsSubsetL:
    """
    For all roots r in L, the product r*s lies in L.

    Wilson proves Ls ⊆ LR = 2B ⊂ L  (where R = L-bar and B is the
    Coxeter-Dickson ring).  We verify the outer inclusion Ls ⊆ L for all
    240 roots.  By linearity over Z, this implies Ls ⊆ L for the full lattice.

    Reference: [Wilson2009] Section 2.
    """

    def test_rs_in_L_for_all_roots(self, all_roots):
        for r in all_roots:
            rs = r * WILSON_S
            assert is_in_L(rs), \
                f"r*s not in L:\n  r = {r}\n  r*s = {rs}"


# ---------------------------------------------------------------------------
# 20: 2L ⊆ Ls
# ---------------------------------------------------------------------------

class TestTwoLSubsetLs:
    """
    For all roots r in L, the element 2r lies in Ls (i.e. Ls contains 2L).

    Wilson proves 2L ⊂ Ls ⊂ L  [Wilson2009 Section 2].  We verify 2r ∈ Ls
    for all 240 roots r via the following equivalence derived from right-
    alternativity of the octonion algebra:

        2r ∈ Ls  iff  (2r) * s^{-1} ∈ L.

    Since s^{-1} = s-bar / N(s) = s-bar / 2, we have:
        (2r) * (s-bar/2) = r * s-bar ∈ L.

    So: 2r ∈ Ls  iff  r * s-bar ∈ L.

    Note: right-alternativity gives  (v * s) * s^{-1} = v * (s * s^{-1}) = v,
    so (r*s)*s^{-1} = r and (2r)*s^{-1} = 2*(r*s^{-1})... wait, we use it
    the other way: if w = u*s then w*s^{-1} = (u*s)*s^{-1} = u*(s*s^{-1}) = u.
    Hence w ∈ Ls iff w*s^{-1} ∈ L.  Taking w = 2r:
        (2r)*s^{-1} = (2r)*(s-bar/2) = r*s-bar.

    Reference: [Wilson2009] Section 2 — "2L ⊂ Ls".
    """

    def test_2r_in_Ls_for_all_roots(self, all_roots):
        s_bar = WILSON_S.conjugate()
        for r in all_roots:
            r_sbar = r * s_bar
            assert is_in_L(r_sbar), \
                f"r*s-bar not in L (so 2r not in Ls):\n  r = {r}\n  r*s-bar = {r_sbar}"
