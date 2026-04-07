"""
test_e8_dixon.py — Verification of Dixon's E8 lattice constructions.

Dixon [Dixon2010] defines two families of 240-element sets inside the octonion
algebra, each forming the minimal-vector shell of an E8 lattice representation:

  Ξ^even = Ξ₀ ∪ Ξ₂  (240 elements, norm² = 1)
  A^odd  = A₁ ∪ A₃  (240 elements, norm² = 2)

Properties verified
-------------------
 1.  Ξ₀ count: 16
 2.  Ξ₂ count: 224
 3.  Ξ^even total: 240
 4.  No duplicates in Ξ^even
 5.  Ξ₀ coordinate structure: single ±1 entry (a basis element ±e_a)
 6.  Ξ₂ coordinate structure: exactly 4 nonzero entries, each ±1/2
 7.  All Ξ^even elements have squared norm 1
 8.  All Ξ₂ elements satisfy the defining condition for some ordering:
       e_a(e_b(e_c*e_d)) = ±e_0
 9.  Exactly 14 valid 4-element subsets arising from Ξ₂
     [Dixon2010 Section 2: 7 Type-A sets and 7 Type-B sets]
10.  Ξ^even spans R^8 (rank 8)
11.  2⟨α,β⟩ ∈ Z for all α,β ∈ Ξ^even (root-system integrality condition)
12.  Root system reflection closure: σ_α(β) = β − 2⟨α,β⟩α ∈ ±Ξ^even
     for all α,β ∈ Ξ^even  [E8 Weyl group acts on the root set]
13.  X-product property [Dixon2010 eq. (5)] for Ξ^even:
       (e_a * X) * (X̄ * e_b) = ±e_c  for all X ∈ Ξ^even, all 64 basis pairs
14.  A₁ count: 112
15.  A₃ count: 128
16.  A^odd total: 240
17.  All A^odd elements have squared norm 2
18.  A₁ = Wilson's type-1 roots (same set)
19.  A₃ ∩ Wilson's type-2 = ∅
20.  A₃ ∪ Wilson's type-2 = all 256 sign patterns at scale ½
21.  X-product property [Dixon2010 eq. (5)] for A^odd (norm 2, X^{-1} = X̄/2):
       (e_a * X) * ((X̄/2) * e_b) = ±e_c  for all X ∈ A^odd, all 64 basis pairs

References
----------
[Dixon2010]  G.M. Dixon, "Integral Octonions, Octonion XY-Product, and the
             Leech Lattice", preprint, 2010.
[Wilson2009] R.A. Wilson, "Octonions and the Leech lattice",
             Journal of Algebra 322 (2009) 2186-2190.
"""

import itertools

import numpy as np
import pytest

from octonions import STANDARD_ALGEBRA
from e8_dixon import (
    dixon_xi0,
    dixon_xi2,
    dixon_xi_even,
    dixon_a1,
    dixon_a3,
    dixon_a_odd,
    is_pm_basis_element,
)
from e8_wilson import wilson_type1_roots, wilson_type2_roots


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def xi0():
    return dixon_xi0()


@pytest.fixture(scope="module")
def xi2():
    return dixon_xi2()


@pytest.fixture(scope="module")
def xi_even():
    return dixon_xi_even()


@pytest.fixture(scope="module")
def a1():
    return dixon_a1()


@pytest.fixture(scope="module")
def a3():
    return dixon_a3()


@pytest.fixture(scope="module")
def a_odd():
    return dixon_a_odd()


@pytest.fixture(scope="module")
def xi_even_gram(xi_even):
    """Full 240×240 Gram matrix G[i][j] = ⟨v_i, v_j⟩ for Ξ^even."""
    R = np.array([v.coords for v in xi_even])  # 240 × 8
    return R @ R.T                              # 240 × 240


@pytest.fixture(scope="module")
def wilson_type1_set():
    """Set of Wilson type-1 coordinate tuples, for comparison."""
    return {tuple(r.coords) for r in wilson_type1_roots()}


@pytest.fixture(scope="module")
def wilson_type2_set():
    """Set of Wilson type-2 coordinate tuples, for comparison."""
    return {tuple(r.coords) for r in wilson_type2_roots()}


# ---------------------------------------------------------------------------
# 1–4: Counts and distinctness
# ---------------------------------------------------------------------------

class TestXiEvenCounts:
    """
    Ξ^even = Ξ₀ ∪ Ξ₂ has exactly 16 + 224 = 240 elements.
    [Dixon2010 Section 2, eq. (4)]
    """

    def test_xi0_count(self, xi0):
        """16 = 2 * 8: two signs for each of the 8 basis elements."""
        assert len(xi0) == 16

    def test_xi2_count(self, xi2):
        """
        224 = 14 * 16: 14 valid 4-element subsets, each contributing 2^4 = 16
        sign-pattern elements.
        """
        assert len(xi2) == 224

    def test_xi_even_count(self, xi_even):
        """Total: 16 + 224 = 240."""
        assert len(xi_even) == 240

    def test_xi_even_no_duplicates(self, xi_even):
        """All 240 coordinate tuples are distinct."""
        coord_set = {tuple(v.coords) for v in xi_even}
        assert len(coord_set) == 240


# ---------------------------------------------------------------------------
# 5–6: Coordinate structure
# ---------------------------------------------------------------------------

class TestXiEvenCoordStructure:

    def test_xi0_single_pm1_entry(self, xi0):
        """
        Ξ₀ = {±e_a : a = 0,...,7}: each element has exactly one nonzero
        coordinate equal to ±1.
        """
        for v in xi0:
            nonzero = np.where(np.abs(v.coords) > 1e-12)[0]
            assert len(nonzero) == 1, \
                f"Expected 1 non-zero coord, got {nonzero} in {v.coords}"
            assert abs(abs(v.coords[nonzero[0]]) - 1.0) < 1e-12, \
                f"Non-zero entry not ±1 in {v.coords}"

    def test_xi2_exactly_four_pm_half_entries(self, xi2):
        """
        Ξ₂ elements are (±e_a ± e_b ± e_c ± e_d)/2: exactly 4 nonzero
        coordinates, each equal to ±1/2.
        """
        for v in xi2:
            nonzero = np.where(np.abs(v.coords) > 1e-12)[0]
            assert len(nonzero) == 4, \
                f"Expected 4 non-zero coords, got {len(nonzero)} in {v.coords}"
            assert np.allclose(np.abs(v.coords[nonzero]), 0.5), \
                f"Non-zero entries not all ±1/2 in {v.coords}"


# ---------------------------------------------------------------------------
# 7: Norms
# ---------------------------------------------------------------------------

class TestXiEvenNorms:
    """
    All Ξ^even elements have squared norm 1.

    Ξ₀: N(±e_a) = 1.
    Ξ₂: N((±e_a ± e_b ± e_c ± e_d)/2) = 4 * (1/2)^2 = 1.

    This is the defining property of Ξ^even as a minimal-vector shell of
    an E8 lattice at scale 1/√2.  Reference: [Dixon2010 Section 2].
    """

    def test_xi0_all_norm_sq_1(self, xi0):
        for v in xi0:
            assert abs(v.norm_sq() - 1.0) < 1e-12, \
                f"xi0 element has norm² = {v.norm_sq()}, expected 1"

    def test_xi2_all_norm_sq_1(self, xi2):
        for v in xi2:
            assert abs(v.norm_sq() - 1.0) < 1e-12, \
                f"xi2 element has norm² = {v.norm_sq()}, expected 1"

    def test_xi_even_all_norm_sq_1(self, xi_even):
        for v in xi_even:
            assert abs(v.norm_sq() - 1.0) < 1e-12, \
                f"xi_even element has norm² = {v.norm_sq()}, expected 1"


# ---------------------------------------------------------------------------
# 8: Ξ₂ condition verification
# ---------------------------------------------------------------------------

class TestXi2Condition:
    """
    Each element of Ξ₂ belongs to a quadruple (a,b,c,d) of distinct indices
    satisfying the defining condition e_a(e_b(e_c*e_d)) = ±e_0.
    [Dixon2010 Section 2, eq. (4)]
    """

    def test_xi2_elements_satisfy_condition(self, xi2):
        """
        For each v ∈ Ξ₂, there exists an ordering (a,b,c,d) of the four
        support indices such that e_a*(e_b*(e_c*e_d)) = ±e_0 and the
        coordinates match ±1/2 at positions a,b,c,d.
        """
        alg = STANDARD_ALGEBRA
        for v in xi2:
            support = np.where(np.abs(v.coords) > 1e-12)[0].tolist()
            assert len(support) == 4
            found = False
            for a, b, c, d in itertools.permutations(support):
                prod = (alg.basis_element(a) *
                        (alg.basis_element(b) *
                         (alg.basis_element(c) * alg.basis_element(d))))
                if (abs(abs(prod.coords[0]) - 1.0) < 1e-9 and
                        np.all(np.abs(prod.coords[1:]) < 1e-9)):
                    found = True
                    break
            assert found, \
                f"No valid ordering found for Ξ₂ element with support {support}"


# ---------------------------------------------------------------------------
# 9: Exactly 14 valid 4-element subsets
# ---------------------------------------------------------------------------

class TestXi2ValidSubsets:
    """
    The condition e_a(e_b(e_c*e_d)) = ±e_0 is satisfied by exactly 14
    unordered 4-element subsets of {0,...,7}:
      Type A (7 sets):  {0} ∪ T for each Fano triple T ⊂ {1,...,7}
      Type B (7 sets):  complement of T within {1,...,7} for each Fano triple T

    Reference: [Dixon2010 Section 2]; described in the module docstring of
    e8_dixon.py.
    """

    def _valid_subsets(self):
        """
        Compute all unordered 4-element subsets of {0,...,7} for which some
        ordering satisfies the condition.
        """
        alg = STANDARD_ALGEBRA
        valid = set()
        indices = list(range(8))
        for quad in itertools.combinations(indices, 4):
            for a, b, c, d in itertools.permutations(quad):
                prod = (alg.basis_element(a) *
                        (alg.basis_element(b) *
                         (alg.basis_element(c) * alg.basis_element(d))))
                if (abs(abs(prod.coords[0]) - 1.0) < 1e-9 and
                        np.all(np.abs(prod.coords[1:]) < 1e-9)):
                    valid.add(frozenset(quad))
                    break
        return valid

    def test_exactly_14_valid_subsets(self):
        """There are exactly 14 valid 4-element subsets."""
        valid = self._valid_subsets()
        assert len(valid) == 14, \
            f"Expected 14 valid subsets, found {len(valid)}: {valid}"

    def test_type_a_subsets_present(self):
        """
        Type-A subsets: {0} ∪ T for each of the 7 standard Fano triples.
        Each contains index 0 plus the three members of a Fano triple.
        """
        from octonions import STANDARD_FANO_TRIPLES
        valid = self._valid_subsets()
        for triple in STANDARD_FANO_TRIPLES:
            expected = frozenset({0} | set(triple))
            assert expected in valid, \
                f"Type-A subset {expected} not found in valid subsets"

    def test_type_b_subsets_present(self):
        """
        Type-B subsets: complement of T in {1,...,7} for each Fano triple T.
        The complement of a Fano triple within {1,...,7} has 7 - 3 = 4 elements.
        """
        from octonions import STANDARD_FANO_TRIPLES
        valid = self._valid_subsets()
        all_imaginary = set(range(1, 8))
        for triple in STANDARD_FANO_TRIPLES:
            complement = frozenset(all_imaginary - set(triple))
            assert complement in valid, \
                f"Type-B subset {complement} not found in valid subsets"


# ---------------------------------------------------------------------------
# 10: Rank
# ---------------------------------------------------------------------------

class TestXiEvenRank:
    """Ξ^even spans all of R^8 (rank 8). [Dixon2010 Section 2]"""

    def test_rank_8(self, xi_even):
        R = np.array([v.coords for v in xi_even])
        assert np.linalg.matrix_rank(R) == 8


# ---------------------------------------------------------------------------
# 11: Root-system integrality (2G integer)
# ---------------------------------------------------------------------------

class TestXiEvenRootSystemIntegrality:
    """
    For a root system with norm² = 1, the Cartan integers are 2⟨α,β⟩/⟨α,α⟩ = 2⟨α,β⟩.
    These must all be integers.  Equivalently, every entry of 2G must be an integer,
    where G is the 240×240 Gram matrix of Ξ^even.

    Inner-product values possible:
      ⟨±e_a, ±e_b⟩  = 0 or ±1   →  2⟨α,β⟩ = 0 or ±2  ∈ Z  ✓
      ⟨±e_a, (½-vec)⟩ = ±½        →  2⟨α,β⟩ = ±1       ∈ Z  ✓
      ⟨(½-vec), (½-vec)⟩ ∈ {-1,-½,0,½,1} → 2⟨α,β⟩ ∈ Z  ✓
    """

    def test_2gram_is_integer(self, xi_even_gram):
        """All entries of 2G are integers (to within floating-point tolerance)."""
        two_gram = 2.0 * xi_even_gram
        assert np.allclose(two_gram, np.round(two_gram), atol=1e-10), \
            "2 * Gram matrix of Ξ^even contains non-integer entries"

    def test_diagonal_is_1(self, xi_even_gram):
        """Diagonal entries G[i][i] = ⟨v,v⟩ = 1 for all v ∈ Ξ^even."""
        assert np.allclose(np.diag(xi_even_gram), 1.0, atol=1e-10)


# ---------------------------------------------------------------------------
# 12: Root system reflection closure
# ---------------------------------------------------------------------------

class TestXiEvenReflectionClosure:
    """
    Ξ^even is closed under the Weyl-group reflections of the root system.

    For α,β ∈ Ξ^even:
        σ_α(β) = β − 2⟨α,β⟩α  ∈  Ξ^even.

    This is equivalent to checking that the root system is well-defined
    (each root system must be closed under its own reflections).
    Reference: standard E8 root system definition.
    """

    def test_reflection_closure(self, xi_even):
        """σ_α(β) ∈ Ξ^even for all α,β ∈ Ξ^even."""
        coord_set = {tuple(np.round(v.coords * 2).astype(int)) for v in xi_even}
        for alpha in xi_even:
            for beta in xi_even:
                ip = float(alpha.coords @ beta.coords)   # ⟨α,β⟩
                n = round(2.0 * ip)                      # 2⟨α,β⟩ ∈ Z
                reflected = beta.coords - n * alpha.coords
                key = tuple(np.round(reflected * 2).astype(int))
                assert key in coord_set, \
                    (f"Reflection σ_α(β) not in Ξ^even:\n"
                     f"  α = {alpha.coords}\n  β = {beta.coords}\n"
                     f"  σ_α(β) = {reflected}")


# ---------------------------------------------------------------------------
# 13: X-product property for Ξ^even
# ---------------------------------------------------------------------------

class TestXiEvenXProduct:
    """
    Dixon's X-product property [Dixon2010 eq. (5)]:

      For all X ∈ Ξ^even and all basis elements e_a, e_b:
          (e_a * X) * (X̄ * e_b) = ±e_c

    for some basis element e_c  (i.e. the result is ±1 times a basis element).

    Since all X ∈ Ξ^even have norm² = 1, the inverse is X^{-1} = X̄/N(X) = X̄.
    So the formula becomes (e_a * X) * (X̄ * e_b) directly.
    Reference: [Dixon2010 Section 2, eq. (5)].
    """

    def test_xproduct_gives_pm_basis_element(self, xi_even):
        """
        For all 240 X ∈ Ξ^even and all 64 pairs (a, b), the X-product
        (e_a * X) * (X̄ * e_b) is a ±basis element.
        """
        alg = STANDARD_ALGEBRA
        for x in xi_even:
            x_conj = x.conjugate()
            for a in range(8):
                ea = alg.basis_element(a)
                for b in range(8):
                    eb = alg.basis_element(b)
                    result = (ea * x) * (x_conj * eb)
                    assert is_pm_basis_element(result), \
                        (f"X-product not ±basis element:\n"
                         f"  X = {x.coords}\n"
                         f"  a = {a}, b = {b}\n"
                         f"  result = {result.coords}")


# ---------------------------------------------------------------------------
# 14–16: A^odd counts
# ---------------------------------------------------------------------------

class TestAOddCounts:
    """
    A^odd = A₁ ∪ A₃ has exactly 112 + 128 = 240 elements.
    [Dixon2010 Section 6]
    """

    def test_a1_count(self, a1):
        """112 = C(8,2) * 4 elements in A₁ = {±e_a ± e_b}."""
        assert len(a1) == 112

    def test_a3_count(self, a3):
        """128 = 2^8 / 2 elements in A₃ (even number of minus signs)."""
        assert len(a3) == 128

    def test_a_odd_count(self, a_odd):
        """Total: 112 + 128 = 240."""
        assert len(a_odd) == 240

    def test_a_odd_no_duplicates(self, a_odd):
        """All 240 coordinate tuples are distinct."""
        coord_set = {tuple(v.coords) for v in a_odd}
        assert len(coord_set) == 240


# ---------------------------------------------------------------------------
# 17: A^odd norms
# ---------------------------------------------------------------------------

class TestAOddNorms:
    """
    All A^odd elements have squared norm 2.

    A₁: N(±e_a ± e_b) = 1² + 1² = 2.
    A₃: N((1/2)(±e_0...±e_7)) = 8*(1/2)² = 2.

    Reference: [Dixon2010 Section 6].
    """

    def test_a1_all_norm_sq_2(self, a1):
        for v in a1:
            assert abs(v.norm_sq() - 2.0) < 1e-12, \
                f"A₁ element has norm² = {v.norm_sq()}, expected 2"

    def test_a3_all_norm_sq_2(self, a3):
        for v in a3:
            assert abs(v.norm_sq() - 2.0) < 1e-12, \
                f"A₃ element has norm² = {v.norm_sq()}, expected 2"

    def test_a_odd_all_norm_sq_2(self, a_odd):
        for v in a_odd:
            assert abs(v.norm_sq() - 2.0) < 1e-12, \
                f"A^odd element has norm² = {v.norm_sq()}, expected 2"


# ---------------------------------------------------------------------------
# 18–20: Comparison with Wilson's roots
# ---------------------------------------------------------------------------

class TestAOddVsWilson:
    """
    Relation between Dixon's A^odd and Wilson's E8 roots.

    A₁ consists of ±e_a ± e_b for distinct a,b — identical to Wilson's 112
    type-1 roots.

    A₃ uses EVEN number of minus signs among the 8 coordinates, whereas
    Wilson's type-2 roots use an ODD number of minus signs.  The two sets are
    therefore disjoint, and together they cover all 2^8 = 256 sign patterns.

    References: [Dixon2010 Section 6]; [Wilson2009 Section 2].
    """

    def test_a1_equals_wilson_type1(self, a1, wilson_type1_set):
        """
        A₁ = Wilson's type-1 roots (same 112-element set).
        """
        a1_set = {tuple(v.coords) for v in a1}
        assert a1_set == wilson_type1_set, \
            "A₁ and Wilson's type-1 roots are not identical"

    def test_a3_disjoint_from_wilson_type2(self, a3, wilson_type2_set):
        """
        A₃ (even minus-count) and Wilson's type-2 (odd minus-count) are disjoint.
        """
        a3_set = {tuple(v.coords) for v in a3}
        overlap = a3_set & wilson_type2_set
        assert len(overlap) == 0, \
            f"A₃ and Wilson's type-2 roots overlap in {len(overlap)} elements"

    def test_a3_union_wilson_type2_covers_all_sign_patterns(self, a3, wilson_type2_set):
        """
        A₃ ∪ Wilson's type-2 = all 256 sign patterns at scale ½.

        Since A₃ uses even minus-count (128 patterns) and Wilson's type-2 uses
        odd minus-count (128 patterns), together they cover all 2^8 = 256 patterns.
        """
        a3_set = {tuple(v.coords) for v in a3}
        union = a3_set | wilson_type2_set
        all_patterns = set()
        for mask in range(256):
            coords = tuple((-0.5 if (mask >> k) & 1 else 0.5) for k in range(8))
            all_patterns.add(coords)
        assert union == all_patterns, \
            f"A₃ ∪ Wilson type-2 does not cover all 256 sign patterns"


# ---------------------------------------------------------------------------
# 21: X-product property for A^odd
# ---------------------------------------------------------------------------

class TestAOddXProduct:
    """
    Dixon's X-product property [Dixon2010 eq. (5)] for A^odd:

      For all X ∈ A^odd and all basis elements e_a, e_b:
          (e_a * X) * (X^{-1} * e_b) = ±e_c

    for some basis element e_c.

    Since all X ∈ A^odd have norm² = 2, the inverse is X^{-1} = X̄/N(X) = X̄/2.
    So the formula becomes (e_a * X) * ((X̄/2) * e_b) = ±e_c.
    Reference: [Dixon2010 Section 2, eq. (5); Section 6].
    """

    def test_xproduct_gives_pm_basis_element(self, a_odd):
        """
        For all 240 X ∈ A^odd and all 64 pairs (a, b), the X-product
        (e_a * X) * ((X̄/2) * e_b) is a ±basis element.
        """
        alg = STANDARD_ALGEBRA
        for x in a_odd:
            # X^{-1} = X̄ / N(X) = X̄ / 2  (since N(X) = 2 for all X in A^odd)
            x_inv = alg.element(x.conjugate().coords / 2.0)
            for a in range(8):
                ea = alg.basis_element(a)
                for b in range(8):
                    eb = alg.basis_element(b)
                    result = (ea * x) * (x_inv * eb)
                    assert is_pm_basis_element(result), \
                        (f"X-product not ±basis element:\n"
                         f"  X = {x.coords}\n"
                         f"  a = {a}, b = {b}\n"
                         f"  result = {result.coords}")
