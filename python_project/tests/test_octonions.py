"""
test_octonions.py — Verification that STANDARD_ALGEBRA and XPRODUCT_ALGEBRA
satisfy all defining properties of the octonion algebra.

Both algebras are exercised by every test via the `alg` fixture, which is
parametrised over [STANDARD_ALGEBRA, XPRODUCT_ALGEBRA].

Properties verified
-------------------
 1.  Identity element
 2.  Imaginary squares: e_i^2 = -e_0
 3.  Anticommutativity of distinct imaginary basis pairs
 4.  Non-commutativity (sanity: the algebra is NOT commutative)
 5.  Non-associativity (sanity: the algebra is NOT associative)
 6.  Composition algebra: N(x*y) = N(x)*N(y)
 7.  Conjugate formula: x * conj(x) = conj(x) * x = N(x) * e_0
 8.  Conjugate is an involution: conj(conj(x)) = x
 9.  Left alternativity: (x*x)*y = x*(x*y)
10.  Right alternativity: x*(y*y) = (x*y)*y
11.  Flexibility: (x*y)*x = x*(y*x)
12.  Total antisymmetry of the associator (over all 512 ordered basis triples)
13.  Moufang identity: z*(x*(z*y)) = ((z*x)*z)*y
14.  Division algebra: every nonzero element has a two-sided inverse
15.  No zero divisors among basis elements
16.  Power-associativity: x^3 and x^4 behave correctly
17.  The two algebras have different multiplication tables
18.  Cross-check: X-product table agrees with Dixon's formula
         A o_X B = (A * l0) * (l0^{-1} * B)  where l0 = (e0+...+e7)/2

Notes on scope
--------------
Properties 9-11 are tested over all 64 ordered pairs of basis elements.
This is a complete check by bilinearity: if the identity holds for all pairs
of basis elements it holds for all real linear combinations.
(More precisely: left alternativity (x*x)*y = x*(x*y) linearises to the
condition that the associator is antisymmetric in its first two arguments,
which is tested exhaustively in test 12 over all 512 basis triples.)

Property 12 (associator antisymmetry) over all 8^3 = 512 basis triples
is the most comprehensive test: by trilinearity it certifies the identity
for ALL real octonion elements, not just basis elements.

References
----------
[Dixon2010]   G.M. Dixon, "Integral Octonions, Octonion XY-Product, and the
              Leech Lattice", preprint, 2010.
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186-2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
[Wikipedia]   https://en.wikipedia.org/wiki/Octonion
[Schafer1966] R.D. Schafer, Introduction to Nonassociative Algebras,
              Academic Press, 1966.
"""

import itertools

import numpy as np
import pytest

from octonions import (
    STANDARD_ALGEBRA,
    XPRODUCT_ALGEBRA,
    Octonion,
    OctonionAlgebra,
    associator,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

ALGEBRAS = [STANDARD_ALGEBRA, XPRODUCT_ALGEBRA]
ALGEBRA_IDS = ["standard", "xproduct"]

ATOL = 1e-10


@pytest.fixture(params=ALGEBRAS, ids=ALGEBRA_IDS)
def alg(request) -> OctonionAlgebra:
    return request.param


# ---------------------------------------------------------------------------
# 1. Identity element
# Property: e_0 * x = x * e_0 = x  for all x.
# Reference: [Wikipedia] Octonion - basic algebra axioms
# ---------------------------------------------------------------------------

class TestIdentity:

    def test_left_identity(self, alg: OctonionAlgebra) -> None:
        e0 = alg.basis_element(0)
        for i in range(8):
            ei = alg.basis_element(i)
            assert (e0 * ei).is_close(ei, ATOL), f"e0 * e{i} != e{i}"

    def test_right_identity(self, alg: OctonionAlgebra) -> None:
        e0 = alg.basis_element(0)
        for i in range(8):
            ei = alg.basis_element(i)
            assert (ei * e0).is_close(ei, ATOL), f"e{i} * e0 != e{i}"


# ---------------------------------------------------------------------------
# 2. Squares of imaginary basis elements
# Property: e_i^2 = -e_0  for i = 1,...,7.
# Reference: [Wikipedia] Octonion - Multiplication
# ---------------------------------------------------------------------------

class TestImaginarySquares:

    def test_squares_equal_minus_identity(self, alg: OctonionAlgebra) -> None:
        e0 = alg.basis_element(0)
        for i in range(1, 8):
            ei = alg.basis_element(i)
            expected = (-1.0) * e0
            assert (ei * ei).is_close(expected, ATOL), f"e{i}^2 != -e0"


# ---------------------------------------------------------------------------
# 3. Anticommutativity of distinct imaginary basis pairs
# Property: e_i * e_j = -(e_j * e_i)  for i != j,  i,j in {1,...,7}.
# The real unit e_0 commutes with everything.
# Reference: [Wikipedia] Octonion - Multiplication
# ---------------------------------------------------------------------------

class TestAnticommutativity:

    def test_distinct_imaginary_units_anticommute(self, alg: OctonionAlgebra) -> None:
        for i, j in itertools.combinations(range(1, 8), 2):
            ei = alg.basis_element(i)
            ej = alg.basis_element(j)
            assert (ei * ej).is_close((-1.0) * (ej * ei), ATOL), \
                f"e{i} * e{j} != -(e{j} * e{i})"

    def test_real_unit_commutes_with_all(self, alg: OctonionAlgebra) -> None:
        e0 = alg.basis_element(0)
        for i in range(8):
            ei = alg.basis_element(i)
            assert (e0 * ei).is_close(ei * e0, ATOL), \
                f"e0 does not commute with e{i}"


# ---------------------------------------------------------------------------
# 4. Non-commutativity (sanity check)
# Reference: [Wikipedia] Octonion - Properties
# ---------------------------------------------------------------------------

class TestNonCommutativity:

    def test_algebra_is_not_commutative(self, alg: OctonionAlgebra) -> None:
        e1 = alg.basis_element(1)
        e2 = alg.basis_element(2)
        assert not (e1 * e2).is_close(e2 * e1, ATOL), \
            "e1 * e2 == e2 * e1: algebra is unexpectedly commutative"


# ---------------------------------------------------------------------------
# 5. Non-associativity (sanity check)
# Reference: [Wikipedia] Octonion - Properties
# ---------------------------------------------------------------------------

class TestNonAssociativity:

    def test_algebra_is_not_associative(self, alg: OctonionAlgebra) -> None:
        e1 = alg.basis_element(1)
        e2 = alg.basis_element(2)
        e3 = alg.basis_element(3)
        lhs = (e1 * e2) * e3
        rhs = e1 * (e2 * e3)
        assert not lhs.is_close(rhs, ATOL), \
            "(e1*e2)*e3 == e1*(e2*e3): algebra is unexpectedly associative"


# ---------------------------------------------------------------------------
# 6. Composition algebra: N(x*y) = N(x) * N(y)
# This is the defining property of a normed (composition) algebra.
# Reference: [Wikipedia] Octonion - Conjugate and norm
# ---------------------------------------------------------------------------

class TestNormMultiplicativity:

    def test_on_all_basis_pairs(self, alg: OctonionAlgebra) -> None:
        for i in range(8):
            for j in range(8):
                ei = alg.basis_element(i)
                ej = alg.basis_element(j)
                lhs = (ei * ej).norm_sq()
                rhs = ei.norm_sq() * ej.norm_sq()
                assert abs(lhs - rhs) < ATOL, \
                    f"N(e{i} * e{j}) != N(e{i}) * N(e{j})"

    def test_on_general_elements(self, alg: OctonionAlgebra) -> None:
        """Spot-check on a few non-basis elements."""
        samples = [
            alg.element([1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            alg.element([0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            alg.element([0.5] * 8),
            alg.element([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]),
            alg.element([0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
        ]
        for x in samples:
            for y in samples:
                lhs = (x * y).norm_sq()
                rhs = x.norm_sq() * y.norm_sq()
                assert abs(lhs - rhs) < ATOL, \
                    "N(x*y) != N(x)*N(y) for general elements"


# ---------------------------------------------------------------------------
# 7. Conjugate formula: x * conj(x) = conj(x) * x = N(x) * e_0
# Reference: [Wikipedia] Octonion - Conjugate and norm
# ---------------------------------------------------------------------------

class TestConjugate:

    def test_x_times_conjugate_is_norm_times_identity(self, alg: OctonionAlgebra) -> None:
        e0 = alg.basis_element(0)
        for i in range(8):
            ei = alg.basis_element(i)
            expected = ei.norm_sq() * e0
            assert (ei * ei.conjugate()).is_close(expected, ATOL), \
                f"e{i} * conj(e{i}) != N(e{i}) * e0"

    def test_conjugate_times_x_is_norm_times_identity(self, alg: OctonionAlgebra) -> None:
        e0 = alg.basis_element(0)
        for i in range(8):
            ei = alg.basis_element(i)
            expected = ei.norm_sq() * e0
            assert (ei.conjugate() * ei).is_close(expected, ATOL), \
                f"conj(e{i}) * e{i} != N(e{i}) * e0"

    # 8. Conjugate is an involution
    def test_double_conjugate_is_identity(self, alg: OctonionAlgebra) -> None:
        for i in range(8):
            ei = alg.basis_element(i)
            assert ei.conjugate().conjugate().is_close(ei, ATOL), \
                f"conj(conj(e{i})) != e{i}"


# ---------------------------------------------------------------------------
# 9. Left alternativity: (x*x)*y = x*(x*y)  for all x, y.
# Reference: [Wikipedia] Octonion - Alternativity; [Schafer1966] Section 3
# ---------------------------------------------------------------------------

class TestLeftAlternativity:

    def test_on_all_basis_pairs(self, alg: OctonionAlgebra) -> None:
        for i in range(8):
            for j in range(8):
                x = alg.basis_element(i)
                y = alg.basis_element(j)
                lhs = (x * x) * y
                rhs = x * (x * y)
                assert lhs.is_close(rhs, ATOL), \
                    f"(e{i}*e{i})*e{j} != e{i}*(e{i}*e{j})"


# ---------------------------------------------------------------------------
# 10. Right alternativity: x*(y*y) = (x*y)*y  for all x, y.
# Reference: [Wikipedia] Octonion - Alternativity; [Schafer1966] Section 3
# ---------------------------------------------------------------------------

class TestRightAlternativity:

    def test_on_all_basis_pairs(self, alg: OctonionAlgebra) -> None:
        for i in range(8):
            for j in range(8):
                x = alg.basis_element(i)
                y = alg.basis_element(j)
                lhs = x * (y * y)
                rhs = (x * y) * y
                assert lhs.is_close(rhs, ATOL), \
                    f"e{i}*(e{j}*e{j}) != (e{i}*e{j})*e{j}"


# ---------------------------------------------------------------------------
# 11. Flexibility: (x*y)*x = x*(y*x)  for all x, y.
# Reference: [Wikipedia] Octonion - Alternativity
# ---------------------------------------------------------------------------

class TestFlexibility:

    def test_on_all_basis_pairs(self, alg: OctonionAlgebra) -> None:
        for i in range(8):
            for j in range(8):
                x = alg.basis_element(i)
                y = alg.basis_element(j)
                lhs = (x * y) * x
                rhs = x * (y * x)
                assert lhs.is_close(rhs, ATOL), \
                    f"(e{i}*e{j})*e{i} != e{i}*(e{j}*e{i})"


# ---------------------------------------------------------------------------
# 12. Total antisymmetry of the associator on all 8^3 = 512 basis triples.
# The associator (x,y,z) = (x*y)*z - x*(y*z) must satisfy:
#   (x,y,z) = -(y,x,z)   [antisymmetry in first two arguments]
#   (x,y,z) = -(x,z,y)   [antisymmetry in last two arguments]
#   (x,y,z) = (y,z,x)    [cyclic symmetry, follows from the two above]
#
# By trilinearity, this exhaustive check on basis elements certifies the
# property for ALL real octonion elements.
#
# Reference: [Wikipedia] Octonion - Associativity;
#            [Schafer1966] Theorem 3.1 (alternativity <=> associator alternating)
# ---------------------------------------------------------------------------

class TestAssociatorAntisymmetry:

    def test_antisymmetry_swap_first_two(self, alg: OctonionAlgebra) -> None:
        """(x,y,z) = -(y,x,z) for all ordered basis triples."""
        basis = alg.basis
        for x in basis:
            for y in basis:
                for z in basis:
                    a_xyz = associator(x, y, z)
                    a_yxz = associator(y, x, z)
                    assert a_xyz.is_close((-1.0) * a_yxz, ATOL)

    def test_antisymmetry_swap_last_two(self, alg: OctonionAlgebra) -> None:
        """(x,y,z) = -(x,z,y) for all ordered basis triples."""
        basis = alg.basis
        for x in basis:
            for y in basis:
                for z in basis:
                    a_xyz = associator(x, y, z)
                    a_xzy = associator(x, z, y)
                    assert a_xyz.is_close((-1.0) * a_xzy, ATOL)

    def test_cyclic_symmetry(self, alg: OctonionAlgebra) -> None:
        """(x,y,z) = (y,z,x) = (z,x,y) for all ordered basis triples."""
        basis = alg.basis
        for x in basis:
            for y in basis:
                for z in basis:
                    a_xyz = associator(x, y, z)
                    a_yzx = associator(y, z, x)
                    assert a_xyz.is_close(a_yzx, ATOL)


# ---------------------------------------------------------------------------
# 13. Moufang identity (one standard form):
#     z * (x * (z * y)) = ((z * x) * z) * y  for all x, y, z.
# Reference: [Wikipedia] Octonion - Moufang identities
# ---------------------------------------------------------------------------

class TestMoufangIdentity:

    def test_left_moufang_on_all_basis_triples(self, alg: OctonionAlgebra) -> None:
        basis = alg.basis
        for z in basis:
            for x in basis:
                for y in basis:
                    lhs = z * (x * (z * y))
                    rhs = ((z * x) * z) * y
                    assert lhs.is_close(rhs, ATOL), \
                        "Left Moufang identity z*(x*(z*y)) = ((z*x)*z)*y failed"


# ---------------------------------------------------------------------------
# 14. Division algebra: every nonzero element has a two-sided inverse.
# x^{-1} = conj(x) / N(x),  satisfying x * x^{-1} = x^{-1} * x = e_0.
# Reference: [Wikipedia] Octonion - Inverse element
# ---------------------------------------------------------------------------

class TestDivisionAlgebra:

    def test_left_inverse_for_all_basis_elements(self, alg: OctonionAlgebra) -> None:
        e0 = alg.basis_element(0)
        for i in range(8):
            ei = alg.basis_element(i)
            assert (ei.inverse() * ei).is_close(e0, ATOL), \
                f"e{i}^(-1) * e{i} != e0"

    def test_right_inverse_for_all_basis_elements(self, alg: OctonionAlgebra) -> None:
        e0 = alg.basis_element(0)
        for i in range(8):
            ei = alg.basis_element(i)
            assert (ei * ei.inverse()).is_close(e0, ATOL), \
                f"e{i} * e{i}^(-1) != e0"


# ---------------------------------------------------------------------------
# 15. No zero divisors among basis elements.
# Follows from N(x*y) = N(x)*N(y) > 0 for nonzero x, y.
# Reference: [Wikipedia] Octonion - Zero divisors
# ---------------------------------------------------------------------------

class TestNoZeroDivisors:

    def test_no_zero_product_among_nonzero_basis_pairs(self, alg: OctonionAlgebra) -> None:
        for i in range(8):
            for j in range(8):
                ei = alg.basis_element(i)
                ej = alg.basis_element(j)
                product_norm = (ei * ej).norm_sq()
                assert product_norm > ATOL, \
                    f"e{i} * e{j} = 0: zero divisors found"


# ---------------------------------------------------------------------------
# 16. Power-associativity: x^m * x^n = x^(m+n).
# For alternative algebras this follows from left/right alternativity.
# Tested here on basis elements for small powers as a direct check.
# Reference: [Wikipedia] Octonion - Power-associativity
# ---------------------------------------------------------------------------

class TestPowerAssociativity:

    def test_x_cubed_is_consistent(self, alg: OctonionAlgebra) -> None:
        """(x*x)*x = x*(x*x) for all basis elements."""
        for i in range(8):
            x = alg.basis_element(i)
            x2 = x * x
            assert (x2 * x).is_close(x * x2, ATOL), \
                f"e{i}: (x^2)*x != x*(x^2)"

    def test_imaginary_fourth_power_is_identity(self, alg: OctonionAlgebra) -> None:
        """(e_i^2)^2 = (-e_0)^2 = e_0 for i = 1,...,7."""
        e0 = alg.basis_element(0)
        for i in range(1, 8):
            x = alg.basis_element(i)
            x2 = x * x        # = -e_0
            x4 = x2 * x2      # = (-e_0)^2 = e_0
            assert x4.is_close(e0, ATOL), f"e{i}^4 != e0"


# ---------------------------------------------------------------------------
# 17. Explicit verification of the stated multiplication rules from the papers.
#
# Standard algebra (Dixon's rule):  e_a * e_{a+1} = e_{a+3}, a=1,...,7 mod 7.
# Reference: [Dixon2010] eq. (3).
#
# Wilson's rule  i_t * i_{t+1} = i_{t+3}  is the same under  i_k <-> e_{k+1}.
# Reference: [Wilson2009] Section 2.
# ---------------------------------------------------------------------------

class TestMultiplicationRulesMatchPapers:

    def _mod7(self, x: int) -> int:
        """Map integer x to {1,...,7} via  ((x-1) % 7) + 1."""
        return ((x - 1) % 7) + 1

    def test_dixon_standard_rule(self) -> None:
        """e_a * e_{a+1} = e_{a+3} for a = 1,...,7 (Dixon's eq. 3)."""
        alg = STANDARD_ALGEBRA
        for a in range(1, 8):
            b = self._mod7(a + 1)
            c = self._mod7(a + 3)
            ea = alg.basis_element(a)
            eb = alg.basis_element(b)
            ec = alg.basis_element(c)
            assert (ea * eb).is_close(ec, ATOL), \
                f"Standard algebra: e{a} * e{b} != e{c}"

    def test_dixon_xproduct_rule(self) -> None:
        """e_a o e_{a+2} = e_{a+3} for a = 1,...,7 (Dixon's X-product rule)."""
        alg = XPRODUCT_ALGEBRA
        for a in range(1, 8):
            b = self._mod7(a + 2)
            c = self._mod7(a + 3)
            ea = alg.basis_element(a)
            eb = alg.basis_element(b)
            ec = alg.basis_element(c)
            assert (ea * eb).is_close(ec, ATOL), \
                f"X-product algebra: e{a} o e{b} != e{c}"


# ---------------------------------------------------------------------------
# 18. The two algebras have different multiplication tables.
# They are isomorphic as abstract algebras but are distinct bilinear maps
# on the same 8-dimensional vector space.
# ---------------------------------------------------------------------------

class TestAlgebrasAreDistinct:

    def test_e1_times_e3_differs_between_algebras(self) -> None:
        """
        Standard algebra: triple (7,1,3) cyclic => e_1 * e_3 = e_7.
        X-product algebra: triple (1,3,4) first  => e_1 * e_3 = e_4.
        """
        # Standard: e_1 * e_3 = e_7
        std_result = (
            STANDARD_ALGEBRA.basis_element(1) * STANDARD_ALGEBRA.basis_element(3)
        )
        assert std_result.is_close(STANDARD_ALGEBRA.basis_element(7), ATOL), \
            "Standard: e1 * e3 should be e7"
        assert not std_result.is_close(STANDARD_ALGEBRA.basis_element(4), ATOL), \
            "Standard: e1 * e3 should not be e4"

        # X-product: e_1 * e_3 = e_4
        xp_result = (
            XPRODUCT_ALGEBRA.basis_element(1) * XPRODUCT_ALGEBRA.basis_element(3)
        )
        assert xp_result.is_close(XPRODUCT_ALGEBRA.basis_element(4), ATOL), \
            "X-product: e1 * e3 should be e4"
        assert not xp_result.is_close(XPRODUCT_ALGEBRA.basis_element(7), ATOL), \
            "X-product: e1 * e3 should not be e7"


# ---------------------------------------------------------------------------
# 19. Cross-check: X-product table agrees with Dixon's defining formula
#     A o_X B = (A * l0) * (l0^{-1} * B),   l0 = (e0 + e1 + ... + e7) / 2.
#
# This test verifies that the XPRODUCT_FANO_TRIPLES in octonions.py are
# the correct combinatorial encoding of Dixon's analytic definition.
# Both sides are computed using STANDARD_ALGEBRA multiplication.
#
# Reference: [Dixon2010] eq. (7), and the multiplication table stated after it.
# ---------------------------------------------------------------------------

class TestXProductConsistency:

    def test_xproduct_formula_matches_table_on_all_basis_pairs(self) -> None:
        """
        For all basis pairs (e_i, e_j), verify that:
          (e_i * l0) * (l0^{-1} * e_j)  [computed in STANDARD_ALGEBRA]
        equals
          e_i o e_j                       [computed in XPRODUCT_ALGEBRA]
        where l0 = (e0 + e1 + ... + e7) / 2.
        """
        # l0 and its inverse in STANDARD_ALGEBRA
        l0 = STANDARD_ALGEBRA.element([0.5] * 8)
        l0_inv = l0.inverse()

        for i in range(8):
            for j in range(8):
                # X-product via Dixon's formula  (using STANDARD_ALGEBRA)
                A = STANDARD_ALGEBRA.basis_element(i)
                B = STANDARD_ALGEBRA.basis_element(j)
                xp_via_formula = (A * l0) * (l0_inv * B)

                # X-product via XPRODUCT_ALGEBRA table
                A_x = XPRODUCT_ALGEBRA.basis_element(i)
                B_x = XPRODUCT_ALGEBRA.basis_element(j)
                xp_via_table = A_x * B_x

                assert np.allclose(
                    xp_via_formula.coords, xp_via_table.coords, atol=ATOL
                ), (
                    f"X-product mismatch for e{i} o e{j}: "
                    f"formula gives {xp_via_formula}, "
                    f"table gives {xp_via_table}"
                )
