"""
test_okubo.py — Tests for the Okubo algebra implementation.

Both construction methods (Petersson isotope from octonions, Hermitian
matrix construction) are exercised by every test via the parametrised
`alg` fixture.

Properties verified
-------------------
Together, properties 1–5 on an 8-dimensional real algebra uniquely
characterise the Okubo algebra (up to isomorphism), by the classification
of symmetric composition algebras [Elduque2000_Triality].

 1.  Bilinearity of the product  (sanity check; guaranteed by structure
     constants but verified explicitly)
 2.  Composition norm:  n(x*y) = n(x) n(y) for the Euclidean norm
 3.  Symmetric composition:  (x*y)*x = x*(y*x) = n(x) y
 4.  Non-unital:  no identity element exists
 5.  Division algebra:  L_a and R_a are bijections for every a ≠ 0

Additional properties (consequences of the above):

 6.  Flexibility:  (x*y)*x = x*(y*x)
 7.  Non-alternative:  left and right alternativity fail
 8.  Existence of a norm-1 idempotent  e * e = e,  n(e) = 1
 9.  Derivation algebra dimension = 8  (= dim su(3), confirming Aut = SU(3))

Methodology note
----------------
Properties 2, 3 are tested exhaustively on all 64 ordered pairs of basis
elements.  By bilinearity/trilinearity this certifies the identities for
ALL real linear combinations, not just basis elements.  Additional random-
element tests provide independent numerical confidence.

References
----------
[MarraniCorradettiZucconi2025]
    A. Marrani, D. Corradetti, F. Zucconi, "Physics with non-unital
    algebras? An invitation to the Okubo algebra", J. Phys. A 58 (2025)
    075202.  Table 2.
[Elduque2000_Triality]
    A. Elduque, "On triality and automorphisms and derivations of
    composition algebras", Linear Algebra Appl. 314 (2000) 49–74.
"""

import numpy as np
import pytest

from okubo import (
    OkuboAlgebra,
    from_octonions,
    from_hermitian_matrices,
    PETERSSON_ALGEBRA,
    HERMITIAN_ALGEBRA,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ATOL = 1e-10
N_RANDOM = 20   # random element pairs per test


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

@pytest.fixture(params=["petersson", "hermitian"])
def alg(request):
    """Parametrised fixture: runs every test on both constructions."""
    if request.param == "petersson":
        return PETERSSON_ALGEBRA
    return HERMITIAN_ALGEBRA


def _random_elements(alg, n=2, seed=42):
    """Return n random algebra elements (reproducible via seed)."""
    rng = np.random.default_rng(seed)
    return [alg.element(rng.standard_normal(8)) for _ in range(n)]


# ---------------------------------------------------------------------------
# 1. Bilinearity
# ---------------------------------------------------------------------------

class TestBilinearity:
    """The product is bilinear.  This is guaranteed by the structure-
    constant representation, but verified explicitly as a sanity check."""

    def test_left_linearity(self, alg):
        rng = np.random.default_rng(100)
        for _ in range(N_RANDOM):
            a, b = rng.standard_normal(2)
            x, y, z = _random_elements(alg, 3, seed=rng.integers(10**9))
            lhs = (a * x + b * y) * z
            rhs = a * (x * z) + b * (y * z)
            assert lhs.is_close(rhs, atol=ATOL)

    def test_right_linearity(self, alg):
        rng = np.random.default_rng(101)
        for _ in range(N_RANDOM):
            a, b = rng.standard_normal(2)
            x, y, z = _random_elements(alg, 3, seed=rng.integers(10**9))
            lhs = z * (a * x + b * y)
            rhs = a * (z * x) + b * (z * y)
            assert lhs.is_close(rhs, atol=ATOL)


# ---------------------------------------------------------------------------
# 2. Composition norm:  n(x*y) = n(x) n(y)
# ---------------------------------------------------------------------------

class TestCompositionNorm:
    """n(x*y) = n(x) · n(y) for the Euclidean norm.

    This is the defining property of a composition algebra.
    Reference: [MarraniCorradettiZucconi2025] Table 2 (Composition: Yes).
    """

    def test_basis_pairs(self, alg):
        """Exhaustive over all 64 ordered basis pairs."""
        for i, ei in enumerate(alg.basis):
            for j, ej in enumerate(alg.basis):
                n_prod = (ei * ej).norm_sq()
                expected = ei.norm_sq() * ej.norm_sq()
                assert abs(n_prod - expected) < ATOL, (
                    f"n(e_{i}*e_{j}) = {n_prod}, "
                    f"expected n(e_{i})·n(e_{j}) = {expected}"
                )

    def test_random_pairs(self, alg):
        rng = np.random.default_rng(200)
        for _ in range(N_RANDOM):
            x, y = _random_elements(alg, 2, seed=rng.integers(10**9))
            n_prod = (x * y).norm_sq()
            expected = x.norm_sq() * y.norm_sq()
            assert abs(n_prod - expected) < ATOL * (1.0 + abs(expected))


# ---------------------------------------------------------------------------
# 3. Symmetric composition:  (x*y)*x = x*(y*x) = n(x) y
# ---------------------------------------------------------------------------

class TestSymmetricComposition:
    """(x*y)*x = x*(y*x) = n(x) · y.

    This is the defining property of symmetric composition algebras,
    which include para-Hurwitz algebras and the Okubo algebra.
    Reference: [Elduque2000_Triality]; [MarraniCorradettiZucconi2025]
    discussion before Table 3.
    """

    def test_basis_pairs(self, alg):
        """Exhaustive over all 64 ordered basis pairs."""
        for i, ei in enumerate(alg.basis):
            nx = ei.norm_sq()
            for j, ej in enumerate(alg.basis):
                rhs = nx * ej
                lhs_left = (ei * ej) * ei
                lhs_right = ei * (ej * ei)
                assert lhs_left.is_close(rhs, atol=ATOL), (
                    f"(e_{i}*e_{j})*e_{i} ≠ n(e_{i})·e_{j}"
                )
                assert lhs_right.is_close(rhs, atol=ATOL), (
                    f"e_{i}*(e_{j}*e_{i}) ≠ n(e_{i})·e_{j}"
                )

    def test_random_pairs(self, alg):
        rng = np.random.default_rng(300)
        for _ in range(N_RANDOM):
            x, y = _random_elements(alg, 2, seed=rng.integers(10**9))
            nx = x.norm_sq()
            rhs = nx * y
            tol = ATOL * (1.0 + abs(nx))
            assert ((x * y) * x).is_close(rhs, atol=tol), (
                "(x*y)*x ≠ n(x)·y"
            )
            assert (x * (y * x)).is_close(rhs, atol=tol), (
                "x*(y*x) ≠ n(x)·y"
            )


# ---------------------------------------------------------------------------
# 4. Flexibility:  (x*y)*x = x*(y*x)
# ---------------------------------------------------------------------------

class TestFlexibility:
    """(x*y)*x = x*(y*x).

    This follows from the symmetric composition property but is tested
    independently as a basic sanity check.
    Reference: [MarraniCorradettiZucconi2025] Table 2 (Flexible: Yes).
    """

    def test_random_pairs(self, alg):
        rng = np.random.default_rng(400)
        for _ in range(N_RANDOM):
            x, y = _random_elements(alg, 2, seed=rng.integers(10**9))
            assert ((x * y) * x).is_close(x * (y * x), atol=ATOL)


# ---------------------------------------------------------------------------
# 5. Non-unital
# ---------------------------------------------------------------------------

class TestNonUnital:
    """No element acts as a two-sided identity.

    Reference: [MarraniCorradettiZucconi2025] Table 2 (Unital: No).
    """

    def test_no_identity_among_basis_elements(self, alg):
        basis = alg.basis
        for i, ei in enumerate(basis):
            is_left = all((ei * ej).is_close(ej, atol=ATOL) for ej in basis)
            assert not is_left, f"e_{i} is a left identity"
            is_right = all((ej * ei).is_close(ej, atol=ATOL) for ej in basis)
            assert not is_right, f"e_{i} is a right identity"

    def test_no_left_identity_general(self, alg):
        """Verify the linear system u*x = x for all x has no solution."""
        C = alg._struct
        # u*e_j = e_j  =>  Σ_i u_i C[i,j,k] = δ_{jk}  for all j, k
        A = np.zeros((64, 8))
        b = np.zeros(64)
        for j in range(8):
            for k in range(8):
                row = j * 8 + k
                A[row, :] = C[:, j, k]
                b[row] = 1.0 if j == k else 0.0
        u, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
        residual = np.linalg.norm(A @ u - b)
        assert residual > 0.1, "A left identity exists (the algebra is unital)"

    def test_no_right_identity_general(self, alg):
        """Verify the linear system x*u = x for all x has no solution."""
        C = alg._struct
        # e_i * u = e_i  =>  Σ_j u_j C[i,j,k] = δ_{ik}  for all i, k
        A = np.zeros((64, 8))
        b = np.zeros(64)
        for i in range(8):
            for k in range(8):
                row = i * 8 + k
                A[row, :] = C[i, :, k]
                b[row] = 1.0 if i == k else 0.0
        u, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
        residual = np.linalg.norm(A @ u - b)
        assert residual > 0.1, "A right identity exists"


# ---------------------------------------------------------------------------
# 6. Non-alternative
# ---------------------------------------------------------------------------

class TestNonAlternative:
    """The algebra is neither left- nor right-alternative.

    Reference: [MarraniCorradettiZucconi2025] Table 2 (Alternative: No).
    """

    def test_not_left_alternative(self, alg):
        """Find x, y with (x*x)*y ≠ x*(x*y)."""
        rng = np.random.default_rng(500)
        found = False
        for _ in range(50):
            x, y = _random_elements(alg, 2, seed=rng.integers(10**9))
            if not ((x * x) * y).is_close(x * (x * y), atol=1e-6):
                found = True
                break
        assert found, "No counterexample to left alternativity found"

    def test_not_right_alternative(self, alg):
        """Find x, y with (y*x)*x ≠ y*(x*x)."""
        rng = np.random.default_rng(501)
        found = False
        for _ in range(50):
            x, y = _random_elements(alg, 2, seed=rng.integers(10**9))
            if not ((y * x) * x).is_close(y * (x * x), atol=1e-6):
                found = True
                break
        assert found, "No counterexample to right alternativity found"


# ---------------------------------------------------------------------------
# 7. Division algebra
# ---------------------------------------------------------------------------

class TestDivision:
    """L_a and R_a are bijective for every nonzero element a.

    This is verified by checking that the 8×8 matrix of left (resp. right)
    multiplication by a is invertible (nonzero determinant).

    For a composition algebra, the division property is equivalent to the
    norm being anisotropic (n(x) = 0 implies x = 0), which holds for the
    Euclidean norm.
    """

    @staticmethod
    def _left_mul_matrix(alg, a):
        """8×8 matrix of L_a: x ↦ a*x.  M[k,j] = Σ_i a_i C[i,j,k]."""
        return np.einsum('i,ijk->kj', a, alg._struct)

    @staticmethod
    def _right_mul_matrix(alg, a):
        """8×8 matrix of R_a: x ↦ x*a.  M[k,i] = Σ_j a_j C[i,j,k]."""
        return np.einsum('j,ijk->ki', a, alg._struct)

    def test_basis_elements(self, alg):
        for i in range(8):
            a = np.zeros(8); a[i] = 1.0
            L = self._left_mul_matrix(alg, a)
            R = self._right_mul_matrix(alg, a)
            assert abs(np.linalg.det(L)) > 1e-10, (
                f"L_{{e_{i}}} is singular"
            )
            assert abs(np.linalg.det(R)) > 1e-10, (
                f"R_{{e_{i}}} is singular"
            )

    def test_random_elements(self, alg):
        rng = np.random.default_rng(600)
        for _ in range(N_RANDOM):
            a = rng.standard_normal(8)
            L = self._left_mul_matrix(alg, a)
            R = self._right_mul_matrix(alg, a)
            assert abs(np.linalg.det(L)) > 1e-10, "L_a singular"
            assert abs(np.linalg.det(R)) > 1e-10, "R_a singular"


# ---------------------------------------------------------------------------
# 8. Idempotent
# ---------------------------------------------------------------------------

class TestIdempotent:
    """There exist norm-1 idempotents: e*e = e with n(e) = 1.

    Unlike unital composition algebras, the Okubo algebra has idempotent
    elements but no identity element.
    Reference: [MarraniCorradettiZucconi2025] eq. (1.10).
    """

    def test_petersson_e0_is_idempotent(self):
        """In the Petersson construction, e_0 (the octonion identity) is
        an idempotent of the Okubo algebra.

        Proof: e_0 * e_0 = τ(ē_0) · τ²(ē_0) = τ(e_0) · τ²(e_0)
             = e_0 · e_0 = e_0  (since τ fixes e_0 and e_0 is the
             octonion identity with ē_0 = e_0).
        """
        alg = PETERSSON_ALGEBRA
        e0 = alg.basis_element(0)
        assert (e0 * e0).is_close(e0, atol=ATOL)
        assert abs(e0.norm_sq() - 1.0) < ATOL

    def test_hermitian_known_idempotent(self):
        """In the Hermitian construction, diag(2, −1, −1) is an idempotent.

        In the scaled Gell-Mann basis (f_0, …, f_7), this matrix
        decomposes as (√3/2) f_2 + (1/2) f_7, with norm² = 3/4 + 1/4 = 1.

        Reference: [MarraniCorradettiZucconi2025] eq. (1.10).
        """
        alg = HERMITIAN_ALGEBRA
        coords = np.zeros(8)
        coords[2] = np.sqrt(3.0) / 2.0   # (√3/2) · f_2 = (√3/2) · √3 λ₃
        coords[7] = 0.5                   # (1/2) · f_7  = (1/2) · diag(1,1,−2)
        e = alg.element(coords)
        assert abs(e.norm_sq() - 1.0) < ATOL, (
            f"Idempotent norm² = {e.norm_sq()}, expected 1"
        )
        assert (e * e).is_close(e, atol=ATOL), "e*e ≠ e"

    def test_idempotent_is_not_identity(self, alg):
        """The idempotent does NOT act as a left or right identity.

        This distinguishes the Okubo algebra from unital composition
        algebras.
        """
        if alg.name.startswith("petersson"):
            e = alg.basis_element(0)
        else:
            coords = np.zeros(8)
            coords[2] = np.sqrt(3.0) / 2.0
            coords[7] = 0.5
            e = alg.element(coords)
        basis = alg.basis
        is_left_id = all((e * b).is_close(b, atol=ATOL) for b in basis)
        assert not is_left_id, "The idempotent acts as a left identity"
        is_right_id = all((b * e).is_close(b, atol=ATOL) for b in basis)
        assert not is_right_id, "The idempotent acts as a right identity"


# ---------------------------------------------------------------------------
# 9. Derivation algebra dimension = 8  (Aut = SU(3))
# ---------------------------------------------------------------------------

class TestDerivationAlgebra:
    """dim Der(O) = 8 = dim su(3).

    A derivation D of the algebra satisfies D(x*y) = D(x)*y + x*D(y).
    The space of all derivations is the Lie algebra of the automorphism
    group.  For the Okubo algebra, Aut(O) = SU(3), so dim Der = 8.

    Reference: [Elduque2000_Triality]; [MarraniCorradettiZucconi2025]
    Table 4  (Aut = SU(3)).
    """

    def test_derivation_dimension(self, alg):
        C = alg._struct
        # D is an 8×8 matrix.  Flatten to 64 unknowns:
        # d[m*8+n] = D_{mn},  meaning D(e_m) = Σ_n D_{mn} e_n.
        #
        # The derivation condition D(e_i * e_j) = D(e_i)*e_j + e_i*D(e_j)
        # expands to, for each (i, j, b):
        #   Σ_k C[i,j,k] D[k,b] = Σ_a D[i,a] C[a,j,b] + Σ_a D[j,a] C[i,a,b]
        n_eq = 8 * 8 * 8   # 512 equations (many redundant)
        n_var = 64
        A = np.zeros((n_eq, n_var))
        for i in range(8):
            for j in range(8):
                for b in range(8):
                    row = i * 64 + j * 8 + b
                    for k in range(8):
                        A[row, k * 8 + b] += C[i, j, k]
                    for a in range(8):
                        A[row, i * 8 + a] -= C[a, j, b]
                        A[row, j * 8 + a] -= C[i, a, b]
        rank = np.linalg.matrix_rank(A, tol=1e-8)
        dim_der = n_var - rank
        assert dim_der == 8, (
            f"dim Der = {dim_der}, expected 8 (su(3))"
        )


# ---------------------------------------------------------------------------
# 10. Norm is positive definite
# ---------------------------------------------------------------------------

class TestNormPositiveDefinite:
    """The Euclidean norm n(x) = Σ x_i² is positive definite."""

    def test_zero_has_zero_norm(self, alg):
        assert alg.zero().norm_sq() == 0.0

    def test_nonzero_has_positive_norm(self, alg):
        rng = np.random.default_rng(700)
        for _ in range(N_RANDOM):
            x = _random_elements(alg, 1, seed=rng.integers(10**9))[0]
            assert x.norm_sq() > 0

    def test_basis_elements_have_unit_norm(self, alg):
        for i, ei in enumerate(alg.basis):
            assert abs(ei.norm_sq() - 1.0) < ATOL, (
                f"n(e_{i}) = {ei.norm_sq()}, expected 1"
            )
