"""
okubo.py — Okubo algebra (pseudo-octonion algebra) implementation.

The Okubo algebra is an 8-dimensional real composition algebra that is
non-unital, non-alternative, and flexible.  It is the unique (up to
isomorphism) 8-dimensional non-unital symmetric composition division
algebra over R.

Two construction methods are provided:

  from_octonions(oct_algebra)
      Petersson isotope construction.  Given an octonion algebra (O, ·, n)
      and an order-3 automorphism τ, the Okubo product is
          x * y = τ(x̄) · τ²(ȳ).
      Reference: [MarraniCorradettiZucconi2025] eq. (1.6)/(1.17);
                 [Elduque2000_Triality].

  from_hermitian_matrices()
      Direct construction on traceless 3×3 Hermitian matrices:
          x * y = μ · xy + μ̄ · yx − (1/3) Tr(xy) · I₃
      where μ = (3 + i√3)/6.
      Reference: [MarraniCorradettiZucconi2025] eq. (1.8).

Key properties (verified in test_okubo.py):

  - Composition algebra:  n(x*y) = n(x) n(y)  (Euclidean norm)
  - Symmetric composition:  (x*y)*x = x*(y*x) = n(x) y
  - Non-unital:  no identity element
  - Non-alternative:  neither left- nor right-alternative
  - Flexible:  (x*y)*x = x*(y*x)
  - Division algebra:  L_a and R_a are bijections for a ≠ 0
  - Automorphism group Aut(O) = SU(3),  dim Der(O) = 8

References
----------
[MarraniCorradettiZucconi2025]
    A. Marrani, D. Corradetti, F. Zucconi, "Physics with non-unital
    algebras? An invitation to the Okubo algebra", J. Phys. A 58 (2025)
    075202.  DOI: 10.1088/1751-8121/adafef

[Elduque2000_Triality]
    A. Elduque, "On triality and automorphisms and derivations of
    composition algebras", Linear Algebra Appl. 314 (2000) 49–74.
    DOI: 10.1016/S0024-3795(00)00105-1
"""

import numpy as np
from typing import List, Optional, Sequence

from octonions import OctonionAlgebra, STANDARD_ALGEBRA


# ---------------------------------------------------------------------------
# OkuboAlgebra
# ---------------------------------------------------------------------------

class OkuboAlgebra:
    """
    An 8-dimensional real algebra encoded by its structure constants.

    The algebra lives on R^8 with basis {e_0, …, e_7}.  The product
    e_i * e_j = Σ_k C[i,j,k] e_k is determined by the structure constants
    tensor C of shape (8, 8, 8).

    The composition norm is the standard Euclidean norm:
    n(x) = x_0² + x_1² + … + x_7².
    """

    def __init__(self, struct_constants: np.ndarray, name: str = "") -> None:
        if struct_constants.shape != (8, 8, 8):
            raise ValueError(
                f"Structure constants must have shape (8, 8, 8), "
                f"got {struct_constants.shape}"
            )
        self._struct = struct_constants.copy()
        self.name = name

    def element(self, coords: Sequence[float]) -> "OkuboElement":
        """Create an element from a length-8 coordinate vector."""
        return OkuboElement(np.asarray(coords, dtype=float), self)

    def basis_element(self, i: int) -> "OkuboElement":
        """Return basis element e_i  (0 ≤ i ≤ 7)."""
        c = np.zeros(8, dtype=float)
        c[i] = 1.0
        return OkuboElement(c, self)

    @property
    def basis(self) -> List["OkuboElement"]:
        """All 8 basis elements [e_0, …, e_7]."""
        return [self.basis_element(i) for i in range(8)]

    def zero(self) -> "OkuboElement":
        """The zero element."""
        return OkuboElement(np.zeros(8, dtype=float), self)

    def _mul_coords(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """Coordinate-level multiplication (used by OkuboElement.__mul__)."""
        return np.einsum('ijk,i,j->k', self._struct, a, b)

    def __repr__(self) -> str:
        return f"OkuboAlgebra(name={self.name!r})"


# ---------------------------------------------------------------------------
# OkuboElement
# ---------------------------------------------------------------------------

class OkuboElement:
    """An element of an OkuboAlgebra."""

    __slots__ = ("coords", "algebra")

    def __init__(self, coords: np.ndarray, algebra: OkuboAlgebra) -> None:
        self.coords: np.ndarray = coords
        self.algebra: OkuboAlgebra = algebra

    def __mul__(self, other: "OkuboElement") -> "OkuboElement":
        return OkuboElement(
            self.algebra._mul_coords(self.coords, other.coords),
            self.algebra,
        )

    def __add__(self, other: "OkuboElement") -> "OkuboElement":
        return OkuboElement(self.coords + other.coords, self.algebra)

    def __sub__(self, other: "OkuboElement") -> "OkuboElement":
        return OkuboElement(self.coords - other.coords, self.algebra)

    def __neg__(self) -> "OkuboElement":
        return OkuboElement(-self.coords, self.algebra)

    def __rmul__(self, scalar: float) -> "OkuboElement":
        return OkuboElement(float(scalar) * self.coords, self.algebra)

    def norm_sq(self) -> float:
        """Squared Euclidean norm  n(x) = x_0² + … + x_7²."""
        return float(np.dot(self.coords, self.coords))

    def is_close(self, other: "OkuboElement", atol: float = 1e-10) -> bool:
        """Return True if all coordinates agree to within atol."""
        return bool(np.allclose(self.coords, other.coords, atol=atol))

    def __repr__(self) -> str:
        parts = [
            f"{c:+g}*e{i}" for i, c in enumerate(self.coords) if abs(c) > 1e-12
        ]
        return "Okubo(" + ("".join(parts) or "0") + ")"


# ---------------------------------------------------------------------------
# Construction 1:  Petersson isotope from octonions
# ---------------------------------------------------------------------------

def _standard_tau_matrix() -> np.ndarray:
    """
    8×8 matrix of the order-3 automorphism τ of the standard octonion
    algebra.

    τ fixes e_0, e_1, e_3, e_7 and applies a 2π/3 rotation to each of
    the planes (e_2, e_5) and (e_4, e_6).

    Reference: [MarraniCorradettiZucconi2025] eq. (1.5).
    """
    T = np.eye(8)
    c, s = -0.5, np.sqrt(3.0) / 2.0   # cos 2π/3, sin 2π/3
    # (e_2, e_5) block
    T[2, 2] = c;  T[2, 5] = s
    T[5, 2] = -s; T[5, 5] = c
    # (e_4, e_6) block
    T[4, 4] = c;  T[4, 6] = s
    T[6, 4] = -s; T[6, 6] = c
    return T


def _validate_order3(T: np.ndarray, atol: float = 1e-12) -> None:
    """Verify T³ = I."""
    T3 = T @ T @ T
    if not np.allclose(T3, np.eye(8), atol=atol):
        raise ValueError("τ does not have order 3 (τ³ ≠ I).")


def _validate_automorphism(oct_algebra: OctonionAlgebra, T: np.ndarray,
                           atol: float = 1e-12) -> None:
    """
    Verify that T is an algebra automorphism of oct_algebra:
    T(e_i · e_j) = T(e_i) · T(e_j) for all basis pairs.
    """
    for i in range(8):
        ei = np.zeros(8); ei[i] = 1.0
        Tei = T @ ei
        for j in range(8):
            ej = np.zeros(8); ej[j] = 1.0
            Tej = T @ ej
            lhs = T @ oct_algebra._mul_coords(ei, ej)
            rhs = oct_algebra._mul_coords(Tei, Tej)
            if not np.allclose(lhs, rhs, atol=atol):
                raise ValueError(
                    f"T is not an automorphism: "
                    f"T(e_{i}·e_{j}) ≠ T(e_{i})·T(e_{j})"
                )


def from_octonions(
    oct_algebra: OctonionAlgebra = STANDARD_ALGEBRA,
    tau_matrix: Optional[np.ndarray] = None,
) -> OkuboAlgebra:
    """
    Construct the Okubo algebra as a Petersson isotope of an octonion
    algebra.

    The Okubo product is:
        x * y = τ(x̄) · τ²(ȳ)
    where · is the octonion product, ¯ is octonion conjugation, and
    τ is an order-3 automorphism.

    The Euclidean norm on R⁸ is the composition norm (inherited from the
    octonion algebra, since τ and conjugation are isometries).

    Parameters
    ----------
    oct_algebra : OctonionAlgebra
        The octonion algebra to isotope (default: STANDARD_ALGEBRA).
    tau_matrix : ndarray of shape (8, 8), optional
        Matrix of an order-3 automorphism.  Defaults to the automorphism
        from [MarraniCorradettiZucconi2025] eq. (1.5), valid for
        STANDARD_ALGEBRA.

    Returns
    -------
    OkuboAlgebra

    References
    ----------
    [MarraniCorradettiZucconi2025] eq. (1.6), (1.17).
    """
    if tau_matrix is None:
        tau_matrix = _standard_tau_matrix()

    _validate_order3(tau_matrix)
    _validate_automorphism(oct_algebra, tau_matrix)

    T = tau_matrix
    T2 = T @ T
    conj = np.diag([1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0])

    struct = np.zeros((8, 8, 8))
    for i in range(8):
        ei = np.zeros(8); ei[i] = 1.0
        tau_conj_ei = T @ (conj @ ei)       # τ(ē_i)
        for j in range(8):
            ej = np.zeros(8); ej[j] = 1.0
            tau2_conj_ej = T2 @ (conj @ ej)  # τ²(ē_j)
            struct[i, j, :] = oct_algebra._mul_coords(
                tau_conj_ei, tau2_conj_ej
            )

    return OkuboAlgebra(struct, name=f"petersson({oct_algebra.name})")


# ---------------------------------------------------------------------------
# Construction 2:  Hermitian matrix construction
# ---------------------------------------------------------------------------

def _gell_mann_basis() -> List[np.ndarray]:
    """
    Orthonormal basis of the 8-dimensional space of traceless 3×3
    Hermitian matrices, under the inner product ⟨x, y⟩ = (1/6) Tr(xy).

    Uses √3 × the standard Gell-Mann matrices λ₁, …, λ₈.
    Index k corresponds to √3 · λ_{k+1}  (k = 0, …, 7).

    Orthonormality:  (1/6) Tr(f_j f_k) = δ_{jk}.
    """
    s3 = np.sqrt(3.0)
    return [
        # f_0 = √3 λ₁
        s3 * np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
        # f_1 = √3 λ₂
        s3 * np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex),
        # f_2 = √3 λ₃
        s3 * np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex),
        # f_3 = √3 λ₄
        s3 * np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
        # f_4 = √3 λ₅
        s3 * np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex),
        # f_5 = √3 λ₆
        s3 * np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
        # f_6 = √3 λ₇
        s3 * np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex),
        # f_7 = √3 λ₈ = diag(1, 1, −2)
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex),
    ]


def _express_in_basis(matrix: np.ndarray,
                      basis: List[np.ndarray]) -> np.ndarray:
    """
    Express a traceless Hermitian 3×3 matrix in the given orthonormal basis.

    Projection formula:  c_k = (1/6) Re Tr(matrix · f_k).
    """
    coords = np.zeros(len(basis))
    for k, fk in enumerate(basis):
        coords[k] = np.real(np.trace(matrix @ fk)) / 6.0
    return coords


def from_hermitian_matrices() -> OkuboAlgebra:
    """
    Construct the Okubo algebra on the space of traceless 3×3 Hermitian
    matrices, with the product

        x * y = μ · xy + μ̄ · yx − (1/3) Tr(xy) · I₃

    where μ = (3 + i√3)/6 and xy is ordinary matrix multiplication.

    An orthonormal scaled Gell-Mann basis is used so that the composition
    norm n(x) = (1/6) Tr(x²) equals the Euclidean norm on coordinates.

    Returns
    -------
    OkuboAlgebra

    References
    ----------
    [MarraniCorradettiZucconi2025] eq. (1.8), (1.10)–(1.11).
    """
    basis = _gell_mann_basis()
    mu = (3.0 + 1j * np.sqrt(3.0)) / 6.0
    mu_bar = mu.conjugate()
    I3 = np.eye(3, dtype=complex)

    struct = np.zeros((8, 8, 8))
    for i in range(8):
        for j in range(8):
            xy = basis[i] @ basis[j]
            yx = basis[j] @ basis[i]
            tr_xy = np.trace(xy)
            product = mu * xy + mu_bar * yx - (tr_xy / 3.0) * I3
            struct[i, j, :] = _express_in_basis(product, basis)

    return OkuboAlgebra(struct, name="hermitian_matrices")


# ---------------------------------------------------------------------------
# Module-level algebra instances
# ---------------------------------------------------------------------------

#: Okubo algebra via Petersson isotope of the standard octonion algebra.
#: Reference: [MarraniCorradettiZucconi2025] eq. (1.5)–(1.6).
PETERSSON_ALGEBRA = from_octonions()

#: Okubo algebra via traceless 3×3 Hermitian matrices.
#: Reference: [MarraniCorradettiZucconi2025] eq. (1.8).
HERMITIAN_ALGEBRA = from_hermitian_matrices()
