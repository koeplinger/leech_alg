"""
octonions.py — Octonion algebra implementation.

Two concrete algebras are provided:

  STANDARD_ALGEBRA
      The product used (up to index relabelling) by both Wilson [Wilson2009]
      and Dixon [Dixon2010] in their respective E8/Leech lattice constructions.
      Fano-plane rule:  e_a * e_{a+1} = e_{a+3},  a = 1,...,7 (mod 7 in {1,...,7}).

  XPRODUCT_ALGEBRA
      Dixon's X-product [Dixon2010], defined as
          A o_X B = (A * l0) * (l0^{-1} * B),   l0 = (e0+e1+...+e7)/2.
      This yields a different but isomorphic octonion algebra on the same vector
      space.  Fano-plane rule:  e_a o e_{a+2} = e_{a+3},  a = 1,...,7 (mod 7).

Index convention (used throughout this file)
--------------------------------------------
Basis elements are labelled e_0, ..., e_7 where e_0 is the real unit (identity).
This matches Dixon's notation directly.

Wilson [Wilson2009] uses the index set PL(7) = {inf} u F_7, with imaginary
units labelled i_0, ..., i_6 (and i_inf for the identity).  The correspondence
to Dixon's labelling is:
    Wilson  i_inf  <->  Dixon  e_0   (identity)
    Wilson  i_k    <->  Dixon  e_{k+1}   for k = 0,...,6.
Under this map Wilson's base rule  i_0 * i_1 = i_3  becomes  e_1 * e_2 = e_4,
which is Dixon's rule at a = 1.  The two papers therefore use the *same*
Fano-plane multiplication table.  See key claim 001 in evidence_and_reasoning/.

Fano-plane triples: convention
-------------------------------
A triple (a, b, c) with a, b, c in {1,...,7} encodes the rule
    e_a * e_b = e_c
together with its two cyclic rotations (positive sign) and three anti-cyclic
rotations (negative sign).  All 42 imaginary-imaginary products are generated
from the 7 stated triples by this rule.

References
----------
[Dixon2010]   G.M. Dixon, "Integral Octonions, Octonion XY-Product, and the
              Leech Lattice", preprint, 2010.
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186-2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
[Wikipedia]   https://en.wikipedia.org/wiki/Octonion
"""

import numpy as np
from typing import List, Sequence, Tuple

# ---------------------------------------------------------------------------
# Fano-plane triple sets
# ---------------------------------------------------------------------------

# Standard product.
# Derived from Dixon's rule: e_a * e_{a+1} = e_{a+3}, a = 1,...,7 (mod 7 in {1,...,7}).
# Equivalently from Wilson's rule: i_t * i_{t+1} = i_{t+3} (mod 7 in {0,...,6}),
# after substituting i_k -> e_{k+1}.
# Reference: [Dixon2010] eq. (3); [Wilson2009] Section 2.
STANDARD_FANO_TRIPLES: Tuple[Tuple[int, int, int], ...] = (
    (1, 2, 4),
    (2, 3, 5),
    (3, 4, 6),
    (4, 5, 7),
    (5, 6, 1),
    (6, 7, 2),
    (7, 1, 3),
)

# X-product (Dixon).
# Derived from Dixon's rule: e_a o_{l0} e_{a+2} = e_{a+3}, a = 1,...,7 (mod 7 in {1,...,7}).
# Reference: [Dixon2010] Section 3, multiplication table stated after eq. (7).
XPRODUCT_FANO_TRIPLES: Tuple[Tuple[int, int, int], ...] = (
    (1, 3, 4),
    (2, 4, 5),
    (3, 5, 6),
    (4, 6, 7),
    (5, 7, 1),
    (6, 1, 2),
    (7, 2, 3),
)


# ---------------------------------------------------------------------------
# OctonionAlgebra
# ---------------------------------------------------------------------------

class OctonionAlgebra:
    """
    An octonion algebra on R^8 determined by 7 Fano-plane triples.

    Basis: e_0 (identity), e_1, ..., e_7 (imaginary units).

    Full multiplication rules:
      - e_0 * e_i = e_i * e_0 = e_i  for all i
      - e_i * e_i = -e_0             for i = 1,...,7
      - For each triple (a, b, c):
            e_a * e_b = +e_c   (and its two cyclic rotations, positive)
            e_b * e_a = -e_c   (and its two anti-cyclic rotations, negative)
    """

    def __init__(
        self,
        fano_triples: Sequence[Tuple[int, int, int]],
        name: str = "",
    ) -> None:
        self.name = name
        self.fano_triples: Tuple[Tuple[int, int, int], ...] = tuple(fano_triples)
        # table[i][j] = (sign, k)  meaning  e_i * e_j = sign * e_k
        self._table: List[List[Tuple[int, int]]] = self._build_table()

    # ------------------------------------------------------------------ build

    def _build_table(self) -> List[List[Tuple[int, int]]]:
        table: List[List[object]] = [[None] * 8 for _ in range(8)]

        # e_0 is the two-sided identity
        for i in range(8):
            table[0][i] = (1, i)
            table[i][0] = (1, i)

        # e_i^2 = -e_0  for i = 1,...,7
        for i in range(1, 8):
            table[i][i] = (-1, 0)

        # Fano-plane triples: cyclic (positive) and anti-cyclic (negative)
        for a, b, c in self.fano_triples:
            table[a][b] = ( 1, c)   # e_a * e_b = +e_c
            table[b][c] = ( 1, a)   # e_b * e_c = +e_a  (cyclic rotation)
            table[c][a] = ( 1, b)   # e_c * e_a = +e_b  (cyclic rotation)
            table[b][a] = (-1, c)   # e_b * e_a = -e_c  (anti-cyclic)
            table[c][b] = (-1, a)   # e_c * e_b = -e_a  (anti-cyclic)
            table[a][c] = (-1, b)   # e_a * e_c = -e_b  (anti-cyclic)

        # Validate completeness: every entry must be assigned
        for i in range(8):
            for j in range(8):
                if table[i][j] is None:
                    raise ValueError(
                        f"Multiplication table entry ({i}, {j}) was not set. "
                        f"The 7 Fano-plane triples must cover all 21 unordered "
                        f"pairs from {{1,...,7}}."
                    )

        return table  # type: ignore[return-value]

    # --------------------------------------------------------------- elements

    def basis_element(self, i: int) -> "Octonion":
        """Return basis element e_i  (0 <= i <= 7)."""
        coords = np.zeros(8, dtype=float)
        coords[i] = 1.0
        return Octonion(coords, self)

    def element(self, coords: Sequence[float]) -> "Octonion":
        """Create an octonion from a length-8 real coordinate list."""
        return Octonion(np.asarray(coords, dtype=float), self)

    @property
    def basis(self) -> List["Octonion"]:
        """All 8 basis elements [e_0, e_1, ..., e_7]."""
        return [self.basis_element(i) for i in range(8)]

    def zero(self) -> "Octonion":
        return Octonion(np.zeros(8, dtype=float), self)

    # --------------------------------------------------------------- multiply

    def _mul_coords(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """Coordinate-level multiplication (used by Octonion.__mul__)."""
        result = np.zeros(8, dtype=float)
        for i in range(8):
            ai = a[i]
            if ai == 0.0:
                continue
            row = self._table[i]
            for j in range(8):
                bj = b[j]
                if bj == 0.0:
                    continue
                sign, k = row[j]
                result[k] += sign * ai * bj
        return result

    def __repr__(self) -> str:
        return f"OctonionAlgebra(name={self.name!r})"


# ---------------------------------------------------------------------------
# Octonion element
# ---------------------------------------------------------------------------

class Octonion:
    """
    An element of an OctonionAlgebra.

    self.coords[i] is the real coefficient of basis element e_i.
    """

    __slots__ = ("coords", "algebra")

    def __init__(self, coords: np.ndarray, algebra: OctonionAlgebra) -> None:
        self.coords: np.ndarray = coords
        self.algebra: OctonionAlgebra = algebra

    # ------------------------------------------------------------ arithmetic

    def __mul__(self, other: "Octonion") -> "Octonion":
        return Octonion(
            self.algebra._mul_coords(self.coords, other.coords),
            self.algebra,
        )

    def __add__(self, other: "Octonion") -> "Octonion":
        return Octonion(self.coords + other.coords, self.algebra)

    def __sub__(self, other: "Octonion") -> "Octonion":
        return Octonion(self.coords - other.coords, self.algebra)

    def __neg__(self) -> "Octonion":
        return Octonion(-self.coords, self.algebra)

    def __rmul__(self, scalar: float) -> "Octonion":
        return Octonion(float(scalar) * self.coords, self.algebra)

    # ---------------------------------------------------- algebraic structure

    def conjugate(self) -> "Octonion":
        """
        Octonion conjugate: x-bar_0 = x_0,  x-bar_i = -x_i  for i >= 1.
        Reference: [Wikipedia] Octonion - Conjugate and norm.
        """
        c = self.coords.copy()
        c[1:] = -c[1:]
        return Octonion(c, self.algebra)

    def norm_sq(self) -> float:
        """
        Squared norm  N(x) = x_0^2 + x_1^2 + ... + x_7^2.
        Equivalently N(x) = Re(x * x-bar) where Re extracts the e_0 coefficient.
        For a composition algebra: N(x*y) = N(x)*N(y).
        Reference: [Wikipedia] Octonion - Conjugate and norm.
        """
        return float(np.dot(self.coords, self.coords))

    def inverse(self) -> "Octonion":
        """
        Multiplicative inverse  x^{-1} = x-bar / N(x),  valid for x != 0.
        Octonions form a normed division algebra, so every nonzero element
        is invertible.
        Reference: [Wikipedia] Octonion - Inverse element.
        """
        n = self.norm_sq()
        if n == 0.0:
            raise ZeroDivisionError("The zero element has no multiplicative inverse.")
        return Octonion(self.conjugate().coords / n, self.algebra)

    # --------------------------------------------------------------- helpers

    def is_close(self, other: "Octonion", atol: float = 1e-10) -> bool:
        """Return True if all coordinates agree to within atol."""
        return bool(np.allclose(self.coords, other.coords, atol=atol))

    def __repr__(self) -> str:
        labels = ["e0", "e1", "e2", "e3", "e4", "e5", "e6", "e7"]
        parts = [
            f"{c:+g}*{lab}"
            for c, lab in zip(self.coords, labels)
            if abs(c) > 1e-12
        ]
        return "Oct(" + ("".join(parts) or "0") + ")"


# ---------------------------------------------------------------------------
# Module-level algebra instances
# ---------------------------------------------------------------------------

#: Standard octonion algebra.
#: Used (up to index relabelling) by both Wilson [Wilson2009] and Dixon [Dixon2010].
STANDARD_ALGEBRA = OctonionAlgebra(STANDARD_FANO_TRIPLES, name="standard")

#: Dixon's X-product algebra.
#: Isomorphic to STANDARD_ALGEBRA as abstract algebras, but the multiplication
#: table on the same vector space is different.
#: Used in Dixon's E8 lattice construction [Dixon2010].
XPRODUCT_ALGEBRA = OctonionAlgebra(XPRODUCT_FANO_TRIPLES, name="xproduct")


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def associator(x: Octonion, y: Octonion, z: Octonion) -> Octonion:
    """
    The associator  (x, y, z) = (x*y)*z - x*(y*z).

    Measures the failure of associativity.  For an alternative algebra the
    associator is totally antisymmetric (alternating in all three arguments):
        (x, y, z) = -(y, x, z) = -(x, z, y)   etc.
    This implies left alternativity, right alternativity, and flexibility.
    Reference: [Wikipedia] Octonion - Associativity;
               R.D. Schafer, Introduction to Nonassociative Algebras,
               Academic Press, 1966, Theorem 3.1.
    """
    return (x * y) * z - x * (y * z)
