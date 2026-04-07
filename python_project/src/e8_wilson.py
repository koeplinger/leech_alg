"""
e8_wilson.py — Wilson's E8 lattice construction.

Wilson [Wilson2009] defines an E8 lattice L inside the real octonion algebra
as the Z-integer span of 240 octonions:

  Type-1 (112 roots):  ±e_a ± e_b  for distinct a, b in {0,...,7}
  Type-2 (128 roots):  (1/2)(±e_0 ± e_1 ± ... ± e_7)  with an ODD number
                        of minus signs.

Index convention
----------------
Wilson labels his octonion basis {i_inf, i_0, ..., i_6} with index set
PL(7) = {inf} u F_7.  This module uses {e_0, ..., e_7}.  The correspondence is:
    Wilson  i_inf  <->  e_0   (real unit / identity)
    Wilson  i_k    <->  e_{k+1}   for k = 0,...,6.

Under this map Wilson's type-1 roots ±i_t ± i_u become ±e_a ± e_b, and his
type-2 roots (1/2)(±1 ± i_0 ± ... ± i_6) with odd # minus signs become
(1/2)(±e_0 ± e_1 ± ... ± e_7) with odd # minus signs.  The two descriptions
are identical in structure.

Key structural element
-----------------------
    s = (1/2)(-1 + i_0 + i_1 + ... + i_6)
      = (1/2)(-e_0 + e_1 + e_2 + ... + e_7)   [under our index map]

Wilson proves: s in L, s-bar in R = L-bar (the conjugate lattice).
Key structural results from [Wilson2009] Section 2:
    N(s) = 2
    Ls = 2B   (B = another copy of E8 / Coxeter-Dickson ring)
    2L ⊂ Ls ⊂ L

Lattice membership criterion
------------------------------
The lattice L equals the D8+ lattice in standard E8 coordinates:

    L = { x in R^8 : all x_i in Z,     sum(x_i) even }
      u { x in R^8 : all x_i in Z+1/2, sum(x_i) odd  }

Verification: every type-1 root ±e_a ± e_b has integer coordinates summing
to 0 or ±2 (both even).  Every type-2 root (1/2)(±e_0...±e_7) with k odd
minus signs has all coordinates in Z+1/2 and sum equal to (8-2k)/2 = 4-k;
since k is odd, 4-k is odd.  Both types satisfy the criterion.

References
----------
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186-2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
"""

import numpy as np
from typing import List

from octonions import STANDARD_ALGEBRA, Octonion


# ---------------------------------------------------------------------------
# Root generators
# ---------------------------------------------------------------------------

def wilson_type1_roots() -> List[Octonion]:
    """
    Return the 112 type-1 roots of L: ±e_a ± e_b for distinct a,b in {0,...,7}.

    Count: C(8,2) * 4 = 28 * 4 = 112.

    Reference: [Wilson2009] Section 2 — "112 octonions ±i_t ± i_u for any
    distinct t, u in PL(7)".  Under our index map i_t -> e_{t+1} (t in F_7)
    and i_inf -> e_0 these become ±e_a ± e_b for distinct a, b in {0,...,7}.
    """
    roots: List[Octonion] = []
    for a in range(8):
        for b in range(a + 1, 8):
            for sa in (+1, -1):
                for sb in (+1, -1):
                    coords = np.zeros(8)
                    coords[a] = float(sa)
                    coords[b] = float(sb)
                    roots.append(STANDARD_ALGEBRA.element(coords))
    return roots


def wilson_type2_roots() -> List[Octonion]:
    """
    Return the 128 type-2 roots of L:
        (1/2)(±e_0 ± e_1 ± ... ± e_7)   with an ODD number of minus signs.

    Count: 2^8 / 2 = 128  (exactly half of all sign patterns have odd minus count).

    Reference: [Wilson2009] Section 2 — "128 octonions (1/2)(±1 ± i_0 ... ± i_6)
    which have an odd number of minus signs".
    """
    roots: List[Octonion] = []
    for mask in range(256):
        # Bit k of mask = 1 means coordinate e_k gets a minus sign.
        num_minus = bin(mask).count("1")
        if num_minus % 2 == 1:  # odd number of minus signs
            signs = np.array(
                [(-1.0 if (mask >> k) & 1 else 1.0) for k in range(8)]
            )
            roots.append(STANDARD_ALGEBRA.element(0.5 * signs))
    return roots


def wilson_e8_roots() -> List[Octonion]:
    """
    Return all 240 roots of Wilson's E8 lattice L.

    The 240 roots are the vectors of minimal norm (norm^2 = 2) in L, and they
    span L over Z.  Reference: [Wilson2009] Section 2.
    """
    return wilson_type1_roots() + wilson_type2_roots()


# ---------------------------------------------------------------------------
# Wilson's special element s
# ---------------------------------------------------------------------------

# s = (1/2)(-1 + i_0 + i_1 + ... + i_6) in Wilson's notation.
# Under our index map: s = (1/2)(-e_0 + e_1 + e_2 + ... + e_7).
# Properties: s in L (it is a type-2 root with exactly 1 minus sign),
# N(s) = 2, and Ls = 2B  where B is the Coxeter-Dickson ring.
# Reference: [Wilson2009] Section 2.
_s_coords = np.array([-0.5] + [0.5] * 7)
WILSON_S: Octonion = STANDARD_ALGEBRA.element(_s_coords)


# ---------------------------------------------------------------------------
# Inner product
# ---------------------------------------------------------------------------

def inner_product(x: Octonion, y: Octonion) -> float:
    """
    Standard real inner product <x, y> = sum_i x_i * y_i.

    For octonions, this equals Re(x * y-bar) where Re extracts the e_0
    coefficient.  The norm satisfies N(x) = <x, x>.
    """
    return float(np.dot(x.coords, y.coords))


# ---------------------------------------------------------------------------
# Lattice membership test
# ---------------------------------------------------------------------------

def is_in_L(x: Octonion, tol: float = 1e-9) -> bool:
    """
    Return True iff x lies in Wilson's E8 lattice L = D8+.

    Criterion:
      - All-integer case:     all x_i in Z  AND  sum(x_i) even.
      - Half-integer case:    all x_i in Z+1/2  AND  sum(x_i) odd integer.
      - All other cases: False.

    Implementation: multiply coords by 2 to obtain integers, then inspect
    parities.  Let c2_i = round(2 * x_i).

      c2_i all even  =>  x_i all integers;  require sum(c2_i) % 4 == 0
                         (equivalent to sum(x_i) divisible by 2).
      c2_i all odd   =>  x_i all half-integers; require sum(c2_i) % 4 == 2
                         (equivalent to sum(x_i) an odd integer; note that
                          sum of 8 odd integers is even, so sum(c2_i) is even,
                          and sum(x_i) = sum(c2_i)/2 is an integer; parity
                          check then becomes sum(c2_i)/2 odd).
      Mixed parities =>  not in L.

    Reference: standard characterisation of E8 = D8+.  Consistency with
    [Wilson2009] Section 2 verified by checking both root types satisfy the
    criterion (see module docstring).
    """
    c2 = np.round(x.coords * 2.0)
    if not np.allclose(x.coords * 2.0, c2, atol=tol):
        # Coords are not multiples of 1/2.
        return False
    c2i = c2.astype(int)
    parities = c2i % 2   # 0 for even, 1 for odd (Python % is non-negative)
    if np.all(parities == 0):
        # All integer coordinates.  Even sum iff sum(c2i) % 4 == 0.
        return int(np.sum(c2i)) % 4 == 0
    if np.all(parities == 1):
        # All half-integer coordinates.  Odd integer sum iff sum(c2i) % 4 == 2.
        return int(np.sum(c2i)) % 4 == 2
    # Mixed: one coord integer, another half-integer — not in L.
    return False
