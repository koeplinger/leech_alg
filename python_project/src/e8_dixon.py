"""
e8_dixon.py — Dixon's E8 lattice constructions.

Dixon [Dixon2010] defines two families of 240-element sets inside the octonion
algebra, each forming the minimal-vector shell of an E8 lattice representation.

Ξ^even (primary, rational coefficients, norm² = 1)
---------------------------------------------------
Ξ₀ = {±e_a : a = 0,...,7}                    16 elements
Ξ₂ = {(±e_a ± e_b ± e_c ± e_d)/2 :           224 elements
       a,b,c,d distinct,
       e_a(e_b(e_c*e_d)) = ±1}
Ξ^even = Ξ₀ ∪ Ξ₂

A^odd (rational coefficients, norm² = 2)
-----------------------------------------
A₁ = {±e_a ± e_b : a,b distinct}              112 elements  (same as Wilson's type-1)
A₃ = {(∑_{a=0}^{7} ±e_a)/2 : even # of +'s}  128 elements  (differ from Wilson's type-2)
A^odd = A₁ ∪ A₃

Dixon defines Ξ₁ and Ξ₃ identically but with irrational normalizations 1/√2 and 1/√8.
In Section 6 of [Dixon2010] he introduces A₁ and A₃ (the rational-coefficient equivalents
used in the Leech lattice construction).  Only A₁, A₃ (rational) are implemented here
since they live in the same vector space as Ξ^even.

Key property [Dixon2010 eq. (5)]
----------------------------------
For all X ∈ Ξ^even ∪ Ξ^odd and all basis pairs (e_a, e_b):
    e_a ◦_X e_b = (e_a * X) * (X^{-1} * e_b) = ±e_c
for some basis element e_c.  We verify this for X ∈ Ξ^even (norm 1, X^{-1} = X̄)
and X ∈ A^odd (norm √2, X^{-1} = X̄/2).

Structure of Ξ₂
-----------------
The condition e_a(e_b(e_c*e_d)) = ±1 is satisfied by exactly 14 unordered 4-element
subsets of {0,...,7}:
  Type A (7 sets): {0} ∪ {Fano triple} for each of the 7 Fano triples.
  Type B (7 sets): complement of each Fano triple within {1,...,7}.
Each valid set contributes 2^4 = 16 elements (all sign choices), giving 14*16 = 224.

References
----------
[Dixon2010]   G.M. Dixon, "Integral Octonions, Octonion XY-Product, and the
              Leech Lattice", preprint, 2010.
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186-2190.
"""

import numpy as np
from typing import List, Set, Tuple

from octonions import STANDARD_ALGEBRA, Octonion


# ---------------------------------------------------------------------------
# Ξ^even
# ---------------------------------------------------------------------------

def dixon_xi0() -> List[Octonion]:
    """
    Return the 16 elements of Ξ₀ = {±e_a : a = 0,...,7}.

    These are the 8 octonion basis elements with both signs.
    All have norm² = 1.
    Reference: [Dixon2010] Section 2, eq. (4).
    """
    result: List[Octonion] = []
    for a in range(8):
        for s in (+1, -1):
            coords = np.zeros(8)
            coords[a] = float(s)
            result.append(STANDARD_ALGEBRA.element(coords))
    return result


def dixon_xi2() -> List[Octonion]:
    """
    Return the 224 elements of Ξ₂.

    Ξ₂ = {(±e_a ± e_b ± e_c ± e_d)/2 : a,b,c,d distinct,
           e_a(e_b(e_c*e_d)) = ±1}

    Implementation: enumerate all ordered 4-tuples (a,b,c,d) of distinct
    indices from {0,...,7}.  For each tuple where the condition holds, add
    all 2^4 = 16 sign-choice elements to the result.  A Python set prevents
    duplicates (different orderings of the same 4-index set yield identical
    16-element families).
    Reference: [Dixon2010] Section 2, eq. (4).
    """
    alg = STANDARD_ALGEBRA
    seen: Set[Tuple[int, ...]] = set()
    result: List[Octonion] = []

    for a in range(8):
        for b in range(8):
            if b == a:
                continue
            for c in range(8):
                if c == a or c == b:
                    continue
                for d in range(8):
                    if d == a or d == b or d == c:
                        continue
                    # Compute e_a * (e_b * (e_c * e_d))
                    prod = (alg.basis_element(a) *
                            (alg.basis_element(b) *
                             (alg.basis_element(c) * alg.basis_element(d))))
                    # Check: result = ±e_0 (the real unit)
                    if (abs(abs(prod.coords[0]) - 1.0) < 1e-9 and
                            np.all(np.abs(prod.coords[1:]) < 1e-9)):
                        # Add all 16 sign-pattern elements (with deduplication)
                        for mask in range(16):
                            coords = np.zeros(8)
                            for k, idx in enumerate((a, b, c, d)):
                                coords[idx] = -0.5 if (mask >> k) & 1 else 0.5
                            key = tuple(np.round(coords * 2).astype(int))
                            if key not in seen:
                                seen.add(key)
                                result.append(alg.element(coords.copy()))
    return result


def dixon_xi_even() -> List[Octonion]:
    """
    Return all 240 elements of Ξ^even = Ξ₀ ∪ Ξ₂.

    Dixon calls these "a representation of the inner shell of an E8 lattice."
    All elements have norm² = 1.
    Reference: [Dixon2010] Section 2.
    """
    return dixon_xi0() + dixon_xi2()


# ---------------------------------------------------------------------------
# A^odd  (rational-coefficient version of Ξ^odd, used in the Leech construction)
# ---------------------------------------------------------------------------

def dixon_a1() -> List[Octonion]:
    """
    Return the 112 elements of A₁ = {±e_a ± e_b : a,b distinct}.

    These are the same as Wilson's type-1 roots.  All have norm² = 2.
    Reference: [Dixon2010] Section 6.
    """
    result: List[Octonion] = []
    for a in range(8):
        for b in range(a + 1, 8):
            for sa in (+1, -1):
                for sb in (+1, -1):
                    coords = np.zeros(8)
                    coords[a] = float(sa)
                    coords[b] = float(sb)
                    result.append(STANDARD_ALGEBRA.element(coords))
    return result


def dixon_a3() -> List[Octonion]:
    """
    Return the 128 elements of A₃.

    A₃ = {(∑_{a=0}^{7} ±e_a)/2 : even number of plus signs}

    'Even number of +'s' means the number of +1 coefficients is even.  With
    8 coordinates this also means even number of -1 coefficients (since 8 is even).
    All elements have norm² = 8*(1/2)² = 2.

    Note: these are DIFFERENT from Wilson's type-2 roots, which require an ODD
    number of minus signs.  A₃ and Wilson's type-2 together cover all 2^8 = 256
    sign patterns.
    Reference: [Dixon2010] Section 6.
    """
    result: List[Octonion] = []
    for mask in range(256):
        # mask bit k = 1 means coordinate k gets a minus sign
        num_minus = bin(mask).count('1')
        # even # of +'s  ↔  even # of -'s  (since 8 - even = even)
        if num_minus % 2 == 0:
            signs = np.array(
                [(-1.0 if (mask >> k) & 1 else 1.0) for k in range(8)]
            )
            result.append(STANDARD_ALGEBRA.element(0.5 * signs))
    return result


def dixon_a_odd() -> List[Octonion]:
    """
    Return all 240 elements of A^odd = A₁ ∪ A₃.

    All elements have norm² = 2.  This includes ℓ₀ = ½(e_0+e_1+...+e_7)
    (all-positive A₃ element), the defining element of XPRODUCT_ALGEBRA.
    Reference: [Dixon2010] Section 6.
    """
    return dixon_a1() + dixon_a3()


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def is_pm_basis_element(x: Octonion, tol: float = 1e-9) -> bool:
    """
    Return True iff x = ±e_k for some k in {0,...,7}.

    Used to verify the X-product closure property eq. (5) of [Dixon2010].
    """
    c = x.coords
    nonzero_mask = np.abs(c) > tol
    if np.sum(nonzero_mask) != 1:
        return False
    return abs(abs(c[nonzero_mask][0]) - 1.0) < tol
