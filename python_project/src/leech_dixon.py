"""
leech_dixon.py — Dixon's octonionic Leech lattice construction.

Dixon [Dixon2010] Section 6 characterises the 196,560 minimal vectors of
the Leech lattice in three explicit families built from the two E8
representations A^odd and A^even, together with the base element ℓ₀.

Key definitions [Dixon2010 Section 6]
--------------------------------------
  ℓ₀ = ½(e₀+e₁+...+e₇)                  N(ℓ₀) = 2  [eq. (6)]
  A^odd  = A₁ ∪ A₃   (240 norm-2 elements)
  A^even = Ξ^even     (240 norm-1 elements)
  e_a    = octonion basis element, a ∈ {0,...,7}

The three families are generated as follows (all three cyclic permutations
of each triple are included):

  Type 1 [eq. (22)]:
    (2P, 0, 0),  P ∈ A^odd
    N_std = N(2P) = 4·N(P) = 4·2 = 8.
    Count: 240 × 3 = 720.

  Type 2 [eq. (21)]:
    (2Q, ±2e_aQ, 0),  Q ∈ A^even, a ∈ {0,...,7}
    N_std = N(2Q) + N(±2e_aQ) = 4·N(Q) + 4·N(e_a)·N(Q) = 4 + 4 = 8.
    Count: 240 × 8 × 2 × 3 = 11,520.

  Type 3 [eq. (20)]:
    (P, ±e_aP, ±(e_aℓ₀)(e_cP)),  P ∈ A^odd, a,c ∈ {0,...,7}, ± independent
    N_std = N(P) + N(e_aP) + N((e_aℓ₀)(e_cP))
          = 2 + N(e_a)·N(P) + N(e_a)·N(ℓ₀)·N(e_c)·N(P)
          = 2 + 2 + 4 = 8.
    Count: 240 × 8 × 8 × 2 × 2 × 3 = 184,320.

  Total: 720 + 11,520 + 184,320 = 196,560.

Relation to Wilson's type-3 [Wilson2009 p. 2189 footnote]
-----------------------------------------------------------
Wilson notes that Dixon's type-3 formula corresponds to his own formula
((λs)j, λk, (λj)k) but with the j-multiplication omitted from the first
component.  Dixon's (P, ±e_aP, ±(e_aℓ₀)(e_cP)) has norm structure (2, 2, 4)
per component; the 3rd cyclic permutation (±(e_aℓ₀)(e_cP), P, ±e_aP) has
structure (4, 2, 2), matching Wilson's type-3 norm structure.

Wilson's type-3 uses s = ½(−e₀+e₁+...+e₇) in the first component of the
third cyclic permutation; Dixon uses ℓ₀ = ½(e₀+...+e₇) in an analogous
role.  Both have norm 2, but they are different elements (s has one minus
sign; ℓ₀ is the all-positive half-integer element).

A^odd versus Wilson's E8 roots
---------------------------------
Dixon's A^odd = A₁ ∪ A₃, where A₃ has an EVEN number of minus signs in the
half-integer representation.  Wilson's E8 roots = A₁ ∪ Wilson-type-2, where
the Wilson type-2 roots have an ODD number of minus signs.  The two sets of
128 half-integer elements are disjoint (they partition all 256 sign patterns).

As a consequence, Dixon's type-1 vectors (2P, 0, 0) for P ∈ A₃ are NOT in
Wilson's lattice Λ (they require 2P ∈ Ls̄ ∩ Ls = 2L, but 2P ∉ 2L because 2P
has coordinates ±1 which are not all even).  See test_leech_dixon.py for the
computational confirmation of this finding.

Integral octonion lattice interpretation
-----------------------------------------
Dixon [Dixon2010] Section 7 asks whether Λ₂₄ can be given its own algebraic
structure: "No idea yet."  The constructions in this module provide the
computational tools for investigating this question.

References
----------
[Dixon2010]   G.M. Dixon, "Integral Octonions, Octonion XY-Product, and the
              Leech Lattice", preprint, 2010.
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186-2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
"""

import numpy as np
from typing import List, Tuple

from octonions import STANDARD_ALGEBRA, Octonion
from e8_dixon import dixon_xi_even, dixon_a_odd, dixon_a1, dixon_a3


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ELLO_0: Octonion = STANDARD_ALGEBRA.element(0.5 * np.ones(8))
"""
Dixon's base element ℓ₀ = ½(e₀+e₁+...+e₇).

N(ℓ₀) = 8 × (½)² = 2.  This is the all-positive element of A₃ (EVEN number
of minus signs = 0).  Dixon uses ℓ₀ to construct the third component of the
type-3 minimal vectors.  Reference: [Dixon2010] Section 6, eq. (6).
"""

# Precomputed at module load time — keeps generator functions O(n) without
# recomputing the same products in different calls.
_a_odd:  List[Octonion] = dixon_a_odd()    # 240 norm-2 elements (A₁ ∪ A₃)
_a_even: List[Octonion] = dixon_xi_even()  # 240 norm-1 elements (Ξ^even)

# Octonion basis elements e_0,...,e_7.
_BASIS: List[Octonion] = [STANDARD_ALGEBRA.basis_element(a) for a in range(8)]

# Precomputed e_a × ℓ₀ for each basis element a (used in type-3 third component).
_ea_l0: List[Octonion] = [ea * ELLO_0 for ea in _BASIS]


# ---------------------------------------------------------------------------
# Minimal vector generation
# ---------------------------------------------------------------------------

def dixon_leech_type1() -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Return the 720 type-1 minimal vectors of Dixon's Leech construction.

    Formula [Dixon2010 eq. (22)]: (2P, 0, 0) and its 3 cyclic permutations,
    where P runs over all 240 elements of A^odd.

    Norm: N_std(2P, 0, 0) = 4·N(P) = 4·2 = 8. ✓
    Count: 240 × 3 = 720.

    Note: A^odd = A₁ ∪ A₃.  Type-1 vectors with P ∈ A₁ satisfy Wilson's
    membership conditions (is_in_leech).  Type-1 vectors with P ∈ A₃ do
    NOT (2P ∉ 2L = Ls̄ ∩ Ls, because 2P has coordinates ±1 which are not
    all even).  See test_leech_dixon.py for computational confirmation.
    Reference: [Dixon2010] Section 6, eq. (22).
    """
    zero = np.zeros(8)
    result: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for P in _a_odd:
        two_P = 2.0 * P.coords
        result.append((two_P.copy(), zero.copy(), zero.copy()))
        result.append((zero.copy(), two_P.copy(), zero.copy()))
        result.append((zero.copy(), zero.copy(), two_P.copy()))
    return result


def dixon_leech_type2() -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Return the 11,520 type-2 minimal vectors of Dixon's Leech construction.

    Formula [Dixon2010 eq. (21)]: (2Q, ±2e_aQ, 0) and its 3 cyclic
    permutations, where Q runs over all 240 elements of A^even and
    e_a ∈ {e₀,...,e₇}.

    Norm: N_std = N(2Q) + N(±2e_aQ) = 4·N(Q) + 4·N(e_a)·N(Q) = 4 + 4 = 8. ✓
    Count: 240 × 8 × 2 × 3 = 11,520.
    Reference: [Dixon2010] Section 6, eq. (21).
    """
    zero = np.zeros(8)
    result: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for Q in _a_even:
        two_Q = 2.0 * Q.coords
        for ea in _BASIS:
            two_ea_Q = 2.0 * (ea * Q).coords   # 2e_aQ
            for sign in (+1, -1):
                comp2 = sign * two_ea_Q          # ±2e_aQ
                # 3 cyclic permutations of (2Q, ±2e_aQ, 0)
                result.append((two_Q.copy(), comp2.copy(), zero.copy()))
                result.append((zero.copy(), two_Q.copy(), comp2.copy()))
                result.append((comp2.copy(), zero.copy(), two_Q.copy()))
    return result


def dixon_leech_type3() -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Return the 184,320 type-3 minimal vectors of Dixon's Leech construction.

    Formula [Dixon2010 eq. (20)]:
        (P, ±e_aP, ±(e_aℓ₀)(e_cP)) and its 3 cyclic permutations,
    where P ∈ A^odd, a,c ∈ {0,...,7}, and the two ± signs are independent.

    Norm: N_std = N(P) + N(e_aP) + N((e_aℓ₀)(e_cP))
               = 2 + N(e_a)·N(P) + N(e_a)·N(ℓ₀)·N(e_c)·N(P)
               = 2 + 2 + 4 = 8. ✓
    Count: 240 × 8 × 8 × 2 × 2 × 3 = 184,320.

    Implementation note: (e_aℓ₀) is precomputed for each a (independent of P
    and c).  For each (P, a), the product e_aP is computed once and reused
    across all values of c.  For each (P, c), the product e_cP is computed
    once and reused across all values of a.
    Reference: [Dixon2010] Section 6, eq. (20).
    """
    alg = STANDARD_ALGEBRA
    result: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for P in _a_odd:
        P_c = P.coords.copy()

        # Precompute e_cP for all c (independent of a).
        ec_P_all: List[np.ndarray] = [(_BASIS[c] * P).coords for c in range(8)]

        for idx_a in range(8):
            ea_P = (_BASIS[idx_a] * P).coords   # e_aP
            ea_l0 = _ea_l0[idx_a]               # e_aℓ₀  (precomputed)

            for idx_c in range(8):
                ec_P = ec_P_all[idx_c]
                # (e_aℓ₀)(e_cP) — note: ea_l0 is an Octonion, ec_P is an array
                prod = (ea_l0 * alg.element(ec_P)).coords

                for s1 in (+1, -1):       # sign on e_aP
                    comp2 = s1 * ea_P
                    for s2 in (+1, -1):   # sign on (e_aℓ₀)(e_cP)
                        comp3 = s2 * prod
                        # 3 cyclic permutations of (P, ±e_aP, ±(e_aℓ₀)(e_cP))
                        result.append((P_c.copy(), comp2.copy(), comp3.copy()))
                        result.append((comp3.copy(), P_c.copy(), comp2.copy()))
                        result.append((comp2.copy(), comp3.copy(), P_c.copy()))
    return result


def dixon_leech_minimal_vectors() -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Return all 196,560 minimal vectors of Dixon's Leech construction.

    Combines the three families from [Dixon2010] Section 6 eqs. (20)–(22).
    Note: this function computes ~200,000 octonion products and may take
    several seconds.  Use the individual dixon_leech_type*() functions when
    only a subset is needed.
    """
    return (dixon_leech_type1()
            + dixon_leech_type2()
            + dixon_leech_type3())
