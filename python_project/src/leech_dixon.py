"""
leech_dixon.py — Dixon's octonionic Leech lattice construction.

Created: 2026-04-07.  Updated: 2026-10-31.

Dixon [Dixon2010] Section 6 characterizes the 196,560 minimal vectors of
a Leech lattice in three explicit families built from the two E8
representations A^odd and A^even, together with the base element ℓ₀.

Key definitions [Dixon2010 Sections 3 and 6]
---------------------------------------------
  ℓ₀ = ½(e₀+e₁+...+e₇)                  N(ℓ₀) = 2  [Section 3, eq. (6)]
  A^odd  = A₁ ∪ A₃   (240 norm-2 elements)
  A^even = Ξ^even     (240 norm-1 elements)
  e_a    = octonion basis element, a ∈ {0,...,7}
All products are the octonion product of Dixon's table (3),
e_a e_{a+1} = e_{a+3}, which is STANDARD_ALGEBRA.

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
    (P, ±e_aP, ±e_a(ℓ₀(e_cP))),  P ∈ A^odd, a,c ∈ {0,...,7}, ± independent
    N_std = N(P) + N(e_aP) + N(e_a(ℓ₀(e_cP)))
          = 2 + N(e_a)·N(P) + N(e_a)·N(ℓ₀)·N(e_c)·N(P)
          = 2 + 2 + 4 = 8.
    Count: 240 × 8 × 8 × 2 × 2 × 3 = 184,320.
    Dixon prints the third component also as (e_aℓ₀e_a)(e_bP), with e_b
    depending on e_a, e_c and P; both forms give the same set of vectors.
    The product nesting matters: (e_aℓ₀)(e_cP) is in general a different
    octonion.

  Total: 720 + 11,520 + 184,320 = 196,560.

Dixon writes "we include all permutations"; the three cyclic shifts used
here already give a set closed under all six orderings of the components.

Relation to Wilson's construction [Wilson2009 Section 3]
--------------------------------------------------------
Dixon derives eq. (20) from Wilson's largest family ((λs)j, λk, (λj)k),
written in Dixon's notation, by conjugating each component, shifting the
components cyclically and replacing P by e_cP [Dixon2010 Section 6].  The
196,560 vectors built here are exactly the componentwise conjugates
(x̄, ȳ, z̄) of Wilson's 196,560 minimal vectors (leech_wilson.py), so they
are the minimal vectors of a copy of the Leech lattice, the image of
Wilson's Λ under componentwise conjugation (test_leech_dixon.py;
verify_dixon_leech_lattice.py, which also checks the lattice directly).

Wilson's bracketed note on Dixon ([Wilson2009] Section 5, p. 2189) concerns
the formula of Dixon's 2005 preprint [Dixon2005], not the 2010 construction
implemented here: "in effect, his formula is equivalent to
(λs, ±λk, ±(λj)k), that is, omitting the (necessary) multiplication by j in
the first coordinate."

A^odd versus Wilson's E8 roots
---------------------------------
Dixon's A^odd = A₁ ∪ A₃, where A₃ has an EVEN number of minus signs in the
half-integer representation.  Wilson's E8 roots = A₁ ∪ Wilson-type-2, where
the Wilson type-2 roots have an ODD number of minus signs.  The two sets of
128 half-integer elements are disjoint (they partition all 256 sign patterns).
Conjugation changes the sign of seven coordinates and so maps the one set
onto the other: the conjugates of A^odd are exactly Wilson's E8 roots.

As a consequence, Dixon's type-1 vectors (2P, 0, 0) for P ∈ A₃ are NOT in
Wilson's lattice Λ (they require 2P ∈ Ls̄ ∩ Ls = 2L, that is, P ∈ L,
but P is a half-integer vector with an even coordinate sum, so P ∉ L; in fact
2P lies in Ls̄ but not in Ls); their conjugates are.  See test_leech_dixon.py
for the computational confirmation.

Integral octonion lattice interpretation
-----------------------------------------
Dixon [Dixon2010] Section 7 asks whether Λ₂₄ can be given its own algebraic
structure: "No idea yet."  The constructions in this module provide the
computational tools for investigating this question.

References
----------
[Dixon2005]   G.M. Dixon, "Integral octonions, octonion multiplication and the
              Leech lattice", preprint, 10 March 2005; reference [12] of
              [Wilson2009].
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
type-3 minimal vectors.  Reference: [Dixon2010] Section 3, eq. (6).
"""

# Precomputed at module load time, so that the generator functions stay O(n)
# and do not recompute the same products in different calls.
_a_odd:  List[Octonion] = dixon_a_odd()    # 240 norm-2 elements (A₁ ∪ A₃)
_a_even: List[Octonion] = dixon_xi_even()  # 240 norm-1 elements (Ξ^even)

# Octonion basis elements e_0,...,e_7.
_BASIS: List[Octonion] = [STANDARD_ALGEBRA.basis_element(a) for a in range(8)]


# ---------------------------------------------------------------------------
# Minimal vector generation
# ---------------------------------------------------------------------------

def dixon_leech_type1() -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Return the 720 type-1 minimal vectors of Dixon's Leech construction.

    Formula [Dixon2010 eq. (22)]: (2P, 0, 0) and its 2 cyclic shifts,
    where P runs over all 240 elements of A^odd.

    Norm: N_std(2P, 0, 0) = 4·N(P) = 4·2 = 8. ✓
    Count: 240 × 3 = 720.

    Note: A^odd = A₁ ∪ A₃.  Type-1 vectors with P ∈ A₁ satisfy Wilson's
    membership conditions (is_in_leech).  Type-1 vectors with P ∈ A₃ do
    NOT (2P ∉ 2L = Ls̄ ∩ Ls, because P ∉ L: P is a half-integer
    vector with an even coordinate sum).  See test_leech_dixon.py for computational confirmation.
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

    Formula [Dixon2010 eq. (21)]: (2Q, ±2e_aQ, 0) and its 2 cyclic
    shifts, where Q runs over all 240 elements of A^even and
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
        (P, ±e_aP, ±e_a(ℓ₀(e_cP))) and its 2 cyclic shifts,
    where P ∈ A^odd, a,c ∈ {0,...,7}, and the two ± signs are independent.

    Norm: N_std = N(P) + N(e_aP) + N(e_a(ℓ₀(e_cP)))
               = 2 + N(e_a)·N(P) + N(e_a)·N(ℓ₀)·N(e_c)·N(P)
               = 2 + 2 + 4 = 8. ✓
    Count: 240 × 8 × 8 × 2 × 2 × 3 = 184,320.

    Implementation note: for each P, the products e_cP and ℓ₀(e_cP) are
    computed once for every c and reused across all values of a; e_aP is
    computed once for each (P, a).
    Reference: [Dixon2010] Section 6, eq. (20).
    """
    result: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for P in _a_odd:
        P_c = P.coords.copy()

        # ℓ₀(e_cP) for all c (independent of a).
        l0_ec_P_all: List[Octonion] = [ELLO_0 * (_BASIS[c] * P) for c in range(8)]

        for idx_a in range(8):
            ea = _BASIS[idx_a]
            ea_P = (ea * P).coords   # e_aP

            for idx_c in range(8):
                prod = (ea * l0_ec_P_all[idx_c]).coords   # e_a(ℓ₀(e_cP))

                for s1 in (+1, -1):       # sign on e_aP
                    comp2 = s1 * ea_P
                    for s2 in (+1, -1):   # sign on e_a(ℓ₀(e_cP))
                        comp3 = s2 * prod
                        # 3 cyclic permutations of (P, ±e_aP, ±e_a(ℓ₀(e_cP)))
                        result.append((P_c.copy(), comp2.copy(), comp3.copy()))
                        result.append((comp3.copy(), P_c.copy(), comp2.copy()))
                        result.append((comp2.copy(), comp3.copy(), P_c.copy()))
    return result


def dixon_leech_minimal_vectors() -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Return all 196,560 minimal vectors of Dixon's Leech construction.

    Combines the three families from [Dixon2010] Section 6 eqs. (20)–(22).
    Note: this function computes 23,040 octonion products and builds
    196,560 vectors (under a second).  Use the individual dixon_leech_type*()
    functions when only a subset is needed.
    """
    return (dixon_leech_type1()
            + dixon_leech_type2()
            + dixon_leech_type3())
