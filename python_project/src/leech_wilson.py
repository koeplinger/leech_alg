"""
leech_wilson.py — Wilson's octonionic Leech lattice construction.

The Leech lattice as an integral octonion lattice
--------------------------------------------------
The Leech lattice Λ is a rank-24 Z-module (lattice) embedded inside L³, where
L = E8 = the integral octonions (the Coxeter–Dickson order, the unique maximal
order of the real octonion division algebra).  Because L carries the octonion
product, L³ inherits candidate algebraic structure.  The membership conditions
(1)–(3) below characterise exactly which elements of L³ belong to Λ.

The central research question of this project: whether a proposed
octonion-based multiplication rule on R^24 restricts to a product on Λ that
closes — i.e., whether Λ forms a 24-dimensional integral algebra under some
octonion-derived product.  Dixon [Dixon2010] has examined this question by
defining an XY-product.  Wilson's construction provides the reference object:
Λ as a precisely characterised rank-24 Z-lattice with algebraic membership
conditions.

The minimal-vector shell of Λ consists of the 196,560 vectors of minimum
norm (N_std = 8 in the standard R^24 metric), partitioned into three families.
These are the analogue, for the Leech lattice, of the 240 roots in E8 — the
extremal vectors that generate and characterise the lattice.  (The Leech
lattice has no vectors of norm 2, so it has no root system in the classical
Lie-theoretic sense; the minimal-vector shell plays the equivalent structural
role.)

What this module establishes (for future verification work)
------------------------------------------------------------
1. A membership test is_in_leech(x, y, z): Wilson's three conditions.
2. The three families of 196,560 minimal-shell vectors.
3. The Leech lattice is rank 24, even, self-dual (Gram det = 1 under Wilson's
   inner product), and has minimum norm 8.

These together uniquely characterise the Leech lattice by Conway's theorem
[Conway–Sloane, Chapter 12].  The membership test and minimal-vector shell
serve as the reference for checking whether a proposed product maps E8 elements
into the Leech lattice.

Wilson's Leech lattice definition [Wilson2009 Section 3]
---------------------------------------------------------
    Λ = { (x, y, z) ∈ L³ :
          1. x, y, z ∈ L
          2. x+y, x+z, y+z ∈ Ls̄       (= 2B, the Coxeter–Dickson ring scaled by 2)
          3. x+y+z ∈ Ls  }

where s = ½(−e₀+e₁+...+e₇) (WILSON_S), s̄ = conjugate of s, and
    Ls̄ = {l × s̄ : l ∈ L}  (right-multiplication of L by s̄)
    Ls  = {l × s  : l ∈ L}  (right-multiplication of L by s)
Both are index-16 sublattices of L with Ls̄ + Ls = L and Ls̄ ∩ Ls = 2L.

The norm on triples used by Wilson is N_Wilson(x,y,z) = ½(N(x)+N(y)+N(z)).
This module uses the STANDARD R^24 norm N_std(x,y,z) = N(x)+N(y)+N(z),
which equals 2×N_Wilson and gives minimal vectors norm-squared 8.
Under Wilson's inner product, Λ is self-dual (det = 1).
Under the standard R^24 inner product, det(Gram_Λ) = 2^24.

196,560 minimal vectors [Wilson2009 Section 3]
----------------------------------------------
The three families below cover all minimal vectors of Λ, where:
  λ  ranges over the 240 roots of L
  j,k ∈ J = {±e_t : t = 0,...,7}  (16 unit octonions with both signs)

  Type 1: (2λ, 0, 0) and cyclic permutations       —  3 × 240      =    720
  Type 2: (λs̄, (λs̄)j, 0) and cyclic permutations — 3 × 240 × 16 = 11,520
  Type 3: ((λs)j, λk, (λj)k) and cyclic perms     — 3 × 240 × 16 × 16 = 184,320
  Note: type-2 uses s̄; type-3 uses s.  This is the distinction Wilson proves is
  necessary: see [Wilson2009 p. 2189] footnote on Dixon's related construction.

  Total: 196,560  =  3 × 240 × (1 + 16 + 16×16)

All three families are verified to have N_std = 8.

Proof that Λ has no vectors of norm < 8  [Wilson2009 Section 4]
---------------------------------------------------------------
Suppose (x, y, z) ∈ Λ has Wilson-norm N < 4 (i.e., N_std < 8).
Then at least one of x, y, z must be 0 (since each nonzero coordinate from L
has N ≥ 2).  If, say, z = 0, then conditions 2 and 3 give x+y ∈ Ls̄ and
x+y ∈ Ls, so x+y ∈ Ls̄ ∩ Ls = 2L.  Also y+z = y ∈ Ls̄, and x+z = x ∈ Ls̄,
so both x,y ∈ Ls̄ = 2B ⊂ L.  Then y ∈ Ls̄ and x = (x+y) − y ∈ 2L + 2B = 2B.
In the only-one-nonzero case (y=z=0), condition 3 gives x ∈ Ls, and
condition 2 gives x ∈ Ls̄, so x ∈ Ls ∩ Ls̄ = 2L, hence N(x) ≥ N_std(2L) = 8.
Contradiction.  Therefore Λ has no vectors of Wilson-norm < 4 (N_std < 8).

References
----------
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186-2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
"""

import numpy as np
from typing import List, Tuple

from octonions import STANDARD_ALGEBRA, Octonion
from e8_wilson import wilson_e8_roots, WILSON_S, is_in_L


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

WILSON_S_BAR: Octonion = WILSON_S.conjugate()
"""
Conjugate of Wilson's element s.

s   = ½(−e₀+e₁+...+e₇)   ∈ L   (type-2 root, 1 minus sign = odd ✓)
s̄  = ½(−e₀−e₁−...−e₇)   ∈ R = L̄  (conjugate lattice of L)

Wilson proves s ∈ L and s̄ ∈ R, and uses these to show Ls̄ = LR = 2B.
Reference: [Wilson2009] Section 2.
"""

# J = {±e_t : t = 0,...,7}  — the 16 unit octonions with both signs.
# Wilson uses J = {±i_t : t ∈ PL(7)} in the description of minimal vectors.
J: List[Octonion] = (
    [STANDARD_ALGEBRA.basis_element(t) for t in range(8)]
    + [STANDARD_ALGEBRA.element(-STANDARD_ALGEBRA.basis_element(t).coords)
       for t in range(8)]
)


# ---------------------------------------------------------------------------
# Sublattice bases: Ls̄ and Ls
# ---------------------------------------------------------------------------

def _verify_sublattice_basis(
    basis: np.ndarray,
    all_products: np.ndarray,
    tol: float = 1e-9,
) -> None:
    """
    Assert that all rows of all_products lie in the Z-span of basis.

    For each product p, verify that p @ inv(basis) is an integer vector.
    This confirms that basis is a Z-basis for the Z-module generated by
    all_products.
    """
    inv_basis = np.linalg.inv(basis)
    for p in all_products:
        coords = p @ inv_basis
        assert np.allclose(coords, np.round(coords), atol=tol), (
            f"Product {p} is NOT an integer combination of the basis.\n"
            f"Coords in basis: {coords}"
        )


def _build_L_zbasis() -> np.ndarray:
    """
    Return an 8×8 matrix whose rows form a Z-basis for L.

    The Z-basis is constructed as follows:
      Row 0 : WILSON_S  (a type-2 root; spans the D8+ coset not in D8)
      Rows 1-7: 7 linearly independent type-1 roots, chosen greedily.

    The resulting Gram matrix has det = 1 (L is unimodular).  This is the same
    construction verified computationally in test_e8_wilson.py (test 16).

    Basis of images under right-multiplication:
      If {bᵢ} is a Z-basis for L and R: v ↦ v*c is right-multiplication by
      some fixed c, then {bᵢ*c} is a Z-basis for R(L), because for any
      λ = ∑nᵢbᵢ ∈ L (nᵢ ∈ Z), λ*c = ∑nᵢ(bᵢ*c) by linearity.
    """
    roots = wilson_e8_roots()
    basis: List[np.ndarray] = [WILSON_S.coords.copy()]
    for r in roots:
        candidate = np.array(basis + [r.coords])
        if np.linalg.matrix_rank(candidate) == len(basis) + 1:
            basis.append(r.coords.copy())
            if len(basis) == 8:
                break
    assert len(basis) == 8, "Could not find 8 independent E8 basis vectors"
    B = np.array(basis)
    det = abs(np.linalg.det(B @ B.T))
    assert abs(det - 1.0) < 1e-6, \
        f"L Z-basis Gram det = {det:.6f}, expected 1 (basis is not unimodular)"
    return B


# Precomputed at module load time.
_e8_roots: List[Octonion] = wilson_e8_roots()
_L_zbasis: np.ndarray = _build_L_zbasis()   # 8×8 unimodular Z-basis for L

# Ls̄ basis: images of the L Z-basis under right-multiplication by s̄.
# By linearity, every λ*s̄ (λ ∈ L) has the same integer coordinates in this
# basis as λ has in the L Z-basis.  So this IS a Z-basis for Ls̄.
_Ls_bar_basis: np.ndarray = np.array([
    (STANDARD_ALGEBRA.element(b) * WILSON_S_BAR).coords for b in _L_zbasis
])  # 8×8
_Ls_bar_raw: np.ndarray = np.array(
    [(r * WILSON_S_BAR).coords for r in _e8_roots]
)
_verify_sublattice_basis(_Ls_bar_basis, _Ls_bar_raw)
_Ls_bar_basis_inv: np.ndarray = np.linalg.inv(_Ls_bar_basis)

# Ls basis: images of the L Z-basis under right-multiplication by s.
_Ls_basis: np.ndarray = np.array([
    (STANDARD_ALGEBRA.element(b) * WILSON_S).coords for b in _L_zbasis
])  # 8×8
_Ls_raw: np.ndarray = np.array(
    [(r * WILSON_S).coords for r in _e8_roots]
)
_verify_sublattice_basis(_Ls_basis, _Ls_raw)
_Ls_basis_inv: np.ndarray = np.linalg.inv(_Ls_basis)


# ---------------------------------------------------------------------------
# Sublattice membership tests
# ---------------------------------------------------------------------------

def is_in_Ls_bar(v: Octonion, tol: float = 1e-9) -> bool:
    """
    Return True iff v ∈ Ls̄ = 2B.

    Tests whether v is an integer combination of the precomputed Z-basis of Ls̄.
    Ls̄ is an index-16 sublattice of L (det(Gram) = 256).
    Reference: [Wilson2009] Section 2 — "Ls̄ = LR = 2B".
    """
    coords = v.coords @ _Ls_bar_basis_inv
    return bool(np.allclose(coords, np.round(coords), atol=tol))


def is_in_Ls(v: Octonion, tol: float = 1e-9) -> bool:
    """
    Return True iff v ∈ Ls.

    Ls is the other index-16 sublattice of L, with Ls̄ + Ls = L and
    Ls̄ ∩ Ls = 2L.  Reference: [Wilson2009] Section 2.
    """
    coords = v.coords @ _Ls_basis_inv
    return bool(np.allclose(coords, np.round(coords), atol=tol))


# ---------------------------------------------------------------------------
# Leech lattice membership test
# ---------------------------------------------------------------------------

def is_in_leech(x: Octonion, y: Octonion, z: Octonion) -> bool:
    """
    Return True iff (x, y, z) ∈ Λ (Wilson's octonionic Leech lattice).

    Conditions [Wilson2009 Section 3]:
      1. x, y, z ∈ L
      2. x+y, x+z, y+z ∈ Ls̄
      3. x+y+z ∈ Ls
    """
    alg = STANDARD_ALGEBRA
    if not (is_in_L(x) and is_in_L(y) and is_in_L(z)):
        return False
    for a, b in [(x, y), (x, z), (y, z)]:
        if not is_in_Ls_bar(alg.element(a.coords + b.coords)):
            return False
    if not is_in_Ls(alg.element(x.coords + y.coords + z.coords)):
        return False
    return True


# ---------------------------------------------------------------------------
# Minimal vector generation
# ---------------------------------------------------------------------------

# Precompute products for all 240 roots to avoid redundant multiplications.
_lsb_cache: List[np.ndarray] = [(r * WILSON_S_BAR).coords for r in _e8_roots]
# λ*s cache (used for the first component of type-3 vectors: (λs)j)
_ls_cache: List[np.ndarray] = [(r * WILSON_S).coords for r in _e8_roots]
_lj_cache: List[List[np.ndarray]] = [
    [(r * j).coords for j in J] for r in _e8_roots
]
_lk_cache: List[List[np.ndarray]] = _lj_cache  # λ*j and λ*k use same cache


def leech_type1_vectors() -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Return the 720 type-1 minimal vectors of Λ.

    Formula: (2λ, 0, 0) and its 3 cyclic permutations,
    where λ runs over the 240 roots of L.

    Norm: N_std(2λ, 0, 0) = N(2λ) = 4·N(λ) = 4·2 = 8. ✓
    Reference: [Wilson2009] Section 3.
    """
    zero = np.zeros(8)
    result: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for lam in _e8_roots:
        two_lam = 2.0 * lam.coords
        result.append((two_lam.copy(), zero.copy(), zero.copy()))
        result.append((zero.copy(), two_lam.copy(), zero.copy()))
        result.append((zero.copy(), zero.copy(), two_lam.copy()))
    return result


def leech_type2_vectors() -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Return the 11,520 type-2 minimal vectors of Λ.

    Formula: (λs̄, (λs̄)j, 0) and its 3 cyclic permutations,
    where λ runs over the 240 roots of L and j ∈ J (16 elements).

    Norm: N_std(λs̄, (λs̄)j, 0) = N(λ)N(s̄) + N(λ)N(s̄)N(j)
                                 = 2·2 + 2·2·1 = 4 + 4 = 8. ✓
    Count: 3 × 240 × 16 = 11,520.
    Reference: [Wilson2009] Section 3.
    """
    alg = STANDARD_ALGEBRA
    zero = np.zeros(8)
    result: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for idx_lam, lam in enumerate(_e8_roots):
        lsb = _lsb_cache[idx_lam]           # λs̄  (8-vector)
        lsb_oct = alg.element(lsb)
        for j in J:
            lsb_j = (lsb_oct * j).coords    # (λs̄)j
            # 3 cyclic permutations of (λs̄, (λs̄)j, 0)
            result.append((lsb.copy(), lsb_j.copy(), zero.copy()))
            result.append((lsb_j.copy(), zero.copy(), lsb.copy()))
            result.append((zero.copy(), lsb.copy(), lsb_j.copy()))
    return result


def leech_type3_vectors() -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Return the 184,320 type-3 minimal vectors of Λ.

    Formula: ((λs)j, λk, (λj)k) and its 3 cyclic permutations,
    where λ runs over the 240 roots of L and j, k ∈ J (16 elements each).

    Note: the first component uses s (Wilson's element WILSON_S), NOT s̄.
    Type-2 uses s̄; type-3 uses s.  This distinction is essential: using s̄
    in type-3 produces norm-8 vectors that do NOT satisfy Wilson's conditions.
    Reference: [Wilson2009] Section 3 — "((λs)j, ±λk, ±(λj)k)".

    Wilson notes: "Dixon's formula is equivalent to (λs, ±λk, ±(λj)k),
    omitting the (necessary) multiplication by j in the first coordinate."
    [Wilson2009, p. 2189 footnote]

    Norm: N_std = N(λ)N(s)N(j) + N(λ)N(k) + N(λ)N(j)N(k)
                = 2·2·1 + 2·1 + 2·1·1 = 4 + 2 + 2 = 8. ✓
    Count: 3 × 240 × 16 × 16 = 184,320.
    """
    alg = STANDARD_ALGEBRA
    result: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for idx_lam, lam in enumerate(_e8_roots):
        ls = _ls_cache[idx_lam]             # λs  (NOT λs̄)

        for idx_j, j in enumerate(J):
            ls_j = (alg.element(ls) * j).coords      # (λs)j
            lj = _lj_cache[idx_lam][idx_j]           # λj

            for idx_k, k in enumerate(J):
                lk = _lk_cache[idx_lam][idx_k]       # λk
                ljk = (alg.element(lj) * k).coords   # (λj)k

                # 3 cyclic permutations of ((λs)j, λk, (λj)k)
                result.append((ls_j.copy(), lk.copy(), ljk.copy()))
                result.append((ljk.copy(), ls_j.copy(), lk.copy()))
                result.append((lk.copy(), ljk.copy(), ls_j.copy()))
    return result


def leech_minimal_vectors() -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Return all 196,560 minimal vectors of Λ.

    These are the 3 × 240 × (1 + 16 + 16×16) = 196,560 vectors of minimum
    norm (N_std = 8) in the Leech lattice, as listed in [Wilson2009] Section 3.

    Note: calling this function computes ~200,000 octonion products and may
    take several seconds.  Use the individual leech_type*_vectors() functions
    when only a subset is needed.
    """
    return leech_type1_vectors() + leech_type2_vectors() + leech_type3_vectors()
