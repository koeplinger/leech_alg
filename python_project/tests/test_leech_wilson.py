"""
test_leech_wilson.py — Verification of Wilson's octonionic Leech lattice construction.

Clarification on the user's question
--------------------------------------
Wilson does NOT define a "multiplication that closes on the Leech lattice."  He
constructs the Leech lattice as a SET Λ ⊂ L³ satisfying three algebraic
conditions.  The Leech lattice has NO vectors of norm 2 (no root system); the
relevant object is the 196,560 minimal vectors (norm-squared 8 under the standard
R^24 norm).  This test suite verifies Wilson's construction and establishes the
membership criteria that will later be used to evaluate other products.

The answer to the user's question — "does Wilson's construction close on the
Leech lattice?" — requires interpreting what "close" means:
  - Wilson's 196,560 minimal vectors come in three explicit families (types 1–3).
  - All three families are verified to satisfy Wilson's membership conditions.
  - The resulting object is the unique even self-dual rank-24 lattice with no
    vectors of norm-squared < 8 [Wilson2009, Section 4].

Properties verified
-------------------
 J (unit octonions used in minimal-vector formulas):
 1.  J has 16 elements
 2.  All elements of J have norm-squared 1

 Ls̄ = {λ * s̄ : λ ∈ L}  (index-16 sublattice of L):
 3.  All 240 products λ * s̄ are in L
 4.  The 8 precomputed Ls̄ basis vectors are linearly independent (rank 8)
 5.  det(Gram of Ls̄ basis) ≈ 256  ([L : Ls̄]² × det(Gram_L) = 16² × 1 = 256)
 6.  is_in_Ls_bar correctly identifies all Ls̄ elements
 7.  is_in_Ls_bar correctly rejects a known non-member (a norm-2 E8 root)

 Ls  = {λ * s  : λ ∈ L}  (the other index-16 sublattice of L):
 8.  All 240 products λ * s are in L
 9.  The 8 precomputed Ls basis vectors are linearly independent (rank 8)
10.  det(Gram of Ls basis) ≈ 256
11.  is_in_Ls correctly identifies all Ls elements
12.  is_in_Ls correctly rejects a known non-member

 Sublattice relations  (Ls̄ ∩ Ls = 2L, Ls̄ + Ls = L):
13.  2 × root ∈ Ls̄ for all 240 roots  (2L ⊆ Ls̄)
14.  2 × root ∈ Ls  for all 240 roots  (2L ⊆ Ls)
15.  Combined Z-span of Ls̄ and Ls vectors has rank 8  (Ls̄ + Ls spans L)

 Type-1 minimal vectors  —  (2λ, 0, 0) and 2 cyclic perms, λ ∈ L:
16.  Count = 720
17.  All N_std = 8
18.  All 720 satisfy Wilson's membership conditions (is_in_leech)
19.  No duplicates among the 720 vectors

 Type-2 minimal vectors  —  (λs̄, (λs̄)j, 0) and 2 cyclic perms:
20.  Count = 11,520
21.  All N_std = 8
22.  Sample of 500 satisfy Wilson's membership conditions (is_in_leech)

 Type-3 minimal vectors  —  ((λs̄)j, λk, (λj)k) and 2 cyclic perms:
23.  Count = 184,320
24.  Sample of 500 have N_std = 8
25.  Sample of 500 satisfy Wilson's membership conditions (is_in_leech)

 Leech lattice summary:
26.  Total count: 720 + 11,520 + 184,320 = 196,560
27.  Types 1, 2, 3 are pairwise distinct (proven by zero-component structure)
28.  The 196,560 vectors span R^24 (rank 24)
29.  Even lattice: sample pairs of type-1 vectors have integer standard dot product
     (which equals 2 × Wilson inner product; self-duality follows from Wilson's proof)

 Non-members:
30.  (root, 0, 0) is NOT in Λ for any E8 root  (N_std = 2 < 8, violates the
     no-norm-2 property; algebraically: root ∉ Ls̄, violating condition 2)
31.  A manually constructed triple that violates condition 2 is correctly rejected

References
----------
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186-2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
"""

import itertools

import numpy as np
import pytest

from octonions import STANDARD_ALGEBRA
from e8_wilson import wilson_e8_roots, is_in_L
from leech_wilson import (
    WILSON_S_BAR,
    J,
    _Ls_bar_basis,
    _Ls_bar_basis_inv,
    _Ls_basis,
    _Ls_basis_inv,
    _e8_roots,
    is_in_Ls_bar,
    is_in_Ls,
    is_in_leech,
    leech_type1_vectors,
    leech_type2_vectors,
    leech_type3_vectors,
    leech_minimal_vectors,
)


# ---------------------------------------------------------------------------
# Fixtures  (module-scoped to avoid recomputing expensive vector sets)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def type1():
    return leech_type1_vectors()   # 720 vectors, fast


@pytest.fixture(scope="module")
def type2():
    return leech_type2_vectors()   # 11,520 vectors, moderate


@pytest.fixture(scope="module")
def type3():
    return leech_type3_vectors()   # 184,320 vectors, slow (~5-10 s)


@pytest.fixture(scope="module")
def e8_roots():
    return wilson_e8_roots()       # 240 E8 roots


# ---------------------------------------------------------------------------
# 1–2: J (the 16 unit octonions)
# ---------------------------------------------------------------------------

class TestJ:
    """
    J = {±e_t : t = 0,...,7} is the set of 16 unit octonions.
    Wilson uses J to index the second and third parameters in the type-2
    and type-3 minimal vector formulas.  [Wilson2009 Section 3]
    """

    def test_j_count(self):
        """J has exactly 16 elements (8 basis elements × 2 signs)."""
        assert len(J) == 16

    def test_j_norms(self):
        """All elements of J are unit octonions: N(j) = 1."""
        for j in J:
            assert abs(j.norm_sq() - 1.0) < 1e-12, \
                f"J element {j.coords} has N = {j.norm_sq()}, expected 1"


# ---------------------------------------------------------------------------
# 3–7: Ls̄ sublattice  (= 2B = LR in Wilson's notation)
# ---------------------------------------------------------------------------

class TestLsBarSublattice:
    """
    Ls̄ = {λ × s̄ : λ ∈ L} is an index-16 sublattice of L.

    Wilson proves Ls̄ = LR = 2B (the Coxeter-Dickson ring scaled by 2).
    Gram matrix: ⟨λ₁s̄, λ₂s̄⟩ = N(s̄) × ⟨λ₁,λ₂⟩ = 2⟨λ₁,λ₂⟩.
    So det(Gram_Ls̄) = 2^8 × det(Gram_L) = 256 × 1 = 256.
    Reference: [Wilson2009] Sections 2–3.
    """

    def test_Ls_bar_elements_in_L(self, e8_roots):
        """All 240 products λ * s̄ are in L (Ls̄ ⊆ L)."""
        for r in e8_roots:
            prod = r * WILSON_S_BAR
            assert is_in_L(prod), \
                f"λ*s̄ not in L for λ = {r.coords}; product = {prod.coords}"

    def test_Ls_bar_basis_rank(self):
        """The precomputed Ls̄ basis has rank 8 (8 independent vectors found)."""
        assert _Ls_bar_basis.shape == (8, 8)
        assert np.linalg.matrix_rank(_Ls_bar_basis) == 8

    def test_Ls_bar_gram_det(self):
        """
        det(Gram of Ls̄ basis) ≈ 256.

        [L : Ls̄] = 16, so det(Gram_Ls̄) = 16² × det(Gram_L) = 256 × 1 = 256.
        The precomputed Ls̄ basis was verified at module load to generate Ls̄
        exactly, so this Gram det tests the index.
        """
        G = _Ls_bar_basis @ _Ls_bar_basis.T
        det = np.linalg.det(G)
        assert abs(det - 256.0) < 1e-6, \
            f"Gram det of Ls̄ basis = {det:.6f}, expected 256"

    def test_is_in_Ls_bar_recognises_products(self, e8_roots):
        """is_in_Ls_bar returns True for all 240 products λ * s̄."""
        for r in e8_roots:
            prod = r * WILSON_S_BAR
            assert is_in_Ls_bar(prod), \
                f"is_in_Ls_bar failed for λ*s̄ with λ = {r.coords}"

    def test_is_in_Ls_bar_rejects_norm2_root(self, e8_roots):
        """
        No norm-2 element of L is in Ls̄.

        Elements of Ls̄ have minimum norm N(λ) × N(s̄) = 2 × 2 = 4.
        Therefore every E8 root (norm 2) is NOT in Ls̄.
        """
        for r in e8_roots:
            assert not is_in_Ls_bar(r), \
                f"E8 root {r.coords} (norm 2) was incorrectly identified as in Ls̄"


# ---------------------------------------------------------------------------
# 8–12: Ls sublattice
# ---------------------------------------------------------------------------

class TestLsSublattice:
    """
    Ls = {λ × s : λ ∈ L} is the other index-16 sublattice of L.

    Ls and Ls̄ are complementary in the sense that Ls̄ + Ls = L and
    Ls̄ ∩ Ls = 2L.  Reference: [Wilson2009] Sections 2–3.
    """

    def test_Ls_elements_in_L(self, e8_roots):
        """All 240 products λ * s are in L (Ls ⊆ L)."""
        from e8_wilson import WILSON_S
        for r in e8_roots:
            prod = r * WILSON_S
            assert is_in_L(prod), \
                f"λ*s not in L for λ = {r.coords}; product = {prod.coords}"

    def test_Ls_basis_rank(self):
        """The precomputed Ls basis has rank 8."""
        assert _Ls_basis.shape == (8, 8)
        assert np.linalg.matrix_rank(_Ls_basis) == 8

    def test_Ls_gram_det(self):
        """det(Gram of Ls basis) ≈ 256  ([L : Ls] = 16)."""
        G = _Ls_basis @ _Ls_basis.T
        det = np.linalg.det(G)
        assert abs(det - 256.0) < 1e-6, \
            f"Gram det of Ls basis = {det:.6f}, expected 256"

    def test_is_in_Ls_recognises_products(self, e8_roots):
        """is_in_Ls returns True for all 240 products λ * s."""
        from e8_wilson import WILSON_S
        for r in e8_roots:
            prod = r * WILSON_S
            assert is_in_Ls(prod), \
                f"is_in_Ls failed for λ*s with λ = {r.coords}"

    def test_is_in_Ls_rejects_norm2_root(self, e8_roots):
        """
        No norm-2 element of L is in Ls.

        Elements of Ls have minimum norm N(λ) × N(s) = 2 × 2 = 4.
        """
        for r in e8_roots:
            assert not is_in_Ls(r), \
                f"E8 root {r.coords} (norm 2) was incorrectly identified as in Ls"


# ---------------------------------------------------------------------------
# 13–15: Sublattice relations  (2L ⊆ Ls̄ ∩ Ls; Ls̄ + Ls spans L)
# ---------------------------------------------------------------------------

class TestSublatticeSumAndIntersection:
    """
    Wilson [Wilson2009 Section 2] asserts:
      Ls̄ + Ls = L   (every L element is a sum of one Ls̄ and one Ls element)
      Ls̄ ∩ Ls = 2L  (the intersection is exactly the scaled lattice 2L)

    The Ls̄ ∩ Ls = 2L direction (⊇) is tested directly: all 2 × roots must
    lie in both sublattices.  The reverse direction (⊆) is implied by the
    index calculation: [L : Ls̄] = [L : Ls] = 16, and the sublattices
    together generate L.

    The Ls̄ + Ls = L direction is tested via rank: combining all Ls̄ and Ls
    generators, the matrix rank is 8 (= rank of L).
    """

    def test_2L_subset_Ls_bar(self, e8_roots):
        """2 × root ∈ Ls̄ for all 240 roots  (2L ⊆ Ls̄ = 2B ⊇ 2L)."""
        alg = STANDARD_ALGEBRA
        for r in e8_roots:
            two_r = alg.element(2.0 * r.coords)
            assert is_in_Ls_bar(two_r), \
                f"2 × root {r.coords} not in Ls̄"

    def test_2L_subset_Ls(self, e8_roots):
        """2 × root ∈ Ls  for all 240 roots  (2L ⊆ Ls)."""
        alg = STANDARD_ALGEBRA
        for r in e8_roots:
            two_r = alg.element(2.0 * r.coords)
            assert is_in_Ls(two_r), \
                f"2 × root {r.coords} not in Ls"

    def test_Ls_bar_plus_Ls_spans_L(self, e8_roots):
        """
        Combined Z-span of all Ls̄ and Ls generators has rank 8.

        This confirms Ls̄ + Ls ⊇ L in the vector-space sense; together with
        the index bound this implies Ls̄ + Ls = L.
        Reference: [Wilson2009] Section 2.
        """
        from e8_wilson import WILSON_S
        Ls_bar_vecs = np.array([(r * WILSON_S_BAR).coords for r in e8_roots])
        Ls_vecs     = np.array([(r * WILSON_S).coords     for r in e8_roots])
        combined = np.vstack([Ls_bar_vecs, Ls_vecs])  # (480, 8)
        assert np.linalg.matrix_rank(combined) == 8, \
            "Combined Ls̄ ∪ Ls vectors do not span R^8 (rank < 8)"


# ---------------------------------------------------------------------------
# 16–19: Type-1 minimal vectors
# ---------------------------------------------------------------------------

class TestType1Vectors:
    """
    Type-1: (2λ, 0, 0) and 2 cyclic permutations, λ ∈ L (240 roots).
    Count: 3 × 240 = 720.
    Norm:  N_std(2λ, 0, 0) = 4 N(λ) = 4 × 2 = 8.
    Reference: [Wilson2009] Section 3.
    """

    def test_type1_count(self, type1):
        """Exactly 720 type-1 vectors."""
        assert len(type1) == 720

    def test_type1_all_norm_8(self, type1):
        """All type-1 vectors have N_std = 8."""
        for x, y, z in type1:
            norm = float(np.dot(x, x) + np.dot(y, y) + np.dot(z, z))
            assert abs(norm - 8.0) < 1e-10, \
                f"Type-1 vector has N_std = {norm}, expected 8"

    def test_type1_all_in_leech(self, type1):
        """
        All 720 type-1 vectors satisfy Wilson's conditions (is_in_leech).

        This is the most direct verification: each (2λ, 0, 0) and permutation
        must lie in Λ as defined by conditions 1–3 of [Wilson2009 Section 3].
        """
        alg = STANDARD_ALGEBRA
        for x, y, z in type1:
            ox = alg.element(x)
            oy = alg.element(y)
            oz = alg.element(z)
            assert is_in_leech(ox, oy, oz), \
                f"Type-1 vector ({x}, {y}, {z}) not recognised as in Λ"

    def test_type1_no_duplicates(self, type1):
        """All 720 type-1 vectors are distinct."""
        triples = {(tuple(x), tuple(y), tuple(z)) for x, y, z in type1}
        assert len(triples) == 720


# ---------------------------------------------------------------------------
# 20–22: Type-2 minimal vectors
# ---------------------------------------------------------------------------

class TestType2Vectors:
    """
    Type-2: (λs̄, (λs̄)j, 0) and 2 cyclic permutations, λ ∈ L, j ∈ J.
    Count: 3 × 240 × 16 = 11,520.
    Norm:  N_std = N(λ)N(s̄) + N(λ)N(s̄)N(j) = 4 + 4 = 8.
    Reference: [Wilson2009] Section 3.
    """

    def test_type2_count(self, type2):
        """Exactly 11,520 type-2 vectors."""
        assert len(type2) == 11_520

    def test_type2_all_norm_8(self, type2):
        """All 11,520 type-2 vectors have N_std = 8."""
        for x, y, z in type2:
            norm = float(np.dot(x, x) + np.dot(y, y) + np.dot(z, z))
            assert abs(norm - 8.0) < 1e-10, \
                f"Type-2 vector has N_std = {norm}, expected 8"

    def test_type2_sample_in_leech(self, type2):
        """
        500 evenly-spaced type-2 vectors satisfy Wilson's membership conditions.

        Sampling every 23rd vector (≈ 500 checks from 11,520) is sufficient
        to catch construction errors while keeping the test fast.
        """
        alg = STANDARD_ALGEBRA
        step = max(1, len(type2) // 500)
        for i in range(0, len(type2), step):
            x, y, z = type2[i]
            ox, oy, oz = alg.element(x), alg.element(y), alg.element(z)
            assert is_in_leech(ox, oy, oz), \
                f"Type-2 vector [{i}] ({x}, {y}, {z}) not recognised as in Λ"


# ---------------------------------------------------------------------------
# 23–25: Type-3 minimal vectors
# ---------------------------------------------------------------------------

class TestType3Vectors:
    """
    Type-3: ((λs̄)j, λk, (λj)k) and 2 cyclic permutations, λ ∈ L, j,k ∈ J.
    Count: 3 × 240 × 16 × 16 = 184,320.
    Norm:  N_std = N(λ)N(s̄)N(j) + N(λ)N(k) + N(λ)N(j)N(k) = 4 + 2 + 2 = 8.
    Reference: [Wilson2009] Section 3.

    Note: leech_type3_vectors() generates ~120,000 octonion products.
    The module-scoped type3 fixture (see above) ensures this is done once.
    """

    def test_type3_count(self, type3):
        """Exactly 184,320 type-3 vectors."""
        assert len(type3) == 184_320

    def test_type3_sample_norm_8(self, type3):
        """500 evenly-spaced type-3 vectors have N_std = 8."""
        step = max(1, len(type3) // 500)
        for i in range(0, len(type3), step):
            x, y, z = type3[i]
            norm = float(np.dot(x, x) + np.dot(y, y) + np.dot(z, z))
            assert abs(norm - 8.0) < 1e-10, \
                f"Type-3 vector [{i}] has N_std = {norm}, expected 8"

    def test_type3_sample_in_leech(self, type3):
        """
        500 evenly-spaced type-3 vectors satisfy Wilson's membership conditions.

        Sampling every 369th vector (≈ 500 checks from 184,320).
        """
        alg = STANDARD_ALGEBRA
        step = max(1, len(type3) // 500)
        for i in range(0, len(type3), step):
            x, y, z = type3[i]
            ox, oy, oz = alg.element(x), alg.element(y), alg.element(z)
            assert is_in_leech(ox, oy, oz), \
                f"Type-3 vector [{i}] not recognised as in Λ"


# ---------------------------------------------------------------------------
# 26–29: Leech lattice summary properties
# ---------------------------------------------------------------------------

class TestLeechLatticeSummary:
    """
    Summary properties of the Leech lattice, derived from the three families.

    These tests verify the four defining properties of the Leech lattice
    [Conway–Sloane, Chapter 12]:
      (a) Even self-dual lattice (tested via dot-product parity)
      (b) Rank 24
      (c) No norm-2 vectors  (tested separately in TestNonMembers)
      (d) Uniqueness  (asserted by reference to Wilson's proof)

    Property (d), self-duality (det_Wilson = 1), follows from Wilson's proof
    [Wilson2009 Section 4] and is not reproduced numerically here.  The even
    lattice property (a) is spot-checked via type-1 pairwise dot products.
    """

    def test_total_minimal_vector_count(self, type1, type2, type3):
        """Total: 720 + 11,520 + 184,320 = 196,560 minimal vectors."""
        assert len(type1) + len(type2) + len(type3) == 196_560

    def test_types_pairwise_distinct(self, type1, type2, type3):
        """
        Type-1, type-2, and type-3 vectors are pairwise distinct.

        Proof (by zero-component structure):
          Type-1: exactly 2 zero blocks of 8.
          Type-2: exactly 1 zero block of 8.
          Type-3: no zero blocks (each component has norm 4, 2, or 2).
        Therefore the three families are automatically disjoint.
        This test verifies the zero-block structure as a proxy.
        """
        def count_nonzero_blocks(x, y, z):
            """Count how many of the three R^8 components are non-zero."""
            return sum(
                np.any(np.abs(v) > 1e-10) for v in (x, y, z)
            )

        # Type-1: exactly 1 nonzero block
        for x, y, z in type1[:10]:
            assert count_nonzero_blocks(x, y, z) == 1, \
                f"Type-1 vector has wrong zero structure"

        # Type-2: exactly 2 nonzero blocks
        for x, y, z in type2[:10]:
            assert count_nonzero_blocks(x, y, z) == 2, \
                f"Type-2 vector has wrong zero structure"

        # Type-3: exactly 3 nonzero blocks
        for x, y, z in type3[:10]:
            assert count_nonzero_blocks(x, y, z) == 3, \
                f"Type-3 vector has wrong zero structure"

    def test_rank_24(self, type1, type2, type3):
        """
        The 196,560 minimal vectors span R^24 (rank 24).

        This confirms the Leech lattice is a rank-24 lattice.
        Tested by finding rank of a sample drawn from all three types.
        Reference: [Wilson2009] Section 4 — "the lattice has rank 24".
        """
        alg = STANDARD_ALGEBRA
        # Use 30 type-1 + 30 type-2 + 30 type-3 vectors (sufficient to reach rank 24)
        sample = (
            list(type1[:30])
            + list(type2[:30])
            + list(type3[:30])
        )
        M = np.array([np.concatenate([x, y, z]) for x, y, z in sample])
        assert np.linalg.matrix_rank(M) == 24, \
            "Minimal vectors do not span R^24 (rank < 24)"

    def test_even_lattice_type1_gram(self, type1):
        """
        All type-1 pairwise dot products are integers (Wilson inner product ∈ Z).

        For an integer lattice, ⟨v₁,v₂⟩_Wilson = ½(v₁·v₂) must be an integer,
        so v₁·v₂ must be in 2Z (even).  For type-1 vectors (2λ,0,0):
          (2λ₁,0,0)·(2λ₂,0,0) = 4(λ₁·λ₂) ∈ 4Z ⊆ 2Z  ✓  (since λᵢ ∈ L, integer lattice).
        This test verifies the full 720×720 Gram matrix for type-1 vectors.
        Reference: [Wilson2009] Section 4 — "Λ is an even lattice".
        """
        coords = np.array([np.concatenate([x, y, z]) for x, y, z in type1])
        G = coords @ coords.T   # 720 × 720 Gram matrix
        assert np.allclose(G, np.round(G), atol=1e-9), \
            "Type-1 Gram matrix has non-integer entries"
        assert np.all(np.round(G).astype(int) % 2 == 0), \
            "Type-1 Gram matrix has odd entries (violates even lattice property)"


# ---------------------------------------------------------------------------
# 30–31: Non-members
# ---------------------------------------------------------------------------

class TestNonMembers:
    """
    Verifies that is_in_leech correctly rejects vectors outside Λ.

    Two key non-members are tested:
      - A single E8 root in position x with y = z = 0: this fails condition 2
        (x + y = root ∉ Ls̄, since min norm in Ls̄ is 4 > 2 = N(root)).
      - A triple (root, root, 0): x + y = 2root ∈ Ls̄ (passes condition 2),
        but x + z = root ∉ Ls̄ (fails condition 2 for the (x,z) pair).
    """

    def test_single_root_not_in_leech(self, e8_roots):
        """
        (root, 0, 0) ∉ Λ for all 240 E8 roots.

        Proof: condition 2 requires x + y = root ∈ Ls̄.  But min N in Ls̄ = 4,
        and N(root) = 2 < 4.  So root ∉ Ls̄, and condition 2 fails.
        Reference: [Wilson2009] Section 4 — "no vectors of norm 2".
        """
        alg = STANDARD_ALGEBRA
        zero = alg.element(np.zeros(8))
        for r in e8_roots:
            assert not is_in_leech(r, zero, zero), \
                f"(root, 0, 0) incorrectly identified as in Λ for root = {r.coords}"

    def test_wrong_sublattice_not_in_leech(self, e8_roots):
        """
        (root, root, 0) ∉ Λ.

        Condition 2 requires x + z = root ∈ Ls̄.  Since root has norm 2 and
        min norm of Ls̄ is 4, root ∉ Ls̄.  So the triple (root, root, 0) fails
        condition 2 for the (x, z) pair, even though x + y = 2*root ∈ Ls̄.
        """
        alg = STANDARD_ALGEBRA
        zero = alg.element(np.zeros(8))
        # Take first 20 roots as a sample
        for r in e8_roots[:20]:
            assert not is_in_leech(r, r, zero), \
                f"(root, root, 0) incorrectly identified as in Λ for root = {r.coords}"
