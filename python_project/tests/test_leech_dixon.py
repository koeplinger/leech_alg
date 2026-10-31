"""
test_leech_dixon.py — Verification of Dixon's octonionic Leech lattice construction.

Created: 2026-04-07.  Updated: 2026-10-31.

This suite verifies Dixon's three families of 196,560 minimal vectors
[Dixon2010 Section 6] and compares them against Wilson's construction
[Wilson2009 Section 3].

Central finding: Dixon's 196,560 vectors are exactly the componentwise
conjugates (x̄, ȳ, z̄) of Wilson's 196,560 minimal vectors, family by family.
Dixon derives the largest family (eq. (20)) from Wilson's in this way
[Dixon2010 Section 6].  Componentwise conjugation is an isometry of R^24, so
Dixon's vectors are the minimal vectors of a copy of the Leech lattice, the
conjugate of Wilson's Λ; their inner products are checked here directly as
well.  The vectors themselves mostly lie outside Wilson's Λ: Dixon's A^odd
uses the half-integer roots with an EVEN number of minus signs (A₃), Wilson's
E8 roots those with an ODD number, and conjugation exchanges the two.

Properties verified
-------------------
 ℓ₀ = ½(e₀+e₁+...+e₇):
 1.  All 8 coordinates of ℓ₀ equal 0.5
 2.  N(ℓ₀) = 2  (ℓ₀ is a norm-2 element like s)
 3.  ℓ₀ ∈ A₃  (all-positive element, 0 minus signs = EVEN = A₃ criterion)

 Type-1 vectors  —  (2P, 0, 0) and 2 cyclic perms, P ∈ A^odd:
 4.  Count = 720
 5.  All N_std = 8
 6.  Exactly 1 nonzero block per triple (first 30 vectors)
 7.  Type-1 vectors with P ∈ A₁ (112 elements): all 336 satisfy is_in_leech
     (A₁ = Wilson's type-1 roots; 2P ∈ 2L, so all conditions are met)
 8.  Type-1 vectors with P ∈ A₃ (128 elements): all 384 FAIL is_in_leech
     (A₃ has even # of minus signs; 2P ∉ 2L = Ls̄ ∩ Ls: 2P lies in Ls̄ but
     not in Ls, violating Wilson's condition 3.  Their conjugates lie in Λ.)
 9.  No duplicates among the 720 vectors

 Type-2 vectors  —  (2Q, ±2e_aQ, 0) and 2 cyclic perms, Q ∈ A^even:
10.  Count = 11,520
11.  All N_std = 8
12.  Exactly 2 nonzero blocks per triple (first 30 vectors)
13.  Their componentwise conjugates are exactly Wilson's 11,520 type-2 vectors

 Type-3 vectors  —  (P, ±e_aP, ±e_a(ℓ₀(e_cP))) and 2 cyclic perms:
14.  Count = 184,320
15.  All N_std = 8
16.  Exactly 3 nonzero blocks per triple (first 30 vectors)
17.  Their componentwise conjugates are exactly Wilson's 184,320 type-3 vectors

 Dixon Leech summary:
18.  Total: 720 + 11,520 + 184,320 = 196,560 minimal vectors (no duplicates)
19.  The 196,560 vectors span R^24 (rank 24; verified using type-1 + type-2)

 Comparison with Wilson's minimal-vector shell, and inner products:
20.  The inner products of 3,000 sampled vectors (and of the pair type-3
     [0], [411]) lie in {0, ±2, ±4, ±8}, ±8 only for v = ±w, as for the
     minimal vectors of the √2-scaled Leech lattice
21.  The componentwise conjugates of all 196,560 vectors are exactly
     Wilson's 196,560 minimal vectors

References
----------
[Dixon2010]   G.M. Dixon, "Integral Octonions, Octonion XY-Product, and the
              Leech Lattice", preprint, 2010.
[Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
              Journal of Algebra 322 (2009) 2186-2190.
              DOI: 10.1016/j.jalgebra.2009.03.021
"""

import numpy as np
import pytest

from octonions import STANDARD_ALGEBRA
from e8_dixon import dixon_a1, dixon_a3, dixon_a_odd, dixon_xi_even
from leech_wilson import (
    is_in_leech,
    leech_type1_vectors,
    leech_type2_vectors,
    leech_type3_vectors,
)
from leech_dixon import (
    ELLO_0,
    dixon_leech_type1,
    dixon_leech_type2,
    dixon_leech_type3,
    dixon_leech_minimal_vectors,
)


# ---------------------------------------------------------------------------
# Fixtures  (module-scoped to avoid recomputing expensive vector sets)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def dix_type1():
    return dixon_leech_type1()    # 720 vectors, fast


@pytest.fixture(scope="module")
def dix_type2():
    return dixon_leech_type2()    # 11,520 vectors, moderate


@pytest.fixture(scope="module")
def dix_type3():
    return dixon_leech_type3()    # 184,320 vectors, the slowest fixture


def to_keyset(vecs):
    """Integer keys (coordinates × 2) of triples, for exact set comparison."""
    return {
        tuple(np.round(np.concatenate([x, y, z]) * 2).astype(int))
        for x, y, z in vecs
    }


def conjugate_keyset(vecs):
    """Integer keys (coordinates × 2) of the componentwise conjugates (x̄, ȳ, z̄)."""
    sign = np.array([1.0] + [-1.0] * 7)
    return {
        tuple(np.round(np.concatenate([sign * x, sign * y, sign * z]) * 2).astype(int))
        for x, y, z in vecs
    }


# ---------------------------------------------------------------------------
# 1–3: ℓ₀ properties
# ---------------------------------------------------------------------------

class TestEllo0:
    """
    ℓ₀ = ½(e₀+e₁+...+e₇) is Dixon's base element [Dixon2010 Section 3, eq. (6)].

    It is the all-positive element of A₃ (even number of minus signs = 0),
    and has norm N(ℓ₀) = 8 × (½)² = 2.  Dixon uses ℓ₀ to build the third
    component of type-3 minimal vectors.
    """

    def test_ello0_coords(self):
        """All 8 coordinates of ℓ₀ equal 0.5."""
        assert np.allclose(ELLO_0.coords, 0.5 * np.ones(8), atol=1e-12), \
            f"ELLO_0 coords = {ELLO_0.coords}, expected all 0.5"

    def test_ello0_norm(self):
        """N(ℓ₀) = 8 × (½)² = 2."""
        assert abs(ELLO_0.norm_sq() - 2.0) < 1e-12, \
            f"N(ℓ₀) = {ELLO_0.norm_sq()}, expected 2"

    def test_ello0_in_a3(self):
        """
        ℓ₀ is the all-positive element of A₃.

        A₃ = {½∑±e_a : even # of minus signs}.  ℓ₀ = ½∑e_a has 0 minus
        signs, which is even.  This test confirms ℓ₀ appears in dixon_a3().
        """
        a3 = dixon_a3()
        keys = {tuple(np.round(v.coords * 2).astype(int)) for v in a3}
        ello0_key = tuple(np.round(ELLO_0.coords * 2).astype(int))
        assert ello0_key in keys, "ℓ₀ is not found among the 128 A₃ elements"


# ---------------------------------------------------------------------------
# 4–9: Type-1 minimal vectors
# ---------------------------------------------------------------------------

class TestType1VectorsDixon:
    """
    Type-1: (2P, 0, 0) and 2 cyclic permutations, P ∈ A^odd = A₁ ∪ A₃.

    A^odd contains 240 norm-2 elements.  However, A₁ and A₃ differ in their
    half-integer component:
      - A₁ (112 elements):  ±e_a ± e_b, integer coordinates.
        These ARE Wilson's type-1 E8 roots.  2P ∈ 2L, so (2P,0,0) ∈ Λ.
      - A₃ (128 elements):  ½∑±e_a, EVEN # of minus signs.
        These are NOT Wilson's E8 roots.  P is a half-integer vector with an
        even coordinate sum, so P ∉ L and 2P ∉ 2L = Ls̄ ∩ Ls; 2P lies in Ls̄
        but not in Ls.  Therefore (2P,0,0) ∉ Λ.

    Wilson's E8 roots use the COMPLEMENTARY half-integer set, with an ODD
    number of minus signs.  Conjugation maps A^odd onto Wilson's E8 roots,
    which is why Dixon's shell is the conjugate of Wilson's (tests 13, 17, 21).
    Reference: [Dixon2010] Section 6, eq. (22).
    """

    def test_type1_count(self, dix_type1):
        """Exactly 720 type-1 vectors (240 × 3 cyclic perms)."""
        assert len(dix_type1) == 720

    def test_type1_all_norm_8(self, dix_type1):
        """All 720 type-1 vectors have N_std = 8."""
        for x, y, z in dix_type1:
            norm = float(np.dot(x, x) + np.dot(y, y) + np.dot(z, z))
            assert abs(norm - 8.0) < 1e-10, \
                f"Type-1 vector has N_std = {norm}, expected 8"

    def test_type1_zero_structure(self, dix_type1):
        """The first 30 type-1 vectors have exactly 1 nonzero block (2 zero blocks)."""
        for x, y, z in dix_type1[:30]:
            nonzero = sum(np.any(np.abs(v) > 1e-10) for v in (x, y, z))
            assert nonzero == 1, \
                f"Type-1 vector has {nonzero} nonzero blocks, expected 1"

    def test_type1_a1_subset_in_leech(self):
        """
        All 336 type-1 vectors with P ∈ A₁ satisfy Wilson's is_in_leech.

        A₁ = {±e_a ± e_b : a ≠ b} = Wilson's type-1 E8 roots.  Then
        2P ∈ 2L (even integer coordinates), and 2L ⊆ Ls̄ ∩ Ls, so all
        three Wilson conditions are met.
        Count: 112 × 3 = 336.
        """
        alg = STANDARD_ALGEBRA
        a1 = dixon_a1()
        zero = np.zeros(8)
        zero_oct = alg.element(zero)
        for P in a1:
            two_P = alg.element(2.0 * P.coords)
            for triple in [
                (two_P, zero_oct, zero_oct),
                (zero_oct, two_P, zero_oct),
                (zero_oct, zero_oct, two_P),
            ]:
                assert is_in_leech(*triple), (
                    f"Type-1 vector with P ∈ A₁ not in Λ: P = {P.coords}"
                )

    def test_type1_a3_subset_not_in_leech(self):
        """
        All 384 type-1 vectors with P ∈ A₃ FAIL Wilson's is_in_leech.

        A₃ = {½∑±e_a : even # of minus signs}.  These are NOT Wilson's E8
        roots (which require ODD # of minus signs).  P ∈ A₃ is a
        half-integer vector with an even coordinate sum, so P ∉ L and
        2P ∉ 2L = Ls̄ ∩ Ls.
        2P lies in Ls̄ but not in Ls, so Wilson's condition 3 (x+y+z = 2P ∈ Ls)
        fails.
        Count: 128 × 3 = 384.
        """
        alg = STANDARD_ALGEBRA
        a3 = dixon_a3()
        zero_oct = alg.element(np.zeros(8))
        fail_count = 0
        for P in a3:
            two_P = alg.element(2.0 * P.coords)
            for triple in [
                (two_P, zero_oct, zero_oct),
                (zero_oct, two_P, zero_oct),
                (zero_oct, zero_oct, two_P),
            ]:
                if not is_in_leech(*triple):
                    fail_count += 1
        assert fail_count == 384, (
            f"Expected all 384 A₃ type-1 vectors to fail is_in_leech, "
            f"but {384 - fail_count} passed unexpectedly"
        )

    def test_type1_no_duplicates(self, dix_type1):
        """All 720 type-1 vectors are distinct."""
        assert len(to_keyset(dix_type1)) == 720


# ---------------------------------------------------------------------------
# 10–13: Type-2 minimal vectors
# ---------------------------------------------------------------------------

class TestType2VectorsDixon:
    """
    Type-2: (2Q, ±2e_aQ, 0) and 2 cyclic perms, Q ∈ A^even, a ∈ {0,...,7}.

    A^even = Ξ^even (240 norm-1 elements: 16 unit basis elements + 224
    four-element half-integer combinations).  N_std = N(2Q)+N(±2e_aQ) = 4+4 = 8.
    Count: 240 × 8 × 2 × 3 = 11,520.
    Reference: [Dixon2010] Section 6, eq. (21).
    """

    def test_type2_count(self, dix_type2):
        """Exactly 11,520 type-2 vectors."""
        assert len(dix_type2) == 11_520

    def test_type2_all_norm_8(self, dix_type2):
        """All 11,520 type-2 vectors have N_std = 8."""
        for x, y, z in dix_type2:
            norm = float(np.dot(x, x) + np.dot(y, y) + np.dot(z, z))
            assert abs(norm - 8.0) < 1e-10, \
                f"Type-2 vector has N_std = {norm}, expected 8"

    def test_type2_zero_structure(self, dix_type2):
        """The first 30 type-2 vectors have exactly 2 nonzero blocks (1 zero block)."""
        for x, y, z in dix_type2[:30]:
            nonzero = sum(np.any(np.abs(v) > 1e-10) for v in (x, y, z))
            assert nonzero == 2, \
                f"Type-2 vector has {nonzero} nonzero blocks, expected 2"

    def test_type2_conjugates_are_wilson_type2(self, dix_type2):
        """
        The componentwise conjugates of the 11,520 type-2 vectors are exactly
        Wilson's 11,520 type-2 vectors.

        Conjugation reverses products: the conjugate of e_aQ is Q̄ē_a, so
        Dixon's left multiplication by e_a becomes Wilson's right
        multiplication by a unit.
        Reference: [Dixon2010] Section 6 eq. (21); [Wilson2009] Section 3.
        """
        assert conjugate_keyset(dix_type2) == to_keyset(leech_type2_vectors())


# ---------------------------------------------------------------------------
# 14–17: Type-3 minimal vectors
# ---------------------------------------------------------------------------

class TestType3VectorsDixon:
    """
    Type-3: (P, ±e_aP, ±e_a(ℓ₀(e_cP))) and 2 cyclic perms,
            P ∈ A^odd, a,c ∈ {0,...,7}, ± signs independent.

    Norm structure (2, 2, 4) per component.  The 3rd cyclic permutation
    has structure (4, 2, 2), matching Wilson's type-3.
    Count: 240 × 8 × 8 × 2 × 2 × 3 = 184,320.
    Reference: [Dixon2010] Section 6, eq. (20).
    """

    def test_type3_count(self, dix_type3):
        """Exactly 184,320 type-3 vectors."""
        assert len(dix_type3) == 184_320

    def test_type3_all_norm_8(self, dix_type3):
        """All 184,320 type-3 vectors have N_std = 8."""
        V = np.array([np.concatenate(t) for t in dix_type3])
        norms = np.einsum("ij,ij->i", V, V)
        assert np.allclose(norms, 8.0, atol=1e-10), \
            f"Type-3 norms range over [{norms.min()}, {norms.max()}], expected 8"

    def test_type3_zero_structure(self, dix_type3):
        """
        The first 30 type-3 vectors have exactly 3 nonzero blocks (no zero
        component).  The three components of (P, ±e_aP, ±e_a(ℓ₀(e_cP))) have
        norms (2, 2, 4) for all P, a, c, so every type-3 vector has three
        nonzero blocks.
        """
        for x, y, z in dix_type3[:30]:
            nonzero = sum(np.any(np.abs(v) > 1e-10) for v in (x, y, z))
            assert nonzero == 3, \
                f"Type-3 vector has {nonzero} nonzero blocks, expected 3"

    def test_type3_conjugates_are_wilson_type3(self, dix_type3):
        """
        The componentwise conjugates of the 184,320 type-3 vectors are exactly
        Wilson's 184,320 type-3 vectors.

        This is the content of Dixon's derivation of eq. (20) from Wilson's
        largest family: conjugate each component, shift cyclically, replace
        P by e_cP.  The nesting e_a(ℓ₀(e_cP)) matters; the octonion
        (e_aℓ₀)(e_cP) differs from it in general, and with it the set fails
        (inner products 1, 3, 5 and 6 occur).
        Reference: [Dixon2010] Section 6 eq. (20); [Wilson2009] Section 3.
        """
        assert conjugate_keyset(dix_type3) == to_keyset(leech_type3_vectors())


# ---------------------------------------------------------------------------
# 18–19: Dixon Leech summary
# ---------------------------------------------------------------------------

class TestDixonLeechSummary:
    """
    Summary properties of Dixon's 196,560 minimal vectors.

    Tests 18–19 verify the two most basic structural properties: the correct
    total count of distinct vectors and that the vectors span all of R^24
    (rank 24).
    """

    def test_total_count(self, dix_type1, dix_type2, dix_type3):
        """Total: 720 + 11,520 + 184,320 = 196,560 distinct minimal vectors."""
        assert len(dix_type1) + len(dix_type2) + len(dix_type3) == 196_560
        keys = to_keyset(dix_type1) | to_keyset(dix_type2) | to_keyset(dix_type3)
        assert len(keys) == 196_560

    def test_rank_24(self, dix_type1, dix_type2):
        """
        The 196,560 Dixon vectors span R^24 (rank 24).

        All 720 type-1 + all 11,520 type-2 = 12,240 vectors are used.
        A small uniform sample alone would fail if all type-3 vectors are
        drawn from the same P (each P generates 768 type-3 entries): the
        first 30 type-3 entries are all from P = A^odd[0], which is not
        enough to span R^24 on its own.  Using the full type-1 and type-2
        sets ensures rank 24 is verified reliably.
        """
        sample = list(dix_type1) + list(dix_type2)   # 720 + 11,520 = 12,240
        M = np.array([np.concatenate([x, y, z]) for x, y, z in sample])
        assert np.linalg.matrix_rank(M) == 24, \
            "Dixon's minimal vectors do not span R^24 (rank < 24)"


# ---------------------------------------------------------------------------
# 20–21: Inner products, and comparison with Wilson's minimal-vector shell
# ---------------------------------------------------------------------------

class TestShellComparison:
    """
    Two checks that the 196,560 vectors are the minimal vectors of a Leech
    lattice at the scale of Wilson's Λ (minimal norm 8).

    Test 20 is direct and does not use Wilson's construction: the inner
    products of minimal vectors of the √2-scaled Leech lattice lie in
    {0, ±2, ±4, ±8}, with ±8 only for v = ±w.  Test 21 is exact: the
    componentwise conjugates are Wilson's 196,560 minimal vectors, and
    componentwise conjugation is an isometry of R^24.
    """

    def test_inner_products_sample(self, dix_type1, dix_type2, dix_type3):
        """
        Inner products of 3,000 sampled vectors, together with the pair
        type-3 [0] and [411], lie in {0, ±2, ±4, ±8}; ±8 occurs only for
        v = ±w.
        """
        V = np.array([np.concatenate(t) for t in
                      list(dix_type1) + list(dix_type2) + list(dix_type3)])
        rng = np.random.default_rng(354)
        start3 = len(dix_type1) + len(dix_type2)
        idx = np.unique(np.concatenate([
            rng.choice(len(V), size=3000, replace=False),
            [start3, start3 + 411],
        ]))
        S = V[idx]
        # Coordinates are half-integers, so S @ S.T is computed exactly in floats.
        G = np.rint(S @ S.T).astype(int)
        assert np.array_equal(S @ S.T, G), "Non-integral inner products"
        assert set(np.unique(G)) <= {-8, -4, -2, 0, 2, 4, 8}, \
            f"Inner products outside {{0, ±2, ±4, ±8}}: {sorted(set(np.unique(G)))}"
        off = G - np.diag(np.diag(G))
        same = np.argwhere(np.abs(off) == 8)
        for i, j in same:
            assert np.allclose(np.abs(S[i]), np.abs(S[j])) and \
                (np.allclose(S[i], S[j]) or np.allclose(S[i], -S[j])), \
                "Inner product ±8 between two vectors that are not ±equal"

    def test_conjugate_shell_equals_wilson_shell(self, dix_type1, dix_type2, dix_type3):
        """
        The componentwise conjugates of all 196,560 Dixon vectors are exactly
        Wilson's 196,560 minimal vectors.

        Hence Dixon's vectors are the minimal vectors of the conjugate of
        Wilson's Λ, a Leech lattice of minimal norm 8 (Wilson's construction
        is verified in test_leech_wilson.py and key claim 004).
        """
        dixon_conj = (
            conjugate_keyset(dix_type1)
            | conjugate_keyset(dix_type2)
            | conjugate_keyset(dix_type3)
        )
        wilson = (
            to_keyset(leech_type1_vectors())
            | to_keyset(leech_type2_vectors())
            | to_keyset(leech_type3_vectors())
        )
        assert len(wilson) == 196_560
        assert dixon_conj == wilson
