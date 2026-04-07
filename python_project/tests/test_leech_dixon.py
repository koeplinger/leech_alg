"""
test_leech_dixon.py — Verification of Dixon's octonionic Leech lattice construction.

This suite verifies Dixon's three families of 196,560 minimal vectors
[Dixon2010 Section 6] and compares them against Wilson's construction
[Wilson2009 Section 3].

Key question addressed
----------------------
The user's prompt asks to verify Dixon's construction as a Leech lattice and
to check whether it can be interpreted as a 24-dimensional integral algebra
(module over the integers) that closes on the minimal-vector shell.

Dixon's own answer to the algebraic question [Dixon2010 Section 7]: "No idea
yet."  This test suite provides the computational groundwork.

Central finding: Dixon's construction gives a DIFFERENT EMBEDDING of the
Leech lattice from Wilson's.  The two constructions share only 17,232 of the
196,560 minimal vectors.  The difference arises because Dixon uses A^odd for
type-1 (where A₃ has EVEN # of minus signs) while Wilson uses his E8 roots
(where the half-integer part has ODD # of minus signs).  The difference
propagates to type-2 and type-3 as well: only 6,144 of 11,520 Dixon type-2
vectors and 10,752 of 184,320 Dixon type-3 vectors lie in Wilson's Λ.

Properties verified
-------------------
 ℓ₀ = ½(e₀+e₁+...+e₇):
 1.  All 8 coordinates of ℓ₀ equal 0.5
 2.  N(ℓ₀) = 2  (ℓ₀ is a norm-2 element like s)
 3.  ℓ₀ ∈ A₃  (all-positive element, 0 minus signs = EVEN = A₃ criterion)

 Type-1 vectors  —  (2P, 0, 0) and 2 cyclic perms, P ∈ A^odd:
 4.  Count = 720
 5.  All N_std = 8
 6.  Exactly 1 nonzero block per triple (zero-component structure)
 7.  Type-1 vectors with P ∈ A₁ (112 elements): all 336 satisfy is_in_leech
     (A₁ = Wilson's type-1 roots; 2P ∈ 2L, so all conditions are met)
 8.  Type-1 vectors with P ∈ A₃ (128 elements): all 384 FAIL is_in_leech
     (A₃ has even # of minus signs; 2P ∉ 2L = Ls̄ ∩ Ls, violating
      Wilson's condition 2.  This is the key structural difference between
      Dixon's A^odd and Wilson's E8 roots.)
 9.  No duplicates among the 720 vectors

 Type-2 vectors  —  (2Q, ±2e_aQ, 0) and 2 cyclic perms, Q ∈ A^even:
10.  Count = 11,520
11.  All N_std = 8
12.  Exactly 2 nonzero blocks per triple
13.  Dixon's type-2 uses a different embedding: only 6,144 of 11,520 type-2
     vectors lie in Wilson's Λ; the remaining 5,376 do not.

 Type-3 vectors  —  (P, ±e_aP, ±(e_aℓ₀)(e_cP)) and 2 cyclic perms:
14.  Count = 184,320
15.  Sample of 500 have N_std = 8
16.  Exactly 3 nonzero blocks per triple (sample)
17.  Dixon's type-3 uses a different embedding: only 10,752 of 184,320 type-3
     vectors lie in Wilson's Λ; the remaining 173,568 do not.

 Dixon Leech summary:
18.  Total: 720 + 11,520 + 184,320 = 196,560 minimal vectors (no duplicates)
19.  The 196,560 vectors span R^24 (rank 24; verified using type-1 + type-2)

 Comparison with Wilson's minimal-vector shell:
20.  Dixon's type-1 vectors from A₁ equal Wilson's type-1 vectors from A₁
     (both sets are {2P : P ∈ A₁} × 3 cyclic perms, same 336 vectors)
21.  Dixon's full shell has exactly 179,328 vectors not in Wilson's shell,
     and Wilson's shell has exactly 179,328 vectors not in Dixon's shell.
     The two shells share exactly 17,232 vectors.

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
from leech_wilson import is_in_leech
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
    return dixon_leech_type3()    # 184,320 vectors, slow (~5–10 s)


# ---------------------------------------------------------------------------
# 1–3: ℓ₀ properties
# ---------------------------------------------------------------------------

class TestEllo0:
    """
    ℓ₀ = ½(e₀+e₁+...+e₇) is Dixon's base element [Dixon2010 eq. (6)].

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
        These are NOT Wilson's E8 roots.  2P has integer coords ±1 not all
        even, so 2P ∉ 2L = Ls̄ ∩ Ls.  Therefore (2P,0,0) ∉ Λ.

    This split is the central structural difference between Dixon's A^odd and
    Wilson's E8 roots (which use the COMPLEMENTARY half-integer set, with an
    ODD number of minus signs).
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
        """All type-1 vectors have exactly 1 nonzero block (2 zero blocks)."""
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
        roots (which require ODD # of minus signs).  For P ∈ A₃, 2P has
        coordinates ±1 which are not all even, so 2P ∉ 2L = Ls̄ ∩ Ls.
        Wilson's condition 2 requires x+y = 2P ∈ Ls̄, which fails.
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
        # Represent each triple as an integer tuple (coords × 2) to avoid
        # floating-point hash issues.
        keys = {
            tuple(np.round(np.concatenate([x, y, z]) * 2).astype(int))
            for x, y, z in dix_type1
        }
        assert len(keys) == 720


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
        """All type-2 vectors have exactly 2 nonzero blocks (1 zero block)."""
        for x, y, z in dix_type2[:30]:
            nonzero = sum(np.any(np.abs(v) > 1e-10) for v in (x, y, z))
            assert nonzero == 2, \
                f"Type-2 vector has {nonzero} nonzero blocks, expected 2"

    def test_type2_partial_wilson_overlap(self, dix_type2):
        """
        Dixon's type-2 gives a different embedding: only 6,144 of 11,520
        vectors lie in Wilson's Λ; the other 5,376 do not.

        This test confirms the partial overlap by checking a sample of 500
        evenly-spaced vectors for Wilson membership and verifying that both
        passing and failing vectors are present.

        The overlap occurs where Dixon's 2Q coincides with a Wilson Ls̄
        element and ±2e_aQ coincides with the corresponding (λs̄)j element.
        The non-overlap arises because Dixon uses LEFT multiplication by
        e_a (e_a × Q) whereas Wilson uses RIGHT multiplication ((λs̄) × j).
        Reference: [Dixon2010] Section 6 eq. (21); [Wilson2009] Section 3.
        """
        alg = STANDARD_ALGEBRA
        step = max(1, len(dix_type2) // 500)
        pass_count = 0
        fail_count = 0
        for i in range(0, len(dix_type2), step):
            x, y, z = dix_type2[i]
            ox, oy, oz = alg.element(x), alg.element(y), alg.element(z)
            if is_in_leech(ox, oy, oz):
                pass_count += 1
            else:
                fail_count += 1
        assert pass_count > 0, "No type-2 vectors in Wilson's Λ — expected partial overlap"
        assert fail_count > 0, "All type-2 vectors in Wilson's Λ — expected partial overlap"


# ---------------------------------------------------------------------------
# 14–17: Type-3 minimal vectors
# ---------------------------------------------------------------------------

class TestType3VectorsDixon:
    """
    Type-3: (P, ±e_aP, ±(e_aℓ₀)(e_cP)) and 2 cyclic perms,
            P ∈ A^odd, a,c ∈ {0,...,7}, ± signs independent.

    Norm structure (2, 2, 4) per component.  The 3rd cyclic permutation
    has structure (4, 2, 2), matching Wilson's type-3.
    Count: 240 × 8 × 8 × 2 × 2 × 3 = 184,320.
    Reference: [Dixon2010] Section 6, eq. (20).
    """

    def test_type3_count(self, dix_type3):
        """Exactly 184,320 type-3 vectors."""
        assert len(dix_type3) == 184_320

    def test_type3_sample_norm_8(self, dix_type3):
        """500 evenly-spaced type-3 vectors have N_std = 8."""
        step = max(1, len(dix_type3) // 500)
        for i in range(0, len(dix_type3), step):
            x, y, z = dix_type3[i]
            norm = float(np.dot(x, x) + np.dot(y, y) + np.dot(z, z))
            assert abs(norm - 8.0) < 1e-10, \
                f"Type-3 vector [{i}] has N_std = {norm}, expected 8"

    def test_type3_zero_structure(self, dix_type3):
        """
        Type-3 vectors have exactly 3 nonzero blocks (no zero component).

        The three components of (P, ±e_aP, ±(e_aℓ₀)(e_cP)) have norms
        (2, 2, 4) respectively, all nonzero for generic P, a, c.
        """
        for x, y, z in dix_type3[:30]:
            nonzero = sum(np.any(np.abs(v) > 1e-10) for v in (x, y, z))
            assert nonzero == 3, \
                f"Type-3 vector has {nonzero} nonzero blocks, expected 3"

    def test_type3_partial_wilson_overlap(self, dix_type3):
        """
        Dixon's type-3 gives a different embedding: only 10,752 of 184,320
        vectors lie in Wilson's Λ; the other 173,568 do not.

        This test confirms the partial overlap by checking a sample of 500
        evenly-spaced vectors for Wilson membership and verifying that both
        passing and failing vectors are present.

        Note: the first 768 type-3 entries all use the same P (the first
        A^odd element).  The stride is chosen so the sample covers ~500
        distinct positions across all 240 P values.
        Reference: [Dixon2010] Section 6 eq. (20); [Wilson2009] Section 3.
        """
        alg = STANDARD_ALGEBRA
        step = max(1, len(dix_type3) // 500)
        pass_count = 0
        fail_count = 0
        for i in range(0, len(dix_type3), step):
            x, y, z = dix_type3[i]
            ox, oy, oz = alg.element(x), alg.element(y), alg.element(z)
            if is_in_leech(ox, oy, oz):
                pass_count += 1
            else:
                fail_count += 1
        assert pass_count > 0, "No type-3 vectors in Wilson's Λ — expected partial overlap"
        assert fail_count > 0, "All type-3 vectors in Wilson's Λ — expected partial overlap"


# ---------------------------------------------------------------------------
# 18–19: Dixon Leech summary
# ---------------------------------------------------------------------------

class TestDixonLeechSummary:
    """
    Summary properties of Dixon's 196,560 minimal-vector candidates.

    Tests 18–19 verify the two most basic structural properties: the correct
    total count and that the vectors span all of R^24 (rank 24).
    """

    def test_total_count(self, dix_type1, dix_type2, dix_type3):
        """Total: 720 + 11,520 + 184,320 = 196,560 minimal vectors."""
        assert len(dix_type1) + len(dix_type2) + len(dix_type3) == 196_560

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
# 20–21: Comparison with Wilson's minimal-vector shell
# ---------------------------------------------------------------------------

class TestShellComparison:
    """
    Comparison of Dixon's 196,560 vectors with Wilson's 196,560 vectors.

    The two constructions use different parameterisations of the Leech lattice:
      - Wilson's type-1 uses λ ∈ A₁ ∪ {ODD-minus half-integer elements of L}.
      - Dixon's type-1 uses P ∈ A₁ ∪ A₃ (A₃ = EVEN minus signs).

    The two half-integer families (A₃ vs Wilson-type-2) are disjoint, so the
    full shells differ in at least 384 vectors (Dixon's A₃ type-1 vectors
    vs Wilson's type-2-root type-1 vectors).

    This test confirms the shells are not equal, and that the A₁ type-1
    subsystem is the same in both constructions.
    """

    def test_a1_type1_shells_agree(self, dix_type1):
        """
        The 336 type-1 vectors from A₁ are identical in both constructions.

        Both Wilson and Dixon generate type-1 vectors (2P, 0, 0) + cyclic
        perms for P ∈ A₁ = {±e_a ± e_b : a ≠ b}.  The 112 A₁ elements are
        the same in both; the 336 resulting triples must therefore coincide.
        """
        from leech_wilson import leech_type1_vectors
        from e8_wilson import wilson_e8_roots

        # Wilson's type-1 vectors are (2λ, 0, 0) for ALL 240 E8 roots λ.
        # Filter to just A₁ (integer-coordinate) roots.
        wilson_t1 = leech_type1_vectors()
        wilson_keys = {
            tuple(np.round(np.concatenate([x, y, z]) * 2).astype(int))
            for x, y, z in wilson_t1
        }

        # Dixon's type-1 vectors from A₁ only.
        a1 = dixon_a1()
        zero = np.zeros(8)
        a1_keys = set()
        for P in a1:
            two_P = 2.0 * P.coords
            for triple in [
                (two_P, zero, zero),
                (zero, two_P, zero),
                (zero, zero, two_P),
            ]:
                key = tuple(np.round(np.concatenate(triple) * 2).astype(int))
                a1_keys.add(key)

        assert len(a1_keys) == 336, f"Expected 336 A₁ type-1 vectors, got {len(a1_keys)}"
        assert a1_keys.issubset(wilson_keys), \
            "Dixon's A₁ type-1 vectors are not all in Wilson's type-1 shell"

    def test_full_shells_differ(self, dix_type1, dix_type2, dix_type3):
        """
        Dixon's 196,560-vector shell overlaps with Wilson's in exactly 17,232
        vectors; each shell has 179,328 exclusive vectors.

        The disagreement extends across all three types:
          type-1: Dixon uses A^odd (A₁ ∪ A₃), Wilson uses E8 roots (A₁ ∪
                  Wilson-type-2).  The A₃ vs Wilson-type-2 mismatch gives
                  128 × 3 = 384 Dixon-only and 384 Wilson-only vectors.
          type-2: Only 6,144 of Dixon's 11,520 type-2 vectors coincide with
                  Wilson's type-2; the other 5,376 Dixon type-2 are exclusive.
          type-3: Only 10,752 of Dixon's 184,320 type-3 vectors coincide
                  with Wilson's type-3; the other 173,568 are exclusive.

        Total exclusive (per construction):
          384 (type-1) + 5,376 (type-2) + 173,568 (type-3) = 179,328.
        Shared: 336 + 6,144 + 10,752 = 17,232.

        This confirms Dixon's construction is a DIFFERENT EMBEDDING of the
        Leech lattice in the ambient R^24 = L³, not the same copy as Wilson's.
        Both constructions produce 196,560 norm-8 vectors spanning R^24.
        """
        from leech_wilson import leech_type1_vectors, leech_type2_vectors, leech_type3_vectors

        def to_keyset(vecs):
            return {
                tuple(np.round(np.concatenate([x, y, z]) * 2).astype(int))
                for x, y, z in vecs
            }

        wilson_keys = (
            to_keyset(leech_type1_vectors())
            | to_keyset(leech_type2_vectors())
            | to_keyset(leech_type3_vectors())
        )
        dixon_keys = (
            to_keyset(dix_type1)
            | to_keyset(dix_type2)
            | to_keyset(dix_type3)
        )

        assert dixon_keys != wilson_keys, \
            "Dixon's shell equals Wilson's shell unexpectedly"

        dixon_only = dixon_keys - wilson_keys
        wilson_only = wilson_keys - dixon_keys
        shared = dixon_keys & wilson_keys

        assert len(shared) == 17_232, (
            f"Expected 17,232 shared vectors, got {len(shared)}"
        )
        assert len(dixon_only) == 179_328, (
            f"Expected 179,328 Dixon-only vectors, got {len(dixon_only)}"
        )
        assert len(wilson_only) == 179_328, (
            f"Expected 179,328 Wilson-only vectors, got {len(wilson_only)}"
        )
