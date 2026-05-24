"""Verify the proposed Discovery 2: V = sigma(Ls_bar)/2L is a totally
isotropic left ideal of L/2L, and W has mixed isotropy.

The quadratic form on L/2L: for v in L (with norm form N(v) = sum v_i^2),
q(v mod 2L) := N(v)/2 mod 2.  Well-defined because L is even (N takes
even values on L) and N(v+2u) = N(v) + 4(v.u + N(u)), which is even
modulo 4, hence mod-2 of N/2 is invariant under v -> v + 2u.

Tests:
  (a) Proposed basis of V spans V = sigma(Ls_bar)/2L.
  (b) V is a left ideal: mu-bar(L/2L, V) <= V.  (Already in
      verify_mod2_quotient.py as Lemma 4.3 mod-2 form.)
  (c) V is totally isotropic: q(v) = 0 for all v in V.
  (d) W's isotropy distribution: how many elements have q = 0 vs 1.
  (e) V (+) W = L/2L (already in verify_mod2_quotient.py).
  (f) Compare V to W's nilpotent ideal N = span(w1, w2, w3): is N = V?
      Is N subset of V? Or are they disjoint?
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from verify_mod2_quotient import (
    MUBAR, mubar, L, sLs, sLsbar, omul, Lcoords, f2span, in_f2span,
)


def f2_span_set(basis):
    out = {(0,)*8}
    for b in basis:
        out |= {tuple(x^y for x,y in zip(t, b)) for t in out}
    return out


def coords_in_basis(v, basis):
    n = len(basis)
    for mask in range(2**n):
        cand = [0]*len(v)
        for k in range(n):
            if (mask >> k) & 1:
                cand = [x^y for x,y in zip(cand, basis[k])]
        if tuple(cand) == tuple(v):
            return tuple((mask >> k) & 1 for k in range(n))
    return None


def L_basis_to_e_basis(coords_in_L):
    """Given F2-coords (c_1,...,c_8) of a coset v + 2L in the L-basis,
    return one integer representative v in R^8 (e_0,...,e_7 coords)."""
    v = [F(0)] * 8
    for k in range(8):
        if coords_in_L[k]:
            for i in range(8):
                v[i] = v[i] + L[k][i]
    return v


def norm_int(v):
    """N(v) = sum v_i^2 as an integer (for v in L; will be even on L)."""
    s = F(0)
    for x in v:
        s = s + x*x
    assert s.denominator == 1, f"non-integer norm: {s}"
    return int(s)


def q_form(coords_in_L):
    """q(v mod 2L) = N(v)/2 mod 2, for v reconstructed from L-basis F2-coords."""
    v = L_basis_to_e_basis(coords_in_L)
    n = norm_int(v)
    assert n % 2 == 0, f"L is even but got N(v) = {n}, v = {v}"
    return (n // 2) % 2


if __name__ == "__main__":
    print("=" * 70)
    print("VERIFICATION OF DISCOVERY 2: V = sigma(Ls_bar)/2L isotropic")
    print("=" * 70)

    # Build V and W from the existing machinery
    V_basis_actual = f2span([tuple(c%2 for c in Lcoords(v)) for v in sLsbar])
    W_basis_actual = f2span([tuple(c%2 for c in Lcoords(v)) for v in sLs])
    V_set = f2_span_set(V_basis_actual)
    W_set = f2_span_set(W_basis_actual)
    print(f"\ndim V = {len(V_basis_actual)}, dim W = {len(W_basis_actual)}")
    print(f"|V| = {len(V_set)}, |W| = {len(W_set)}  (each should be 16)")

    # ----------------------------------------------------------------
    # (a) Proposed basis of V
    # ----------------------------------------------------------------
    proposed_V = [
        (0,1,1,0,0,0,0,0),   # v0
        (0,0,1,0,1,0,1,0),   # v1
        (0,0,0,1,0,0,1,1),   # v2
        (0,0,0,0,1,1,0,1),   # v3
    ]
    print("\n(a) Proposed basis of V and membership in actual V:")
    for i, v in enumerate(proposed_V):
        print(f"    v{i} = {v}: in V = {v in V_set}")
    prop_V_span = f2_span_set(proposed_V)
    print(f"    span of proposed (4 vectors) has size {len(prop_V_span)} (expect 16)")
    print(f"    proposed spans exactly V: {prop_V_span == V_set}")

    # ----------------------------------------------------------------
    # (b) V is a left ideal: mu-bar(L/2L, V) <= V
    # ----------------------------------------------------------------
    print("\n(b) V is a left ideal: mu-bar(L/2L, V) <= V")
    e = [tuple(1 if k==i else 0 for k in range(8)) for i in range(8)]
    full_LbarL = f2_span_set(e)  # all of L/2L as F2 vectors
    is_left_ideal = all(mubar(a, v) in V_set for a in full_LbarL for v in V_set)
    print(f"    V is a left ideal: {is_left_ideal}")
    # Also: is V a right ideal? two-sided?
    is_right_ideal = all(mubar(v, a) in V_set for v in V_set for a in full_LbarL)
    print(f"    V is a right ideal: {is_right_ideal}")

    # ----------------------------------------------------------------
    # (c) V is totally isotropic under q(v) = N(v)/2 mod 2
    # ----------------------------------------------------------------
    print("\n(c) V totally isotropic under q(v) = N(v)/2 mod 2:")
    V_iso = {q_form(v) for v in V_set}
    print(f"    {{q(v) : v in V}} = {V_iso}  (expect {{0}})")
    V_count = {0: 0, 1: 0}
    for v in V_set:
        V_count[q_form(v)] += 1
    print(f"    distribution in V: q=0 -> {V_count[0]}, q=1 -> {V_count[1]}")

    # ----------------------------------------------------------------
    # (d) W's isotropy distribution
    # ----------------------------------------------------------------
    print("\n(d) W's isotropy distribution:")
    W_count = {0: 0, 1: 0}
    for w in W_set:
        W_count[q_form(w)] += 1
    print(f"    distribution in W: q=0 -> {W_count[0]}, q=1 -> {W_count[1]}")
    print(f"    (proposed: q=0 -> 8, q=1 -> 8)")

    # ----------------------------------------------------------------
    # (e) Full L/2L isotropy (sanity check on the +-type form)
    # ----------------------------------------------------------------
    full_count = {0: 0, 1: 0}
    for x in full_LbarL:
        full_count[q_form(x)] += 1
    print(f"\n(e) L/2L isotropy distribution:")
    print(f"    q=0 -> {full_count[0]}, q=1 -> {full_count[1]}")
    print(f"    (plus-type form on F2^8: q=0 -> 136, q=1 -> 120;")
    print(f"     minus-type form on F2^8: q=0 -> 120, q=1 -> 136)")

    # ----------------------------------------------------------------
    # (f) Compare V (left ideal) to N (nilpotent ideal of W)
    # ----------------------------------------------------------------
    print("\n(f) Compare V (left ideal of L/2L) to W's nilpotent backbone")
    print("    N = span(w1, w2, w3):")
    w1 = (0,1,0,1,1,0,1,0)
    w2 = (0,0,0,1,0,1,1,1)
    w3 = (0,0,0,0,1,0,1,1)
    N_set = f2_span_set([w1, w2, w3])
    print(f"    |N| = {len(N_set)}  (expect 8 = 2^3)")
    print(f"    N subset of V?      {N_set <= V_set}")
    print(f"    V subset of N?      {V_set <= N_set}")
    print(f"    |N intersect V|     {len(N_set & V_set)}")
    print(f"    |N intersect W|     {len(N_set & W_set)}  (expect 8 -- N is in W)")
    print(f"    N totally isotropic?  {all(q_form(n) == 0 for n in N_set)}")
    # The user's question: is V the same/similar ideal as N?
    # V is in V (the "other half" of L/2L from W), N is in W.
    # They live in disjoint summands (modulo the zero vector).
    print(f"\n    Summary: V lives in V (V ∩ W = {{0}}), N lives in W.")
    print(f"    V is a 4-dim left ideal of (L/2L, mu-bar).")
    print(f"    N is a 3-dim nilpotent ideal of W (the subalgebra).")
    print(f"    They are structurally distinct: V is the left-ideal half")
    print(f"    of the L/2L = V (+) W decomposition; N is the radical-like")
    print(f"    piece *inside* the subalgebra half W.")
