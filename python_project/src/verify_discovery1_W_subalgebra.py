"""Verify the proposed Discovery 1: explicit subalgebra structure of
W = sigma(Ls)/2L inside the F2-octonion-algebra L/2L.

Tests:
  (a) Is the F2-octonion-algebra L/2L commutative? (Expected: yes.)
  (b) Does the proposed basis (w0, w1, w2, w3) really span W?
  (c) What is the actual multiplication table of W in this basis?
  (d) Does that match the proposed table?
  (e) Structural features: 3-dim nilpotent ideal? universal zero
      divisor? identity not idempotent in W?

Uses MUBAR from verify_mod2_quotient.py (structure constants of the
F2-octonion-algebra L/2L in the paper's L-basis).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_mod2_quotient import MUBAR, mubar, sLs, Lcoords, f2span, in_f2span


def f2_span_set(basis):
    """All 2^k linear combinations of the basis vectors."""
    out = {(0,)*8}
    for b in basis:
        out |= {tuple(x^y for x,y in zip(t, b)) for t in out}
    return out


def coords_in_basis(v, basis):
    """If v is in span(basis), return F2-coordinates (c0,...,c_{k-1}); else None."""
    n = len(basis)
    for mask in range(2**n):
        cand = [0]*len(v)
        for k in range(n):
            if (mask >> k) & 1:
                cand = [x^y for x,y in zip(cand, basis[k])]
        if tuple(cand) == tuple(v):
            return tuple((mask >> k) & 1 for k in range(n))
    return None


if __name__ == "__main__":
    print("=" * 70)
    print("VERIFICATION OF DISCOVERY 1: W = sigma(Ls)/2L subalgebra of L/2L")
    print("=" * 70)

    # ---------------------------------------------------------------
    # (a) Is L/2L commutative as an F2-algebra?
    # ---------------------------------------------------------------
    print("\n(a) Is the F2-octonion-algebra L/2L commutative?")
    e = [tuple(1 if k==i else 0 for k in range(8)) for i in range(8)]
    commutes = all(mubar(e[i], e[j]) == mubar(e[j], e[i])
                   for i in range(8) for j in range(8))
    print(f"    L/2L commutative on F2-basis (all 64 pairs): {commutes}")

    # ---------------------------------------------------------------
    # (b) Compute actual W = sigma(Ls)/2L; check proposed basis lies in W
    # ---------------------------------------------------------------
    print("\n(b) Compare proposed basis to actual W = sigma(Ls)/2L")
    W_gens = [tuple(c % 2 for c in Lcoords(v)) for v in sLs]
    W_basis = f2span(W_gens)
    W_set = f2_span_set(W_basis)
    print(f"    dim W = {len(W_basis)}  (expect 4)")
    print(f"    W (as 16-element F2-subspace of L/2L) computed.")

    proposed = {
        'w0': (1,1,1,0,0,0,0,0),
        'w1': (0,1,0,1,1,0,1,0),
        'w2': (0,0,0,1,0,1,1,1),
        'w3': (0,0,0,0,1,0,1,1),
    }
    print(f"\n    Proposed basis vectors and membership in W:")
    for name, v in proposed.items():
        is_in = v in W_set
        print(f"      {name} = {v}: in W = {is_in}")

    prop_list = [proposed['w0'], proposed['w1'], proposed['w2'], proposed['w3']]
    prop_span = f2_span_set(prop_list)
    print(f"\n    Proposed basis spans a {len(prop_span)}-element subspace "
          f"(expect 16 for a 4-dim span).")
    matches = prop_span == W_set
    print(f"    Proposed basis spans exactly W: {matches}")

    if not matches:
        print(f"\n    Extra in proposed (not in W): {len(prop_span - W_set)}")
        print(f"    Missing from proposed (in W but not span): {len(W_set - prop_span)}")
        # If proposed basis is not in W, abort.
        if not all(v in W_set for v in prop_list):
            print("\n    STOP: proposed basis vectors are not all in W.")
            sys.exit(1)

    # ---------------------------------------------------------------
    # (c) Compute the actual multiplication table of W in {w0, w1, w2, w3}
    # ---------------------------------------------------------------
    print("\n(c) Actual multiplication table of W in basis (w0, w1, w2, w3):")
    labels = ['w0', 'w1', 'w2', 'w3']
    actual = {}
    print(f"\n    {'·':>5} | {'w0':>14} {'w1':>14} {'w2':>14} {'w3':>14}")
    print(f"    {'-'*5} + {'-'*14} {'-'*14} {'-'*14} {'-'*14}")
    for i in range(4):
        row = []
        for j in range(4):
            prod = mubar(prop_list[i], prop_list[j])
            c = coords_in_basis(prod, prop_list)
            actual[(i,j)] = c
            if c is None:
                row.append("NOT_IN_W")
            else:
                terms = [labels[k] for k in range(4) if c[k]]
                row.append("+".join(terms) if terms else "0")
        print(f"    {labels[i]:>5} | " + " ".join(f"{r:>14}" for r in row))

    # ---------------------------------------------------------------
    # (d) Compare to proposed table
    # ---------------------------------------------------------------
    proposed_table = {
        (0,0): (1,1,0,0),  # w0+w1
        (0,1): (0,0,0,0),  # 0
        (0,2): (0,0,1,0),  # w2
        (0,3): (0,0,0,0),  # 0
        (1,0): (0,1,1,1),  # w1+w2+w3
        (1,1): (0,0,0,0),
        (1,2): (0,0,0,0),
        (1,3): (0,0,1,0),  # w2
        (2,0): (0,0,0,0),
        (2,1): (0,0,0,0),
        (2,2): (0,0,0,0),
        (2,3): (0,0,0,0),
        (3,0): (0,0,0,1),  # w3
        (3,1): (0,0,1,0),  # w2
        (3,2): (0,0,0,0),
        (3,3): (0,0,0,0),
    }
    print("\n(d) Comparison of proposed table vs actual:")
    mismatches = 0
    for i in range(4):
        for j in range(4):
            if actual[(i,j)] != proposed_table[(i,j)]:
                mismatches += 1
                print(f"    MISMATCH at ({labels[i]}, {labels[j]}): "
                      f"actual={actual[(i,j)]}, proposed={proposed_table[(i,j)]}")
    print(f"    Total mismatches: {mismatches} of 16")

    # ---------------------------------------------------------------
    # (e) Check symmetry of the actual table (must be symmetric)
    # ---------------------------------------------------------------
    print("\n(e) Symmetry of actual multiplication table on W:")
    sym = all(actual[(i,j)] == actual[(j,i)] for i in range(4) for j in range(4))
    print(f"    Symmetric (commutative): {sym}")
    if not sym:
        print(f"    (But L/2L itself is commutative from (a) -- contradiction.)")

    # ---------------------------------------------------------------
    # (f) If actual matches, check structural claims
    # ---------------------------------------------------------------
    print("\n(f) Structural features of W (using ACTUAL table):")
    # Universal zero divisor: w2 multiplied by anything in W gives 0?
    w2 = prop_list[2]
    w2_universal_zd = all(mubar(w2, w) == (0,)*8 for w in W_set)
    print(f"    w2 is a universal zero divisor (w2 . everything in W = 0): "
          f"{w2_universal_zd}")

    # Identity of L/2L: e_0 (since e_0 is the octonion identity).
    # Is e_0 in W? If so, is w0 the e_0 in W's basis decomposition?
    e0 = (1,0,0,0,0,0,0,0)
    e0_in_W = e0 in W_set
    print(f"    e_0 (identity of L/2L) in W: {e0_in_W}")
    if e0_in_W:
        e0_coords = coords_in_basis(e0, prop_list)
        e0_terms = [labels[k] for k in range(4) if e0_coords[k]]
        print(f"    e_0 in W basis = {'+'.join(e0_terms) if e0_terms else '0'}")

    # Nilpotent ideal N = span{w1, w2, w3}? Check it's an ideal first
    N_basis = [prop_list[1], prop_list[2], prop_list[3]]
    N_set = f2_span_set(N_basis)
    # N is an ideal iff for all w in W, n in N: w*n in N
    N_is_ideal = all(mubar(w, n) in N_set and mubar(n, w) in N_set
                     for w in W_set for n in N_set)
    print(f"    N = span(w1, w2, w3) is an ideal of W: {N_is_ideal}")

    # N^2 = span of all n1.n2 for n1, n2 in N
    N2_gens = [mubar(n1, n2) for n1 in N_set for n2 in N_set]
    N2 = f2_span_set(f2span([g for g in N2_gens if g != (0,)*8]) or [(0,)*8])
    print(f"    |N^2| = {len(N2)}  (proposed: span{{w2}} -> 2 elements)")

    # N^3 = N^2 . N
    N3_gens = [mubar(a, b) for a in N2 for b in N_set] + \
              [mubar(a, b) for a in N_set for b in N2]
    N3 = f2_span_set(f2span([g for g in N3_gens if g != (0,)*8]) or [(0,)*8])
    print(f"    |N^3| = {len(N3)}  (proposed: 0 -> 1 element)")

    # Number of non-zero entries in the table (counting both directions)
    nonzero = sum(1 for i in range(4) for j in range(4) if any(actual[(i,j)]))
    print(f"    Non-zero entries in 4x4 table: {nonzero}  "
          f"(proposed: 6)")
