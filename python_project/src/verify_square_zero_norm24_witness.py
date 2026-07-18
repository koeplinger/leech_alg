"""A square-zero vector of Lambda of norm 24: the norm-12 stratum is not all
of the square-zero locus.

Context (paper Section 5.4).  u = (x, y, z) in Lambda satisfies u star u = 0
iff x, y, z lie in (L sbar) ^ Im(O), have equal norms, and sum to zero
(verify_square_zero_classification.py).  The 4,032 hexagonal triples of the
126 norm-4 roots are exactly the NORM-12 stratum of that locus.  This script
shows the locus does not stop there: it constructs a square-zero vector of
Lambda of norm 24.

Construction.  Take two hexagonal triples (x1,y1,z1) and (x2,y2,z2) of E_7
roots that are BLOCKWISE ORTHOGONAL (<x1,x2> = <y1,y2> = <z1,z2> = 0), and
sum them blockwise:
    u = (x1+x2, y1+y2, z1+z2).
Each block is purely imaginary and lies in (L sbar) ^ Im(O) (a lattice
intersected with a subspace is closed under addition); the blocks have equal
norm 4 + 4 = 8 (orthogonality); and they sum to zero.  By the
characterization, u star u = 0, and N(u) = 3 * 8 = 24.

Both the membership u in Lambda and u star u = 0 are then re-checked
directly, with the same implementations used by the classification script.

Consistency: 24 lies in 12Z, so the proved constraint that every square-zero
norm of Lambda is a multiple of 12 (N(u) = 3 N(x) with N(x) even, and every
Lambda-norm in 4Z) is respected; the witness shows the multiples above 12 are
realized.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fractions import Fraction as F
from itertools import combinations, product as iproduct

from verify_mod2_quotient import omul, s, sbar
from verify_square_zero_classification import star, is_zero, inLambda


def ip(u, v):
    return sum(a * b for a, b in zip(u, v))


def N(v):
    return sum(c * c for c in v)


def e7_roots():
    """The 126 norm-4 vectors of (L sbar) ^ Im(O), as images lam*sbar of the
    E_8 roots lam with <lam, s> = 0."""
    h = F(1, 2)
    roots = []
    for i, j in combinations(range(8), 2):
        for si in (F(1), F(-1)):
            for sj in (F(1), F(-1)):
                v = [F(0)] * 8
                v[i] = si
                v[j] = sj
                roots.append(v)
    for sg in iproduct([h, -h], repeat=8):
        if sum(1 for t in sg if t < 0) % 2 == 1:
            roots.append(list(sg))
    assert len(roots) == 240
    V = [tuple(omul(lam, sbar)) for lam in roots if ip(lam, s) == 0]
    assert len(V) == 126
    return V


def main():
    print("=" * 68)
    print("Square-zero vector of Lambda at norm 24 (stratum-not-census witness)")
    print("=" * 68)

    V = e7_roots()
    Vset = set(V)

    # first hexagonal triple
    t1 = None
    for x in V:
        for y in V:
            if y == x:
                continue
            z = tuple(-(a + b) for a, b in zip(x, y))
            if z in Vset:
                t1 = (x, y, z)
                break
        if t1:
            break
    x1, y1, z1 = t1

    # second triple, blockwise orthogonal to the first
    t2 = None
    for x in V:
        if ip(x, x1) != 0:
            continue
        for y in V:
            if ip(y, y1) != 0:
                continue
            z = tuple(-(a + b) for a, b in zip(x, y))
            if z in Vset and ip(z, z1) == 0:
                t2 = (x, y, z)
                break
        if t2:
            break
    assert t2 is not None, "no blockwise-orthogonal second triple found"
    x2, y2, z2 = t2

    X = [a + b for a, b in zip(x1, x2)]
    Y = [a + b for a, b in zip(y1, y2)]
    Z = [a + b for a, b in zip(z1, z2)]
    u = (X, Y, Z)

    print(f"\n  block norms          : {N(X)}, {N(Y)}, {N(Z)}   (expect 8, 8, 8)")
    print(f"  purely imaginary     : {all(v[0] == 0 for v in (X, Y, Z))}")
    print(f"  blocks sum to zero   : {[a + b + c for a, b, c in zip(X, Y, Z)] == [F(0)] * 8}")
    in_lambda = inLambda(u)
    sq_zero = is_zero(star(u, u))
    norm = N(X) + N(Y) + N(Z)
    print(f"  u in Lambda          : {in_lambda}")
    print(f"  N(u)                 : {norm}   (expect 24)")
    print(f"  u star u == 0        : {sq_zero}")

    ok = (in_lambda and sq_zero and norm == 24
          and N(X) == N(Y) == N(Z) == 8)
    print()
    print("ALL PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
