"""
00_selftest.py -- sanity tests of the primitives in leechlib.py.

* octonion table: e1e2=e4 etc. as given by the Fano triples, e_i^2 = -e0,
  norm multiplicativity N(xy)=N(x)N(y), alternativity (xx)y = x(xy),
  y(xx) = (yx)x, and non-associativity, on random integer octonions;
* the exact-integer HNF / SNF / determinant / F2-rank routines against
  sympy (independent implementation) on random integer matrices.

Run:  python3 00_selftest.py > 00_selftest.out.txt
"""
import random
from math import gcd
import sympy
from sympy.matrices.normalforms import smith_normal_form
from leechlib import (MUL, FANO, omul, norm, hnf, hnf_index, snf_diagonal,
                      det_bareiss, f2_rank)

random.seed(20261013)

# --- octonion table -------------------------------------------------------
for (a, b, c) in FANO:
    for (p, q, r) in ((a, b, c), (b, c, a), (c, a, b)):
        assert MUL[p][q] == (1, r) and MUL[q][p] == (-1, r)
for i in range(1, 8):
    assert MUL[i][i] == (-1, 0)
print("Fano table: e_a e_b = +e_c on cyclic orders of the 7 triples, -e_c anticyclic, e_i^2 = -e0: OK")
print("e1e2 =", MUL[1][2], " e2e1 =", MUL[2][1], " e1e3 =", MUL[1][3], " (sign, index)")

def rnd():
    return [random.randint(-5, 5) for _ in range(8)]

for _ in range(2000):
    x, y, z = rnd(), rnd(), rnd()
    assert norm(omul(x, y)) == norm(x) * norm(y)
    assert omul(omul(x, x), y) == omul(x, omul(x, y))
    assert omul(omul(y, x), x) == omul(y, omul(x, x))
nonassoc = sum(omul(omul(x, y), z) != omul(x, omul(y, z))
               for x, y, z in ((rnd(), rnd(), rnd()) for _ in range(200)))
print("norm multiplicativity and left/right alternativity on 2000 random triples: OK")
print("non-associative triples among 200 random ones:", nonassoc)

# --- lattice routines vs sympy ---------------------------------------------
SIZES = [(random.randint(2, 7), None) for _ in range(60)] + [(24, 40)] * 6
for trial, (n, mm) in enumerate(SIZES):
    m = mm if mm else random.randint(n, n + 6)
    if mm:   # sparse-ish 24-column cases, closer to the real 576 x 24 use
        M = [[random.choice([0, 0, 0, 1, -1, 2]) for _ in range(n)] for _ in range(m)]
    else:
        M = [[random.randint(-6, 6) for _ in range(n)] for _ in range(m)]
    idx, H = hnf_index(M, n)
    S = sympy.Matrix(M)
    rank = S.rank()
    if rank < n:
        assert idx is None and len(H) == rank
        continue
    # index of the row lattice = gcd of the n x n minors = product of the SNF
    snf = smith_normal_form(S, domain=sympy.ZZ)
    d_sym = [abs(snf[i, i]) for i in range(n)]
    prod = 1
    for d in d_sym:
        prod *= d
    assert idx == prod, (idx, prod)
    # SNF of the square HNF vs sympy SNF of the original
    d_mine = snf_diagonal(H)
    assert sorted(d_mine) == sorted(d_sym), (d_mine, d_sym)
    # the HNF rows generate the same lattice: each original row is an integer
    # combination of H (solve the triangular system exactly)
    for r in M:
        v = list(r)
        for i, h in enumerate(H):
            piv_col = next(k for k in range(n) if h[k] != 0)
            assert v[piv_col] % h[piv_col] == 0
            q = v[piv_col] // h[piv_col]
            v = [a - q * b for a, b in zip(v, h)]
        assert all(a == 0 for a in v)
    # determinant
    sq = M[:n]
    assert det_bareiss(sq) == sympy.Matrix(sq).det()
    # F2 rank vs sympy's DomainMatrix rank over GF(2)
    from sympy.polys.matrices import DomainMatrix
    from sympy import GF
    dm = DomainMatrix([[GF(2)(a % 2) for a in row] for row in M], (m, n), GF(2))
    assert f2_rank(M) == dm.rank()
print("hnf / hnf_index / snf_diagonal / det_bareiss / f2_rank agree with sympy on 66 random matrices (60 small, 6 of size 40 x 24): OK")

# --- HNF-free invariants used by 06_certificates.py, against sympy SNF ------
from leechlib import rank_mod_p, log2_order_mod4_span
for trial in range(40):
    n = random.randint(2, 6)
    m = random.randint(n, n + 8)
    M = [[random.choice([0, 1, 2, 4, -2, 3, 8]) for _ in range(n)] for _ in range(m)]
    S = sympy.Matrix(M)
    if S.rank() < n:
        continue
    snf = smith_normal_form(S, domain=sympy.ZZ)
    d = [abs(snf[i, i]) for i in range(n)]
    for p in (2, 3, 5, 7):
        assert rank_mod_p(M, p) == n - sum(1 for x in d if x % p == 0)
    # |A/4A| = prod gcd(d_i, 4) ; |M| = 4^n / |A/4A|
    a4 = 1
    for x in d:
        a4 *= gcd(x, 4)
    assert 2 ** log2_order_mod4_span(M) * a4 == 4 ** n, (M, d)
print("rank_mod_p (p = 2,3,5,7) and log2_order_mod4_span agree with the sympy SNF on random full-rank matrices: OK")
