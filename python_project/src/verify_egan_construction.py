"""The four links of the Baez-Egan verification chain (paper Section 6,
the sentence "We have likewise independently reproduced the full
Baez--Egan verification chain") that no other script of the repository
checks, in exact arithmetic.  The other two links have their own
scripts: the lower bound 244,035,421 (verify_egan_lower_bound.py) and
Egan's count 17,280 = 270 x 64 (egan_baez_count.py).

  (1) Egan's simple roots of E8 and the two sublattices built from them
      (Baez, "Integral octonions", part 9, Egan's text, with its explicit
      rows r_i, w_{i,1}, w_{i,2}): the Gram matrix of r_1, ..., r_8 is a
      simply laced Cartan matrix of determinant 1, that is, of type E8;
      w_{i,1} + w_{i,2} = r_i; w_{i,1}.r_i = w_{i,2}.r_i = 1; and the
      Gram matrices of the w_{i,1} and of the w_{i,2} are twice that of
      the r_i.
  (2) The lemma that these meet in twice that E8 lattice: with L0, L1,
      L2 the lattices spanned by the r_i, the w_{i,1} and the w_{i,2},
      2L0 lies in L1 and in L2, [L0 : L1] = [L0 : L2] = 16, and L1 + L2 =
      L0 (since r_i = w_{i,1} + w_{i,2}); counting in L0/2L0 then gives
      L1 cap L2 = 2L0.
  (3) The Leech-lattice theorem for that pair: Lambda(L1, L2) := {(a, b,
      c) in L0^3 : a+b, a+c, b+c in L1; a+b+c in L2} contains (2L0)^3,
      and counting its classes in (L0/2L0)^3 gives index 2^12 in L0^3;
      the 24 vectors (a, 0, -a), (b, -b, 0) (a, b in the basis of L1)
      and (-w, w, w) (w in the basis of L2) lie in it with determinant
      2^12, so they are a Z-basis.  Half its Gram matrix is integral, even and
      of determinant 1, and it contains none of the 720 vectors of norm 2
      and none of the 179,280 vectors of norm 4 of L0^3.  So
      Lambda(L1, L2)/sqrt2 is an even unimodular lattice of rank 24
      without roots: the Leech lattice, by Conway's characterization.
  (4) The doubled-Jordan-product closure (Baez, part 10, Theorem 4) for a
      Leech lattice at the standard scale, in the paper's conventions:
      L_L := a Lambda, with a = (e1+e2)/2 applied on the left in each
      block and Lambda Wilson's lattice (paper Definition 2.1), with the
      Z-basis B of verify_consecutive_twists_exact.py (checked here: B
      lies in Lambda with determinant 2^12 = [L^3 : Lambda] in L^3
      coordinates).  Since a L = MO, the maximal order of Section 2.2 (its
      definition (2.2)), L_L lies in MO^3, and its Gram matrix is half
      that of Lambda.  For
      the 27-dimensional lattice Z^3 + L_L of Hermitian matrices
          ( alpha   X      Y   )
          ( X*      beta   Z   )      alpha, beta, gamma in Z,
          ( Y*      Z*     gamma)     (X, Y, Z) in L_L
      (Baez's placement in part 10), 2(PQ + QP) lies in Z^3 + L_L for
      every pair of its 27 basis elements, while PQ + QP does not (a
      witness is printed).

Egan's data are in the coordinates of part 9: L0 is the E8 lattice
whose half-integer vectors have even coordinate sum, not the paper's L.
The rows below are those of the part-9 post, transcribed from the page.

Run: python3 verify_egan_construction.py  ->  ALL PASS (exit 0).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from itertools import product as iproduct
from sympy import Matrix
from verify_consecutive_twists_exact import B, l3, omul, inLambda, add

h = F(1, 2)
def fr(*c): return [F(x) for x in c]

# ---------------------------------------------------------------- Egan's rows (part 9)
r = [fr(-1, -1, 0, 0, 0, 0, 0, 0), fr(0, 1, -1, 0, 0, 0, 0, 0),
     fr(0, 0, 1, -1, 0, 0, 0, 0), fr(0, 0, 0, 1, -1, 0, 0, 0),
     fr(0, 0, 0, 0, 1, -1, 0, 0), fr(0, 0, 0, 0, 0, 1, -1, 0),
     fr(0, 0, 0, 0, 0, 0, 1, -1), fr(h, -h, -h, -h, -h, h, h, h)]
w1 = [fr(0, -1, 1, 1, 1, 0, 0, 0), fr(-1, 0, -1, 0, -1, 1, 0, 0),
      fr(-h, h, h, -h, h, -3 * h, h, h), fr(1, -1, 0, 0, -1, 0, 0, -1),
      fr(0, 1, -1, 1, 1, 0, 0, 0), fr(-h, h, 3 * h, -h, -h, h, -h, -h),
      fr(h, -h, -h, -h, h, h, 3 * h, h), fr(0, -1, 0, -1, 0, 0, -1, 1)]
w2 = [fr(-1, 0, -1, -1, -1, 0, 0, 0), fr(1, 1, 0, 0, 1, -1, 0, 0),
      fr(h, -h, h, -h, -h, 3 * h, -h, -h), fr(-1, 1, 0, 1, 0, 0, 0, 1),
      fr(0, -1, 1, -1, 0, -1, 0, 0), fr(h, -h, -3 * h, h, h, h, -h, h),
      fr(-h, h, h, h, -h, -h, -h, -3 * h), fr(h, h, -h, h, -h, h, 3 * h, -h)]

def dot(u, v): return sum(a * b for a, b in zip(u, v))
def gram(vs): return [[dot(u, v) for v in vs] for u in vs]
def sub(u, v): return [a - b for a, b in zip(u, v)]
def neg(u): return [-a for a in u]

checks = []
def check(name, ok):
    checks.append(bool(ok))
    print(("  ok    " if ok else "  FAIL  ") + name)


class Lat:
    """Membership in the lattice spanned by 8 rows (doubled coordinates
    are integers): v is in it iff 2v . adj(2W) = 0 mod det(2W)."""
    def __init__(self, rows):
        W2 = Matrix([[int(2 * x) for x in row] for row in rows])
        self.det = int(W2.det())
        self.adj = [[int(x) for x in row] for row in W2.adjugate().tolist()]
        assert self.det != 0
    def __contains__(self, v):
        v2 = [int(2 * x) for x in v]
        if any(2 * x != y for x, y in zip(v, v2)):
            return False
        d = self.det
        return all(sum(v2[i] * self.adj[i][j] for i in range(8)) % d == 0
                   for j in range(8))


def link1():
    print("(1) Egan's simple roots and the two sublattices built from them")
    C = gram(r)
    simply_laced = all(C[i][i] == 2 for i in range(8)) and all(
        C[i][j] in (0, -1) for i in range(8) for j in range(8) if i != j)
    edges = sum(1 for i in range(8) for j in range(i + 1, 8) if C[i][j] == -1)
    check("Gram(r) is a simply laced Cartan matrix with 7 edges, det %s"
          % Matrix(C).det(), simply_laced and edges == 7 and Matrix(C).det() == 1)
    check("w_{i,1} + w_{i,2} = r_i for all i",
          all(add(w1[i], w2[i]) == r[i] for i in range(8)))
    check("w_{i,1}.r_i = w_{i,2}.r_i = 1 for all i",
          all(dot(w1[i], r[i]) == 1 and dot(w2[i], r[i]) == 1 for i in range(8)))
    check("Gram(w_{.,1}) = Gram(w_{.,2}) = 2 Gram(r)",
          gram(w1) == [[2 * x for x in row] for row in C] == gram(w2))


def link2(L0, L1, L2):
    print("(2) the lemma L1 cap L2 = 2 L0")
    check("w_{i,j} lie in L0", all(v in L0 for v in w1 + w2))
    check("2 L0 lies in L1 and in L2",
          all([2 * x for x in v] in L1 and [2 * x for x in v] in L2 for v in r))
    i1 = abs(Matrix(w1).det() / Matrix(r).det())
    i2 = abs(Matrix(w2).det() / Matrix(r).det())
    check("[L0 : L1] = %s, [L0 : L2] = %s" % (i1, i2), i1 == 16 and i2 == 16)
    # |L1/2L0| = |L2/2L0| = 256/16 = 16 and L1 + L2 = L0 (r_i = w_{i,1} + w_{i,2}),
    # so |(L1 cap L2)/2L0| = 16 * 16 / 256 = 1.
    check("hence L1 cap L2 = 2 L0", i1 * i2 == 256)


def roots_and_norm4(L0):
    """The 240 roots and 2160 vectors of norm 4 of L0, by search over
    entries in {-2, ..., 2} and {+-1/2, +-3/2}."""
    out2, out4 = [], []
    for vals in (range(-2, 3), (-3 * h, -h, h, 3 * h)):
        for t in iproduct([F(x) for x in vals], repeat=8):
            n = dot(t, t)
            if n in (2, 4) and list(t) in L0:
                (out2 if n == 2 else out4).append(list(t))
    return out2, out4


def link3(L0, L1, L2):
    print("(3) the Leech-lattice theorem for that pair")
    Z = [F(0)] * 8
    gens = ([(a, Z, neg(a)) for a in w1] + [(b, neg(b), Z) for b in w1]
            + [(neg(w), w, w) for w in w2])
    def inLam(u):
        x, y, z = u
        return (add(x, y) in L1 and add(x, z) in L1 and add(y, z) in L1
                and add(x, y, z) in L2)
    check("the 24 generators lie in Lambda(L1, L2)", all(inLam(g) for g in gens))
    Rinv = Matrix(r).inv()
    def coords(v): return [int(x) for x in (Matrix([v]) * Rinv)]
    # Lambda(L1, L2) contains (2L0)^3 (2L0 lies in L1 and in L2), so its index
    # in L0^3 is 2^24 / |Lambda/(2L0)^3|; count its classes in r-coordinates mod 2.
    def mask(v): return sum(1 << i for i, x in enumerate(coords(v)) if x % 2)
    def span2(vs):
        S = {0}
        for m in map(mask, vs):
            S |= {t ^ m for t in S}
        return S
    A1, A2 = span2(w1), span2(w2)
    n = sum(1 for x in range(256) for y in range(256) if x ^ y in A1
            for z in range(256) if x ^ z in A1 and y ^ z in A1 and x ^ y ^ z in A2)
    check("[L0^3 : Lambda(L1, L2)] = 2^24/%d = 2^12" % n, n == 2 ** 12)
    G = Matrix([coords(x) + coords(y) + coords(z) for x, y, z in gens])
    check("det of the generators in L0^3 coordinates = %s = 2^12, so they are a Z-basis"
          % abs(G.det()), abs(G.det()) == 2 ** 12)
    full = [x + y + z for x, y, z in gens]
    Gr = [[dot(u, v) / 2 for v in full] for u in full]
    integral = all(x.denominator == 1 for row in Gr for x in row)
    even = all(Gr[i][i] % 2 == 0 for i in range(24))
    d = Matrix([[int(x) for x in row] for row in Gr]).det() if integral else None
    check("half the Gram matrix is integral, even, det %s" % d,
          integral and even and d == 1)
    n2, n4 = roots_and_norm4(L0)
    check("L0 has 240 roots and 2160 vectors of norm 4", len(n2) == 240 and len(n4) == 2160)
    short2 = [tuple(v if k == c else Z for k in range(3)) for v in n2 for c in range(3)]
    short4 = ([tuple(v if k == c else Z for k in range(3)) for v in n4 for c in range(3)]
              + [tuple(u if k == c1 else (v if k == c2 else Z) for k in range(3))
                 for (c1, c2) in ((0, 1), (0, 2), (1, 2)) for u in n2 for v in n2])
    check("L0^3 has 720 vectors of norm 2 and 179,280 of norm 4",
          len(short2) == 720 and len(short4) == 179280)
    check("none of them lies in Lambda(L1, L2)",
          not any(inLam(u) for u in short2) and not any(inLam(u) for u in short4))


# ---------------------------------------------------------------- (4) Jordan closure
a = fr(0, h, h, 0, 0, 0, 0, 0)                     # (e1 + e2)/2
abar = [a[0]] + [-x for x in a[1:]]

def conj(x): return [x[0]] + [-t for t in x[1:]]
def scal(c): return [F(c)] + [F(0)] * 7

def herm(alpha, beta, gamma, X, Y, Z):
    return [[scal(alpha), X, Y], [conj(X), scal(beta), Z], [conj(Y), conj(Z), scal(gamma)]]

def mmul(P, Q):
    return [[add(*[omul(P[i][k], Q[k][j]) for k in range(3)]) for j in range(3)]
            for i in range(3)]

def jordan(P, Q, c):
    PQ, QP = mmul(P, Q), mmul(Q, P)
    return [[[c * (x + y) for x, y in zip(PQ[i][j], QP[i][j])] for j in range(3)]
            for i in range(3)]

def in_ZLL(M):
    """M in Z^3 + L_L: Hermitian, integral real diagonal, (X, Y, Z) in a Lambda."""
    for i in range(3):
        d = M[i][i]
        if any(t != 0 for t in d[1:]) or d[0].denominator != 1:
            return False
    if not (M[1][0] == conj(M[0][1]) and M[2][0] == conj(M[0][2])
            and M[2][1] == conj(M[1][2])):
        return False
    X, Y, Z = M[0][1], M[0][2], M[1][2]
    back = tuple([2 * t for t in omul(abar, w)] for w in (X, Y, Z))   # a^{-1} = 2 abar
    return inLambda(back)


def link4():
    print("(4) the doubled-Jordan-product closure at the standard scale")
    dB = abs(Matrix([l3(b) for b in B]).det())
    check("the 24 vectors B lie in Lambda with |det| %s = [L^3 : Lambda] = 2^12: a Z-basis"
          % dB, dB == 2 ** 12 and all(inLambda(b) for b in B))
    LL = [tuple(omul(a, x) for x in b) for b in B]
    Gl = [[sum(dot(p, q) for p, q in zip(u, v)) for v in LL] for u in LL]
    Gb = [[sum(dot(p, q) for p, q in zip(u, v)) for v in B] for u in B]
    check("Gram(a Lambda) = Gram(Lambda)/2 (the standard scale)",
          Gl == [[x / 2 for x in row] for row in Gb])
    check("2 abar . (a x) = x on the basis of Lambda",
          all(tuple([2 * t for t in omul(abar, w)] for w in u) == tuple(b)
              for u, b in zip(LL, B)))
    Z8 = [F(0)] * 8
    basis = ([herm(1, 0, 0, Z8, Z8, Z8), herm(0, 1, 0, Z8, Z8, Z8), herm(0, 0, 1, Z8, Z8, Z8)]
             + [herm(0, 0, 0, *u) for u in LL])
    pairs = [(i, j) for i in range(27) for j in range(i, 27)]
    check("2(PQ + QP) in Z^3 + L_L for all %d pairs of basis elements" % len(pairs),
          all(in_ZLL(jordan(basis[i], basis[j], 2)) for i, j in pairs))
    fails = [(i, j) for i, j in pairs if not in_ZLL(jordan(basis[i], basis[j], 1))]
    check("PQ + QP leaves Z^3 + L_L (%d of %d pairs; first witness %s)"
          % (len(fails), len(pairs), fails[:1]), len(fails) > 0)


if __name__ == "__main__":
    L0, L1, L2 = Lat(r), Lat(w1), Lat(w2)
    link1()
    link2(L0, L1, L2)
    link3(L0, L1, L2)
    link4()
    print("ALL PASS" if all(checks) else "FAIL -- INVESTIGATE")
    sys.exit(0 if all(checks) else 1)
