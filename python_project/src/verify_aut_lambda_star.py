"""Task B, part 2: Aut(Lambda, +, star) -- exact order, structure, generators.

    Aut(Lambda, +, star) := { g : Lambda -> Lambda  Z-linear bijection
                              with  g(u star v) = g(u) star g(v) }.

The ambient is GL_24(Z), NOT Co_0; whether every star-automorphism is an
isometry of Lambda is one of the things settled here (it is -- see step 6).

INPUT FACTS (from verify_star_algebra_structure.py, all proved there):
  R^24 = D (+) T with D = {(a,a,a)}, T = {sum zero}, two-sided ideals with
  D star T = T star D = 0, Ann(D) = T, Ann(T) = D.  Every element of
  Aut(R^24, star) preserves D and T.  NOTE the route to that last statement
  (12 July 2026): it is NOT "Ann is intrinsic, so g permutes D and T" -- that
  is circular, since g(Ann(D)) = Ann(g(D)) gives g(T) = T only once g(D) = D
  is known.  It is: D and T are SIMPLE and annihilate each other, so pi_T . g
  restricted to D is an algebra map onto an ideal of T of dimension <= 8 < 16,
  hence 0; so g(D) <= D, and only then does g(T) = T follow.  See the
  CORRECTION block in verify_star_algebra_structure.py.  Consequently
        Aut(R^24, star) = { h_{A,B} rho_pi : A,B in Aut(O,.s), pi in S_3 },
        h_{A,B}(u)_alpha = B u_alpha + (A - B)(u_0+u_1+u_2)/3,
        rho_pi(u)_alpha  = u_{pi^{-1}(alpha)},
  with h_{A,B} rho_pi = rho_pi h_{A,B}, i.e. (G_2 x G_2) x S_3.

WHAT THIS SCRIPT COMPUTES (exact arithmetic throughout):

  (a) Lambda ^ D = {(a,a,a) : a in Ls}   and   Lambda ^ T = {(p,q,r) in
      (Ls_bar)^3 : p+q+r = 0}.                       [lattice computation]
      Hence  h_{A,B} rho_pi in Aut(Lambda,+,star)  =>  A(Ls) = Ls  and
      B(Ls_bar) = Ls_bar.  Both stabilisers are finite, so the group is
      FINITE.

  (b) COMPLETE enumerations of the finite groups
          Aut(O,.s) ^ Stab(L)      -- order  1344
          Aut(O,.s) ^ Stab(Ls)     -- order   168   (candidates for A)
          Aut(O,.s) ^ Stab(Ls_bar) -- order 12096   (candidates for B)
      Completeness: an automorphism of (O,.s) is orthogonal and is determined
      by the images of any three elements that generate O as an algebra; it
      permutes the minimal vectors of the stabilised lattice; so the search
      over image-triples of three generating minimal vectors, pruned by the
      Gram condition, is exhaustive.

  (c) K := {(A,B) : h_{A,B}(Lambda) = Lambda}, by TWO independent methods:
        method 1 -- glue: Lambda_0 := (Lambda^D) (+) (Lambda^T) has index
          3^8 = 6561 in Lambda; h_{A,B} preserves Lambda_0 automatically, so
          the condition is that the induced F_3-linear map preserves the
          8-dimensional glue code Lambda/Lambda_0 <= F_3^24.
        method 2 -- integrality: the matrix of h_{A,B} in a Z-basis of Lambda
          must be integral; test that directly, modulo the common
          denominator, with no reference to Lambda_0 or to F_3.
      Both give |K| = 6, and in every one of the 6 pairs A = B.

  (d) THE ANSWER.  Aut(Lambda,+,star) = K x S_3, of ORDER 36, and since
      A = B throughout, every element is
              (blockwise A) o (block permutation),   A in C_6.
      All 36 are re-verified one by one: they preserve Lambda and satisfy
      g(B_i star B_j) = g(B_i) star g(B_j) on all 576 basis pairs.

  (e) Co_0.  Every A in Aut(O,.s) is orthogonal (an algebra automorphism of
      the real octonions preserves the positive-definite norm form, because
      N(x) 1 = x x_bar and conjugation is determined by the algebra).  D and T
      are orthogonal complements in R^24 and h_{A,B} acts as A on D and
      blockwise-B on T, so h_{A,B} rho_pi is in O(24) -- for ALL A,B, not
      merely the lattice-preserving ones.  Hence
              Aut(Lambda,+,star) <= O(24) ^ GL(Lambda) = Co_0 .
      Confirmed here elementwise on the Gram matrix of Lambda.

A GAP file with the 36 matrices in the Lambda-basis is written for the
StructureDescription (see gap_project/aut_lambda_star.g).
"""
import os, sys, time
from fractions import Fraction as F
from itertools import combinations, product as iproduct
from sympy import Matrix, Rational, eye

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_consecutive_twists_exact import (
    omul, Lb, Lsb, Lsbarb, inL, inLs, inLsbar, inLambda, B as LAMBDA_BASIS)

# ---------------------------------------------------------------- octonions
SIGMA = (0,2,1,3,4,5,6,7)                     # sigma = (1 2)
def sig(v):
    r = [F(0)]*8
    for j in range(8): r[SIGMA[j]] = v[j]
    return r
def odot(a,b):  return sig(omul(sig(a), sig(b)))        # a .s b
def ip(a,b):    return sum(x*y for x,y in zip(a,b))
E  = [[F(1) if k==t else F(0) for k in range(8)] for t in range(8)]
Z8 = [F(0)]*8
H  = F(1,2)

# ------------------------------------------- minimal vectors of L, Ls, Lsbar
def _norm2_vectors():                          # the 240 roots of L (norm 2)
    out = []
    for i,j in combinations(range(8),2):
        for si in (1,-1):
            for sj in (1,-1):
                v = [F(0)]*8; v[i] = F(si); v[j] = F(sj); out.append(v)
    for sg in iproduct((H,-H), repeat=8):
        if sum(1 for x in sg if x < 0) % 2 == 1: out.append(list(sg))
    return [v for v in out if inL(v)]
def _norm4_vectors():                          # the 2160 norm-4 vectors of L
    out = []
    for pos in combinations(range(8),4):
        for sgn in iproduct((1,-1), repeat=4):
            v = [F(0)]*8
            for p,s0 in zip(pos,sgn): v[p] = F(s0)
            out.append(v)
    for i in range(8):
        for s0 in (2,-2):
            v = [F(0)]*8; v[i] = F(s0); out.append(v)
    for i in range(8):
        for s3 in (F(3,2), F(-3,2)):
            for rest in iproduct((H,-H), repeat=7):
                v = [F(0)]*8; v[i] = s3; k = 0
                for j in range(8):
                    if j != i: v[j] = rest[k]; k += 1
                out.append(v)
    return [v for v in out if inL(v) and ip(v,v) == 4]

MIN_L     = _norm2_vectors()
_N4       = _norm4_vectors()
MIN_LS    = [v for v in _N4 if inLs(v)]
MIN_LSBAR = [v for v in _N4 if inLsbar(v)]

# ------------------- complete enumeration of Aut(O,.s) ^ Stab(Gamma) --------
def _mat(vs):  return Matrix([[Rational(x.numerator, x.denominator) for x in v]
                              for v in vs]).T
def _mul8(X,Y): return [[sum(X[r][k]*Y[k][c] for k in range(8)) for c in range(8)]
                        for r in range(8)]
def _app(M,v):  return [sum(M[r][k]*v[k] for k in range(8)) for r in range(8)]

def octonion_autos_stabilising(minvecs, in_gamma, gamma_basis, start=0):
    """ALL A with A(x .s y) = Ax .s Ay and A(Gamma) = Gamma.

    Complete: A is orthogonal and fixes e_0, so it permutes minvecs; and A is
    determined by the images of m1,m2,m3 below, which generate (O,.s) as an
    algebra.  The three loops with the Gram filters therefore cover every
    possible A.  `start` selects the first generator m1: a different start
    gives a different search, and the answer must not change (used as an
    independent completeness check in verify_aut_octonion_crosscheck.py)."""
    m1 = minvecs[start]; cand = None
    for m2 in minvecs:
        S = [E[0], m1, m2, odot(m1,m2)]
        if _mat(S).rank() == 4:
            for m3 in minvecs:
                S2 = S + [m3, odot(m1,m3), odot(m2,m3), odot(odot(m1,m2), m3)]
                if _mat(S2).rank() == 8:
                    cand = (m1,m2,m3,S2); break
            if cand: break
    m1,m2,m3,VB = cand
    Vinv = _mat(VB).inv()
    Vi = [[F(int(Vinv[r,c].p), int(Vinv[r,c].q)) for c in range(8)] for r in range(8)]
    g  = [[ip(x,y) for y in (m1,m2,m3)] for x in (m1,m2,m3)]
    g0 = [ip(E[0],x) for x in (m1,m2,m3)]
    out = []
    for n1 in [n for n in minvecs if ip(E[0],n)==g0[0] and ip(n,n)==g[0][0]]:
        for n2 in [n for n in minvecs if ip(E[0],n)==g0[1] and ip(n,n)==g[1][1]
                                       and ip(n1,n)==g[0][1]]:
            for n3 in [n for n in minvecs if ip(E[0],n)==g0[2] and ip(n,n)==g[2][2]
                                           and ip(n1,n)==g[0][2] and ip(n2,n)==g[1][2]]:
                NB = [E[0], n1, n2, odot(n1,n2), n3, odot(n1,n3), odot(n2,n3),
                      odot(odot(n1,n2), n3)]
                Nm = [[NB[c][r] for c in range(8)] for r in range(8)]
                A  = _mul8(Nm, Vi)
                if not all(in_gamma(_app(A, bv)) for bv in gamma_basis): continue
                if not all(_app(A, odot(E[a],E[b])) == odot(_app(A,E[a]), _app(A,E[b]))
                           for a in range(8) for b in range(8)): continue
                out.append(A)
    return out

# ------------------------------------------------------- Lambda machinery
def flat(u):  return [c for blk in u for c in blk]
def Rn(x):    return Rational(x.numerator, x.denominator)
def blockwise(A):
    return lambda u: tuple(_app(A, blk) for blk in u)
def h_AB(A, Bm):
    def g(u):
        s  = [F(1,3)*(u[0][k]+u[1][k]+u[2][k]) for k in range(8)]
        As, Bs = _app(A, s), _app(Bm, s)
        return tuple([_app(Bm, u[al])[k] + As[k] - Bs[k] for k in range(8)]
                     for al in range(3))
    return g
def rho(perm):
    inv = [0]*3
    for a in range(3): inv[perm[a]] = a
    return lambda u: tuple(u[inv[a]] for a in range(3))
PERMS = [(0,1,2),(1,2,0),(2,0,1),(1,0,2),(0,2,1),(2,1,0)]

def star(u,v):
    x,y,z = u; xp,yp,zp = v
    def a3(p,q,r): return [p[k]+q[k]+r[k] for k in range(8)]
    return (a3(odot(x,xp), odot(z,yp), odot(y,zp)),
            a3(odot(y,yp), odot(x,zp), odot(z,xp)),
            a3(odot(z,zp), odot(y,xp), odot(x,yp)))

def rref3(vs, n):
    vs = [[int(x) % 3 for x in v] for v in vs]; piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, len(vs)) if vs[i][c] % 3), None)
        if pr is None: continue
        vs[r], vs[pr] = vs[pr], vs[r]
        iv = pow(vs[r][c], -1, 3); vs[r] = [(x*iv) % 3 for x in vs[r]]
        for i in range(len(vs)):
            if i != r and vs[i][c] % 3:
                f = vs[i][c] % 3
                vs[i] = [(vs[i][k]-f*vs[r][k]) % 3 for k in range(n)]
        piv.append(c); r += 1
        if r == len(vs): break
    return vs[:r], piv

def main():
    t0 = time.time()
    print("="*72); print("Aut(Lambda, +, star)"); print("="*72)

    # ---------------- (a) Lambda ^ D and Lambda ^ T
    LD = [(a,a,a) for a in Lsb]                                   # a in Ls
    LT = ([(p, Z8, [-c for c in p]) for p in Lsbarb] +
          [(Z8, q, [-c for c in q]) for q in Lsbarb])             # p,q in Lsbar
    print("(a) Lambda ^ D  = {(a,a,a) : a in Ls}      basis in Lambda:",
          all(inLambda(u) for u in LD))
    print("    Lambda ^ T  = {(p,q,r) in Lsbar^3, sum 0}  basis in Lambda:",
          all(inLambda(u) for u in LT))
    # Lambda ^ D exactly: (a,a,a) in Lambda  <=>  a in L, 2a in Lsbar, 3a in Ls.
    # 2L <= Ls <= L, so the three conditions are constant on the 2^8 classes of
    # L/2L; test all 256 of them.  Exactly the 16 classes of Ls survive.
    print("    2L <= Ls:", all(inLs([2*c for c in v]) for v in Lb),
          "  2L <= Lsbar:", all(inLsbar([2*c for c in v]) for v in Lb),
          "  (so the three conditions are constant on the classes of L/2L)")
    hits = 0
    for cvec in iproduct((0,1), repeat=8):
        a = [F(0)]*8
        for k,ck in enumerate(cvec):
            if ck: a = [a[i]+Lb[k][i] for i in range(8)]
        if (inL(a) and inLsbar([2*c for c in a]) and inLs([3*c for c in a])):
            hits += 1
            if not inLs(a): print("    *** class outside Ls satisfies the conditions")
    print("    of the 256 classes of L/2L, the ones with (a,a,a) in Lambda:", hits,
          "= |Ls/2L| => Lambda ^ D = {(a,a,a) : a in Ls} EXACTLY")
    BLam = Matrix([[Rn(c) for c in flat(u)] for u in LAMBDA_BASIS]).T
    B0   = Matrix([[Rn(c) for c in flat(u)] for u in LD+LT]).T
    idx  = abs(B0.det()/BLam.det())
    print("    [Lambda : (Lambda^D) (+) (Lambda^T)] =", idx, "= 3^8")

    # ---------------- (b) the three octonion-automorphism stabilisers
    SL = octonion_autos_stabilising(MIN_L,     inL,      Lb)
    SA = octonion_autos_stabilising(MIN_LS,    inLs,     Lsb)
    SB = octonion_autos_stabilising(MIN_LSBAR, inLsbar,  Lsbarb)
    print("(b) |Aut(O,.s) ^ Stab(L)|      =", len(SL))
    print("    |Aut(O,.s) ^ Stab(Ls)|     =", len(SA), "  <- candidates for A")
    print("    |Aut(O,.s) ^ Stab(Lsbar)|  =", len(SB), "  <- candidates for B")
    def kk(A): return tuple(tuple(r) for r in A)
    sL, sA, sB = {kk(x) for x in SL}, {kk(x) for x in SA}, {kk(x) for x in SB}
    print("    |Stab(L) ^ Stab(Ls) ^ Stab(Lsbar)| =", len(sL & sA & sB),
          "  (= the A with blockwise A preserving Lambda)")

    # ---------------- (c) method 1: the mod-3 glue code
    import numpy as np
    Cm = B0.inv()*BLam
    Gcols = [[int(3*Cm[r,c]) % 3 for r in range(24)] for c in range(24)]
    Cb, piv = rref3(Gcols, 24)
    G = np.array(Cb, dtype=np.int64).T                       # 24 x 8
    free = [c for c in range(24) if c not in piv]
    Hs = []
    for fc in free:
        v = [0]*24; v[fc] = 1
        for i,pc in enumerate(piv): v[pc] = (-Cb[i][fc]) % 3
        Hs.append(v)
    Hm = np.array(Hs, dtype=np.int64)                        # 16 x 24, Hm G = 0 mod 3
    assert G.shape == (24,8) and Hm.shape == (16,24)
    assert not (Hm.dot(G) % 3).any()
    MLs   = _mat(Lsb);    MLsI  = MLs.inv()
    MLsb  = _mat(Lsbarb); MLsbI = MLsb.inv()
    def m8(A): return Matrix([[Rn(A[r][c]) for c in range(8)] for r in range(8)])
    def mod3(M): return np.array([[int(M[r,c]) % 3 for c in range(8)]
                                  for r in range(8)], dtype=np.int64)
    H1, H2 = Hm[:, :8], Hm[:, 8:]
    G1, G2 = G[:8, :],  G[8:, :]
    UA = {}
    for i,A in enumerate(SA):
        M = mod3(MLsI*m8(A)*MLs)
        UA.setdefault(((H1.dot(M).dot(G1)) % 3).tobytes(), []).append(i)
    VB = {}
    for j,Bm in enumerate(SB):
        M = mod3(MLsbI*m8(Bm)*MLsb)
        MM = np.zeros((16,16), dtype=np.int64); MM[:8,:8] = M; MM[8:,8:] = M
        VB.setdefault(((-H2.dot(MM).dot(G2)) % 3).tobytes(), []).append(j)
    K1 = [(i,j) for k,ia in UA.items() if k in VB for i in ia for j in VB[k]]
    print("(c) method 1 (mod-3 glue code):  |K| =", len(K1))

    # ---------------- (c) method 2: integrality of the matrix of h in the Lambda basis
    Kinv = BLam.inv()
    dd = 1
    for e in Kinv: dd = dd*Rational(e).q // __import__('math').gcd(dd, Rational(e).q)
    Ki = np.array([[int(Kinv[r,c]*dd) for c in range(24)] for r in range(24)],
                  dtype=object)
    D2 = np.array([[int(2*BLam[r,c]) for c in range(24)] for r in range(24)],
                  dtype=object)
    m  = 24*dd
    J  = np.ones((3,3), dtype=object); I3 = np.eye(3, dtype=object)
    def kron(M3, A8):
        A = np.array([[F(A8[r][c]) for c in range(8)] for r in range(8)], dtype=object)
        out = np.zeros((24,24), dtype=object)
        for a in range(3):
            for b in range(3):
                if M3[a][b] == 0: continue
                out[8*a:8*a+8, 8*b:8*b+8] = M3[a][b]*A
        return out
    def to_int(M):
        return np.array([[int(x) for x in row] for row in M], dtype=object)
    FA, GB = {}, {}
    for i,A in enumerate(SA):                      # 12h = 2J(x)A2 + (6I-2J)(x)B2
        M = to_int(2*kron([[1,1,1]]*3, [[2*x for x in r] for r in A]))
        FA.setdefault(((Ki.dot(M).dot(D2)) % m).tobytes(), []).append(i)
    for j,Bm in enumerate(SB):
        C3 = [[6-2 if a==b else -2 for b in range(3)] for a in range(3)]
        M = to_int(kron(C3, [[2*x for x in r] for r in Bm]))
        GB.setdefault(((-Ki.dot(M).dot(D2)) % m).tobytes(), []).append(j)
    K2 = [(i,j) for k,ia in FA.items() if k in GB for i in ia for j in GB[k]]
    print("    method 2 (integrality in the Lambda basis, mod", m, "):  |K| =", len(K2))
    print("    the two methods agree on the same pair set:",
          sorted(K1) == sorted(K2))
    print("    in every surviving pair A = B:", all(SA[i] == SB[j] for i,j in K1))

    # ---------------- (d) the 36 elements
    K0 = [SA[i] for i,j in K1]
    els = [(A, p) for A in K0 for p in PERMS]
    def gmap(A,p):
        b, r = blockwise(A), rho(p)
        return lambda u: b(r(u))
    okL = okS = okG = 0
    def ipv(u,v): return sum(a*b for bu,bv in zip(u,v) for a,b in zip(bu,bv))
    for (A,p) in els:
        g = gmap(A,p)
        if all(inLambda(g(b)) for b in LAMBDA_BASIS): okL += 1
        if all(flat(g(star(bi,bj))) == flat(star(g(bi),g(bj)))
               for bi in LAMBDA_BASIS for bj in LAMBDA_BASIS): okS += 1
        if all(ipv(g(bi),g(bj)) == ipv(bi,bj)
               for bi in LAMBDA_BASIS for bj in LAMBDA_BASIS): okG += 1
    sigs = {tuple(tuple(flat(gmap(A,p)(b))) for b in LAMBDA_BASIS) for (A,p) in els}
    print("(d) |Aut(Lambda,+,star)| = |K| * |S_3| =", len(K1), "*", 6, "=", len(els))
    print("    all", len(els), "preserve Lambda:", okL == len(els),
          "  preserve star:", okS == len(els))
    print("    all", len(els), "are distinct linear maps:", len(sigs) == len(els))
    print("(e) all", len(els), "preserve the Leech inner product => inside Co_0:",
          okG == len(els))

    # element orders of K0 (a group of order 6)
    def order(A):
        Y, n = A, 1
        while kk(Y) != kk([[F(1) if r==c else F(0) for c in range(8)] for r in range(8)]):
            Y = _mul8(Y, A); n += 1
        return n
    print("    element orders in K (the octonion part):", sorted(order(A) for A in K0),
          "=> K = C_6")
    print("    K is abelian:", all(kk(_mul8(X,Y)) == kk(_mul8(Y,X)) for X in K0 for Y in K0))

    # ---------------- GAP output: 24x24 integer matrices in the Lambda basis
    BinvL = BLam.inv()
    def mat24(g):
        Mg = Matrix([[Rn(c) for c in flat(g(b))] for b in LAMBDA_BASIS]).T
        return BinvL*Mg
    gens = []
    Agen = max(K0, key=order)                       # an order-6 octonion automorphism
    for g in [gmap(Agen,(0,1,2)), gmap(K0[0],(1,2,0)), gmap(K0[0],(1,0,2))]:
        gens.append(mat24(g))
    for Mg in gens:
        assert all(Rational(e).q == 1 for e in Mg)
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "gap_project", "aut_lambda_star_gens.g")
    with open(os.path.abspath(path), "w") as f:
        f.write("# generators of Aut(Lambda,+,star) as 24x24 integer matrices\n")
        f.write("# in the Z-basis of Lambda used by verify_aut_lambda_star.py\n")
        f.write("gens := [\n")
        for Mg in gens:
            f.write("[" + ",".join("[" + ",".join(str(int(Mg[r,c])) for c in range(24))
                                   + "]" for r in range(24)) + "],\n")
        f.write("];\n")
    print("    GAP generators written to gap_project/aut_lambda_star_gens.g")
    print("    (blockwise C_6 generator; block 3-cycle; block transposition)")
    print("\nruntime %.0f s" % (time.time()-t0))


if __name__ == "__main__":
    main()
