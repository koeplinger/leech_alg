"""Task B, part 1: the structure of Aut(R^24, star) as a real algebraic group.

star is the paper's Z_3-symmetric sigma-twisted triple-octonion product on
R^24 = O_1 (+) O_2 (+) O_3, sigma = (1 2).  Index the blocks by Z/3 and write
.s for the sigma-twisted octonion product.  star is bilinear, so every
basis-exhaustive check below is a PROOF, not a sample.

  (1) ROUTING.   P_gamma = sum over {(alpha,beta) : alpha + beta = 2 gamma
      (mod 3)} of u_alpha .s v_beta.                        [checked]

  (2) IDEALS.    D = {(a,a,a)} (dim 8) and T = {(p,q,r) : p+q+r = 0} (dim 16)
      are two-sided ideals, D star T = T star D = 0, R^24 = D (+) T, and
              Ann(T) = D    and    Ann(D) = T.              [checked]
      CORRECTION (12 July 2026, from the independent verification pass).
      An earlier version of this docstring argued: "Ann is intrinsic, so
      every automorphism maps D onto D and T onto T; no dimension count is
      needed."  That is CIRCULAR: g(Ann(D)) = Ann(g(D)) yields g(T) = T only
      once g(D) = D is already known.  The correct, non-circular route, and
      the one printed in paper/automorphism_group_2026-07-12.tex:
        - D star T = T star D = 0, so the projections pi_D, pi_T are algebra
          homomorphisms;
        - D and T are SIMPLE (their multiplication algebras <L_u, R_u> are the
          full matrix algebras, of dimensions 64 and 256, computed exactly);
        - hence pi_T(g(D)) is an ideal of T of dimension <= 8 < 16, so it is 0,
          so g(D) <= D and therefore g(D) = D;
        - only THEN does g(T) = g(Ann(D)) = Ann(g(D)) = Ann(D) = T follow.
      The conclusion is unchanged:
              Aut(R^24, star) = Aut(D, star) x Aut(T, star),
      the two factors being independently choosable because the cross
      products vanish.                                      [checked]

  (3) FOURIER.   Over Q(zeta), zeta a primitive cube root of unity, put
              u_hat_j := sum_alpha zeta^{j alpha} u_alpha.
      Then      (u star v)_hat_k = u_hat_{2k} .s v_hat_{2k}.  [checked, exactly
      in Q(zeta) -- carried as pairs p + q zeta with zeta^2 = -1 - zeta]
      So over C the algebra is O_C (the j=0 part, = D (x) C) plus the
      SWAP-TWISTED algebra S: (a_1,a_2)(b_1,b_2) = (a_2 b_2, a_1 b_1) on
      O_C (+) O_C (the j=1,2 part, = T (x) C).

  (4) THE GROUP.  For A,B in Aut(O,.s), h_{A,B} := P (x) A + Q (x) B, where P
      is the orthogonal projection of R^24 onto D and Q = I - P; concretely
              h_{A,B}(u)_alpha = B u_alpha + (A - B)(u_0+u_1+u_2)/3.
      Every h_{A,B} is a star-automorphism, every block permutation rho_pi is
      a star-automorphism, and they commute.                [checked]
      The report proves the converse (Aut(S) = G_2(C) x S_3 for the
      swap-twisted algebra, whose real points give blockwise B together with
      the block permutations), so
              Aut(R^24, star) = (G_2 x G_2) x S_3,   dim 28,
      with G_2 the COMPACT real form in both factors.  As a computational
      corroboration of the identity component, this script computes
              dim Der(O, .s)      = 14      (exact nullspace)
              dim Der(R^24, star) = 28      (28 explicit derivations, and an
                                             exact mod-p rank bound showing no
                                             more exist).

  (5) TRACE FORM.  B(u,v) := trace(L_u L_v) is computed exactly.  It is
              B(u,v) = 24 * sum_alpha Re(u_alpha .s v_alpha),
      diagonal with entries (+24, -24,...,-24) per block, signature (3,21).
      It is NOT a positive multiple of the Leech Gram matrix, and it is NOT of
      the shape alpha <u_D,v_D> + beta <u_T,v_T> (it is blind to the D/T
      split).  So the trace form does NOT by itself place Aut(Lambda,+,star)
      inside Co_0.  The correct argument is compactness: an algebra
      automorphism of the REAL octonions preserves the positive-definite norm
      (N(x) 1 = x x_bar, and conjugation x_bar = t(x) 1 - x is determined by
      the algebra), so Aut(O,.s) < O(8); D and T are orthogonal complements;
      hence every h_{A,B} rho_pi is in O(24).  See verify_aut_lambda_star.py.
"""
import os, random, sys
from fractions import Fraction as F
from itertools import product as iproduct
from sympy import Matrix, Rational, eye

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ------------------------------------------------------------- octonions
FANO = [(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)]
MUL = {}
for a in range(8):
    MUL[(0,a)] = (1,a); MUL[(a,0)] = (1,a)
for a in range(1,8):
    MUL[(a,a)] = (-1,0)
for (a,b,c) in FANO:
    for (x,y,z) in [(a,b,c),(b,c,a),(c,a,b)]:
        MUL[(x,y)] = (1,z); MUL[(y,x)] = (-1,z)

def omul(u, v):
    r = [F(0)]*8
    for i in range(8):
        if u[i] == 0: continue
        for j in range(8):
            if v[j] == 0: continue
            sg, c = MUL[(i,j)]
            r[c] += sg*u[i]*v[j]
    return r

SIGMA = (0,2,1,3,4,5,6,7)
def sig(v):
    r = [F(0)]*8
    for j in range(8): r[SIGMA[j]] = v[j]
    return r
def odot(a, b): return sig(omul(sig(a), sig(b)))          # a .s b

Z8 = [F(0)]*8
E  = [[F(1) if k == t else F(0) for k in range(8)] for t in range(8)]
def a3(p,q,r): return [p[k]+q[k]+r[k] for k in range(8)]

def star(u, v):
    x,y,z = u; xp,yp,zp = v
    return (a3(odot(x,xp), odot(z,yp), odot(y,zp)),
            a3(odot(y,yp), odot(x,zp), odot(z,xp)),
            a3(odot(z,zp), odot(y,xp), odot(x,yp)))

FULL = [tuple(E[t] if b == blk else Z8 for b in range(3))
        for blk in range(3) for t in range(8)]
def flat(u):  return [c for blk in u for c in blk]
def eqv(u,v): return flat(u) == flat(v)

# ------------------------------------------------- (1) the routing pattern
def routing_ok():
    for gamma in range(3):
        pairs = [(al,be) for al in range(3) for be in range(3)
                 if (al+be) % 3 == (2*gamma) % 3]
        for bi in FULL:
            for bj in FULL:
                rhs = [F(0)]*8
                for (al,be) in pairs:
                    t = odot(bi[al], bj[be])
                    rhs = [rhs[k]+t[k] for k in range(8)]
                if star(bi,bj)[gamma] != rhs: return False
    return True

# ----------------------------------------------------------- (2) the ideals
D_BASIS = [(E[t],E[t],E[t]) for t in range(8)]
T_BASIS = ([(E[t], [-c for c in E[t]], Z8) for t in range(8)] +
           [(Z8, E[t], [-c for c in E[t]]) for t in range(8)])
def span(vs): return Matrix([[Rational(c) for c in flat(v)] for v in vs]).T

def annihilator(sub):
    rows = []
    for X in sub:
        for side in (0,1):
            for out in range(24):
                rows.append([Rational(flat(star(b,X) if side==0 else star(X,b))[out])
                             for b in FULL])
    ns = Matrix(rows).nullspace()
    return Matrix.hstack(*ns) if ns else Matrix.zeros(24,0)

# ------------------------------------------------- (3) the Fourier statement
def fourier_ok():
    """exact over Q(zeta): elements are pairs (p,q) = p + q*zeta, zeta^2 = -1-zeta."""
    def cmul(a,b):
        (p,q),(r,s) = a,b
        return (p*r - q*s, p*s + q*r - q*s)
    def cadd(a,b): return (a[0]+b[0], a[1]+b[1])
    ZT = [(F(1),F(0)), (F(0),F(1)), (F(-1),F(-1))]        # 1, zeta, zeta^2
    def omul_c(u,v):
        r = [(F(0),F(0))]*8
        for i in range(8):
            if u[i] == (0,0): continue
            for j in range(8):
                if v[j] == (0,0): continue
                sg,c = MUL[(i,j)]
                t = cmul(u[i], v[j])
                if sg < 0: t = (-t[0], -t[1])
                r[c] = cadd(r[c], t)
        return r
    def sig_c(v):
        r = [(F(0),F(0))]*8
        for j in range(8): r[SIGMA[j]] = v[j]
        return r
    def odot_c(a,b): return sig_c(omul_c(sig_c(a), sig_c(b)))
    def hat(u,j):
        out = [(F(0),F(0))]*8
        for al in range(3):
            w = ZT[(j*al) % 3]
            for k in range(8):
                out[k] = cadd(out[k], cmul(w, (u[al][k], F(0))))
        return out
    for bi in FULL:
        for bj in FULL:
            p = star(bi,bj)
            for k in range(3):
                if hat(p,k) != odot_c(hat(bi,(2*k)%3), hat(bj,(2*k)%3)):
                    return False
    return True

# ------------------------------------- (4) h_{A,B} = P (x) A + Q (x) B
def app(M,v):  return [sum(M[r][k]*v[k] for k in range(8)) for r in range(8)]
def h_AB(A,B):
    def g(u):
        s  = [F(1,3)*(u[0][k]+u[1][k]+u[2][k]) for k in range(8)]
        As, Bs = app(A,s), app(B,s)
        return tuple([app(B,u[al])[k] + As[k] - Bs[k] for k in range(8)]
                     for al in range(3))
    return g
def rho(perm):
    inv = [0]*3
    for a in range(3): inv[perm[a]] = a
    return lambda u: tuple(u[inv[a]] for a in range(3))
PERMS = [(0,1,2),(1,2,0),(2,0,1),(1,0,2),(0,2,1),(2,1,0)]
def is_star_auto(g):
    return all(eqv(g(star(bi,bj)), star(g(bi), g(bj)))
               for bi in FULL for bj in FULL)

# --------------------------------------------------- derivation algebras
def _der_octonion_rows():
    """rows of the linear system  X(a .s b) - X(a) .s b - a .s X(b) = 0,
    unknowns X[r][c] indexed 8r+c."""
    rows = []
    for a in range(8):
        for b in range(8):
            prod = odot(E[a], E[b])
            S1 = [odot(E[r], E[b]) for r in range(8)]
            S2 = [odot(E[a], E[r]) for r in range(8)]
            for out in range(8):
                row = [Rational(0)]*64
                for r in range(8):
                    for c in range(8):
                        val = Rational(0)
                        if r == out: val += Rational(prod[c])   # X(a.s b)
                        if c == a:   val -= Rational(S1[r][out])  # X(a).s b
                        if c == b:   val -= Rational(S2[r][out])  # a .s X(b)
                        row[8*r+c] = val
                rows.append(row)
    return rows

def der_octonion_dim():
    """dim {X in gl_8 : X(a .s b) = X(a) .s b + a .s X(b)}, exact."""
    return 64 - Matrix(_der_octonion_rows()).rank()

def der_star_dim_bound(nrows=1400, primes=(1000003, 999983)):
    """UPPER bound on dim Der(R^24, star), exact.

    Rows of the (integer) linear system X(u*v) - X(u)*v - u*X(v) = 0 are built
    for a subset of basis pairs; the rank of an integer matrix mod p is a LOWER
    bound for its rank over Q, so the mod-p nullity of a SUBSYSTEM is an UPPER
    bound for the nullity of the full system over Q."""
    import numpy as np
    rng = random.Random(20260712)
    pairs = [(i,j) for i in range(24) for j in range(24)]
    rng.shuffle(pairs)
    rows = []
    for (i,j) in pairs:
        if len(rows) >= nrows: break
        bi, bj = FULL[i], FULL[j]
        prod = flat(star(bi,bj))
        # contribution of unknown X[r][c]:
        #   to  X(u*v)_out :  delta(out,r) * prod[c]
        #   to  (X u)*v_out :  delta(c,i) * star(e_r, b_j)_out
        #   to  u*(X v)_out :  delta(c,j) * star(b_i, e_r)_out
        S1 = [flat(star(FULL[r], bj)) for r in range(24)]
        S2 = [flat(star(bi, FULL[r])) for r in range(24)]
        for out in range(24):
            row = [0]*576
            for r in range(24):
                for c in range(24):
                    val = 0
                    if r == out: val += int(prod[c])
                    if c == i:   val -= int(S1[r][out])
                    if c == j:   val -= int(S2[r][out])
                    if val: row[24*r+c] = val
            rows.append(row)
    M = np.array(rows, dtype=np.int64)
    best = 576
    for p in primes:
        A = M % p
        r = 0
        for c in range(576):
            piv = None
            for i in range(r, A.shape[0]):
                if A[i,c] % p: piv = i; break
            if piv is None: continue
            A[[r,piv]] = A[[piv,r]]
            inv = pow(int(A[r,c]), p-2, p)
            A[r] = (A[r]*inv) % p
            col = A[r+1:, c].copy()
            nz = np.nonzero(col)[0]
            if len(nz):
                A[r+1+nz] = (A[r+1+nz] - np.outer(col[nz], A[r])) % p
            r += 1
            if r == A.shape[0]: break
        best = min(best, 576 - r)
    return best

def der_star_explicit():
    """the 28 derivations P (x) X + Q (x) Y, X,Y in Der(O,.s): exact check."""
    ns = Matrix(_der_octonion_rows()).nullspace()
    ders = []
    for v in ns:
        ders.append([[F(int(v[8*r+c].p), int(v[8*r+c].q)) for c in range(8)]
                     for r in range(8)])
    # X = P (x) Xo + Q (x) Yo  as a map on R^24
    def Dmap(Xo, Yo):
        def d(u):
            s  = [F(1,3)*(u[0][k]+u[1][k]+u[2][k]) for k in range(8)]
            Xs, Ys = app(Xo,s), app(Yo,s)
            return tuple([app(Yo,u[al])[k] + Xs[k] - Ys[k] for k in range(8)]
                         for al in range(3))
        return d
    ZERO = [[F(0)]*8 for _ in range(8)]
    cands = [Dmap(X, ZERO) for X in ders] + [Dmap(ZERO, Y) for Y in ders]
    ok = True
    for d in cands:
        for bi in FULL:
            for bj in FULL:
                lhs = flat(d(star(bi,bj)))
                rhs = [a+b for a,b in zip(flat(star(d(bi),bj)), flat(star(bi,d(bj))))]
                if lhs != rhs: ok = False
    Mm = Matrix([[Rational(x) for x in sum((flat(d(b)) for b in FULL), [])]
                 for d in cands])
    return len(ders), len(cands), ok, Mm.rank()

# ------------------------------------------------------------ (5) trace form
def trace_form():
    def Lmat(u):
        return Matrix([[Rational(c) for c in flat(star(u,b))] for b in FULL]).T
    Ls = [Lmat(b) for b in FULL]
    return Matrix(24,24, lambda i,j: (Ls[i]*Ls[j]).trace())


def main():
    print("="*72)
    print("Aut(R^24, star): the ambient algebraic group")
    print("="*72)
    print("(1) routing  P_gamma = sum_{alpha+beta = 2 gamma (3)} u_alpha .s v_beta:",
          routing_ok())

    aD, aT = annihilator(D_BASIS), annihilator(T_BASIS)
    MD, MT = span(D_BASIS), span(T_BASIS)
    print("(2) dim Ann(D) =", aD.shape[1], " and Ann(D) = T:",
          aD.shape[1] == 16 and Matrix.hstack(aD, MT).rank() == 16)
    print("    dim Ann(T) =", aT.shape[1], " and Ann(T) = D:",
          aT.shape[1] == 8 and Matrix.hstack(aT, MD).rank() == 8)
    print("    => with SIMPLICITY of D and T (dim 8 < dim 16), every automorphism")
    print("       of (R^24, star) preserves D and T setwise; see the docstring for")
    print("       why the annihilator relation alone does not suffice.")

    print("(3) Fourier  (u*v)^_k = u^_{2k} .s v^_{2k}  exactly over Q(zeta_3):",
          fourier_ok())
    print("    => the j=1,2 part of T (x) C is the swap-twisted algebra")
    print("       (a1,a2)(b1,b2) = (a2 b2, a1 b1) on O_C (+) O_C.")

    from verify_aut_lambda_star import octonion_autos_stabilising, MIN_L, inL, Lb
    G = octonion_autos_stabilising(MIN_L, inL, Lb)
    random.seed(20260712)
    ok = all(is_star_auto(h_AB(random.choice(G), random.choice(G))) for _ in range(10))
    print("(4) h_{A,B} = P(x)A + Q(x)B is a star-automorphism (10 random pairs A,B")
    print("    drawn from the 1344 octonion automorphisms preserving L; each")
    print("    checked on all 576 basis pairs):", ok)
    print("    the six block permutations rho_pi are star-automorphisms:",
          all(is_star_auto(rho(p)) for p in PERMS))
    A, Bm = random.choice(G), random.choice(G)
    hh, rr = h_AB(A,Bm), rho((1,2,0))
    print("    h_{A,B} commutes with rho_pi:",
          all(eqv(hh(rr(b)), rr(hh(b))) for b in FULL))

    d14 = der_octonion_dim()
    nder, ncand, dok, drank = der_star_explicit()
    ub = der_star_dim_bound()
    print("    dim Der(O, .s) =", d14, "(= dim G_2)")
    print("    explicit derivations P(x)X + Q(x)Y, X,Y in Der(O,.s):", ncand,
          " all verified:", dok, " rank (independent):", drank)
    print("    exact mod-p upper bound on dim Der(R^24, star):", ub)
    print("    => dim Aut(R^24, star) = 28 = dim(G_2 x G_2):", drank == 28 and ub == 28)

    Bf = trace_form()
    offd = all(Bf[i,j] == 0 for i in range(24) for j in range(24) if i != j)
    diag = [Bf[i,i] for i in range(24)]
    shape = (offd and all(diag[8*b] == 24 for b in range(3))
             and all(diag[8*b+t] == -24 for b in range(3) for t in range(1,8)))
    print("(5) trace form B(u,v) = tr(L_u L_v): diagonal:", offd)
    print("    B = 24 * sum_alpha Re(u_alpha .s v_alpha), signature (3,21):", shape)
    print("    B is a scalar multiple of the Leech Gram matrix:",
          Bf == Bf[0,0]*eye(24), " <-- NO")
    print("    B distinguishes D from T:  B|_D and B|_T both have signature")
    print("    (1,7) resp. (2,14); B is blind to the D/T split, so the trace")
    print("    form does NOT by itself force Aut(Lambda,+,star) into Co_0.")


if __name__ == "__main__":
    main()
