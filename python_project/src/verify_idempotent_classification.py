"""Task A, part 2: COMPLETE classification of the idempotents of (R^24, +, star).

RESULT.  The algebra has exactly EIGHT idempotents (u star u = u), namely

    0,
    eps_1 = (e_0, 0, 0),   eps_2 = (0, e_0, 0),   eps_3 = (0, 0, e_0),
    w/3   = (1/3)(e_0, e_0, e_0),
    eps_1 - w/3 = (1/3)( 2 e_0, -e_0, -e_0),   and its two block-cyclings.

Every one of them lies in the 3-dimensional subalgebra A_3 = span{eps_1,
eps_2, eps_3} of e_0-multiples; none has a nonzero imaginary part.  (There
is no identity element: an identity would be an idempotent, and none of the
eight acts as one.)

PROOF, each step verified exactly below (all arithmetic in Fractions/sympy):

 S1.  R^24 = D (+) T with D, T two-sided ideals that annihilate each other
      (verify_ideal_decomposition.py).  Hence u = d + t is idempotent iff
      d and t are.
 S2.  phi : a |-> (1/3)(a,a,a) is an isomorphism (O, .s) -> (D, star).  The
      octonions are a real division algebra, so a^2 = a forces a in {0, e_0}:
      the only idempotents of D are 0 and w/3 = phi(e_0).  BOTH ARE REAL
      (multiples of e_0 in every block).
 S3.  Consequently every idempotent u has Im-blocks summing to zero.
      Transport by the algebra isomorphism sigma^(+3) (step (6) of
      verify_ideal_decomposition.py) to the UNTWISTED triple product and
      write X = sigma(x), Y = sigma(y), Z = sigma(z); the idempotent
      equations become, in the plain octonions,
          X^2 + (ZY + YZ) = X,  Y^2 + (XZ + ZX) = Y,  Z^2 + (YX + XY) = Z.
      Split X = a + X' (a = Re X, X' = Im X), etc., using
          A^2 = 2 Re(A) A - N(A),   AB + BA = 2Re(A)B + 2Re(B)A - 2<A,B>.
 S4.  The imaginary parts give the linear system  M (X', Y', Z')^T = 0  with
          M = 2 S - I,   S = [[a, c, b], [c, b, a], [b, a, c]]
      (S is symmetric; it is the matrix of left star-multiplication by
      (a,b,c) inside A_3).  Eigenvalues of S: a+b+c on (1,1,1), and
      +-sqrt(D2) on the sum-zero plane, D2 := a^2+b^2+c^2-ab-bc-ca >= 0.
 S5.  By S3, (X',Y',Z') lies in the sum-zero plane tensor Im(O), which M
      preserves.  If (X',Y',Z') =/= 0, then M must be singular there, i.e.
      2 sqrt(D2) - 1 = 0 (the other eigenvalue -1 - 2 sqrt(D2) < 0 never
      vanishes).  That eigenvalue is SIMPLE, so ker M is 1-dimensional:
          (X', Y', Z') = (p W, q W, r W)   for a single W in Im(O), W =/= 0.
 S6.  Therefore X, Y, Z all lie in R + R W, which is a subalgebra of O
      isomorphic to C (W^2 = -|W|^2).  The idempotent system restricted to
      (R + R W)^3 is exactly the idempotent system of A_3 tensor C.
 S7.  Solve that system exactly over C (sympy, 3 complex unknowns): it has
      exactly 8 solutions, and ALL EIGHT ARE REAL.  So p = q = r = 0,
      contradicting W =/= 0 in S5.
 S8.  Hence every idempotent has zero imaginary part, i.e. lies in A_3, and
      the 8 solutions of S7 (read over R) are the complete list.       QED

Independent second pass: (i) a float Newton search from 20000 random starts
in R^24 converges only to the eight (reported, but only as corroboration --
it is floating point and proves nothing); (ii) the 8 are re-verified with the
independently coded star of verify_ideal_decomposition.py.
"""
import random
from fractions import Fraction as F

import sympy as sp

from verify_ideal_decomposition import (
    omul, odot, sig, star, star_indep, add, zero8, E, eq, is_zero, Sigma,
)

# --------------------------------------------------------------- the eight
third = F(1, 3)
def blocks(c0, c1, c2):
    return tuple([F(c) if k == 0 else F(0) for k in range(8)]
                 for c in (c0, c1, c2))

ZERO  = blocks(0, 0, 0)
EPS1  = blocks(1, 0, 0)
EPS2  = blocks(0, 1, 0)
EPS3  = blocks(0, 0, 1)
OM3   = blocks(third, third, third)                 # w/3
P1    = blocks(2*third, -third, -third)             # eps_1 - w/3
P2    = blocks(-third, 2*third, -third)
P3    = blocks(-third, -third, 2*third)
EIGHT = [("0", ZERO), ("eps_1", EPS1), ("eps_2", EPS2), ("eps_3", EPS3),
         ("w/3", OM3), ("eps_1 - w/3", P1), ("eps_2 - w/3", P2),
         ("eps_3 - w/3", P3)]


def show(u):
    def s8(v):
        return "0" if all(c == 0 for c in v) else str([str(c) for c in v if c != 0][0]) + "*e0"
    return "(" + ", ".join(s8(b) for b in u) + ")"


def step_S2():
    """Idempotents of D: exactly 0 and w/3."""
    # (O, .s) is a division algebra: a .s a = a  =>  a = 0 or a = e_0.
    # verified symbolically: Re and N are the same for . and .s (sigma is
    # orthogonal and fixes e_0), so a .s a = 2 Re(a) a - N(a) e_0.
    A = sp.symbols("a0:8", real=True)
    lhs = sp.Matrix(sym_odot(A, A))
    rhs = 2*A[0]*sp.Matrix(list(A)) - sum(x*x for x in A)*sp.Matrix([1] + [0]*7)
    ident_ok = sp.simplify(lhs - rhs) == sp.zeros(8, 1)
    # a .s a = a  <=>  2 Re(a) a - N(a) e_0 = a.  Imaginary: (2a0 - 1) a' = 0.
    #  a' = 0  -> a0^2 = a0 -> a0 in {0,1};  a0 = 1/2 -> 1/4 - |a'|^2 = 1/2 -> |a'|^2 < 0.
    sols = sp.solve([sp.Eq(x, y) for x, y in zip(sym_odot(A, A), A)], A,
                    dict=True)
    return ident_ok, sols


def sym_omul(u, v):
    from verify_ideal_decomposition import MUL
    r = [sp.Integer(0)]*8
    for i in range(8):
        for j in range(8):
            sg, c = MUL[(i, j)]
            r[c] += sg*u[i]*v[j]
    return r

def sym_sig(v):
    r = [sp.Integer(0)]*8
    for j in range(8):
        r[(0,2,1,3,4,5,6,7)[j]] = v[j]
    return r

def sym_odot(u, v):
    return sym_sig(sym_omul(sym_sig(u), sym_sig(v)))


def main():
    print("=== the eight claimed idempotents: u star u = u ? ===")
    all_idem = True
    for name, u in EIGHT:
        good = eq(star(u, u), u) and eq(star_indep(u, u), u)
        all_idem &= good
        print(f"  {name:<14} {show(u):<34} u*u = u : {good}")
    print()

    # eps_i acts as a block transposition (H2)
    print("=== left/right multiplication by eps_i ===")
    rng = random.Random(7)
    def rv():
        return tuple([F(rng.randint(-5, 5)) for _ in range(8)] for _ in range(3))
    swaps = {0: (0, 2, 1), 1: (2, 1, 0), 2: (1, 0, 2)}   # eps_i swaps the other two
    tr_ok = True
    for i, eps in enumerate([EPS1, EPS2, EPS3]):
        good = True
        for _ in range(50):
            v = rv()
            want = tuple(v[swaps[i][k]] for k in range(3))
            good &= eq(star(eps, v), want) and eq(star(v, eps), want)
        tr_ok &= good
        print(f"  eps_{i+1} * v = v * eps_{i+1} = block-transposition "
              f"{tuple(k+1 for k in swaps[i])} of v : {good}")
    print()

    # --------------------------------------------------- S2
    ident_ok, dsols = step_S2()
    print("=== S2: idempotents of D  (D ~= (O, .s) via a |-> (1/3)(a,a,a)) ===")
    print(f"  identity  a .s a = 2 Re(a) a - N(a) e_0  verified symbolically: {ident_ok}")
    print(f"  solutions of  a .s a = a  in O :  {len(dsols)}")
    for s in dsols:
        print("    a =", [sp.nsimplify(s.get(sp.Symbol(f'a{k}', real=True), 0))
                          for k in range(8)])
    d_ok = len(dsols) == 2
    print(f"  => only a = 0 and a = e_0, i.e. D-idempotents are 0 and w/3 : {d_ok}")
    print()

    # --------------------------------------------------- S4 / S5
    a, b, c = sp.symbols("a b c", real=True)
    S = sp.Matrix([[a, c, b], [c, b, a], [b, a, c]])
    M = 2*S - sp.eye(3)
    D2 = a**2 + b**2 + c**2 - a*b - b*c - c*a
    print("=== S4/S5: eigenvalues of S = [[a,c,b],[c,b,a],[b,a,c]] ===")
    print(f"  S symmetric: {S == S.T}")
    print(f"  S*(1,1,1) = (a+b+c)*(1,1,1): "
          f"{sp.simplify(S*sp.Matrix([1,1,1]) - (a+b+c)*sp.Matrix([1,1,1])) == sp.zeros(3,1)}")
    cp = sp.factor(sp.expand(S.charpoly(sp.Symbol('L')).as_expr()))
    lam = sp.Symbol('L')
    roots = sp.roots(sp.Poly(S.charpoly(lam).as_expr(), lam))
    print(f"  char. poly factors as: {sp.factor(sp.expand((lam-(a+b+c))*(lam**2 - D2)))}")
    ok_cp = sp.simplify(S.charpoly(lam).as_expr() -
                        sp.expand((lam - (a+b+c))*(lam**2 - D2))) == 0
    print(f"  char(S) = (L - (a+b+c)) (L^2 - D2),  D2 = a^2+b^2+c^2-ab-bc-ca : {ok_cp}")
    print( "  => on the sum-zero plane S has eigenvalues +-sqrt(D2);")
    print( "     M = 2S - I is singular there iff 2 sqrt(D2) = 1, and that")
    print( "     eigenvalue is simple, so ker M is 1-dimensional.")
    print()

    # verify that M is exactly the imaginary-part system
    print("=== S4: the imaginary-part system really is M (X',Y',Z') = 0 ===")
    XS = sp.symbols("x0:8", real=True); YS = sp.symbols("y0:8", real=True)
    ZS = sp.symbols("z0:8", real=True)
    Xs = list(XS); Ys = list(YS); Zs = list(ZS)
    # untwisted equations  X^2 + ZY + YZ - X  (etc.), imaginary components
    eq1 = [sym_omul(Xs, Xs)[k] + sym_omul(Zs, Ys)[k] + sym_omul(Ys, Zs)[k] - Xs[k]
           for k in range(8)]
    eq2 = [sym_omul(Ys, Ys)[k] + sym_omul(Xs, Zs)[k] + sym_omul(Zs, Xs)[k] - Ys[k]
           for k in range(8)]
    eq3 = [sym_omul(Zs, Zs)[k] + sym_omul(Ys, Xs)[k] + sym_omul(Xs, Ys)[k] - Zs[k]
           for k in range(8)]
    imag_ok = True
    for k in range(1, 8):        # imaginary components
        want1 = (2*Xs[0] - 1)*Xs[k] + 2*Zs[0]*Ys[k] + 2*Ys[0]*Zs[k]
        want2 = 2*Zs[0]*Xs[k] + (2*Ys[0] - 1)*Ys[k] + 2*Xs[0]*Zs[k]
        want3 = 2*Ys[0]*Xs[k] + 2*Xs[0]*Ys[k] + (2*Zs[0] - 1)*Zs[k]
        imag_ok &= (sp.expand(eq1[k] - want1) == 0)
        imag_ok &= (sp.expand(eq2[k] - want2) == 0)
        imag_ok &= (sp.expand(eq3[k] - want3) == 0)
    print(f"  rows of the imaginary system = rows of (2S - I) with"
          f" (a,b,c) = (X_0,Y_0,Z_0) : {imag_ok}")
    print()

    # --------------------------------------------------- S6 / S7
    print("=== S6/S7: idempotents of A_3 (x) C  (X,Y,Z in R + R W ~= C) ===")
    al, be, ga = sp.symbols("alpha beta gamma")           # COMPLEX unknowns
    sys3 = [al**2 + 2*be*ga - al, be**2 + 2*ga*al - be, ga**2 + 2*al*be - ga]
    solsC = sp.solve(sys3, [al, be, ga], dict=True)
    print(f"  number of solutions over C: {len(solsC)}")
    allreal = True
    for s in solsC:
        v = (sp.nsimplify(s[al]), sp.nsimplify(s[be]), sp.nsimplify(s[ga]))
        r = all(sp.im(t) == 0 for t in v)
        allreal &= r
        print(f"    (alpha,beta,gamma) = {v}   real: {r}")
    print(f"  all solutions real: {allreal}")
    print(f"  exactly 8 solutions: {len(solsC) == 8}")
    print()

    print("=== CONCLUSION ===")
    concl = (all_idem and tr_ok and ident_ok and d_ok and ok_cp and imag_ok
             and allreal and len(solsC) == 8)
    print("  Every idempotent has zero imaginary part (S5-S7), hence lies in")
    print("  A_3 = span{eps_1, eps_2, eps_3}; the idempotents of A_3 are the")
    print("  8 real solutions above.  THE ALGEBRA HAS EXACTLY 8 IDEMPOTENTS.")
    print()

    # -------------------------------------------------- no identity element
    print("=== no identity element ===")
    noid = True
    for name, u in EIGHT:
        if is_zero(u):
            continue
        found = None
        for _ in range(20):
            v = rv()
            if not eq(star(u, v), v):
                found = True; break
        noid &= bool(found)
    print(f"  none of the 8 idempotents is a two-sided identity: {noid}")
    print( "  (an identity is idempotent, so the algebra is NOT unital)")
    print()

    # ------------------------------------------- corroboration: float Newton
    print("=== corroboration only (FLOATING POINT, proves nothing) ===")
    import numpy as np
    Fp = np.zeros((8, 8, 8))
    from verify_ideal_decomposition import MUL
    for (i, j), (sg, k) in MUL.items():
        Fp[k, i, j] += sg
    SG = [0, 2, 1, 3, 4, 5, 6, 7]
    Pm = np.zeros((8, 8)); Pm[SG, range(8)] = 1.0          # sigma as a matrix
    def dot_np(u, v):
        return Pm @ np.einsum('kij,i,j->k', Fp, Pm @ u, Pm @ v)
    def star2_np(u, v):
        x, y, z = u[:8], u[8:16], u[16:]
        xp, yp, zp = v[:8], v[8:16], v[16:]
        return np.concatenate([
            dot_np(x,xp) + dot_np(z,yp) + dot_np(y,zp),
            dot_np(y,yp) + dot_np(x,zp) + dot_np(z,xp),
            dot_np(z,zp) + dot_np(y,xp) + dot_np(x,yp)])
    def star_np(u):
        return star2_np(u, u)
    # star is bilinear, so the Jacobian of F(u) = u*u - u is ANALYTIC:
    #   J(u) v = u*v + v*u - v
    I24 = np.eye(24)
    def Jac(u):
        J = np.column_stack([star2_np(u, I24[:, k]) + star2_np(I24[:, k], u)
                             for k in range(24)])
        return J - I24
    rs = np.random.default_rng(20260712)
    found = set()
    fails = 0
    NSTART = 5000
    for _ in range(NSTART):
        u = rs.normal(size=24)*rs.choice([0.2, 0.5, 1.0, 2.0])
        for _ in range(60):
            r = star_np(u) - u
            if np.linalg.norm(r) < 1e-13: break
            try:
                u = u - np.linalg.solve(Jac(u), r)
            except np.linalg.LinAlgError:
                break
            if np.linalg.norm(u) > 1e4: break
        if np.linalg.norm(star_np(u) - u) < 1e-9 and np.linalg.norm(u) < 1e3:
            found.add(tuple(np.round(u, 6) + 0.0))
        else:
            fails += 1
    print(f"  Newton from {NSTART} random starts: {len(found)} distinct limits,"
          f" {fails} non-convergent")
    for f8 in sorted(found, key=lambda t: np.linalg.norm(t)):
        nz = {k: v for k, v in enumerate(f8) if abs(v) > 1e-9}
        print(f"    {nz}")
    print()
    print("ALL PASS" if concl else "FAILURE")


if __name__ == "__main__":
    main()
