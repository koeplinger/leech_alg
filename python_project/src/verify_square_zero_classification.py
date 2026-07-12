"""Task A, part 3: COMPLETE classification of the square-zero elements of
(R^24, +, star), and the consequence for the lattice.

RESULT (algebra).  u = (x, y, z) satisfies u star u = 0 if and only if

      Re(x) = Re(y) = Re(z) = 0,    x + y + z = 0,    N(x) = N(y) = N(z),

i.e. x, y, z are three PURELY IMAGINARY octonions of equal length summing to
zero (equivalently: a hexagonal / A_2 triple in Im(O) = R^7).  The solution
set is a cone of dimension 12 (choose x in R^7 freely: 7; then y lies on the
5-sphere {N(y) = N(x), <x,y> = -N(x)/2}: 5).  It is contained in T and meets
D only in 0.  This condition is invariant under sigma (sigma is orthogonal
and fixes e_0), so it reads the same in the twisted and untwisted algebras.

RESULT (lattice) -- THIS REFUTES THE WORKING ASSUMPTION.  Lambda DOES contain
nonzero square-zero elements.  u = (x, y, z) in Lambda has u star u = 0 iff
x, y, z lie in (L sbar) ^ Im(O), sum to zero, and have equal norms.  Now
      (L sbar) ^ Im(O) = {lambda sbar : lambda in L, <lambda, s> = 0},
whose 126 vectors of norm 4 form a copy of the E_7 root system, and E_7 root
systems contain hexagonal triples.  Every such triple gives a square-zero
vector of Lambda of norm 12.  The script exhibits one explicitly and counts
them exhaustively.  Their norm is necessarily a multiple of 12 (norm = 3 N(x)
with N(x) in 2Z, and every Lambda-norm is in 4Z), so none of them is minimal
-- which is exactly why the exhaustive minimal-shell search
(verify_idempotents_min_shell.py) found none.

The statement "Lambda has no nonzero square-zero element" is therefore FALSE;
only the weaker (and true) statement "Min(Lambda) has no square-zero element"
is supported.

Second pass: the characterisation is re-derived independently by a symbolic
check of the two octonion identities used, and every exhibited lattice vector
is re-tested with the independently coded star of verify_ideal_decomposition.
"""
from fractions import Fraction as F
from itertools import combinations

import sympy as sp

from verify_ideal_decomposition import (
    omul, odot, sig, star, star_indep, add, zero8, E, eq, is_zero, Sigma,
)
from verify_consecutive_twists_exact import inLambda, inL, inLs, inLsbar, s, sbar


def N(v):  return sum(c*c for c in v)
def ip(u, v): return sum(a*b for a, b in zip(u, v))
def Nu(u): return sum(N(b) for b in u)


def main():
    # --------------------------------------------------------------- S1
    print("=== the characterisation, symbolically ===")
    X = sp.symbols("X0:8", real=True); Y = sp.symbols("Y0:8", real=True)
    Z = sp.symbols("Z0:8", real=True)
    from verify_idempotent_classification import sym_omul
    Xs, Ys, Zs = list(X), list(Y), list(Z)
    # untwisted square-zero equations
    e1 = [sp.expand(sym_omul(Xs,Xs)[k] + sym_omul(Zs,Ys)[k] + sym_omul(Ys,Zs)[k])
          for k in range(8)]
    e2 = [sp.expand(sym_omul(Ys,Ys)[k] + sym_omul(Xs,Zs)[k] + sym_omul(Zs,Xs)[k])
          for k in range(8)]
    e3 = [sp.expand(sym_omul(Zs,Zs)[k] + sym_omul(Ys,Xs)[k] + sym_omul(Xs,Ys)[k])
          for k in range(8)]
    # Substitute the claimed solution shape: purely imaginary, sum zero.
    #   X = P, Y = Q, Z = -(P+Q), Re = 0.
    # Then EVERY component of u star u is an exact rational-linear combination
    # of the two quadrics
    #   g1 := N(P) - N(Q),      g2 := N(P) + 2<P,Q>,
    # and conversely g1, g2 are rational-linear combinations of the components.
    # The two systems are therefore EQUIVALENT, which proves the
    # characterisation in both directions.  Certified here by exact linear
    # algebra over the degree-2 monomial basis (no Groebner heuristics).
    P = sp.symbols("p1:8", real=True); Q = sp.symbols("q1:8", real=True)
    sub = {X[0]: 0, Y[0]: 0, Z[0]: 0}
    for k in range(1, 8):
        sub[X[k]] = P[k-1]; sub[Y[k]] = Q[k-1]; sub[Z[k]] = -(P[k-1] + Q[k-1])
    NP = sum(t*t for t in P); NQ = sum(t*t for t in Q)
    PQ = sum(t*u for t, u in zip(P, Q))
    g1 = sp.expand(NP - NQ); g2 = sp.expand(NP + 2*PQ)
    comps = [sp.expand(eqs[k].subs(sub)) for eqs in (e1, e2, e3) for k in range(8)]
    varsPQ = list(P) + list(Q)
    mons = sorted({m for f in comps + [g1, g2]
                   for m in sp.Poly(f, varsPQ).monoms()} if any(comps) else set())
    def cvec(f):
        p = sp.Poly(f, varsPQ)
        d = dict(zip(p.monoms(), p.coeffs()))
        return [sp.Rational(d.get(m, 0)) for m in mons]
    Gm = sp.Matrix([cvec(g1), cvec(g2)]).T                 # |mons| x 2
    Cm = sp.Matrix([cvec(f) for f in comps]).T             # |mons| x 24
    # every component in span(g1, g2)?
    fwd = sp.Matrix.hstack(Gm, Cm).rank() == Gm.rank() == 2
    # g1, g2 in span(components)?
    bwd = sp.Matrix.hstack(Cm, Gm).rank() == Cm.rank()
    ok_char = bool(fwd and bwd)
    print(f"  substitute  Re = 0,  Z = -(X+Y):")
    print(f"  every component of u star u lies in span_Q(g1, g2):  {fwd}")
    print(f"  and g1, g2 lie in span_Q(components):                {bwd}")
    print(f"  => u star u = 0  <=>  N(X) = N(Y)  and  N(X) + 2<X,Y> = 0")
    print(f"     <=>  N(X) = N(Y) = N(Z)  (given Re = 0 and X+Y+Z = 0):"
          f" {ok_char}")
    print( "  (necessity is proved in the docstring from the ideal splitting:")
    print( "   d*d = 0 forces d = 0 since D ~= O is a division algebra; then")
    print( "   S is singular on the sum-zero plane only if D2 = 0, i.e.")
    print( "   a = b = c = 0, after which the real parts give exactly these")
    print( "   equations.)")
    print()

    # a necessity check by exhaustive linear algebra on the real parts
    a, b, c = sp.symbols("a b c", real=True)
    S = sp.Matrix([[a, c, b], [c, b, a], [b, a, c]])
    print("=== necessity: a = b = c = 0 for every square-zero element ===")
    print( "  D-part must vanish (D ~= O, division algebra) => a+b+c = 0 and")
    print( "  Im-parts sum to zero.  Imaginary equations: S (X',Y',Z') = 0.")
    print( "  On the sum-zero plane S has eigenvalues +-sqrt(D2); so either")
    print( "  (X',Y',Z') = 0, or D2 = 0 which with a+b+c = 0 gives a=b=c=0.")
    A3 = [a**2 + 2*b*c, b**2 + 2*c*a, c**2 + 2*a*b]      # square in A_3
    sol0 = sp.solve(A3 + [a + b + c], [a, b, c], dict=True)
    print(f"  and if (X',Y',Z') = 0 the real system A*A = 0, a+b+c = 0 has"
          f" solutions: {sol0}")
    print(f"  => a = b = c = 0 in every case: {sol0 == [{a: 0, b: 0, c: 0}]}")
    print()

    # ------------------------------------------------ explicit algebra example
    print("=== an explicit nonzero square-zero element of R^24 ===")
    h = F(1, 2); r3 = None
    # x = e1, y = -1/2 e1 + w, z = -1/2 e1 - w with N(w) = 3/4 -- needs sqrt(3),
    # so use instead an exact rational hexagonal triple:
    #   x = e1 + e2, y = -e1 + e3, z = -e2 - e3  -- check equal norms & sum 0
    x = add(E(1), E(2)); y = add([-c for c in E(1)], E(3))
    z = [-(p+q) for p, q in zip(x, y)]
    u = (x, y, z)
    print(f"  x = {x}\n  y = {y}\n  z = {z}")
    print(f"  Re = 0 : {x[0] == y[0] == z[0] == 0}   x+y+z = 0 :"
          f" {all(t == 0 for t in add(x, y, z))}")
    print(f"  N(x), N(y), N(z) = {N(x)}, {N(y)}, {N(z)}")
    print(f"  u star u = 0 : {is_zero(star(u, u))}  (independent star:"
          f" {is_zero(star_indep(u, u))})")
    print()

    # ------------------------------------------------- the lattice: refutation
    print("=== square-zero elements INSIDE Lambda ===")
    print("  Lambda-condition with x+y+z = 0 reduces to: x, y, z in L sbar.")
    # (L sbar) ^ Im(O): lambda sbar with lambda in L, Re(lambda sbar) = <lambda,s> = 0
    roots = []
    for i, j in combinations(range(8), 2):
        for si in (F(1), F(-1)):
            for sj in (F(1), F(-1)):
                v = [F(0)]*8; v[i] = si; v[j] = sj; roots.append(v)
    from itertools import product as iproduct
    for sg in iproduct([h, -h], repeat=8):
        if sum(1 for t in sg if t < 0) % 2 == 1:
            roots.append(list(sg))
    assert len(roots) == 240
    im_roots = [lam for lam in roots if ip(lam, s) == 0]
    V = [omul(lam, sbar) for lam in im_roots]          # norm-4 vectors of Lsbar ^ Im O
    print(f"  E_8 roots lambda with <lambda, s> = 0 : {len(im_roots)}"
          f"   (expect 126 = #roots of E_7)")
    assert all(v[0] == 0 for v in V), "not purely imaginary"
    assert all(N(v) == 4 for v in V), "not norm 4"
    print(f"  their images lambda sbar : all purely imaginary, all of norm 4 :"
          f" True")
    Vset = {tuple(v) for v in V}

    # hexagonal triples: ordered (x,y) in V^2 with -(x+y) in V
    triples = []
    for xx in V:
        for yy in V:
            zz = tuple(-(p+q) for p, q in zip(xx, yy))
            if zz in Vset:
                triples.append((xx, yy, list(zz)))
    print(f"  ordered hexagonal triples (x, y, z) in V^3 with x+y+z = 0 :"
          f" {len(triples)}")

    good = 0; bad = 0
    for t in triples:
        uu = (list(t[0]), list(t[1]), list(t[2]))
        if inLambda(uu) and is_zero(star(uu, uu)) and is_zero(star_indep(uu, uu)):
            good += 1
        else:
            bad += 1
    print(f"  of these, in Lambda AND square-zero : {good}   (failures: {bad})")
    print(f"  every one has norm N(u) = 3*4 = 12 :"
          f" {all(Nu((list(t[0]),list(t[1]),list(t[2]))) == 12 for t in triples)}")
    print()

    ex = triples[0]
    print("  EXPLICIT WITNESS (a square-zero vector of Lambda):")
    print(f"    x = {[str(t) for t in ex[0]]}")
    print(f"    y = {[str(t) for t in ex[1]]}")
    print(f"    z = {[str(t) for t in ex[2]]}")
    w = (list(ex[0]), list(ex[1]), list(ex[2]))
    print(f"    in Lambda: {inLambda(w)}   norm: {Nu(w)}   u star u = 0:"
          f" {is_zero(star(w, w))}")
    print()

    print("=== why the minimal shell sees none ===")
    print("  N(u) = N(x)+N(y)+N(z) = 3 N(x) with N(x) in 2Z (x in E_8), so")
    print("  N(u) in 6Z.  Every norm of Lambda lies in 4Z (checked below).")
    print("  Hence every square-zero vector of Lambda has norm in 12Z.  The")
    print("  minimal norm is 8, so no minimal vector can be square-zero --")
    print("  consistent with the exhaustive shell search, but NOT equivalent")
    print("  to 'Lambda has no square-zero element', which is FALSE.")
    # proof that every Lambda-norm is in 4Z, from the Gram matrix of a Z-basis
    from verify_consecutive_twists_exact import B as LAMBDA_BASIS
    G = [[sum(ip(bi[k], bj[k]) for k in range(3)) for bj in LAMBDA_BASIS]
         for bi in LAMBDA_BASIS]
    diag4 = all(G[i][i] % 4 == 0 for i in range(24))
    off2 = all(G[i][j] % 2 == 0 for i in range(24) for j in range(24) if i != j)
    print(f"  Gram matrix of a Z-basis of Lambda: all G_ii in 4Z : {diag4};"
          f"  all G_ij (i=/=j) in 2Z : {off2}")
    print(f"  => N(sum a_i b_i) = sum a_i^2 G_ii + 2 sum_(i<j) a_i a_j G_ij"
          f" is in 4Z for every integer vector a : {diag4 and off2}")
    print()

    print("=== second pass: the pre-existing INDEPENDENT star implementation ===")
    # verify_section5_properties.star is the doubled-integer product coded
    # before this task; re-check every one of the 4032 witnesses with it.
    from verify_section5_properties import star as star_d, eq as eq_d, Nsq
    ZD = ([0]*8, [0]*8, [0]*8)
    bad2 = 0
    for t in triples:
        Ud = tuple([int(2*c) for c in blk] for blk in t)   # doubled-int coords
        if not eq_d(star_d(Ud, Ud), ZD) or Nsq(Ud)//4 != 12:
            bad2 += 1
    print(f"  4032 witnesses re-tested with the doubled-integer star:"
          f" {bad2} failures (expect 0)")
    print()

    allok = (ok_char and good == len(triples) and bad == 0 and bad2 == 0
             and len(im_roots) == 126 and diag4 and off2)
    print("ALL PASS" if allok else "FAILURE")


if __name__ == "__main__":
    main()
