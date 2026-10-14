"""
f7_structure_checks.py -- machine checks of every step of the two short proofs
used in the review, the "trivial explanation" tests for S = Lambda, the
basis dependence of failure counts, and the right-multiple convention variant.

Conventions: leech_final.py.

PROOF 1 (closure criterion).  Mod 2L write u = (x,y,z) in Lambda_g as
  w = x+y+z (in W), a = x+y, b = x+z (in V); then x = w+a+b, y = w+b, z = w+a.
  For ANY bilinear block product, in characteristic 2:
     P1+P2+P3 = w w',   P1+P2 = a a' + a b' + b a',   P1+P3 = a b' + b a' + b b'.
  Checked here by formal (symbolic) expansion over Z, reduced mod 2.
PROOF 2 (index formula).  R^24 = D + T, D = {(a,a,a)}, T = {p+q+r = 0}.
  (i)  D*T = T*D = 0 and (a,a,a)*(b,b,b) = 3(ab,ab,ab) for ANY bilinear block
       product (formal check);
  (ii) on T = O (x) T_Z the product is (x(x)alpha)*(y(x)beta) = xy (x) (alpha o beta)
       with the routing product o on T_Z = {alpha in Z^3 : sum = 0}, and
       [T_Z : <T_Z o T_Z>] = 3;
  (iii) Lambda cap D = Delta(Ls), Lambda cap T = Lsbar (x) T_Z,
       [Lambda : (Lambda cap D) + (Lambda cap T)] = 3^8 (a 2-adic unit);
  (iv) (2a,0,0)*(2b,0,0) = (4ab,0,0) and <L._g L> = L give (4L)^3 <= S_g.
  Hence [Lambda:S_g] = [Ls : <Ls._g Ls>] * [Lsbar : <Lsbar._g Lsbar>]^2.
TRIVIAL-EXPLANATION TESTS for S = Lambda under s_r:
  * no u in Lambda has a unimodular left or right multiplication: for every
    twist, det(L_u) = N(x+y+z)^4 * (integer), so 16 | det(L_u) or det = 0
    (checked on all basis vectors, 600 minimal vectors and 300 random vectors,
    for sigma and for s_r);
  * the span index is a ring-isomorphism invariant, so (Lambda,*_sigma) and
    (Lambda,*_r) are not isomorphic rings.
IDEAL PROPERTIES of V_g, W_g for sigma, the 3-cycle and s_r: the paper's
  Lemma 4.3 analogue (V_g a left ideal) holds for sigma but fails for the
  3-cycle and for s_r, both of which nevertheless close.
BASIS DEPENDENCE: failing-product counts of non-closing twists in five
  different Z-bases of Lambda (this folder's and the four agents'), versus the
  basis-independent obstruction dimension.
RIGHT-MULTIPLE VARIANT: Lambda_R built from sbar.l and s.l.
"""
import re
import random
from fractions import Fraction
from itertools import product, permutations
from collections import Counter
import numpy as np
from leech_final import (Setup, C_STD, T_STD, perm_map, perm_inverse, reflection_map, gamma_of,
                         twisted_constants, triple_products, lambda_coords, f2_rank, f2_span_basis,
                         z4_structure, hnf, hnf_index, int_det, bil, Lcoords, from_Lcoords, omul, dmul,
                         in_L, in_Lambda_def, S2, SBAR2, LB, wilson_lattice_basis, lattice_inverse2,
                         h_index, mod2_ok)

import os as _os
# the four independent implementations whose certificates and bases are re-checked here
REPRO = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "independent") + _os.sep
ST = Setup()
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


# ------------------------------------------------------------------ formal algebra
# An element is a dict {symbol: integer coefficient}; a product is a dict
# {(left symbol, right symbol): coefficient}.  Primes mark the right factor.
def lin(**kw):
    return {k: v for k, v in kw.items() if v}


def add(*xs):
    r = {}
    for x in xs:
        for k, v in x.items():
            r[k] = r.get(k, 0) + v
    return {k: v for k, v in r.items() if v}


def pr(x, y):
    r = {}
    for a, u in x.items():
        for b, v in y.items():
            r[(a, b)] = r.get((a, b), 0) + u * v
    return {k: v for k, v in r.items() if v}


def triple(u, v):
    x, y, z = u
    X, Y, Zp = v
    return (add(pr(x, X), pr(z, Y), pr(y, Zp)),
            add(pr(y, Y), pr(x, Zp), pr(z, X)),
            add(pr(z, Zp), pr(y, X), pr(x, Y)))


def mod2(d):
    return {k: v % 2 for k, v in d.items() if v % 2}


say("PROOF 1: closure criterion identities (formal, mod 2)")
w, a, b = lin(w=1), lin(a=1), lin(b=1)
u = (add(w, a, b), add(w, b), add(w, a))
v = (add(lin(W=1), lin(A=1), lin(B=1)), add(lin(W=1), lin(B=1)), add(lin(W=1), lin(A=1)))
P1, P2, P3 = triple(u, v)
ok1 = mod2(add(P1, P2, P3)) == {("w", "W"): 1}
ok2 = mod2(add(P1, P2)) == mod2(add(pr(a, lin(A=1)), pr(a, lin(B=1)), pr(b, lin(A=1))))
ok3 = mod2(add(P1, P3)) == mod2(add(pr(a, lin(B=1)), pr(b, lin(A=1)), pr(b, lin(B=1))))
say("  P1+P2+P3 = ww':", ok1, "  P1+P2 = aa'+ab'+ba':", ok2, "  P1+P3 = ab'+ba'+bb':", ok3)
assert ok1 and ok2 and ok3
say("  => closure iff W.W <= W and V.V <= V (take b = b' = 0, resp. a = b = 0, for necessity;")
say("     test vectors (v,0,v) and (w,w,w) lie in Lambda_g).")

say("")
say("PROOF 2: index formula")
# (i) D*T = T*D = 0, D*D
d1 = (lin(a=1), lin(a=1), lin(a=1))
t1 = (lin(p=1), lin(q=1), lin(p=-1, q=-1))
dd = triple(d1, (lin(B=1),) * 3)
dt = triple(d1, (lin(P=1), lin(Q=1), lin(P=-1, Q=-1)))
td = triple(t1, (lin(B=1),) * 3)
say("  (i)   D*T = 0:", all(not c for c in dt), "  T*D = 0:", all(not c for c in td),
    "  D*D = 3(ab,ab,ab):", all(c == {("a", "B"): 3} for c in dd))
assert all(not c for c in dt) and all(not c for c in td)


# (ii) routing product on T_Z
def route(al, be):
    return (al[0] * be[0] + al[2] * be[1] + al[1] * be[2],
            al[1] * be[1] + al[0] * be[2] + al[2] * be[0],
            al[2] * be[2] + al[1] * be[0] + al[0] * be[1])


tb = [(1, -1, 0), (0, 1, -1)]
prods = [route(x, y) for x in tb for y in tb]
assert all(sum(p) == 0 for p in prods)
co = [[p[0], -p[2]] for p in prods]        # coordinates in basis (1,-1,0),(0,1,-1): p = c1 t1 + c2 t2
idxT, _ = hnf_index(co, 2)
say("  (ii)  routing products of T_Z basis:", prods, "-> [T_Z : <T_Z o T_Z>] =", idxT)
assert idxT == 3
# (iii) Lambda cap D and Lambda cap T
D_gens = [list(wv) * 3 for wv in ST.Ls_L]
T_gens = [list(x) + [-t for t in x] + [0] * 8 for x in ST.Lsbar_L] + [[0] * 8 + list(x) + [-t for t in x] for x in ST.Lsbar_L]
allg = np.array(D_gens + T_gens, dtype=np.int64)
c, ok = lambda_coords(allg, ST.N)
assert ok.all()
iDT, _ = hnf_index(c.tolist(), 24)
say("  (iii) Delta(Ls) and Lsbar (x) T_Z lie in Lambda; [Lambda : Delta(Ls) + Lsbar(x)T_Z] =", iDT, "= 3^8:", iDT == 3 ** 8)
assert iDT == 3 ** 8
# maximality: Lambda cap D is exactly Delta(Ls) -- (a,a,a) in Lambda iff 3a in Ls iff a in Ls (2a in 2L <= Ls)
cnt = 0
for m in range(256):
    vv = [(m >> k) & 1 for k in range(8)]
    P = np.array([vv * 3], dtype=np.int64)
    inLam = mod2_ok(P, ST.cond)[0]
    inLs = f2_rank(np.vstack([ST.W0, np.array(vv)])) == 4
    cnt += (inLam == inLs)
say("        (a,a,a) in Lambda <=> a in Ls, on all 256 classes mod 2L:", cnt == 256)
say("  (iv)  (4L)^3 <= S_g needs <L._g L> = L, checked for every twist in f2/f4 (and [L:<L.L>] = 1 in f1).")

# ------------------------------------------------------------------ trivial explanation tests
say("")
say("TRIVIAL-EXPLANATION TESTS for S = Lambda under s_r")
rng = random.Random(5)
Be = ST.basis_e()


def mult_matrices(uL, Cg):
    """Matrices of v -> u*v and v -> v*u on Lambda (Lambda-coordinates, rows = images of b_q)."""
    U = np.array([uL], dtype=np.int64)
    Bm = ST.B
    def prodsUV(A, Bv):
        X, Y, Z = A[:, 0:8], A[:, 8:16], A[:, 16:24]
        Xp, Yp, Zp = Bv[:, 0:8], Bv[:, 8:16], Bv[:, 16:24]
        P1 = bil(X, Xp, Cg) + bil(Z, Yp, Cg) + bil(Y, Zp, Cg)
        P2 = bil(Y, Yp, Cg) + bil(X, Zp, Cg) + bil(Z, Xp, Cg)
        P3 = bil(Z, Zp, Cg) + bil(Y, Xp, Cg) + bil(X, Yp, Cg)
        return np.concatenate([P1, P2, P3], axis=2).reshape(-1, 24)
    Lm, ok1 = lambda_coords(prodsUV(U, Bm), ST.N)
    Rm, ok2 = lambda_coords(prodsUV(Bm, U), ST.N)
    assert ok1.all() and ok2.all()
    return Lm, Rm


# sample vectors of Lambda: basis, minimal vectors (types 1-3), random combinations
J = [[2 * s_ * int(i == t) for i in range(8)] for t in range(8) for s_ in (1, -1)]
roots = [list(X) for X in product(range(-2, 3), repeat=8) if sum(t * t for t in X) == 8 and in_L(X)]
samples = [list(map(int, r)) for r in ST.B]
for _ in range(200):
    lam = rng.choice(roots)
    j, k = rng.choice(J), rng.choice(J)
    t1v = [2 * t for t in lam] + [0] * 16
    ls = [p // 2 for p in omul(lam, SBAR2)]
    t2v = ls + dmul(ls, j) + [0] * 8
    t3v = dmul(dmul(lam, S2), j) + dmul(lam, k) + dmul(dmul(lam, j), k)
    for vv in (t1v, t2v, t3v):
        assert in_Lambda_def(vv)
        samples.append(sum((Lcoords(vv[8 * q:8 * q + 8]) for q in range(3)), []))
for _ in range(300):
    coef = [rng.randint(-2, 2) for _ in range(24)]
    samples.append([int(x) for x in np.array(coef) @ ST.B])
for name, fmap in (("sigma", perm_map((2, 1, 3, 4, 5, 6, 7))), ("s_r", reflection_map((1, 1, 1, 1, 1, 1, 1, -1)))):
    Cg = twisted_constants(gamma_of(fmap))
    stats = Counter()
    for uL in samples:
        s_e = [sum(from_Lcoords(uL[8 * q:8 * q + 8])[i] for q in range(3)) for i in range(8)]
        Nsum = Fraction(sum(t * t for t in s_e), 4)
        assert Nsum.denominator == 1 and Nsum % 2 == 0
        Lm, Rm = mult_matrices(uL, Cg)
        for M in (Lm, Rm):
            dt_ = int_det(M.tolist())
            if Nsum == 0:
                stats["x+y+z = 0, det = 0"] += (dt_ == 0)
                stats["x+y+z = 0"] += 1
            else:
                stats["N^4 | det"] += (dt_ % int(Nsum) ** 4 == 0)
                stats["nonzero cases"] += 1
            stats["unimodular"] += (abs(dt_) == 1)
    say("  %-5s: %d vectors x (left, right): %s" % (name, len(samples), dict(stats)))
    assert stats["unimodular"] == 0 and stats["N^4 | det"] == stats["nonzero cases"]
say("  => no element of Lambda acts unimodularly on either side (16 | det), so S = Lambda is not")
say("     produced by a single unit-like element; it needs genuinely many products.")

# ------------------------------------------------------------------ ideal properties
say("")
say("IDEAL PROPERTIES of V_g = g(Lsbar)/2L and W_g = g(Ls)/2L in L/2L (standard product mod 2)")


def ideal_props(Ub):
    U = np.array(Ub, dtype=np.int64) % 2
    E = np.eye(8, dtype=np.int64)
    r0 = f2_rank(U)
    inU = lambda P: f2_rank(np.vstack([U, P % 2])) == r0
    sub = inU(bil(U, U, C_STD).reshape(-1, 8))
    lft = inU(bil(E, U, C_STD).reshape(-1, 8))
    rgt = inU(bil(U, E, C_STD).reshape(-1, 8))
    return "subalgebra %s, left ideal %s, right ideal %s" % (sub, lft, rgt)


for name, fmap in (("sigma", perm_map((2, 1, 3, 4, 5, 6, 7))), ("3-cycle (5 6 7)", perm_map((1, 2, 3, 4, 6, 7, 5))),
                   ("s_r", reflection_map((1, 1, 1, 1, 1, 1, 1, -1)))):
    G = gamma_of(fmap)
    Vg = f2_span_basis((np.array(ST.V0) @ G) % 2)
    Wg = f2_span_basis((np.array(ST.W0) @ G) % 2)
    say("  %-16s V_g: %s" % (name, ideal_props(Vg)))
    say("  %-16s W_g: %s" % ("", ideal_props(Wg)))
say("  => for s_r the analogue of the paper's Lemma 4.3 (L . g(Lsbar) <= g(Lsbar)) FAILS, yet *_g closes:")
say("     closure needs only the two subalgebra properties (PROOF 1).")

# ------------------------------------------------------------------ basis dependence
say("")
say("BASIS DEPENDENCE of failing-product counts (non-closing twists)")


def load_e_rows(path, pat=None, scale=1, frac=False):
    rows = []
    txt = open(path).read()
    if pat:
        for m in re.finditer(pat, txt):
            vals = [Fraction(t.strip()) for t in m.group(1).split(",")]
            rows.append([int(scale * t) for t in vals])
    return rows


bases = {"this folder": Be}
bases["py"] = [[int(t) for t in ln.split("|L3|")[0].split()] for ln in open(REPRO + "py/lambda_basis.txt")
               if ln.strip() and not ln.startswith("#")]
gtxt = open(REPRO + "gap/cert_index1_reflection.txt").read()
gm = re.search(r"BASIS2\s*=\s*\[(.*?)\n\s*\]", gtxt, re.S)
bases["gap"] = [[int(x) for x in r.split(",")] for r in re.findall(r"\[([^\[\]]*)\]", gm.group(1))]
bases["structure"] = [[int(x) for x in ln.split()] for ln in open(REPRO + "structure/lambda_basis_ecoords_doubled.txt")
                      if ln.strip() and not ln.startswith("#")]
bases["refute"] = load_e_rows(REPRO + "refute/c3_witness.gp.out", r"b_\d+ = \[([^\]]*)\]", scale=2)
for k, bb in bases.items():
    assert len(bb) == 24 and all(in_Lambda_def(r) for r in bb), k
    Lc = [sum((Lcoords(r[8 * q:8 * q + 8]) for q in range(3)), []) for r in bb]
    c, ok = lambda_coords(np.array(Lc, dtype=np.int64), ST.N)
    assert abs(int_det(c.tolist())) == 1, k
    bases[k] = np.array(Lc, dtype=np.int64)
for name, pi in (("identity", (1, 2, 3, 4, 5, 6, 7)), ("inverse 3-cycle e5->e7->e6", (1, 2, 3, 4, 7, 5, 6))):
    Cg = twisted_constants(gamma_of(perm_map(pi)))
    row = []
    dims = set()
    for k, Bm in bases.items():
        P = triple_products(Bm, Cg)
        nf = int((~mod2_ok(P, ST.cond)).sum())
        dims.add(f2_rank((ST.cond @ (P.T % 2)) % 2))
        row.append("%s %d" % (k, nf))
    say("  %-28s failing basis products: %s;  obstruction dimension in every basis: %s" % (name, ", ".join(row), dims))

# ------------------------------------------------------------------ right-multiple variant
say("")
say("RIGHT-MULTIPLE VARIANT: Lambda_R from sbar.l and s.l (NOT the paper's convention)")
LsR = [Lcoords(dmul(S2, LB[k])) for k in range(8)]
LsbR = [Lcoords([p // 2 for p in omul(SBAR2, LB[k])]) for k in range(8)]
VR, WR = f2_span_basis(LsbR), f2_span_basis(LsR)
BR, condR = wilson_lattice_basis(VR, WR)
NR = lattice_inverse2(BR)
sameLam = hnf(BR.tolist(), 24) == hnf(ST.B.tolist(), 24)
say("  Lambda_R == Lambda:", sameLam, "; det(L^3-coords) =", abs(int_det(BR.tolist())))


def analyze_R(fmap):
    Cg = twisted_constants(gamma_of(fmap))
    P = triple_products(BR, Cg)
    c, ok = lambda_coords(P, NR)
    if not ok.all():
        return False, None
    a, b = z4_structure(c)
    return True, 2 ** (48 - 2 * a - b)


say("  sigma = (1 2):", analyze_R(perm_map((2, 1, 3, 4, 5, 6, 7))))
say("  s_r, r = (e0+...+e6-e7)/2:", analyze_R(reflection_map((1, 1, 1, 1, 1, 1, 1, -1))))
rres = Counter()
for X in product(range(-2, 3), repeat=8):
    if sum(t * t for t in X) == 8 and in_L(X) and next(t for t in X if t) > 0:
        rres[analyze_R(reflection_map(X))] += 1
say("  120 reflections on Lambda_R: (closes, index) counts:", dict(rres))
pres = Counter()
for pi in permutations(range(1, 8)):
    sgn = (-1) ** sum(1 for i in range(7) for j in range(i + 1, 7) if pi[i] > pi[j])
    cl, idx = analyze_R(perm_map(pi))
    if cl:
        pres[(sgn, idx)] += 1
say("  5040 permutations (convention A) on Lambda_R: closing (sign, index) counts:", dict(pres))
open("f7_structure_checks.out.txt", "w").write("\n".join(out) + "\n")
