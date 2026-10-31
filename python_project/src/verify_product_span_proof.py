"""The proof of the span of the triple product (paper Section 5.6,
Proposition 5.4, label prop:span), step by step, in exact arithmetic.  The
paper cites this script at Proposition 5.4, with verify_product_span_index.py
and verify_product_span_structure.py.  Notation as in the paper: L = D8+ with the
octonion product of (2.1), s and sbar of Section 2.3, sigma = (1 2), x.y
the sigma-twisted product (Definition 3.2), star the triple product
(Definition 3.3), Lambda Wilson's lattice (Definition 2.1), MO =
(1/2)(e1+e2)L the maximal order of (2.2), p = e2 - e4 + e5 - e6 of (4.1).
As in the paper's proof, juxtaposition is the octonion product of (2.1)
in (1) and (2), and the twisted product x.y in (4) and (5).

Created: 2026-10-23.  Updated: 2026-10-28.

  (1) The products (yp)z with y, z in MO span MO, and so do the products
      y(pz), by hand: v := -(1/2)(e1+e2)s = (1/2)(e0+e2+e3-e7) lies in MO,
      pv + vp = 2 v_0 p - 2<p, v> e0 = p - e0, and the alternative laws
      give y = (yp)e0 - (yp)v - ((yv)p)e0 and y = e0(py) - v(py) -
      e0(p(vy)) for every y in MO (checked on a Z-basis b_1, ..., b_8 of
      MO; both sides are linear in y).  As a cross-check, the 64 products
      (b_i p) b_j span MO (index 1), and so do the b_i (p b_j); the
      conjugate of (yp)z is -zbar (p ybar), and MO is closed under
      conjugation.  For contrast: p has norm 4 and is not in 2 MO; the
      norm is integral on MO, so every product (ap)b, of norm
      4 N(a) N(b), has norm divisible by 4; and the products ap and pb
      (a, b in MO) span a sublattice of index 64.
  (2) The sublattice facts, on basis pairs: Lsbar.Lsbar lies in 2 Lsbar;
      Lsbar.Ls and Ls.Lsbar lie in Lsbar; Ls.Ls lies in Ls; Ls.Ls spans
      Ls, and Lsbar.Ls spans Lsbar.  Also (4.1) again: sigma(Lsbar) =
      2 MO and sigma(Ls) = p MO.
  (3) The vectors alpha(a) = (a, 0, -a), beta(a) = (a, -a, 0) (a in a
      Z-basis of Lsbar) and gamma(w) = (-w, w, w) (w in a Z-basis of Ls)
      lie in Lambda with determinant 2^12 = [L^3 : Lambda]: a Z-basis.
  (4) The nine kinds of products of these vectors, on all 576 pairs,
      writing (k1, k2, k3) c for (k1 c, k2 c, k3 c) and xy for x.y:
          alpha(a) * alpha(a') = (1, -2, 1) aa'
          beta(b)  * beta(b')  = (1, 1, -2) bb'
          alpha(a) * beta(b)   = (2, -1, -1) ab,   beta(b) * alpha(a) = (2, -1, -1) ba
          alpha(a) * gamma(w)  = (-2, 2, 0) aw,    gamma(w) * alpha(a) = (-2, 2, 0) wa
          beta(b)  * gamma(w)  = (-2, 0, 2) bw,    gamma(w) * beta(b)  = (-2, 0, 2) wb
          gamma(w) * gamma(w') = (3, -1, -1) ww'
  (5) For c in Lsbar, (c, -c, 0) and (c, 0, -c) lie in Lambda, hence
      every (k1, k2, k3) c with integers k1 + k2 + k3 = 0; the first eight
      kinds lie in 2 Lambda; gamma(w) * gamma(w') = (4c, 0, 0) - (c, c, c) with c =
      ww', where (2c, 0, 0) and (c, c, c) lie in Lambda.
  (6) With Delta = {(w, w, w) : w in Ls}: the span S of the 576 products
      contains Delta and 2 Lambda and lies in Delta + 2 Lambda, so S =
      Delta + 2 Lambda; [Lambda : S] = 2^16 and [Delta + 2 Lambda :
      2 Lambda] = 2^8, that is, Delta cap 2 Lambda = 2 Delta.
  (7) The cyclic shift (x, y, z) -> (z, x, y) maps Lambda to itself and
      commutes with star (on the basis of (3)); and star agrees, on the same
      basis pairs, with the triple product that verify_consecutive_twists_exact.py
      builds for sigma = (1 2).

Run: python3 verify_product_span_proof.py  ->  ALL PASS (exit 0).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from sympy import Matrix
from verify_consecutive_twists_exact import (omul, Lb, Lsb, Lsbarb, Lcoords,
                                             inL, inLs, inLsbar, inLambda,
                                             l3, twisted, detB)

h = F(1, 2)
def vec(*c): return [F(x) for x in c]
def add(*vs): return [sum(t) for t in zip(*vs)]
def sc(k, v): return [k * t for t in v]
def conj(x): return [x[0]] + [-t for t in x[1:]]
def norm(x): return sum(t * t for t in x)
def sig(v): return [v[0], v[2], v[1]] + v[3:]           # sigma = (1 2)
def dot(x, y): return sig(omul(sig(x), sig(y)))         # x . y = x ._sigma y
def star(u, v):                                         # Definition 3.3
    x, y, z = u; xp, yp, zp = v
    return (add(dot(x, xp), dot(z, yp), dot(y, zp)),
            add(dot(y, yp), dot(x, zp), dot(z, xp)),
            add(dot(z, zp), dot(y, xp), dot(x, yp)))
def shape(k, c): return tuple(sc(F(t), c) for t in k)   # (k1, k2, k3) c
def blocksum(u): return add(*u)
def half(u): return tuple(sc(h, b) for b in u)
def dbl(u): return tuple(sc(2, b) for b in u)
def cyc(u): return (u[2], u[0], u[1])

checks = []
def check(name, ok):
    checks.append(bool(ok))
    print(("  ok    " if ok else "  FAIL  ") + name)

# ---------------------------------------------------------------- integer spans
def egcd(a, b):
    if b == 0: return (abs(a), (1 if a >= 0 else -1), 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)
def span(vectors, dim):
    """Echelon basis {pivot column: row} of the Z-span of integer vectors."""
    piv = {}
    for v in vectors:
        v = list(v)
        for c in range(dim):
            if v[c] == 0: continue
            if c not in piv:
                piv[c] = v; break
            r = piv[c]; a, b = r[c], v[c]
            g, x, y = egcd(a, b)
            piv[c] = [x * ri + y * vi for ri, vi in zip(r, v)]
            v = [(b // g) * ri - (a // g) * vi for ri, vi in zip(r, v)]
    return piv
def member(piv, v):
    v = list(v)
    for c in range(len(v)):
        if v[c] == 0: continue
        if c not in piv or v[c] % piv[c][c]: return False
        q = v[c] // piv[c][c]
        v = [vi - q * ri for vi, ri in zip(v, piv[c])]
    return True
def index(piv, dim):                     # index in Z^dim (None if not full rank)
    if len(piv) < dim: return None
    out = 1
    for c in piv: out *= abs(piv[c][c])
    return out
def same(p1, p2, dim):
    return (all(member(p1, r) for r in p2.values())
            and all(member(p2, r) for r in p1.values()))

# ---------------------------------------------------------------- coordinates
a0 = vec(0, h, h, 0, 0, 0, 0, 0)                        # (e1 + e2)/2
MOb = [omul(a0, b) for b in Lb]                          # Z-basis of MO = a0 L
MOinv = Matrix([[MOb[k][r] for k in range(8)] for r in range(8)]).inv()
def ocoords(v):                                          # integer coords in MO
    c = list(MOinv * Matrix(v))
    assert all(x.q == 1 for x in c), "not in MO"
    return [int(x) for x in c]
def lcoords(v):                                          # integer coords in L
    c = Lcoords(v)
    assert all(x.q == 1 for x in c), "not in L"
    return [int(x) for x in c]
def inMO(v): return all(x.q == 1 for x in (MOinv * Matrix(v)))

p = vec(0, 0, 1, 0, -1, 1, -1, 0)
e0 = vec(1, 0, 0, 0, 0, 0, 0, 0)
e1, e2 = vec(0, 1, 0, 0, 0, 0, 0, 0), vec(0, 0, 1, 0, 0, 0, 0, 0)
from verify_consecutive_twists_exact import s as s_
v = vec(h, 0, h, h, 0, 0, 0, -h)                         # (e0 + e2 + e3 - e7)/2

# ---------------------------------------------------------------- (1)
print("(1) the products (yp)z span MO")
check("e0 and p lie in MO; p has norm 4 and p/2 is not in MO",
      inMO(e0) and inMO(p) and norm(p) == 4 and not inMO(sc(h, p)))
check("v = -(1/2)(e1+e2)s = (1/2)(e0+e2+e3-e7), with s in L, lies in MO",
      v == sc(-h, omul(add(e1, e2), s_)) and inL(s_) and inMO(v))
check("pv + vp = 2 v_0 p - 2<p, v> e0 = p - e0",
      add(omul(p, v), omul(v, p))
      == add(sc(2 * v[0], p), sc(-2 * sum(x * y for x, y in zip(p, v)), e0))
      == add(p, sc(-1, e0)))
check("y = (yp)e0 - (yp)v - ((yv)p)e0 on the basis of MO",
      all(y == add(omul(omul(y, p), e0), sc(-1, omul(omul(y, p), v)),
                   sc(-1, omul(omul(omul(y, v), p), e0))) for y in MOb))
check("y = e0(py) - v(py) - e0(p(vy)) on the basis of MO",
      all(y == add(omul(e0, omul(p, y)), sc(-1, omul(v, omul(p, y))),
                   sc(-1, omul(e0, omul(p, omul(v, y))))) for y in MOb))
check("MO is closed under conjugation (on its basis)", all(inMO(conj(b)) for b in MOb))
P1 = [omul(omul(bi, p), bj) for bi in MOb for bj in MOb]
P2 = [omul(bi, omul(p, bj)) for bi in MOb for bj in MOb]
check("the 64 products (b_i p) b_j span MO (index 1)",
      index(span([ocoords(v) for v in P1], 8), 8) == 1)
check("the 64 products b_i (p b_j) span MO (index 1)",
      index(span([ocoords(v) for v in P2], 8), 8) == 1)
check("conj((b_i p) b_j) = -conj(b_j) (p conj(b_i)) on basis pairs",
      all(conj(omul(omul(bi, p), bj)) == sc(-1, omul(conj(bj), omul(p, conj(bi))))
          for bi in MOb for bj in MOb))
check("the norm is integral on MO (basis norms and 2<b_i, b_j> are integers)",
      all(norm(bi).denominator == 1 for bi in MOb)
      and all((2 * sum(x * y for x, y in zip(bi, bj))).denominator == 1
              for bi in MOb for bj in MOb))
check("norm((b_i p) b_j) = 4 norm(b_i) norm(b_j)",
      all(norm(omul(omul(bi, p), bj)) == 4 * norm(bi) * norm(bj)
          for bi in MOb for bj in MOb))
i64 = index(span([ocoords(omul(b, p)) for b in MOb] + [ocoords(omul(p, b)) for b in MOb], 8), 8)
check("for contrast: MO p + p MO has index 64 in MO (got %s)" % i64, i64 == 64)

# ---------------------------------------------------------------- (2)
print("(2) the sublattice facts")
LsP = span([lcoords(w) for w in Lsb], 8)
LsbarP = span([lcoords(a) for a in Lsbarb], 8)
check("(4.1): sigma(Lsbar) = 2 MO",
      same(span([ocoords(sig(a)) for a in Lsbarb], 8), span([ocoords(sc(2, b)) for b in MOb], 8), 8))
check("(4.1): sigma(Ls) = p MO",
      same(span([ocoords(sig(w)) for w in Lsb], 8), span([ocoords(omul(p, b)) for b in MOb], 8), 8))
check("Lsbar . Lsbar lies in 2 Lsbar",
      all(inLsbar(sc(h, dot(a, b))) for a in Lsbarb for b in Lsbarb))
check("Lsbar . Ls and Ls . Lsbar lie in Lsbar",
      all(inLsbar(dot(a, w)) and inLsbar(dot(w, a)) for a in Lsbarb for w in Lsb))
check("Ls . Ls lies in Ls", all(inLs(dot(w, v)) for w in Lsb for v in Lsb))
check("Ls . Ls spans Ls",
      same(span([lcoords(dot(w, v)) for w in Lsb for v in Lsb], 8), LsP, 8))
check("Lsbar . Ls spans Lsbar",
      same(span([lcoords(dot(a, w)) for a in Lsbarb for w in Lsb], 8), LsbarP, 8))

# ---------------------------------------------------------------- (3)
print("(3) alpha, beta, gamma: a Z-basis of Lambda")
Z = [F(0)] * 8
def alpha(a): return (a, Z, sc(-1, a))
def beta(a): return (a, sc(-1, a), Z)
def gamma(w): return (sc(-1, w), w, w)
gens = ([('alpha', a, alpha(a)) for a in Lsbarb] + [('beta', a, beta(a)) for a in Lsbarb]
        + [('gamma', w, gamma(w)) for w in Lsb])
check("the 24 vectors lie in Lambda", all(inLambda(g) for _, _, g in gens))
dG = abs(Matrix([l3(g) for _, _, g in gens]).det())
check("their determinant in L^3 is 2^12 = [L^3 : Lambda] (got %s)" % dG,
      dG == 2 ** 12 == detB)

# ---------------------------------------------------------------- (4), (5)
print("(4) the nine kinds of products, (5) where they lie")
SHAPES = {('alpha', 'alpha'): (1, -2, 1), ('beta', 'beta'): (1, 1, -2),
          ('alpha', 'beta'): (2, -1, -1), ('beta', 'alpha'): (2, -1, -1),
          ('alpha', 'gamma'): (-2, 2, 0), ('gamma', 'alpha'): (-2, 2, 0),
          ('beta', 'gamma'): (-2, 0, 2), ('gamma', 'beta'): (-2, 0, 2),
          ('gamma', 'gamma'): (3, -1, -1)}
check("(c, -c, 0) and (c, 0, -c) lie in Lambda for c in the basis of Lsbar",
      all(inLambda((c, sc(-1, c), Z)) and inLambda((c, Z, sc(-1, c))) for c in Lsbarb))
prods, ok_shape, ok_where = [], True, True
for k1, f1, g1 in gens:
    for k2, f2, g2 in gens:
        u = star(g1, g2)
        prods.append(u)
        c = dot(f1, f2)
        ok_shape &= (u == shape(SHAPES[(k1, k2)], c))
        if (k1, k2) != ('gamma', 'gamma'):
            ok_where &= inLambda(half(u))
        else:
            ok_where &= (tuple(add(x, c) for x in u) == (sc(4, c), Z, Z)
                         and inLambda((sc(2, c), Z, Z)) and inLambda((c, c, c)))
check("all 576 products have the stated shapes", ok_shape)
check("the first eight kinds lie in 2 Lambda; gamma*gamma = (4c,0,0) - (c,c,c) "
      "with (2c,0,0), (c,c,c) in Lambda", ok_where)

# ---------------------------------------------------------------- (6)
print("(6) S = Delta + 2 Lambda")
S = span([l3(u) for u in prods], 24)
Delta = [(w, w, w) for w in Lsb]
twoLam = [dbl(g) for _, _, g in gens]
X = span([l3(u) for u in Delta + twoLam], 24)
check("Delta lies in S", all(member(S, l3(u)) for u in Delta))
check("2 Lambda lies in S", all(member(S, l3(u)) for u in twoLam))
check("S lies in Delta + 2 Lambda", all(member(X, r) for r in S.values()))
iS, iX = index(S, 24), index(X, 24)
check("[Lambda : S] = 2^16 (got %s)" % (iS // detB), iS == iX and iS == detB * 2 ** 16)
check("[Delta + 2 Lambda : 2 Lambda] = 2^8, so Delta cap 2 Lambda = 2 Delta",
      (2 ** 24 * detB) // iX == 2 ** 8)

# ---------------------------------------------------------------- (7)
print("(7) the cyclic shift of the blocks")
check("it maps the basis of Lambda into Lambda", all(inLambda(cyc(g)) for _, _, g in gens))
check("it commutes with star on basis pairs",
      all(cyc(star(g1, g2)) == star(cyc(g1), cyc(g2)) for _, _, g1 in gens for _, _, g2 in gens))
shared = twisted((0, 2, 1, 3, 4, 5, 6, 7))
check("star agrees with the shared twisted product for (1 2) on basis pairs",
      all(shared(g1, g2) == star(g1, g2) for _, _, g1 in gens for _, _, g2 in gens))

print()
print("ALL PASS" if all(checks) else "FAIL -- INVESTIGATE")
sys.exit(0 if all(checks) else 1)
