"""
s2_criterion.py -- machine check of every step in the proofs of
  (C5)   *_g closes on Lambda  <=>  V_g = g(Lsbar)/2L and W_g = g(Ls)/2L are
         subalgebras of A = L/2L,
  (IDX)  [Lambda : S_g] = h(W_g) * h(V_g)^2,  h(U) := [U~ : Z-span(U~ . U~)],
for g in O(L).  Conventions: leech_core.py.

Reduction used throughout.  For g in O(L) put G = g x g x g.  Then
G(u *_g v) = G(u) * G(v) with * the UNTWISTED triple product, and
G(Lambda) = Lambda_{V,W} := {(x,y,z) in L^3 : x+y, x+z, y+z in V~; x+y+z in W~}
with V~ = g(Lsbar), W~ = g(Ls).  So it suffices to treat the untwisted product
on Lambda_{V,W} for an arbitrary complementary pair of Lagrangians (V, W).

PART 1 (facts the proof needs) -- checked in s1_basics.py:
   L.L in L;  2L in Lsbar and in Ls;  Lsbar + Ls = L;  Lsbar cap Ls = 2L.
   Descent: (x + 2c) y = xy + 2(cy) with cy in L, so the product is defined on A.
PART 2 (sufficiency): formal identities, valid for ANY bilinear product, checked
   here symbolically in the free bilinear algebra:
     with x = a+w, y = b+w, z = c+w, x' = a'+w', y' = b'+w', z' = c'+w',
     a+b+c = 0 = a'+b'+c' (mod 2):
       P1+P2 = a a' + b b',  P2+P3 = b b' + c c',  P1+P3 = a a' + c c'  (mod 2)
       P1+P2+P3 = (x+y+z)(x'+y'+z')   (exactly, over Z)  = w w'  (mod 2).
PART 3 (necessity): explicit elements u, v of Lambda_{V,W} whose product
   detects V.V and W.W; demonstrated on Wilson's own pair (W0 not closed).
PART 4 (2-adic splitting): Lambda' contains Delta(W~) + (V~ (x) T) with index
   3^8, Delta(a)*Delta(b) = Delta(3ab), Delta * T = T * Delta = 0,
   (x(x)alpha)*(y(x)beta) = xy (x) (alpha o beta), [T_Z : T_Z o T_Z] = 3.
   Consequences (proved in the review text): closure <=> V~.V~ in V~ and
   W~.W~ in W~, and [Lambda':S'] = [W~ : W~W~] [V~ : V~V~]^2 (checked on all
   4,704 closed pairs in s3_lagrangians.py).
PART 5: for the paper's sigma, S = 2 Lambda + Delta(Ls) exactly.
Output: s2_criterion.out
"""
from collections import Counter
from fractions import Fraction as Fr
import itertools
import numpy as np
import leech_core as C

lines = []


def log(*a):
    s_ = " ".join(str(t) for t in a)
    print(s_)
    lines.append(s_)


# ---------------------------------------------------------------------------
# PART 2: formal identities in the free bilinear algebra
# An element of the "first factor" is a dict {generator: coeff}; a product of a
# first-factor element with a second-factor element is a dict {(g1, g2): coeff}.
# ---------------------------------------------------------------------------
def lin(**kw):
    return {k: v for k, v in kw.items() if v}


def add(*xs):
    out = Counter()
    for x in xs:
        for k, v in x.items():
            out[k] += v
    return {k: v for k, v in out.items() if v}


def bprod(x, y):
    out = Counter()
    for g1, c1 in x.items():
        for g2, c2 in y.items():
            out[(g1, g2)] += c1 * c2
    return {k: v for k, v in out.items() if v}


def mod2(p):
    return {k: v % 2 for k, v in p.items() if v % 2}


def triple(u, v):
    (x, y, z), (xp, yp, zp) = u, v
    P1 = add(bprod(x, xp), bprod(z, yp), bprod(y, zp))
    P2 = add(bprod(y, yp), bprod(x, zp), bprod(z, xp))
    P3 = add(bprod(z, zp), bprod(y, xp), bprod(x, yp))
    return P1, P2, P3


# generic u = (x,y,z), v = (x',y',z') over Z: exact identity P1+P2+P3 = (x+y+z)(x'+y'+z')
X, Y, Z = lin(x=1), lin(y=1), lin(z=1)
Xp, Yp, Zp = lin(xp=1), lin(yp=1), lin(zp=1)
P1, P2, P3 = triple((X, Y, Z), (Xp, Yp, Zp))
log("PART 2  exact identity P1+P2+P3 == (x+y+z)(x'+y'+z') in the free bilinear algebra:",
    add(P1, P2, P3) == bprod(add(X, Y, Z), add(Xp, Yp, Zp)))
log("        Z3-equivariance (x,y,z)->(y,z,x) maps (P1,P2,P3) to (P2,P3,P1):",
    [mod2(p) for p in triple((Y, Z, X), (Yp, Zp, Xp))] ==
    [mod2(p) for p in (lambda q: (q[1], q[2], q[0]))(triple((X, Y, Z), (Xp, Yp, Zp)))])
# mod-2 normal form: x = a+w, y = b+w, z = (a+b)+w   (c = a+b mod 2)
x = lin(a=1, w=1)
y = lin(b=1, w=1)
z = lin(a=1, b=1, w=1)
xp = lin(ap=1, wp=1)
yp = lin(bp=1, wp=1)
zp = lin(ap=1, bp=1, wp=1)
P1, P2, P3 = triple((x, y, z), (xp, yp, zp))
c, cp = lin(a=1, b=1), lin(ap=1, bp=1)
aa, bb, cc = bprod(lin(a=1), lin(ap=1)), bprod(lin(b=1), lin(bp=1)), bprod(c, cp)
ww = bprod(lin(w=1), lin(wp=1))
log("        P1+P2 == a a' + b b'  (mod 2):", mod2(add(P1, P2)) == mod2(add(aa, bb)))
log("        P2+P3 == b b' + c c'  (mod 2):", mod2(add(P2, P3)) == mod2(add(bb, cc)))
log("        P1+P3 == a a' + c c'  (mod 2):", mod2(add(P1, P3)) == mod2(add(aa, cc)))
log("        P1+P2+P3 == w w'      (mod 2):", mod2(add(P1, P2, P3)) == mod2(ww))
ab_ = bprod(lin(a=1), lin(bp=1))
ba_ = bprod(lin(b=1), lin(ap=1))
log("        P1 == a a' + a b' + b a' + w w'  (mod 2):", mod2(P1) == mod2(add(aa, ab_, ba_, ww)))
log("   => if V.V in V and W.W in W then P1+P2, P2+P3, P1+P3 lie in V and P1+P2+P3 in W,")
log("      i.e. u*v satisfies Wilson's conditions for (V, W); with L.L in L this is closure.")

# ---------------------------------------------------------------------------
# PART 3: necessity, made concrete on Wilson's own pair
# ---------------------------------------------------------------------------
s, sb = C.S_VEC, C.SBAR_VEC
H_sb = C.hnf([C.to_L_int(C.omul(C.L_BASIS[k], sb)) for k in range(8)])
H_s = C.hnf([C.to_L_int(C.omul(C.L_BASIS[k], s)) for k in range(8)])
V0 = C.span_bits([C.to_bits(v) for v in H_sb])
W0 = C.span_bits([C.to_bits(v) for v in H_s])
TAB = C.mul_table_mod2(C.T_STD)
log("\nPART 3  necessity.  For a~, a~' in V~:  u = (a~,0,a~), v = (a~',0,a~') lie in Lambda_{V,W}")
log("        and u*v = (a~a~', 2a~a~', a~a~'), so P1+P2 = 3 a~a~' : closure forces V.V in V.")
log("        For w~, w~' in W~:  u = (w~,w~,w~), v = (w~',w~',w~') lie in Lambda_{V,W}")
log("        and u*v = 3(w~w~', w~w~', w~w~'), so P1+P2+P3 = 9 w~w~' : closure forces W.W in W.")
# concrete: paper's Remark 4.5 elements u1, u2 in Ls with u1.u2 not in Ls
u1 = [-1, 0, 0, 0, 1, 0, 1, 1]
u2 = [-1, 1, 0, 0, 0, 1, 0, 1]
c1, c2 = C.to_L_int(u1), C.to_L_int(u2)
p = C.omul(u1, u2)
log("        Wilson's pair: u1 =", u1, "u2 =", u2, "in Ls:", C.to_bits(c1) in W0, C.to_bits(c2) in W0)
log("        u1.u2 =", p, "in Ls:", C.to_bits(C.to_L_int(p)) in W0, "| in Lsbar:", C.to_bits(C.to_L_int(p)) in V0)
U0 = np.array(C.leech_basis_from_pair(V0, W0), dtype=np.int64)
u = np.array(c1 + c1 + c1, dtype=np.int64)[None, :]
v = np.array(c2 + c2 + c2, dtype=np.int64)[None, :]
in_lam = C.in_pair_lattice_mod2(np.vstack([u, v]), V0, W0)
prod = np.concatenate([sum(np.einsum('i,j,ijk->k', u[0, 8 * a:8 * a + 8], v[0, 8 * b:8 * b + 8], C.T_STD)
                           for (a, b) in pairs) for pairs in
                       ([(0, 0), (2, 1), (1, 2)], [(1, 1), (0, 2), (2, 0)], [(2, 2), (1, 0), (0, 1)])])
log("        u = (u1,u1,u1), v = (u2,u2,u2) in Lambda:", bool(in_lam.all()),
    "| u*v in Lambda:", bool(C.in_pair_lattice_mod2(prod[None, :], V0, W0)[0]))
log("        V0 subalgebra:", all(TAB[a, b] in V0 for a in V0 for b in V0),
    "| W0 subalgebra:", all(TAB[a, b] in W0 for a in W0 for b in W0))

# ---------------------------------------------------------------------------
# PART 4: the 2-adic splitting Lambda'_2 = Delta(W~)_2 (+) (V~ (x) T)_2
# ---------------------------------------------------------------------------
log("\nPART 4  2-adic splitting (shown for Wilson's Lambda; the argument is the same for every pair)")
Vb, Wb = C.preimage_basis(V0), C.preimage_basis(W0)
z8 = [0] * 8
Delta = [w + w + w for w in Wb]                                   # Delta(W~)
Tpart = [a + z8 + [-t for t in a] for a in Vb] + [z8 + a + [-t for t in a] for a in Vb]   # V~ (x) T
sub = C.hnf(Delta + Tpart)
idx_sub = C.hnf_det(sub, 24) // C.hnf_det(U0.tolist(), 24)
log("        [Lambda : Delta(Ls) + (Lsbar (x) T)] =", idx_sub, "= 3^%d" % round(np.log(idx_sub) / np.log(3)))
inside = C.in_pair_lattice_mod2(np.array(Delta + Tpart), V0, W0).all()
_, ok3 = C.coords_in_basis(3 * U0, np.array(sub, dtype=np.int64))
log("        Delta(Ls) and Lsbar(x)T lie in Lambda:", bool(inside), "| 3*Lambda inside their sum:", bool(ok3.all()))
# products: Delta(a)*Delta(b) = Delta(3ab), Delta * T = 0 = T * Delta
def tp(u, v, T=C.T_STD):
    u = np.asarray(u, dtype=np.int64)
    v = np.asarray(v, dtype=np.int64)
    blk = lambda w, i: w[8 * i:8 * i + 8]
    bp = lambda i, j: np.einsum('i,j,ijk->k', blk(u, i), blk(v, j), T)
    return np.concatenate([bp(0, 0) + bp(2, 1) + bp(1, 2), bp(1, 1) + bp(0, 2) + bp(2, 0), bp(2, 2) + bp(1, 0) + bp(0, 1)])
okDD = all(np.array_equal(tp(d1, d2), np.concatenate([3 * C.lmul(d1[:8], d2[:8])] * 3)) for d1 in Delta for d2 in Delta)
okDT = all(not tp(d, t).any() and not tp(t, d).any() for d in Delta for t in Tpart)
log("        Delta(a)*Delta(b) == Delta(3ab) on all basis pairs:", okDD, "| Delta*T == T*Delta == 0:", okDT)
# the 3-dim routing algebra C = R^3: e_i e_i = e_i, e_i e_j = e_k ({i,j,k} = {1,2,3})
def route(al, be):
    p, q, r = al
    pp, qq, rr = be
    return (p * pp + r * qq + q * rr, q * qq + p * rr + r * pp, r * rr + q * pp + p * qq)
TZ = [(1, -1, 0), (0, 1, -1)]
prods = [route(a, b) for a in TZ for b in TZ]
idxN = C.hnf_det([[t[0], t[1]] for t in prods], 2) // C.hnf_det([[t[0], t[1]] for t in TZ], 2)
log("        routing algebra on T_Z = {p+q+r = 0}: products", prods, "-> [T_Z : span(T_Z o T_Z)] =", idxN)
okVT = True
for a in Vb:
    for b in Vb:
        for al in TZ:
            for be in TZ:
                u_ = [x_ * al[0] for x_ in a] + [x_ * al[1] for x_ in a] + [x_ * al[2] for x_ in a]
                v_ = [x_ * be[0] for x_ in b] + [x_ * be[1] for x_ in b] + [x_ * be[2] for x_ in b]
                ab = C.lmul(a, b)
                g = route(al, be)
                okVT &= np.array_equal(tp(u_, v_), np.concatenate([ab * g[0], ab * g[1], ab * g[2]]))
log("        (x (x) alpha)*(y (x) beta) == xy (x) (alpha o beta) on all basis elements x, y of V~ and alpha, beta of T_Z:", okVT)
log("   => over Z_2 (3 a unit): Lambda'_2 = Delta(W~)_2 (+) (V~ (x) T_Z)_2 as rings, with")
log("      S'_2 = Delta(span W~W~)_2 (+) (span(V~V~) (x) T_Z)_2, hence closure <=> V~V~ in V~ and W~W~ in W~,")
log("      and [Lambda':S'] = [W~:W~W~] [V~:V~V~]^2 (odd primes: both sides equal [L:span LL]^3_odd = 1).")

# ---------------------------------------------------------------------------
# PART 5: sigma.  S_sigma = 2 Lambda + Delta(Ls)
# ---------------------------------------------------------------------------
sig = C.perm_matrix_e({1: 2, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7})
Tsig = C.tensor_int(C.structure_tensor_L(sig))
P = C.triple_products(U0, Tsig)
Cc, ok = C.coords_in_basis(P, U0)
S_hnf = C.hnf(P.tolist())
cand = C.hnf([[2 * t for t in row] for row in U0.tolist()] + [w + w + w for w in H_s])
log("\nPART 5  sigma = (1 2): closes:", bool(ok.all()), "| [Lambda:S] =", C.span_index(Cc)[0],
    "| S == 2*Lambda + Delta(Ls) (equal HNF in L^3-coords):", S_hnf == cand)

with open('s2_criterion.out', 'w') as f:
    f.write("\n".join(lines) + "\n")
