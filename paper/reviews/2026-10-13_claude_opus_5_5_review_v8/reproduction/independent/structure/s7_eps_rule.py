"""
s7_eps_rule.py -- the rule that decides the span index, in both frames.

Conventions: leech_core.py.  eps := class of 2e0 in A = L/2L.  For a twist g,
1_g := g^{-1}(e0) is the identity of the twisted algebra (R^8, ._g), and
eps in V_g = g(Lsbar)/2L  <=>  2*1_g = g^{-1}(2e0) in Lsbar  (2L in Lsbar),
eps in W_g = g(Ls)/2L     <=>  2*1_g in Ls.

RULE (all g in O(L) = W(E8) whose twisted product closes on Lambda):
     [Lambda : S_g] = 2^(2k) if 2*1_g in Lsbar,  2^k if 2*1_g in Ls,  1 otherwise,
     k = 5 if det g = +1,  k = 8 if det g = -1.
Status: index = h(W_g) h(V_g)^2 is proved (s2_criterion.py, 2-adic splitting);
the values of h on the 142 subalgebra Lagrangians are an exhaustive finite
computation (s3); the inequality h(U) >= 2 when eps in U is proved by the
lemma below and machine-checked here on all 30 Lagrangians containing eps.

LEMMA.  Let U be a Lagrangian of (A, q) with eps in U, U~ its preimage in L.
 (i)  U is orthogonal to eps, and b(eps, u) = 2 u_0 mod 2, so every element
      of U~ has integer coordinates: U~ lies in D8 (Gravesian integers).
 (ii) For x, y in Z^8: Re(xy) = x0 y0 - sum_{i>0} x_i y_i = <x,y> mod 2, and
      <x,y> = b(x,y) = 0 mod 2 on U~ (U totally singular).  So every element
      of Z-span(U~ . U~) has EVEN real part.
 (iii) 2s = (-1,1,...,1) lies in 2L, inside U~, and has real part -1.
 Hence Z-span(U~ U~) != U~, i.e. h(U) >= 2.  In particular every twist that
 fixes e0 (all permutation and signed-permutation twists) has
 2*1_g = 2e0 = s.sbar in Lsbar, so its index is at least h(V_g)^2 >= 4
 (in fact 2^10 or 2^16): index 1 is impossible without moving the identity.
Output: s7_eps_rule.out
"""
import json
from fractions import Fraction as Fr
import leech_core as C

lines = []


def log(*a):
    s_ = " ".join(str(t) for t in a)
    print(s_)
    lines.append(s_)


D3 = json.load(open('s3_lagrangians.json'))
LAG = [frozenset(U) for U in D3['lagrangians']]
P = D3['props']
EPS = C.to_bits(C.to_L_int([2, 0, 0, 0, 0, 0, 0, 0]))
twos = C.to_L_int([2 * t for t in C.S_VEC])     # 2s in L-coords

with_eps = [i for i, U in enumerate(LAG) if EPS in U]
log("Lagrangians containing eps:", len(with_eps), "| by family:",
    {f: sum(1 for i in with_eps if P[i]['family'] == f) for f in ('same', 'other')},
    "| subalgebras among them:", sum(P[i]['subalgebra'] for i in with_eps))
ok_i = ok_ii = ok_iii = True
for i in with_eps:
    Ub = C.preimage_basis(LAG[i])
    E = [C.from_L(c) for c in Ub]                          # e-coordinates of a Z-basis of U~
    ok_i &= all(all(Fr(t).denominator == 1 for t in e) for e in E)
    ok_ii &= all(C.omul(a, b)[0] % 2 == 0 for a in E for b in E)
    ok_iii &= C.to_bits(twos) in LAG[i] and C.from_L(twos)[0] == -1
log("Lemma (i)  U~ has integer coordinates, all 30:", ok_i)
log("Lemma (ii) Re(xy) even on all 64 basis products, all 30:", ok_ii)
log("Lemma (iii) 2s in U~ with real part -1, all 30:", ok_iii)
sub = [i for i, U in enumerate(LAG) if P[i]['subalgebra']]
from collections import Counter
log("h(U) over the 142 subalgebra Lagrangians by (family, eps in U):",
    dict(Counter((P[i]['family'], EPS in LAG[i], P[i]['h']) for i in sub)))


# twisted-frame form of the rule on the C3 witness ---------------------------
def in_lat_by_inverse(t, right):
    """t in L.right  <=>  t . right^{-1} in L (inverse property);  right^{-1} = conj(right)/N(right)."""
    inv = [x / C.norm(right) for x in C.conj(right)]
    return C.in_L_e(C.omul(t, inv))


r = [Fr(1, 2)] * 7 + [Fr(-1, 2)]
G = C.reflection_matrix_e(r)
two1 = C.matvec(G, [2, 0, 0, 0, 0, 0, 0, 0])               # g^{-1}(2e0) = s_r(2e0)
log("\nC3 witness r = (e0+...+e6-e7)/2: 2*1_g = s_r(2e0) =", [str(t) for t in two1],
    "| in L:", C.in_L_e(two1), "| in Lsbar:", in_lat_by_inverse(two1, C.SBAR_VEC),
    "| in Ls:", in_lat_by_inverse(two1, C.S_VEC), "| det g =", int(C.mat_det(G)))
log("   => rule predicts index 1; direct computation (s4, s5, s6) gives 1.")
sig = C.perm_matrix_e({1: 2, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7})
log("sigma: 2*1_g = 2e0 = s.sbar:", C.omul(C.S_VEC, C.SBAR_VEC) == [2, 0, 0, 0, 0, 0, 0, 0],
    "| in Lsbar:", in_lat_by_inverse([2, 0, 0, 0, 0, 0, 0, 0], C.SBAR_VEC),
    "| det g = -1 => rule predicts 2^16")

# all 120 reflections in the twisted frame ----------------------------------
roots = C.roots_of_L()
reps = [rr for rr in roots if next(t for t in rr if t != 0) > 0]
res = {}
for line in open('s4_reflections.tsv'):
    if line.startswith('root'):
        continue
    f = line.rstrip('\n').split('\t')
    res[f[0]] = (f[2] == 'True', f[3])
tab = Counter()
ok = 0
for rr in reps:
    G = C.reflection_matrix_e(rr)
    t = C.matvec(G, [2, 0, 0, 0, 0, 0, 0, 0])
    loc = 'Lsbar' if in_lat_by_inverse(t, C.SBAR_VEC) else ('Ls' if in_lat_by_inverse(t, C.S_VEC) else 'neither')
    closed, lg = res[",".join(str(x) for x in rr)]
    tab[(loc, 'closes' if closed else 'no close', lg)] += 1
    if closed:
        pred = {'Lsbar': '16', 'Ls': '8', 'neither': '0'}[loc]
        ok += (pred == lg)
log("\n120 reflections: (location of 2*1_g = s_r(2e0), closure, log2 index) : count")
for k_ in sorted(tab, key=str):
    log("   ", k_, ":", tab[k_])
log("rule (det = -1, k = 8) matches the directly computed index on", ok, "of", sum(v for k_, v in tab.items() if k_[1] == 'closes'),
    "closing reflections")
with open('s7_eps_rule.out', 'w') as f:
    f.write("\n".join(lines) + "\n")
