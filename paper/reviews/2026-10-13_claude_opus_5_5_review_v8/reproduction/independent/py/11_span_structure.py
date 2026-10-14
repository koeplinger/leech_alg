r"""
11_span_structure.py -- an explicit description of the span S for the
index-2^16 class (sigma, every closing odd permutation, the 28 integral root
reflections), addressing the paper's remark that identifying S/2Lambda
structurally is open (Section 5.5).

Candidate (derived in the report from the mod-2 analysis):
    S_16 := { (x, y, z) : x, y, z in L.s ;  x - y, x - z in 2 (L.sbar) }
          = Z-span of (w, w, w) [w in L.s], (0, 2v, 0), (0, 0, 2v) [v in L.sbar].
Checks: S_16 <= Lambda, [Lambda : S_16] = 2^16, and S_g == S_16 (equal HNFs in
Lambda-coordinates) for every closing twist g of index 2^16 in the census
(permutations and reflections).  Also reports, for comparison, whether the
index-2^10 twists all share one span lattice.
Run:  python3 11_span_structure.py > 11_span_structure.out.txt
"""
import json
from fractions import Fraction as F
from importlib import import_module

import numpy as np
import leechlib as LL

census_mod = import_module('03_twist_census')     # reuses products_and_coords()
D = census_mod.D
BLam = D['BLam']
N2 = D['N2']


def to_lambda_coords(rows_L3):
    P = np.array(rows_L3, dtype=np.int64)
    ok, c2 = LL.lambda_coordinates(P, N2)
    assert ok.all()
    return [list(map(int, r)) for r in (c2 // 2)]


Z8 = [0] * 8
gens = [w + w + w for w in D['Wrows']]
gens += [Z8 + [2 * a for a in v] + Z8 for v in D['Vrows']]
gens += [Z8 + Z8 + [2 * a for a in v] for v in D['Vrows']]
S16_coords = to_lambda_coords(gens)            # asserts S_16 <= Lambda
idx16, H16 = LL.hnf_index(S16_coords, 24)
print(f"S_16 = span{{(w,w,w), (0,2v,0), (0,0,2v) : w in L.s, v in L.sbar}}: contained in Lambda; "
      f"[Lambda : S_16] = {idx16} = 2^{idx16.bit_length() - 1}")
# the same lattice described as {(x,y,z) in (L.s)^3 : x-y, x-z in 2 L.sbar}: check both inclusions
# on generators: (x, x+2b, x+2c) with x in L.s, b, c in L.sbar is exactly what the generators span.

census = json.load(open('twist_results.json'))
classes = {}
for kind, items in (('perm', census['permutations']), ('refl', census['reflections'])):
    for it in items:
        if not it['closes']:
            continue
        spec = ('perm', tuple(it['pi'])) if kind == 'perm' else ('refl', tuple(it['root_doubled_e_coords']))
        G, Cg, P, ok, c2 = census_mod.products_and_coords(spec)
        rows = [list(map(int, r)) for r in np.unique(c2 // 2, axis=0) if np.any(r)]
        idx, H = LL.hnf_index(rows, 24)
        key = tuple(tuple(h) for h in H)
        classes.setdefault(idx, {}).setdefault(key, []).append(spec)

for idx in sorted(classes):
    lat = classes[idx]
    n = sum(len(v) for v in lat.values())
    print(f"index 2^{idx.bit_length() - 1}: {n} closing twists, {len(lat)} distinct span lattice(s) S_g")
    if idx == idx16:
        same = list(lat.keys()) == [tuple(tuple(h) for h in H16)]
        print(f"   all of them equal S_16: {same}")

# Independent (HNF-free) check for sigma with Fractions from the definitions:
# every product b_i *_sigma b_j satisfies the defining conditions of S_16
#   x, y, z in L.s   (a in L.s <=> (a.sbar)/2 in L)
#   (x - y)/2, (x - z)/2 in L.sbar   (a in L.sbar <=> (a.s)/2 in L)
# so S_sigma <= S_16; both have index 2^16 in Lambda (F2-rank / 06), hence equal.
B = []
for line in open('lambda_basis.txt'):
    if not line.startswith('#'):
        B.append([F(int(a), 2) for a in line.split('|L3|')[0].split()])
pi = (2, 1, 3, 4, 5, 6, 7)
mul = lambda x, y: LL.perm_apply(pi, LL.omul(LL.perm_apply(pi, x), LL.perm_apply(pi, y)))   # sigma = sigma^{-1}
n_ok = 0
for i in range(24):
    for j in range(24):
        u, v = B[i], B[j]
        x, y, z = u[0:8], u[8:16], u[16:24]
        xp, yp, zp = v[0:8], v[8:16], v[16:24]
        P1 = LL.vadd(mul(x, xp), mul(z, yp), mul(y, zp))
        P2 = LL.vadd(mul(y, yp), mul(x, zp), mul(z, xp))
        P3 = LL.vadd(mul(z, zp), mul(y, xp), mul(x, yp))
        ok = (LL.in_L(P1) and LL.in_L(P2) and LL.in_L(P3)
              and LL.in_L_s(P1) and LL.in_L_s(P2) and LL.in_L_s(P3)
              and LL.in_L_sbar(LL.vscale(F(1, 2), LL.vsub(P1, P2)))
              and LL.in_L_sbar(LL.vscale(F(1, 2), LL.vsub(P1, P3))))
        n_ok += ok
print(f"sigma, Fractions from the definitions: {n_ok}/576 products b_i *_sigma b_j satisfy "
      f"x,y,z in L.s and (x-y)/2, (x-z)/2 in L.sbar  ->  S_sigma <= S_16, and equal indices give S_sigma = S_16")
