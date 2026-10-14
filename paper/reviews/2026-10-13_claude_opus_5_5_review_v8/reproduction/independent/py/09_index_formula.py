r"""
09_index_formula.py -- an F2-level formula for the span index, tested on all
4,704 closing polarizations (V, W) of 07_polarizations.py.

Notation: Abar = L/2L with the reduced product; for a subspace X of Abar
  d(X) = 4 - dim span(X X)                  ("square defect")
  f(X) = dim span(X + Abar X + X Abar) - 4  (how far X is from a two-sided ideal)
  m(X) = 4 - f(X)                           ("ideal defect")
Proven (report, and checked in 07):
  [Lambda : S + (2L)^3] = 2^( d(W) + 2 d(V) ).
Tested here (observation, not a proof):
  [(2L)^3 : S /\ (2L)^3] = 2^( m(W) + 2 m(V) ),
i.e.  log2 [Lambda : S] = d(W) + 2 d(V) + m(W) + 2 m(V).
Heuristic reason: products with one factor in (2L)^3 give 2*(a x', a z', a y')
etc., whose images in (2L)^3 / 2Lambda = Abar/V + Abar/V + Abar/W are
(Abar V + V Abar)/V twice and (Abar W + W Abar)/W once.
Run:  python3 09_index_formula.py > 09_index_formula.out.txt
"""
import json
import numpy as np
import leechlib as LL

D = LL.build_all()
CM2 = (D['C'] % 2).astype(np.int64)


def vec(m):
    return [(m >> i) & 1 for i in range(8)]


def prod(a, b):
    return list(np.einsum('i,j,ijk->k', np.array(vec(a)), np.array(vec(b)), CM2) % 2)


pol = json.load(open('polarization_results.json'))
lag = pol['lagrangians']
inv = {}
for k, s in enumerate(lag):
    X = [vec(m) for m in s]
    XX = [prod(a, b) for a in s for b in s]
    AX = [prod(a, b) for a in range(256) for b in s] + [prod(b, a) for a in range(256) for b in s]
    inv[k] = {'d': 4 - LL.f2_rank(XX), 'm': 4 - (LL.f2_rank(X + AX) - 4)}

ok = 0
table = {}
closing = [p for p in pol['pairs'] if p['closes']]
for p in closing:
    dV, mV = inv[p['a']]['d'], inv[p['a']]['m']
    dW, mW = inv[p['b']]['d'], inv[p['b']]['m']
    pred1 = dW + 2 * dV
    pred2 = mW + 2 * mV
    e = p['index'].bit_length() - 1
    e1 = p['idx_S_2L3'].bit_length() - 1
    good = (e1 == pred1) and (e - e1 == pred2)
    ok += good
    key = (e, dV, mV, dW, mW, good)
    table[key] = table.get(key, 0) + 1
print("log2 index | d(V) m(V) | d(W) m(W) | formula d(W)+2d(V)+m(W)+2m(V) holds | number of ordered pairs")
for k, v in sorted(table.items()):
    print(f"   {k[0]:2d}      |  {k[1]}    {k[2]}   |  {k[3]}    {k[4]}   |  {k[5]}  |  {v}")
print(f"\nformula holds on {ok} of {len(closing)} closing polarizations")
print("\nper-Lagrangian invariants among the 142 subalgebras: (d, m) -> count")
h = {}
for k, s in enumerate(lag):
    if pol['lagrangian_info'][k]['sub']:
        key = (inv[k]['d'], inv[k]['m'])
        h[key] = h.get(key, 0) + 1
print("  ", dict(sorted(h.items())))
