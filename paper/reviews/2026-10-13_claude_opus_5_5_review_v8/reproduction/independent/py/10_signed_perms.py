r"""
10_signed_perms.py -- claim C4: twists by signed permutations of the imaginary
units (e0 fixed).

CONVENTION: a signed permutation g = (pi, eps) acts by g(e0) = e0,
g(e_a) = eps_a e_{pi(a)} (a = 1..7, eps_a = +-1), i.e. (g x)_{pi(a)} = eps_a x_a.
Twisted product x ._g y = g^{-1}(g(x).g(y)), triple product as in leechlib.

Steps
 1. All 2^7 * 5040 = 645,120 signed permutations; the twisted multiplication
    table e_a ._g e_b of each; number of distinct tables and the size of each
    fibre (= the signed permutations that are automorphisms of the standard
    table, if the fibres are cosets of that group).
 2. g(L) = L iff the number of minus signs is even (one sign change moves the
    half-integer vectors of odd coordinate sum to even sum).  For odd sign
    count, g = kappa o g'' with kappa = conjugation (7 sign changes) and
    g'' of even sign count, and x ._g y = y ._g'' x (the opposite product);
    for the Z3-routed triple product this gives u *_g v = v *_g'' u, so
    closure and span index of g and g'' coincide.  We do NOT rely on this:
    all distinct tables are computed directly.
 3. For every distinct table: structure constants in the L-basis computed
    with Fractions directly from g (asserted integral, i.e. L ._g L <= L),
    576 products, closure, exact HNF index, F2-rank.
Outputs: signed_perm_results.json, 10_signed_perms.out.txt.
Run:  python3 10_signed_perms.py > 10_signed_perms.out.txt
"""
import json
import time
from fractions import Fraction as F
from itertools import permutations, product
from multiprocessing import Pool

import numpy as np
import leechlib as LL

D = LL.build_all()
BL, LC = D['BL'], D['LC']
MUL = LL.MUL


def table_key(pi, eps):
    full = (0,) + tuple(pi)
    e = (1,) + tuple(eps)
    pinv = [0] * 8
    for a in range(8):
        pinv[full[a]] = a
    key = []
    for a in range(1, 8):
        for b in range(1, 8):
            s, cp = MUL[full[a]][full[b]]
            c = pinv[cp]
            key.append(s * e[a] * e[b] * e[c] * (c + 1))   # signed index, c+1 to keep the sign of e0
    return tuple(key)


def g_apply(pi, eps, x):
    y = [None] * 8
    y[0] = x[0]
    for a in range(1, 8):
        y[pi[a - 1]] = eps[a - 1] * x[a]
    return y


def g_inv_apply(pi, eps, y):
    x = [None] * 8
    x[0] = y[0]
    for a in range(1, 8):
        x[a] = eps[a - 1] * y[pi[a - 1]]
    return x


def analyse(job):
    pi, eps = job
    Cg = np.zeros((8, 8, 8), dtype=np.int64)
    for i in range(8):
        gi = g_apply(pi, eps, BL[i])
        for j in range(8):
            gj = g_apply(pi, eps, BL[j])
            Cg[i, j, :] = LC.to_L(g_inv_apply(pi, eps, LL.omul(gi, gj)))   # asserts integrality
    P = LL.all_basis_products(D['BLam_np'], Cg).reshape(576, 24)
    ok, c2 = LL.lambda_coordinates(P, D['N2'])
    out = {'pi': list(pi), 'eps': list(eps), 'n_minus': sum(1 for t in eps if t < 0),
           'closes': bool(ok.all()), 'n_fail': int((~ok).sum())}
    if out['closes']:
        rows = [list(map(int, r)) for r in np.unique(c2 // 2, axis=0) if np.any(r)]
        rows.sort(key=lambda r: (sum(abs(x) for x in r), r))
        idx, H = LL.hnf_index(rows, 24)
        out.update({'index': idx, 'f2rank': LL.f2_rank(rows)})
    return out


def main():
    t0 = time.time()
    fibres = {}
    for pi in permutations(range(1, 8)):
        for eps in product((1, -1), repeat=7):
            k = table_key(pi, eps)
            if k not in fibres:
                fibres[k] = [(pi, eps), 0, set()]
            fibres[k][1] += 1
            fibres[k][2].add(sum(1 for t in eps if t < 0) % 2)
    sizes = {}
    for v in fibres.values():
        sizes[v[1]] = sizes.get(v[1], 0) + 1
    print(f"signed permutations: {5040 * 128}; distinct twisted tables: {len(fibres)}; fibre sizes: {sizes}")
    print(f"each fibre has a single parity of the number of sign changes: {all(len(v[2]) == 1 for v in fibres.values())}")
    print(f"tables from even / odd sign count: {sum(1 for v in fibres.values() if 0 in v[2])} / "
          f"{sum(1 for v in fibres.values() if 1 in v[2])}")
    reps = [v[0] for v in fibres.values()]
    with Pool(14) as pool:
        res = pool.map(analyse, reps, chunksize=8)
    for par in (0, 1):
        sub = [r for r in res if r['n_minus'] % 2 == par]
        cl = [r for r in sub if r['closes']]
        h = {}
        for r in cl:
            h[r['index']] = h.get(r['index'], 0) + 1
        print(f"{'even' if par == 0 else 'odd '} sign count: {len(sub)} tables, {len(cl)} close; indices "
              f"{ {f'2^{k.bit_length() - 1}': v for k, v in sorted(h.items())} }")
    cl = [r for r in res if r['closes']]
    h = {}
    for r in cl:
        h[r['index']] = h.get(r['index'], 0) + 1
    print(f"all {len(res)} tables: {len(cl)} close; indices { {f'2^{k.bit_length() - 1}': v for k, v in sorted(h.items())} }")
    print(f"2 Lambda <= S for every closing table (index = 2^(24 - F2-rank)): "
          f"{all(r['index'] == 2 ** (24 - r['f2rank']) for r in cl)}")
    # pairing g <-> kappa o g (kappa = conjugation) gives identical results
    bykey = {table_key(tuple(r['pi']), tuple(r['eps'])): r for r in res}
    pair_ok = True
    for r in res:
        k2 = table_key(tuple(r['pi']), tuple(-t for t in r['eps']))
        r2 = bykey[k2]
        pair_ok &= (r2['closes'] == r['closes']) and (r2.get('index') == r.get('index'))
    print(f"g and conj o g (opposite twisted product) have identical closure and index for all tables: {pair_ok}")
    json.dump(res, open('signed_perm_results.json', 'w'), indent=0)
    print(f"total time {time.time() - t0:.1f} s")


if __name__ == '__main__':
    main()
