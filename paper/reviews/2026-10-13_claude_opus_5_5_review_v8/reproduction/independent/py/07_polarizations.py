r"""
07_polarizations.py -- the closure criterion and the span index over ALL
twists g in Aut(L) = W(E8), organized by the polarization (V_g, W_g).

Background (proved in the report; checked here exhaustively):
  * u *_g v = g^{-1}(g(u) * g(v)) blockwise, so *_g closes on Lambda iff the
    UNTWISTED product * closes on Lambda_g := g(Lambda) =
        {(x,y,z) in L^3 : x+y, x+z, y+z in g(L.sbar); x+y+z in g(L.s)},
    and S_g = g^{-1}(span of Lambda_g * Lambda_g), so the span index depends
    only on the pair of F2-subspaces V_g = g(L.sbar)/2L, W_g = g(L.s)/2L of
    the F2-algebra L/2L (product reduced from L.L <= L).
  * W(E8) maps onto O+(8,2) = O(L/2L, q), which is transitive on ordered pairs
    of complementary Lagrangians, so every such pair (V, W) is (V_g, W_g) for
    some g in Aut(L).

This script
  1. enumerates the Lagrangians (maximal totally singular subspaces for
     q(x) = N(x)/2 mod 2) of L/2L as the orbit of L.sbar/2L under the 120
     transvections x -> x + B(x,r) r (r the 120 nonsingular vectors);
  2. classifies them: subalgebra / left ideal / right ideal / two-sided ideal;
  3. enumerates ordered complementary pairs (V, W), V /\ W = 0;
  4. for EVERY ordered complementary pair builds Lambda(V,W) (Gram det 2^24 checked),
     tests closure of the untwisted product * on its 576 basis products, and
     compares with "V and W are subalgebras" (claim C5);
  5. for every closing pair computes the span index exactly (HNF), the F2-rank
     of S mod 2Lambda, [Lambda : S + (2L)^3], dim span(VV), dim span(WW);
  6. cross-checks against the census (03): for each permutation and reflection
     twist, the index looked up from (V_g, W_g) equals the index computed
     directly for the twisted product.
Outputs: polarization_results.json, 07_polarizations.out.txt (stdout).
Run:  python3 07_polarizations.py > 07_polarizations.out.txt
"""
import json
import time
from fractions import Fraction as F
from multiprocessing import Pool

import numpy as np
import leechlib as LL

D = LL.build_all()
BL, LC = D['BL'], D['LC']
C = D['C']
CM2 = (C % 2).astype(np.int64)
GL = [[int(LL.ip(a, b)) for b in BL] for a in BL]          # Gram of the L-basis


def bits(v):
    return sum((int(a) & 1) << i for i, a in enumerate(v))


def vec(m):
    return [(m >> i) & 1 for i in range(8)]


# products in L/2L for all 256 x 256 pairs (bitmask table)
PROD = [[0] * 256 for _ in range(256)]
for a in range(256):
    va = np.array(vec(a), dtype=np.int64)
    Ka = np.einsum('i,ijk->jk', va, CM2) % 2               # 8 x 8 : row j -> a.e_j
    for b in range(256):
        vb = np.array(vec(b), dtype=np.int64)
        PROD[a][b] = bits((vb @ Ka) % 2)


def q(m):
    c = vec(m)
    return (sum(c[i] * GL[i][j] * c[j] for i in range(8) for j in range(8)) // 2) % 2


def Bf(a, b):
    ca, cb = vec(a), vec(b)
    return sum(ca[i] * GL[i][j] * cb[j] for i in range(8) for j in range(8)) % 2


def span(gens):
    s = {0}
    for g in gens:
        s |= {x ^ g for x in s}
    return frozenset(s)


def main():
    t0 = time.time()
    V0 = span(bits(r) for r in D['Vrows'])
    W0 = span(bits(r) for r in D['Wrows'])
    nonsing = [m for m in range(1, 256) if q(m) == 1]
    sing = [m for m in range(1, 256) if q(m) == 0]
    print(f"L/2L: {len(sing)} nonzero singular and {len(nonsing)} nonsingular vectors (expected 135 / 120)")
    # the 120 nonsingular vectors are exactly the classes of the 240 roots
    root_classes = {bits(LC.to_L(r)) for r in LL.roots_of_L()}
    print(f"classes of the 240 roots mod 2L: {len(root_classes)}, equal to the nonsingular set: {root_classes == set(nonsing)}")
    trans = []
    for r in nonsing:
        trans.append([x ^ r if Bf(x, r) else x for x in range(256)])
    seen = {V0}
    frontier = [V0]
    while frontier:
        nxt = []
        for Ls in frontier:
            for t in trans:
                im = frozenset(t[x] for x in Ls)
                if im not in seen:
                    seen.add(im)
                    nxt.append(im)
        frontier = nxt
    lag = sorted(seen, key=lambda s: sorted(s))
    assert all(len(s) == 16 and all(q(x) == 0 for x in s) for s in lag)
    print(f"Lagrangians (orbit of L.sbar/2L under the 120 transvections): {len(lag)}  (theory: 2*3*5*9 = 270)")
    assert W0 in seen

    def is_sub(s):
        return all(PROD[a][b] in s for a in s for b in s)

    def is_left(s):
        return all(PROD[a][b] in s for a in range(256) for b in s)

    def is_right(s):
        return all(PROD[b][a] in s for a in range(256) for b in s)

    info = []
    for s in lag:
        sub, lft, rgt = is_sub(s), is_left(s), is_right(s)
        vv = LL.f2_rank([vec(PROD[a][b]) for a in s for b in s])
        info.append({'sub': sub, 'left': lft, 'right': rgt, 'dimVV': vv})
    print(f"  subalgebras: {sum(i['sub'] for i in info)}; left ideals: {sum(i['left'] for i in info)}; "
          f"right ideals: {sum(i['right'] for i in info)}; two-sided ideals: {sum(i['left'] and i['right'] for i in info)}")
    hist = {}
    for i in info:
        if i['sub']:
            k = (i['dimVV'], i['left'], i['right'])
            hist[k] = hist.get(k, 0) + 1
    print("  subalgebras by (dim span(VV), left ideal, right ideal):", dict(sorted(hist.items())))
    i0, iw = lag.index(V0), lag.index(W0)
    print(f"  L.sbar/2L: {info[i0]};  L.s/2L: {info[iw]}")

    pairs = [(a, b) for a in range(len(lag)) for b in range(len(lag)) if lag[a] & lag[b] == {0}]
    print(f"ordered complementary pairs (V, W): {len(pairs)}  (each Lagrangian has "
          f"{len(pairs) // len(lag)} complements; theory 2^6 = 64)")
    # orbit of (L.sbar/2L, L.s/2L) under the group generated by the 120 transvections
    # (= image of W(E8) = Aut(L) in O+(8,2)); if it is all ordered complementary pairs,
    # every polarization is (V_g, W_g) for some g in Aut(L) (a product of root reflections)
    pos = {s: k for k, s in enumerate(lag)}
    tperm = [[pos[frozenset(t[x] for x in s)] for s in lag] for t in trans]
    orbit = {(i0, iw)}
    frontier = [(i0, iw)]
    while frontier:
        nxt = []
        for (a, b) in frontier:
            for tp in tperm:
                im = (tp[a], tp[b])
                if im not in orbit:
                    orbit.add(im)
                    nxt.append(im)
        frontier = nxt
    print(f"orbit of the Wilson pair (L.sbar/2L, L.s/2L) under the reflection group: {len(orbit)} ordered pairs; "
          f"equal to the set of all ordered complementary pairs: {orbit == set(pairs)}")
    both = [(a, b) for a, b in pairs if info[a]['sub'] and info[b]['sub']]
    print(f"ordered complementary pairs with V and W both subalgebras: {len(both)}")

    jobs = [(a, b, sorted(lag[a]), sorted(lag[b])) for a, b in pairs]
    with Pool(14) as pool:
        res = pool.map(analyse_pair, jobs, chunksize=32)
    crit_ok = all(r['closes'] == (info[r['a']]['sub'] and info[r['b']]['sub']) for r in res)
    print(f"closure of * on Lambda(V,W) tested on all {len(res)} ordered pairs (576 products each); "
          f"closes <=> (V and W subalgebras) on every pair: {crit_ok}")
    print(f"  all {len(res)} lattices Lambda(V,W) have Gram det 2^24 and even Gram/2 diagonal: "
          f"{all(r['gram_ok'] for r in res)}")
    closing = [r for r in res if r['closes']]
    print(f"closing pairs: {len(closing)}")
    h = {}
    for r in closing:
        k = (r['index'], r['f2rank'], r['dimVV'], r['dimWW'], r['idx_S_2L3'], info[r['a']]['left'],
             info[r['a']]['right'], info[r['b']]['left'], info[r['b']]['right'])
        h[k] = h.get(k, 0) + 1
    print("  (index, F2-rank of S mod 2Lambda, dim VV, dim WW, [Lambda:S+(2L)^3], V left ideal, V right ideal,"
          " W left ideal, W right ideal) : number of ordered pairs")
    for k, v in sorted(h.items(), key=lambda kv: (kv[0][0], str(kv[0]))):
        print(f"    index {k[0]:>6d} = 2^{k[0].bit_length() - 1:<2d} F2-rank {k[1]:2d}  dimVV {k[2]}  dimWW {k[3]}  "
              f"[Lam:S+(2L)^3] {k[4]:>4d}  V(l,r) {k[5]!s:5s} {k[6]!s:5s}  W(l,r) {k[7]!s:5s} {k[8]!s:5s} : {v}")
    ih = {}
    for r in closing:
        ih[r['index']] = ih.get(r['index'], 0) + 1
    print("  span index distribution over the closing polarizations:",
          {f"2^{k.bit_length() - 1}": v for k, v in sorted(ih.items())})
    print(f"  every index a power of 2 with 2Lambda <= S: {all(r['two_Lambda_in_S'] for r in closing)}")
    pred = all(r['idx_S_2L3'] == 2 ** ((4 - r['dimWW']) + 2 * (4 - r['dimVV'])) for r in closing)
    print(f"  [Lambda : S + (2L)^3] = 2^((4 - dim WW) + 2(4 - dim VV)) on every closing pair: {pred}")

    # cross-check with the direct twisted census (03)
    lookup = {(lag[r['a']], lag[r['b']]): r for r in res}
    census = json.load(open('twist_results.json'))
    n = agree = 0
    for kind, items in (('perm', census['permutations']), ('refl', census['reflections'])):
        for it in items:
            if kind == 'perm':
                gf = (lambda pi: (lambda x: LL.perm_apply(pi, x)))(tuple(it['pi']))
            else:
                rr = [F(a, 2) for a in it['root_doubled_e_coords']]
                gf = (lambda rv: (lambda x: LL.refl_apply(rv, x)))(rr)
            G = np.array(LL.twist_matrix(gf, BL, LC), dtype=np.int64)
            Vg = span(bits(r) for r in (np.array(D['Vrows']) @ G) % 2)
            Wg = span(bits(r) for r in (np.array(D['Wrows']) @ G) % 2)
            pr = lookup[(Vg, Wg)]
            n += 1
            agree += (pr['closes'] == it['closes']) and (not it['closes'] or pr['index'] == it['span_index'])
    print(f"cross-check: for {agree} of {n} census twists (5,040 permutations + 120 reflections) the closure"
          f" and index of the polarization (V_g, W_g) equal the directly computed ones")
    json.dump({'lagrangians': [sorted(s) for s in lag], 'lagrangian_info': info,
               'pairs': res}, open('polarization_results.json', 'w'))
    print(f"total time {time.time() - t0:.1f} s")


def analyse_pair(job):
    a, b, Vs, Ws = job
    Vb, _ = LL.f2_rowspace_basis([vec(m) for m in Vs], 8)
    Wb, _ = LL.f2_rowspace_basis([vec(m) for m in Ws], 8)
    Z8 = [0] * 8
    gens = [w + w + w for w in Wb] + [v + Z8 + v for v in Vb] + [v + v + Z8 for v in Vb]
    gens += [[2 * int(i == j) for j in range(24)] for i in range(24)]
    BLam = LL.hnf(gens, 24)
    assert len(BLam) == 24
    detB = 1
    for i, h in enumerate(BLam):
        detB *= h[i]
    assert detB == 2 ** 12
    # Gram check in e-coordinates
    Be = [LC.to_e(r[0:8]) + LC.to_e(r[8:16]) + LC.to_e(r[16:24]) for r in BLam]
    Gm = [[LL.ip(x, y) for y in Be] for x in Be]
    gram_ok = (all(g.denominator == 1 and g % 2 == 0 for row in Gm for g in row)
               and LL.det_bareiss([[int(g) for g in row] for row in Gm]) == 2 ** 24
               and all((Gm[i][i] // 2) % 2 == 0 for i in range(24)))
    inv, den = LL.int_inverse_scaled([[F(x) for x in r] for r in BLam])
    assert all((2 * x) % den == 0 for row in inv for x in row)
    N2 = np.array([[2 * x // den for x in row] for row in inv], dtype=np.int64)
    Bnp = np.array(BLam, dtype=np.int64)
    P = LL.all_basis_products(Bnp, C).reshape(576, 24)
    ok, c2 = LL.lambda_coordinates(P, N2)
    out = {'a': a, 'b': b, 'closes': bool(ok.all()), 'gram_ok': gram_ok}
    Vset, Wset = set(Vs), set(Ws)
    out['dimVV'] = LL.f2_rank([vec(PROD[x][y]) for x in Vs for y in Vs])
    out['dimWW'] = LL.f2_rank([vec(PROD[x][y]) for x in Ws for y in Ws])
    if out['closes']:
        coords = c2 // 2
        rows = [list(map(int, r)) for r in np.unique(coords, axis=0) if np.any(r)]
        rows.sort(key=lambda r: (sum(abs(x) for x in r), r))
        r2 = LL.f2_rank(rows)
        idx, H = LL.hnf_index(rows, 24)
        idx2, _ = LL.hnf_index(rows + [list(map(int, r)) for r in N2], 24)
        out.update({'index': idx, 'f2rank': r2, 'two_Lambda_in_S': idx == 2 ** (24 - r2), 'idx_S_2L3': idx2})
    return out


if __name__ == '__main__':
    main()
