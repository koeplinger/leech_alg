r"""
03_twist_census.py -- closure and span index of the twisted triple product
                      *_g on Wilson's Leech lattice, for
   * sigma = (1 2)                                            (claim C1)
   * all 5,040 permutations pi of the imaginary units e1..e7
     (the identity is computed as a reference; the census is over the
      5,039 non-identity permutations)                        (claim C2)
   * all 120 root reflections s_r(x) = x - <x,r> r, one per pair +-r
                                                              (claim C3)

CONVENTIONS (see leechlib.py for the full list)
  * twisted product  x ._g y := g^{-1}( g(x) . g(y) )
  * permutation pi stored as (pi(1),...,pi(7)) with pi(e_a) = e_{pi(a)}
    (so (pi(x))_{pi(a)} = x_a); the paper's x ._pi y = pi^{-1}(pi(x).pi(y))
    is the case g = pi.  Cycle notation (a b c) means pi(a)=b, pi(b)=c, pi(c)=a.
  * the reflection s_r is an involution, so g = g^{-1} and no convention issue.
  * the representative of a pair +-r is the one whose first nonzero
    e-coordinate is positive.

METHOD, per twist g
  1. G = matrix of g in the Z-basis of L (asserted integral and unimodular,
     i.e. g(L) = L).  Twisted structure constants C_g[i,j,:] = L-coordinates
     of g^{-1}(g(l_i).g(l_j)) = (G_i (x) G_j : C) G^{-1}, all integers.
  2. The 576 products b_i *_g b_j of the Lambda basis, in L^3-coordinates.
  3. Closure, two independent tests that must agree product by product:
       (a) Lambda-coordinates c = P . BLam^{-1} are integral
           (computed as P @ (2 BLam^{-1}), all entries even);
       (b) the 12 F2 conditions Q (x+y, x+z, y+z in L.sbar/2L, x+y+z in L.s/2L).
  4. If closed: the 576 x 24 integer matrix of Lambda-coordinates generates
     S_g; exact HNF (pure Python ints) gives [Lambda : S_g] as the product of
     the pivots; SNF of the HNF gives the elementary divisors.  Also:
     r = rank over F2 of that matrix mod 2 (= dim (S_g + 2 Lambda)/2 Lambda),
     so [Lambda : S_g + 2 Lambda] = 2^(24-r) is a LOWER bound for the index,
     and 2 Lambda <= S_g  <=>  [Lambda : S_g] = 2^(24-r).
     Also [Lambda : S_g + (2L)^3].
  5. Mod-2 algebra data (closure criterion of claim C5): V_g = g(L.sbar)/2L and
     W_g = g(L.s)/2L as subspaces of the F2-algebra L/2L (product reduced from
     L.L <= L): subalgebra?  ideal?  dimension of span(V_g V_g), span(W_g W_g).

OUTPUTS
  twist_results.json        per-twist results (5,040 permutations, 120 reflections)
  03_twist_census.out.txt   human-readable summary (this script's stdout)
  pari_input_matrices.gp    576 x 24 Lambda-coordinate matrices of the products
                            for the PARI/GP cross-check set (04_pari_crosscheck.gp)
  hnf_<name>.txt            HNF of S_g for sigma and the two witnesses
Run:  python3 03_twist_census.py > 03_twist_census.out.txt
"""
import json
import random
import time
from fractions import Fraction as F
from itertools import permutations
from multiprocessing import Pool

import numpy as np
import leechlib as LL

D = LL.build_all()
QNP = np.array(D['Q'], dtype=np.int64)            # 16 x 24 (rank 12)
CMOD2 = (D['C'] % 2).astype(np.int64)             # product of L/2L
VROWS = np.array(D['Vrows'], dtype=np.int64)
WROWS = np.array(D['Wrows'], dtype=np.int64)
UNIT8 = np.eye(8, dtype=np.int64)


def cycles_str(pi):
    seen, out = set(), []
    for a in range(1, 8):
        if a in seen or pi[a - 1] == a:
            seen.add(a)
            continue
        cyc, b = [], a
        while b not in seen:
            seen.add(b)
            cyc.append(b)
            b = pi[b - 1]
        out.append('(' + ' '.join(map(str, cyc)) + ')')
    return ''.join(out) if out else '()'


def f2_basis(rows):
    basis, _ = LL.f2_rowspace_basis([list(map(int, r)) for r in rows], 8)
    return [np.array(b, dtype=np.int64) for b in basis]


def mod2_prod(a, b):
    return np.einsum('i,j,ijk->k', a, b, CMOD2) % 2


def contains(basis, vecs):
    """Are all vecs in the F2 span of basis?"""
    k = LL.f2_rank(basis)
    return all(LL.f2_rank(list(basis) + [v]) == k for v in vecs)


def mod2_data(G):
    Gm = np.array(G, dtype=np.int64)
    Vb = f2_basis((VROWS @ Gm) % 2)               # V_g = g(L.sbar)/2L
    Wb = f2_basis((WROWS @ Gm) % 2)               # W_g = g(L.s)/2L
    assert len(Vb) == 4 and len(Wb) == 4
    VV = [mod2_prod(a, b) for a in Vb for b in Vb]
    WW = [mod2_prod(a, b) for a in Wb for b in Wb]
    LV = [mod2_prod(e, b) for e in UNIT8 for b in Vb]
    VLp = [mod2_prod(b, e) for e in UNIT8 for b in Vb]
    LW = [mod2_prod(e, b) for e in UNIT8 for b in Wb]
    WLp = [mod2_prod(b, e) for e in UNIT8 for b in Wb]
    return {
        'V_g_subalgebra': contains(Vb, VV),
        'W_g_subalgebra': contains(Wb, WW),
        'V_g_left_ideal': contains(Vb, LV),
        'V_g_right_ideal': contains(Vb, VLp),
        'W_g_left_ideal': contains(Wb, LW),
        'W_g_right_ideal': contains(Wb, WLp),
        'dim_span_VgVg': LL.f2_rank(VV),
        'dim_span_WgWg': LL.f2_rank(WW),
    }


def gfun_of(spec):
    kind, data = spec
    if kind == 'perm':
        return lambda x: LL.perm_apply(data, x)
    r = [F(a, 2) for a in data]                   # data = doubled root
    return lambda x: LL.refl_apply(r, x)


def products_and_coords(spec):
    G = LL.twist_matrix(gfun_of(spec), D['BL'], D['LC'])
    Cg = LL.twisted_structure_constants(D['C'], G)
    P = LL.all_basis_products(D['BLam_np'], Cg).reshape(576, 24)   # row 24*i+j = b_i * b_j
    ok, coords2 = LL.lambda_coordinates(P, D['N2'])
    okQ = np.all((P @ QNP.T) % 2 == 0, axis=1)
    assert np.array_equal(ok, okQ), "closure tests (a) and (b) disagree"
    return G, Cg, P, ok, coords2


def analyse(spec):
    kind, data = spec
    G, Cg, P, ok, coords2 = products_and_coords(spec)
    res = {'kind': 'permutation' if kind == 'perm' else 'reflection'}
    if kind == 'perm':
        res['pi'] = list(data)
        res['cycles'] = cycles_str(data)
        ct = LL.cycle_type(data)
        res['cycle_type'] = LL.CYCLE_TYPE_NAMES[ct]
        res['parity'] = 'even' if LL.parity(data) == 0 else 'odd'
    else:
        res['root_doubled_e_coords'] = list(data)
        res['root'] = '(' + ', '.join(str(F(a, 2)) for a in data) + ')'
        res['root_type'] = 'integral (+-e_i +- e_j)' if all(a % 2 == 0 for a in data) else 'half-integral'
    res['closes'] = bool(ok.all())
    res['n_products_not_in_Lambda'] = int((~ok).sum())
    res.update(mod2_data(G))
    if res['closes']:
        coords = coords2 // 2
        rows = np.unique(coords, axis=0)
        rows = [list(map(int, r)) for r in rows if np.any(r)]
        rows.sort(key=lambda r: (sum(abs(a) for a in r), r))
        r2 = LL.f2_rank(rows)
        idx, H = LL.hnf_index(rows, 24)
        assert idx is not None, "S_g not of full rank"
        snf = LL.snf_diagonal(H)
        prod = 1
        for d in snf:
            prod *= d
        assert prod == idx
        # S_g + (2L)^3: (2L)^3 is spanned, in Lambda-coordinates, by the rows of N2
        idx2L, _ = LL.hnf_index(rows + [list(map(int, r)) for r in D['N2']], 24)
        res['span_index'] = idx
        res['log2_span_index'] = (idx.bit_length() - 1) if idx & (idx - 1) == 0 else None
        res['elementary_divisors_nontrivial'] = [d for d in snf if d != 1]
        res['f2_rank_S_mod_2Lambda'] = r2
        res['index_S_plus_2Lambda'] = 2 ** (24 - r2)
        res['two_Lambda_in_S'] = (idx == 2 ** (24 - r2))
        res['index_S_plus_2L3'] = idx2L
        res['n_distinct_nonzero_product_rows'] = len(rows)
    return res


def roots_reps():
    reps = []
    for r in LL.roots_of_L():
        first = next(a for a in r if a != 0)
        if first > 0:
            reps.append([int(2 * a) for a in r])
    assert len(reps) == 120
    return reps


def main():
    t0 = time.time()
    perms = [tuple(p) for p in permutations(range(1, 8))]
    roots = roots_reps()
    specs = [('perm', p) for p in perms] + [('refl', tuple(r)) for r in roots]
    with Pool(14) as pool:
        results = pool.map(analyse, specs, chunksize=16)
    out = {'conventions': {
        'twisted_product': 'x ._g y := g^{-1}( g(x) . g(y) )',
        'permutation': 'pi stored as [pi(1),...,pi(7)], pi(e_a) = e_{pi(a)}, pi fixes e0; cycle (a b c): a->b->c->a',
        'reflection': 's_r(x) = x - <x,r> r, r a root of L (norm 2); representative of +-r has first nonzero coordinate > 0',
        'triple_product': 'P1 = x.x\' + z.y\' + y.z\', P2 = y.y\' + x.z\' + z.x\', P3 = z.z\' + y.x\' + x.y\' (all ._g)',
        'span_index': '[Lambda : Z-span of the 576 basis-pair products]; exact HNF',
        'f2_rank_S_mod_2Lambda': 'dim_F2 (S + 2 Lambda)/2 Lambda',
        'index_S_plus_2L3': '[Lambda : S + (2L)^3]',
        'V_g, W_g': 'g(L.sbar)/2L and g(L.s)/2L inside the F2-algebra L/2L',
    }, 'permutations': [], 'reflections': []}
    for spec, res in zip(specs, results):
        (out['permutations'] if spec[0] == 'perm' else out['reflections']).append(res)
    with open('twist_results.json', 'w') as fh:
        json.dump(out, fh, indent=1)

    P = out['permutations']
    R = out['reflections']
    byperm = {tuple(r['pi']): r for r in P}

    print("=" * 76)
    print("C1: sigma = (1 2), i.e. pi = (2,1,3,4,5,6,7)")
    print("=" * 76)
    sg = byperm[(2, 1, 3, 4, 5, 6, 7)]
    for k in ('closes', 'n_products_not_in_Lambda', 'span_index', 'log2_span_index',
              'elementary_divisors_nontrivial', 'f2_rank_S_mod_2Lambda', 'two_Lambda_in_S',
              'index_S_plus_2L3', 'V_g_subalgebra', 'W_g_subalgebra', 'V_g_left_ideal',
              'V_g_right_ideal', 'dim_span_VgVg', 'dim_span_WgWg'):
        print(f"  {k:32s} {sg[k]}")
    ident = byperm[(1, 2, 3, 4, 5, 6, 7)]
    print(f"  [reference] identity (untwisted): closes = {ident['closes']}, "
          f"products not in Lambda = {ident['n_products_not_in_Lambda']}/576, "
          f"V subalgebra {ident['V_g_subalgebra']}, W subalgebra {ident['W_g_subalgebra']}")

    print()
    print("=" * 76)
    print("C2: all 5,039 non-identity permutations of e1..e7")
    print("=" * 76)
    paper = {'2-cycle': (21, 21), '3-cycle': (70, 35), '(2,2)': (105, 42), '4-cycle': (210, 0),
             '(3,2)': (420, 126), '5-cycle': (504, 252), '(2,2,2)': (105, 21), '(3,3)': (280, 112),
             '(4,2)': (630, 294), '6-cycle': (840, 210), '(3,2,2)': (210, 105), '(4,3)': (420, 168),
             '(5,2)': (504, 42), '7-cycle': (720, 336)}
    print(f"  {'cycle type':10s} {'parity':6s} {'total':>6s} {'close':>6s} {'paper':>11s}  span indices of closing twists")
    tot = cl = 0
    all_match = True
    for ctn, (pt, pc) in paper.items():
        rows = [r for r in P if r['cycle_type'] == ctn]
        c = [r for r in rows if r['closes']]
        hist = {}
        for r in c:
            hist[r['span_index']] = hist.get(r['span_index'], 0) + 1
        hs = ', '.join(f"2^{k.bit_length()-1}: {v}" if k & (k - 1) == 0 else f"{k}: {v}"
                       for k, v in sorted(hist.items()))
        match = (len(rows), len(c)) == (pt, pc)
        all_match &= match
        print(f"  {ctn:10s} {rows[0]['parity']:6s} {len(rows):6d} {len(c):6d} {str((pt, pc)):>11s}  {hs}"
              f"{'' if match else '   <-- MISMATCH'}")
        tot += len(rows)
        cl += len(c)
    print(f"  {'total':10s} {'':6s} {tot:6d} {cl:6d}   all cycle types match the paper's table: {all_match}")
    for par in ('even', 'odd'):
        c = [r for r in P if r['closes'] and r['parity'] == par and r['cycles'] != '()']
        hist = {}
        for r in c:
            key = (r['span_index'], tuple(r['elementary_divisors_nontrivial']), r['two_Lambda_in_S'])
            hist[key] = hist.get(key, 0) + 1
        print(f"  closing {par} permutations: {len(c)};  (index, elementary divisors != 1, 2Lambda<=S) -> count:")
        for (k, e, t), v in sorted(hist.items()):
            print(f"      index {k} = 2^{k.bit_length()-1}, divisors {{2^1 x {len(e)}}} "
                  f"{'(all = 2)' if all(x == 2 for x in e) else str(e)}, 2Lambda<=S: {t}  ->  {v}")
    wit = (1, 2, 3, 4, 6, 7, 5)          # e5 -> e6 -> e7 -> e5
    witinv = LL.perm_inverse(wit)
    for name, p in (('C2 witness pi = (5 6 7), i.e. pi(5)=6, pi(6)=7, pi(7)=5', wit),
                    ('its inverse (5 7 6)', witinv)):
        r = byperm[p]
        print(f"  {name}: pi = {p}, closes = {r['closes']}, index = {r.get('span_index')}"
              f" (2^{r.get('log2_span_index')}), F2-rank = {r.get('f2_rank_S_mod_2Lambda')}, "
              f"2Lambda<=S: {r.get('two_Lambda_in_S')}")
    # closure criterion check (C5) on all permutations
    crit = all(r['closes'] == (r['V_g_subalgebra'] and r['W_g_subalgebra']) for r in P)
    print(f"  closure <=> (V_g and W_g subalgebras of L/2L) on all 5,040 permutations: {crit}")
    three = [r for r in P if r['cycle_type'] == '3-cycle']
    sets = {}
    for r in three:
        supp = tuple(sorted(a for a in range(1, 8) if r['pi'][a - 1] != a))
        sets.setdefault(supp, []).append(r['closes'])
    print(f"  3-cycles: number of 3-subsets where exactly one of the two orientations closes: "
          f"{sum(1 for v in sets.values() if sorted(v) == [False, True])} of {len(sets)}")

    print()
    print("=" * 76)
    print("C3: the 120 root reflections s_r")
    print("=" * 76)
    hist = {}
    for r in R:
        key = ('closes' if r['closes'] else 'does not close', r.get('span_index'), r['root_type'])
        hist[key] = hist.get(key, 0) + 1
    for k, v in sorted(hist.items(), key=lambda kv: (kv[0][0], kv[0][1] or 0, kv[0][2])):
        print(f"  {k[0]:15s} index {str(k[1]):>6s}  {k[2]:26s} : {v}")
    for target, label in ((1, 'index 1'), (2 ** 16, 'index 2^16')):
        lst = [r for r in R if r['closes'] and r['span_index'] == target]
        print(f"  roots r with {label} ({len(lst)}):")
        for r in lst:
            print(f"     r = {r['root']:44s} F2-rank {r['f2_rank_S_mod_2Lambda']:2d}, "
                  f"divisors {r['elementary_divisors_nontrivial'] if len(r['elementary_divisors_nontrivial']) < 3 else '2^1 x ' + str(len(r['elementary_divisors_nontrivial']))}, "
                  f"dim VgVg {r['dim_span_VgVg']}, WgWg {r['dim_span_WgWg']}, [Lam:S+(2L)^3] {r['index_S_plus_2L3']}")
    other = sorted({r.get('span_index') for r in R if r['closes']} - {1, 2 ** 16})
    print(f"  other indices among closing reflections: {other}")
    wroot = (1, 1, 1, 1, 1, 1, 1, -1)
    wr = next(r for r in R if tuple(r['root_doubled_e_coords']) == wroot)
    print(f"  C3 witness r = (e0+...+e6-e7)/2: closes = {wr['closes']}, index = {wr.get('span_index')}, "
          f"F2-rank = {wr.get('f2_rank_S_mod_2Lambda')}")
    s12 = next(r for r in R if tuple(r['root_doubled_e_coords']) == (0, 2, -2, 0, 0, 0, 0, 0))
    print(f"  reflection in e1-e2 (= sigma): closes = {s12['closes']}, index = {s12.get('span_index')}")
    crit = all(r['closes'] == (r['V_g_subalgebra'] and r['W_g_subalgebra']) for r in R)
    print(f"  closure <=> (V_g and W_g subalgebras of L/2L) on all 120 reflections: {crit}")

    print()
    print("=" * 76)
    print("Mod-2 structure of the closing twists: [Lambda : S + (2L)^3] vs dims of span(VgVg), span(WgWg)")
    print("=" * 76)
    pred_ok = True
    tab = {}
    for r in P + R:
        if not r['closes']:
            continue
        pred = 2 ** ((4 - r['dim_span_WgWg']) + 2 * (4 - r['dim_span_VgVg']))
        pred_ok &= (pred == r['index_S_plus_2L3'])
        key = (r['kind'], r.get('parity', r.get('root_type')), r['span_index'], r['f2_rank_S_mod_2Lambda'],
               r['dim_span_VgVg'], r['dim_span_WgWg'], r['index_S_plus_2L3'],
               r['V_g_left_ideal'] and r['V_g_right_ideal'], r['W_g_left_ideal'] and r['W_g_right_ideal'])
        tab[key] = tab.get(key, 0) + 1
    print("  prediction [Lambda : S+(2L)^3] = 2^((4 - dim WgWg) + 2 (4 - dim VgVg)) holds for all closing twists:", pred_ok)
    print("  (kind, parity/root type, index, F2-rank, dim VgVg, dim WgWg, [Lam:S+(2L)^3], Vg 2-sided ideal, Wg 2-sided ideal) : count")
    for k, v in sorted(tab.items(), key=lambda kv: str(kv[0])):
        print("   ", k, ":", v)

    # ------------------------------------------------------------------
    # export matrices for the PARI/GP cross-check
    # ------------------------------------------------------------------
    random.seed(13102026)
    closing_perms = [r for r in P if r['closes'] and tuple(r['pi']) not in
                     ((2, 1, 3, 4, 5, 6, 7), wit, witinv)]
    rnd = random.sample(closing_perms, 50)
    # (the inverse (5 7 6) of the C2 witness does not close, so it is not exported)
    chosen = [('sigma (1 2)', ('perm', (2, 1, 3, 4, 5, 6, 7))),
              ('C2 witness (5 6 7)', ('perm', wit)),
              ('C3 witness r=(1,1,1,1,1,1,1,-1)/2', ('refl', wroot))]
    for r in R:
        if r['closes'] and tuple(r['root_doubled_e_coords']) != wroot:
            chosen.append((f"reflection r={r['root']}", ('refl', tuple(r['root_doubled_e_coords']))))
    for r in rnd:
        chosen.append((f"perm {r['cycles']} pi={tuple(r['pi'])}", ('perm', tuple(r['pi']))))
    with open('pari_input_matrices.gp', 'w') as fh:
        fh.write('\\\\ 576 x 24 matrices: row 24*(i-1)+j = Lambda-coordinates of b_i *_g b_j\n')
        fh.write(f'NAMES = [{", ".join(chr(34) + n + chr(34) for n, _ in chosen)}];\n')
        fh.write(f'MATS = vector({len(chosen)});\n')
        for k, (name, spec) in enumerate(chosen, 1):
            G, Cg, Pm, ok, coords2 = products_and_coords(spec)
            assert ok.all()
            coords = coords2 // 2
            fh.write(f'MATS[{k}] = [' + ';'.join(','.join(str(int(a)) for a in row) for row in coords) + '];\n')
    with open('pari_crosscheck_names.json', 'w') as fh:
        json.dump([[n, s[0], list(s[1])] for n, s in chosen], fh, indent=1)
    print()
    print(f"exported {len(chosen)} product matrices for the PARI/GP cross-check "
          f"(sigma, C2 witness, C3 witness, the other 55 closing reflections, "
          f"50 random closing permutations, seed 13102026)")

    # HNFs of S for the three named twists
    for fname, spec in (('hnf_sigma.txt', ('perm', (2, 1, 3, 4, 5, 6, 7))),
                        ('hnf_C2_witness_567.txt', ('perm', wit)),
                        ('hnf_C3_witness_root.txt', ('refl', wroot))):
        G, Cg, Pm, ok, coords2 = products_and_coords(spec)
        coords = coords2 // 2
        rows = [list(map(int, r)) for r in coords]
        idx, H = LL.hnf_index(rows, 24)
        with open(fname, 'w') as fh:
            fh.write(f'# HNF (rows) of S_g in Lambda-coordinates for {spec}; index = product of pivots = {idx}\n')
            for h in H:
                fh.write(' '.join(f'{a:3d}' for a in h) + '\n')
    print(f"total time {time.time() - t0:.1f} s")


if __name__ == '__main__':
    main()
