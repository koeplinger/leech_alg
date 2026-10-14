"""
s4_twists.py -- direct computation of closure and span index for the twist
families of claims C2-C4, on Wilson's Lambda with the TWISTED product, and a
comparison with the pair data of s3_lagrangians.py.

Conventions: leech_core.py.  For a twist g (8x8 rational matrix on column
e-coordinates):   x ._g y = g^{-1}(g(x).g(y)),  P1 = x.x'+z.y'+y.z', ...
  * permutations:  convention A: g = g_pi, g_pi(e_a) = e_{pi(a)}  (the paper's
    x ._pi y = pi^{-1}(pi(x).pi(y)));  convention B: g = g_pi^{-1}.
  * root reflections: g = s_r = I - r r^T, one root r per pair +-r.
  * signed permutations of e1..e7 (e0 fixed): g(e_a) = eps_a e_{pi(a)}.
Direct method for each twist: integral twisted structure tensor on L, the 576
products of the Z-basis of Lambda (s1), exact coordinates in that basis
(closure <=> all integral), span index = product of the HNF diagonal of the
576 x 24 coordinate matrix (leech_core.span_index; exact).
For g in O(L) the pair (V_g, W_g) = (g(Lsbar)/2L, g(Ls)/2L) is looked up in
s3_lagrangians.json, and we check: closure <=> V_g, W_g subalgebras; index ==
h(W_g) h(V_g)^2; index == rule(eps) where eps = class of 2e0 in L/2L:
    eps in V_g -> 2^(2k),  eps in W_g -> 2^k,  else 1,
    k = 5 if det g = +1 (V_g in the family of Lsbar/2L), k = 8 if det g = -1.
Also: W(E8) is transitive on the 17,280 ordered complementary pairs (BFS with
the 120 reflections), and explicit twists realizing every index value.
Outputs: s4_twists.out, s4_permutations.tsv, s4_reflections.tsv,
s4_signed_tables.tsv, s4_index_realizers.json
"""
import itertools
import json
import time
from collections import Counter, deque
from fractions import Fraction as Fr
from multiprocessing import Pool
import numpy as np
import leech_core as C

D3 = json.load(open('s3_lagrangians.json'))
LAG = [frozenset(U) for U in D3['lagrangians']]
LAGIDX = {U: i for i, U in enumerate(LAG)}
PROPS = D3['props']
s, sb = C.S_VEC, C.SBAR_VEC
H_sb = C.hnf([C.to_L_int(C.omul(C.L_BASIS[k], sb)) for k in range(8)])
H_s = C.hnf([C.to_L_int(C.omul(C.L_BASIS[k], s)) for k in range(8)])
V0 = C.span_bits([C.to_bits(v) for v in H_sb])
W0 = C.span_bits([C.to_bits(v) for v in H_s])
U0 = np.array(C.leech_basis_from_pair(V0, W0), dtype=np.int64)
EPS = C.to_bits(C.to_L_int([2, 0, 0, 0, 0, 0, 0, 0]))     # class of 2e0 in A


def rank_mod2(M):
    rows = [int("".join(str(int(t) & 1) for t in row), 2) for row in M]
    rows = [r for r in rows if r]
    rank = 0
    while rows:
        piv = max(rows)
        top = piv.bit_length() - 1
        rows = [r ^ piv if (r >> top) & 1 else r for r in rows if r != piv]
        rows = [r for r in rows if r]
        rank += 1
    return rank


def direct(G):
    """Closure and span index of the g-twisted triple product on Lambda."""
    T = C.tensor_int(C.structure_tensor_L(G))
    P = C.triple_products(U0, T)
    Cc, ok = C.coords_in_basis(P, U0)
    res = dict(closed=bool(ok.all()), n_fail=int((~ok).sum()))
    # consistency of the exact solve with the mod-2 membership test
    assert np.array_equal(ok, C.in_pair_lattice_mod2(P, V0, W0))
    if res['closed']:
        idx, piv = C.span_index(Cc)
        r2 = rank_mod2(Cc)
        res.update(index=idx, log2=idx.bit_length() - 1, rank2=r2, twoLam_in_S=(idx == 2 ** (24 - r2)))
    return res


def pair_data(G):
    """(V_g, W_g) and their properties, for g in O(L)."""
    GL = C.L_matrix(G)
    detg = int(C.mat_det(GL))
    img = lambda S: C.span_bits([C.to_bits(C.vecmat(C.from_bits(v), GL)) for v in C.basis_bits(S)])
    V, W = img(V0), img(W0)
    iV, iW = LAGIDX[V], LAGIDX[W]
    pv, pw = PROPS[iV], PROPS[iW]
    k = 5 if detg == 1 else 8
    rule = 2 ** (2 * k) if EPS in V else (2 ** k if EPS in W else 1)
    return dict(det=detg, iV=iV, iW=iW, V_sub=pv['subalgebra'], W_sub=pw['subalgebra'],
                V_left=pv['left_ideal'], V_right=pv['right_ideal'], W_left=pw['left_ideal'],
                W_right=pw['right_ideal'], V_2s=pv['two_sided'], W_2s=pw['two_sided'],
                dimVV=pv['dimUU'], dimWW=pw['dimUU'], hV=pv['h'], hW=pw['h'],
                family=pv['family'], eps_in_V=EPS in V, eps_in_W=EPS in W, rule_index=rule,
                formula_index=(pw['h'] * pv['h'] ** 2) if (pv['subalgebra'] and pw['subalgebra']) else None)


def full(G):
    r = direct(G)
    r.update(pair_data(G))
    r['crit_ok'] = (r['closed'] == (r['V_sub'] and r['W_sub']))
    if r['closed']:
        r['formula_ok'] = (r['index'] == r['formula_index'])
        r['rule_ok'] = (r['index'] == r['rule_index'])
    return r


def cycle_type(pi):
    seen, ct = set(), []
    for a in range(1, 8):
        if a in seen:
            continue
        n, b = 0, a
        while b not in seen:
            seen.add(b)
            b = pi[b]
            n += 1
        if n > 1:
            ct.append(n)
    return tuple(sorted(ct, reverse=True))


def parity(pi):
    return sum(n - 1 for n in cycle_type(pi)) % 2


def perm_job(p):
    pi = {a: p[a - 1] for a in range(1, 8)}
    inv = {pi[a]: a for a in pi}
    rA = full(C.perm_matrix_e(pi))
    rB = full(C.perm_matrix_e(inv))
    return pi, rA, rB


def refl_job(r):
    return full(C.reflection_matrix_e(r))


CT_NAME = {(): 'identity', (2,): '2-cycle', (3,): '3-cycle', (2, 2): '(2,2)', (4,): '4-cycle',
           (3, 2): '(3,2)', (5,): '5-cycle', (2, 2, 2): '(2,2,2)', (3, 3): '(3,3)', (4, 2): '(4,2)',
           (6,): '6-cycle', (3, 2, 2): '(3,2,2)', (4, 3): '(4,3)', (5, 2): '(5,2)', (7,): '7-cycle'}

if __name__ == '__main__':
    t0 = time.time()
    lines = []

    def log(*a):
        s_ = " ".join(str(t) for t in a)
        print(s_, flush=True)
        lines.append(s_)

    pool = Pool(14)
    # ------------------------------------------------------------------ perms
    perms = list(itertools.permutations(range(1, 8)))
    R = pool.map(perm_job, perms, chunksize=40)
    log("PERMUTATIONS (%d incl. identity), %.0fs" % (len(R), time.time() - t0))
    for conv, key in (('A (g = pi, the paper)', 1), ('B (g = pi^-1)', 2)):
        rs = [(x[0], x[key]) for x in R]
        tab = Counter()
        for pi, r in rs:
            ct = cycle_type(pi)
            tab[(ct, 'n')] += 1
            if r['closed']:
                tab[(ct, 'c')] += 1
        log("  convention %s: closing permutations = %d (identity closes: %s)" %
            (conv, sum(r['closed'] for _, r in rs), rs[0][1]['closed']))
        log("    cycle type : closing / total :", ", ".join(
            "%s %d/%d" % (CT_NAME[ct], tab[(ct, 'c')], tab[(ct, 'n')])
            for ct in sorted({cycle_type(p) for p, _ in rs}, key=lambda c: (sum(c), c)) if ct))
        idx = Counter((parity(pi), r['log2']) for pi, r in rs if r['closed'])
        log("    (parity, log2 index) of closing permutations:", dict(sorted(idx.items())))
        log("    criterion agrees with direct closure on", sum(r['crit_ok'] for _, r in rs), "of", len(rs),
            "| formula h(W)h(V)^2 ok on", sum(r['formula_ok'] for _, r in rs if r['closed']),
            "| eps-rule ok on", sum(r['rule_ok'] for _, r in rs if r['closed']), "| 2Lambda in S on",
            sum(r['twoLam_in_S'] for _, r in rs if r['closed']), "of", sum(r['closed'] for _, r in rs))
        log("    eps = [2e0] lies in V_g for every permutation:", all(r['eps_in_V'] for _, r in rs),
            "| det g = sign(pi) for all:", all(r['det'] == (-1) ** parity(pi) for pi, r in rs))
        prop = Counter((parity(pi), r['closed'], r['V_sub'], r['V_left'], r['V_right'], r['W_sub'],
                        r['W_left'], r['W_right']) for pi, r in rs)
        log("    (parity, closes, V sub, V left, V right, W sub, W left, W right) : count")
        for k_ in sorted(prop, key=str):
            log("       ", k_, ":", prop[k_])
        log("    distinct (V_g, W_g) pairs among all permutations: %d; among closing ones: %d" %
            (len({(r['iV'], r['iW']) for _, r in rs}), len({(r['iV'], r['iW']) for _, r in rs if r['closed']})))
    # witness of C2
    w = next(x for x in R if x[0] == {1: 1, 2: 2, 3: 3, 4: 4, 5: 6, 6: 7, 7: 5})
    log("  C2 witness pi(5)=6, pi(6)=7, pi(7)=5:  convention A closes %s index %s | convention B closes %s index %s" %
        (w[1]['closed'], w[1].get('index'), w[2]['closed'], w[2].get('index')))
    sg = next(x for x in R if x[0] == {1: 2, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7})
    log("  sigma = (1 2): closes %s index %s, V_sigma two-sided ideal %s, dim V.V %s, h(V) %s, h(W) %s" %
        (sg[1]['closed'], sg[1]['index'], sg[1]['V_2s'], sg[1]['dimVV'], sg[1]['hV'], sg[1]['hW']))
    with open('s4_permutations.tsv', 'w') as f:
        f.write("# convention A: g = g_pi (paper), B: g = g_pi^-1; V/W columns refer to convention A\n")
        f.write("pi(1..7)\tcycle_type\tparity\tA_closed\tA_log2index\tB_closed\tB_log2index\tA_V_sub\tA_V_left"
                "\tA_V_right\tA_V_2sided\tA_W_sub\tA_W_left\tA_W_right\tA_W_2sided\tA_dimVV\tA_dimWW\tA_hV\tA_hW"
                "\tA_eps_in_V\tA_eps_in_W\n")
        for pi, rA, rB in R:
            f.write("\t".join(map(str, ["".join(str(pi[a]) for a in range(1, 8)), CT_NAME[cycle_type(pi)], parity(pi),
                                         rA['closed'], rA.get('log2', ''), rB['closed'], rB.get('log2', ''),
                                         rA['V_sub'], rA['V_left'], rA['V_right'], rA['V_2s'],
                                         rA['W_sub'], rA['W_left'], rA['W_right'], rA['W_2s'],
                                         rA['dimVV'], rA['dimWW'], rA['hV'], rA['hW'],
                                         rA['eps_in_V'], rA['eps_in_W']])) + "\n")

    # ------------------------------------------------------------ reflections
    roots = C.roots_of_L()
    reps = []
    for r in roots:
        first = next(t for t in r if t != 0)
        if first > 0:
            reps.append(r)
    assert len(reps) == 120
    RR = pool.map(refl_job, reps)

    def rtype(r):
        if all(t.denominator == 1 for t in r):
            return 'int, e0 not involved' if r[0] == 0 else 'int, +-e0+-ej'
        return 'half-integer'
    log("\nROOT REFLECTIONS (120, one per +-r), %.0fs" % (time.time() - t0))
    log("  closing:", sum(x['closed'] for x in RR), "| log2 index distribution:",
        dict(Counter(x.get('log2', 'no close') for x in RR)))
    tab = Counter((rtype(r), x.get('log2', 'no close')) for r, x in zip(reps, RR))
    for k_ in sorted(tab, key=str):
        log("    ", k_, ":", tab[k_])
    tab = Counter((x['eps_in_V'], x['eps_in_W'], x.get('log2', 'no close')) for x in RR)
    log("  (eps in V_g, eps in W_g, log2 index):", dict(sorted(tab.items(), key=str)))
    prop = Counter((x['closed'], x.get('log2'), x['V_sub'], x['V_left'], x['V_right'], x['W_sub'], x['W_left'],
                    x['W_right']) for x in RR)
    log("  (closes, log2 index, V sub, V left, V right, W sub, W left, W right) : count")
    for k_ in sorted(prop, key=str):
        log("       ", k_, ":", prop[k_])
    log("  criterion ok:", sum(x['crit_ok'] for x in RR), "/ 120 | formula ok:",
        sum(x['formula_ok'] for x in RR if x['closed']), "| eps-rule ok:", sum(x['rule_ok'] for x in RR if x['closed']),
        "| 2Lambda in S:", sum(x['twoLam_in_S'] for x in RR if x['closed']), "of", sum(x['closed'] for x in RR))
    w = C.reflection_matrix_e([Fr(1, 2)] * 7 + [Fr(-1, 2)])
    wr = full(w)
    unit = C.matvec(w, [1, 0, 0, 0, 0, 0, 0, 0])
    log("  C3 witness r = (e0+...+e6-e7)/2: closes %s, index %s, V_g subalgebra %s (two-sided %s, dim VV %s),"
        " W_g subalgebra %s (dim WW %s), eps in V_g %s, eps in W_g %s" %
        (wr['closed'], wr.get('index'), wr['V_sub'], wr['V_2s'], wr['dimVV'], wr['W_sub'], wr['dimWW'],
         wr['eps_in_V'], wr['eps_in_W']))
    log("    identity element of the twisted algebra: g^{-1}(e0) = s_r(e0) =", [str(t) for t in unit],
        "| in L:", C.in_L_e(unit), "| 2*it in L:", C.in_L_e([2 * t for t in unit]))
    with open('s4_reflections.tsv', 'w') as f:
        f.write("root_r(e-coords)\ttype\tclosed\tlog2_index\tV_sub\tV_left\tV_right\tV_2sided\tW_sub\tW_left"
                "\tW_right\tW_2sided\tdimVV\tdimWW\thV\thW\teps_in_V\teps_in_W\n")
        for r, x in zip(reps, RR):
            f.write("\t".join(map(str, [",".join(str(t) for t in r), rtype(r), x['closed'], x.get('log2', ''),
                                         x['V_sub'], x['V_left'], x['V_right'], x['V_2s'],
                                         x['W_sub'], x['W_left'], x['W_right'], x['W_2s'], x['dimVV'], x['dimWW'],
                                         x['hV'], x['hW'], x['eps_in_V'], x['eps_in_W']])) + "\n")

    # -------------------------------------------------- signed permutations
    log("\nSIGNED PERMUTATIONS of e1..e7 (e0 fixed): 2^7 * 7! = %d" % (128 * 5040))
    eps_all = np.array([[1] + [1 - 2 * ((m >> i) & 1) for i in range(7)] for m in range(128)], dtype=np.int64)
    tables = {}
    std_key = None
    aut = []
    for p in perms:
        pi = [0] + list(p)
        inv = [0] * 8
        for a in range(8):
            inv[pi[a]] = a
        # e_a ._g e_b = eps_a eps_b sgn(pi a, pi b) eps_d e_d,  with e_{pi a} e_{pi b} = sgn e_{pi d}
        idx = np.zeros((8, 8), dtype=np.int64)
        base = np.zeros((8, 8), dtype=np.int64)
        for a in range(8):
            for b in range(8):
                sg_, c = C.MULT[pi[a]][pi[b]]
                idx[a, b] = inv[c]
                base[a, b] = sg_
        signs = (eps_all[:, :, None] * eps_all[:, None, :] * base[None] *
                 eps_all[:, idx.reshape(-1)].reshape(128, 8, 8))
        keys = signs * 8 + idx[None]           # encodes (sign, index) per entry
        for m in range(128):
            kk = keys[m].tobytes()
            if kk not in tables:
                tables[kk] = (tuple(pi), tuple(int(t) for t in eps_all[m]))
            if tuple(pi) == tuple(range(8)) and m == 0:
                std_key = kk
            if kk == std_key:
                aut.append((tuple(pi), tuple(int(t) for t in eps_all[m])))
    # automorphisms found only after std_key is set (identity is the first permutation, m = 0 first)
    log("  distinct twisted tables:", len(tables))
    chi = Counter(int(np.prod(e)) for _, e in aut)
    log("  signed permutations preserving the standard table (automorphisms):", len(aut),
        "| product of signs on them:", dict(chi))
    jobs = []
    for kk, (pi, e) in tables.items():
        jobs.append((pi, e))

    def chi_of(e):
        return int(np.prod(e))
    res_tables = []
    for pi, e in jobs:
        G = C.signed_perm_matrix_e({a: pi[a] for a in range(1, 8)}, list(e))
        res_tables.append((pi, e, G))
    TR = pool.map(direct, [G for _, _, G in res_tables])
    n_even = sum(1 for pi, e, _ in res_tables if chi_of(e) == 1)
    log("  representatives with even number of sign changes (g in W(D8), g(L) = L): %d; odd: %d" %
        (n_even, len(res_tables) - n_even))
    closing = [(pi, e, r) for (pi, e, _), r in zip(res_tables, TR) if r['closed']]
    log("  closing tables:", len(closing), "| log2 index distribution:", dict(Counter(r['log2'] for _, _, r in closing)))
    log("  by sign parity:", dict(Counter((chi_of(e), r['closed'], r.get('log2')) for (pi, e, _), r in zip(res_tables, TR))))
    # opposite relation: for odd chi, the table of g is the opposite of the table of g' = conj o g,
    # and u *_g v = v *_{g'} u; check closure/index agree with g'
    agree = 0
    n_odd = 0
    pdat = Counter()
    for (pi, e, G), r in zip(res_tables, TR):
        if chi_of(e) == 1:
            pd = pair_data(G)
            if r['closed']:
                pdat[(pd['det'], pd['eps_in_V'], r['log2'], r['index'] == pd['rule_index'])] += 1
            continue
        n_odd += 1
        e2 = [1] + [-t for t in e[1:]]
        G2 = C.signed_perm_matrix_e({a: pi[a] for a in range(1, 8)}, e2)
        r2 = direct(G2)
        agree += (r2['closed'] == r['closed'] and r2.get('index') == r.get('index'))
    log("  odd-sign tables whose closure and index equal those of g' = conj o g:", agree, "of", n_odd)
    log("  even-sign closing tables by (det g, eps in V_g, log2 index, eps-rule ok):", dict(pdat))
    with open('s4_signed_tables.tsv', 'w') as f:
        f.write("pi(1..7)\teps(0..7)\tsign_product\tclosed\tlog2_index\n")
        for (pi, e, _), r in zip(res_tables, TR):
            f.write("\t".join(map(str, ["".join(map(str, pi[1:])), ",".join(map(str, e)), chi_of(e),
                                         r['closed'], r.get('log2', '')])) + "\n")

    # ---------------------------------------- W(E8) orbit of (V0, W0) by BFS
    acts = []
    for r in reps:
        GL = C.L_matrix(C.reflection_matrix_e(r))
        perm = [C.to_bits(C.vecmat(C.from_bits(v), GL)) for v in range(256)]
        acts.append(perm)
    start = (LAGIDX[V0], LAGIDX[W0])
    parent = {start: None}
    dq = deque([start])
    while dq:
        cur = dq.popleft()
        for ri, perm in enumerate(acts):
            nV = LAGIDX[frozenset(perm[v] for v in LAG[cur[0]])]
            nW = LAGIDX[frozenset(perm[v] for v in LAG[cur[1]])]
            nxt = (nV, nW)
            if nxt not in parent:
                parent[nxt] = (cur, ri)
                dq.append(nxt)
    log("\nW(E8) ORBIT of (Lsbar/2L, Ls/2L) under the 120 reflections: %d ordered pairs (all complementary pairs: 17280)"
        % len(parent))

    def word(state):
        w_ = []
        while parent[state] is not None:
            state, ri = parent[state]
            w_.append(ri)
        return w_[::-1]      # apply reflections in this order: g = s_{w[-1]} ... s_{w[0]}
    pairres = {(p['iV'], p['iW']): p for p in D3['pairs']}
    realizers = {}
    for target in (0, 5, 8, 10, 16):
        cands = [st for st in parent if pairres[st]['closed'] and pairres[st]['log2'] == target]
        st = min(cands, key=lambda s_: (len(word(s_)), s_))
        wd = word(st)
        G = [[Fr(int(i == j)) for j in range(8)] for i in range(8)]
        for ri in wd:
            G = C.matmul(C.reflection_matrix_e(reps[ri]), G)
        r = full(G)
        realizers[target] = dict(reflections=[[str(t) for t in reps[ri]] for ri in wd], closed=r['closed'],
                                 log2=r.get('log2'), det=r['det'], eps_in_V=r['eps_in_V'], eps_in_W=r['eps_in_W'])
        log("  index 2^%d realized by g = product of %d reflection(s) %s: direct closes %s, log2 index %s, det %d,"
            " eps in V %s, eps in W %s" % (target, len(wd), [[str(t) for t in reps[ri]] for ri in wd], r['closed'],
                                           r.get('log2'), r['det'], r['eps_in_V'], r['eps_in_W']))
    json.dump(realizers, open('s4_index_realizers.json', 'w'), indent=1)
    log("total time %.0fs" % (time.time() - t0))
    with open('s4_twists.out', 'w') as f:
        f.write("\n".join(lines) + "\n")
