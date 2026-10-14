"""
s3_lagrangians.py -- the 270 Lagrangians of (A = L/2L, q) and all 17,280
ordered complementary pairs (V, W), with closure and span index computed
DIRECTLY for every pair.

Conventions: leech_core.py.  A = L/2L with the product induced by the
standard octonion product (L.L in L), q(v) = N(v)/2 mod 2.

For an ordered complementary pair (V, W) of Lagrangians let V~, W~ be the
full preimages in L (they contain 2L) and
    Lambda_{V,W} = {(x,y,z) in L^3 : x+y, x+z, y+z in V~ ; x+y+z in W~}.
Wilson's Lambda is Lambda_{V0,W0} with V0 = Lsbar/2L, W0 = Ls/2L, and for a
twist g in O(L) the map G = g x g x g carries (Lambda, *_g) isomorphically
onto (Lambda_{g V0, g W0}, *), * = UNTWISTED triple product.  So closure and
span index of every twist g in O(L) are functions of the pair (gV0, gW0).

For each Lagrangian U we record: family (relative to V0), subalgebra,
left ideal (A.U in U), right ideal (U.A in U), two-sided ideal,
dim_F2 span(U.U), and for subalgebras h(U) := [U~ : Z-span(U~ . U~)].

For each of the 17,280 pairs we compute directly (576 basis products of a
Z-basis of Lambda_{V,W}, untwisted triple product):
  closure; if closed: span index [Lambda_{V,W} : S] (exact HNF), the F_2-rank
  of S in Lambda/2Lambda, whether 2Lambda in S;
and compare with
  criterion  C5:  closure <=> V and W are subalgebras of A;
  formula:        [Lambda_{V,W} : S] = h(W) * h(V)^2.
Outputs: s3_lagrangians.out, s3_lagrangians.json, s3_pairs.tsv
"""
import json
import time
from collections import Counter
from multiprocessing import Pool
import numpy as np
import leech_core as C

TAB = C.mul_table_mod2(C.T_STD)
s, sb = C.S_VEC, C.SBAR_VEC
H_sb = C.hnf([C.to_L_int(C.omul(C.L_BASIS[k], sb)) for k in range(8)])
H_s = C.hnf([C.to_L_int(C.omul(C.L_BASIS[k], s)) for k in range(8)])
V0 = C.span_bits([C.to_bits(v) for v in H_sb])
W0 = C.span_bits([C.to_bits(v) for v in H_s])
ALL = range(256)


# --- enumerate the Lagrangians (maximal totally singular subspaces) ---------
def lagrangians():
    Q = [v for v in range(1, 256) if C.q_form(v) == 0]
    B = {(u, v): C.b_form(u, v) for u in Q for v in Q}
    level = {frozenset({0})}
    counts = [1]
    for k in range(4):
        nxt = set()
        for S in level:
            basis = C.basis_bits(S)
            for v in Q:
                if v in S:
                    continue
                if all(B[(v, b)] == 0 for b in basis):
                    nxt.add(frozenset(S | {x ^ v for x in S}))
        level = nxt
        counts.append(len(level))
    return Q, counts, sorted(level, key=lambda S: sorted(S))


def props(U):
    Ul = list(U)
    sub = all(TAB[u, v] in U for u in Ul for v in Ul)
    left = all(TAB[a, u] in U for a in ALL for u in Ul)
    right = all(TAB[u, a] in U for a in ALL for u in Ul)
    UU = C.span_bits([int(TAB[u, v]) for u in Ul for v in Ul])
    fam = 'same' if C.dim_bits(U & V0) % 2 == 0 else 'other'
    return dict(subalgebra=sub, left_ideal=left, right_ideal=right,
                two_sided=left and right, dimUU=C.dim_bits(UU), family=fam,
                meetV0=C.dim_bits(U & V0), meetW0=C.dim_bits(U & W0))


def h_index(U):
    """[U~ : Z-span(U~ . U~)] (exact), and whether span(U~U~) is inside U~."""
    Ub = C.preimage_basis(U)
    prods = [[int(t) for t in C.lmul(a, b)] for a in Ub for b in Ub]
    H = C.hnf(prods)
    if len(H) < 8:
        return None, False, None
    inside = all(C.to_bits(p) in U for p in prods)
    idx_in_L = 1
    for i in range(8):
        idx_in_L *= H[i][i]
    return (idx_in_L, inside, H)


def pair_job(args):
    iV, iW, V, W = args
    U = np.array(C.leech_basis_from_pair(V, W), dtype=np.int64)
    P = C.triple_products(U, C.T_STD)
    member = C.in_pair_lattice_mod2(P, V, W)
    closed = bool(member.all())
    res = dict(iV=iV, iW=iW, closed=closed, n_fail=int((~member).sum()))
    if closed:
        Cc, ok = C.coords_in_basis(P, U)
        assert ok.all()
        idx, piv = C.span_index(Cc)
        # F_2-rank of the image of S in Lambda/2Lambda
        rank2 = rank_mod2(Cc)
        res.update(index=idx, log2=idx.bit_length() - 1, hnf_diagonal=piv, rank2=rank2,
                   twoLambda_in_S=(idx == 2 ** (24 - rank2)))
    return res


def rank_mod2(M):
    rows = [int("".join(str(int(t) & 1) for t in row), 2) for row in M]
    rank = 0
    rows = [r for r in rows if r]
    while rows:
        piv = max(rows)
        top = piv.bit_length() - 1
        rows = [r ^ piv if (r >> top) & 1 else r for r in rows if r != piv]
        rows = [r for r in rows if r]
        rank += 1
    return rank


if __name__ == '__main__':
    t0 = time.time()
    lines = []

    def log(*a):
        s_ = " ".join(str(t) for t in a)
        print(s_, flush=True)
        lines.append(s_)

    Q, counts, LAG = lagrangians()
    log("singular nonzero vectors of q:", len(Q))
    log("totally singular subspaces by dimension 0..4:", counts)
    log("Lagrangians (4-dim totally singular):", len(LAG))
    idxV0, idxW0 = LAG.index(V0), LAG.index(W0)
    log("V0 = Lsbar/2L is Lagrangian #%d, W0 = Ls/2L is Lagrangian #%d" % (idxV0, idxW0))

    P = [props(U) for U in LAG]
    for i, U in enumerate(LAG):
        hL, inside, H = h_index(U)
        P[i]['h_in_L'] = hL
        P[i]['UU_inside'] = inside
        P[i]['h'] = (hL // 16) if inside else None          # [U~ : span U~U~] = [L : span]/[L : U~]
        assert inside == P[i]['subalgebra']                   # integer-level vs mod-2 closure agree
    fam = Counter(p['family'] for p in P)
    log("families (relative to V0):", dict(fam))
    log("subalgebras:", sum(p['subalgebra'] for p in P),
        "| left ideals:", sum(p['left_ideal'] for p in P),
        "| right ideals:", sum(p['right_ideal'] for p in P),
        "| two-sided ideals:", sum(p['two_sided'] for p in P))
    log("subalgebras by family:", dict(Counter(p['family'] for p in P if p['subalgebra'])))
    tab = Counter((p['family'], p['subalgebra'], p['left_ideal'], p['right_ideal'], p['dimUU'], p['h'])
                  for p in P)
    log("\nclasses of Lagrangians: (family, subalgebra, left ideal, right ideal, dim span(UU), h(U)) : count")
    for k in sorted(tab, key=lambda k: (k[0], not k[1], str(k))):
        log("   ", k, ":", tab[k])
    log("\nV0 properties:", P[idxV0])
    log("W0 properties:", P[idxW0])

    # complementary pairs ----------------------------------------------------
    pairs = [(i, j) for i in range(len(LAG)) for j in range(len(LAG)) if not (LAG[i] & LAG[j]) - {0}]
    log("\nordered complementary pairs:", len(pairs))
    per = Counter(i for i, j in pairs)
    log("complements per Lagrangian:", set(per.values()))
    both = [(i, j) for (i, j) in pairs if P[i]['subalgebra'] and P[j]['subalgebra']]
    log("ordered complementary pairs with both members subalgebras:", len(both))
    fam_pairs = Counter(P[i]['family'] for i, j in pairs)
    log("pairs by family of V:", dict(fam_pairs), "(V and W always in the same family:",
        all(P[i]['family'] == P[j]['family'] for i, j in pairs), ")")

    jobs = [(i, j, LAG[i], LAG[j]) for (i, j) in pairs]
    with Pool(14) as pool:
        R = pool.map(pair_job, jobs, chunksize=64)
    log("direct closure computed for all %d pairs (%.0fs)" % (len(R), time.time() - t0))

    closed = [r for r in R if r['closed']]
    crit = [P[r['iV']]['subalgebra'] and P[r['iW']]['subalgebra'] for r in R]
    agree = sum(1 for r, c in zip(R, crit) if r['closed'] == c)
    log("closed pairs:", len(closed), "| criterion (V, W subalgebras) agrees with direct closure on",
        agree, "of", len(R), "pairs")
    log("non-closed pairs where exactly V fails / exactly W fails / both fail:",
        sum(1 for r in R if not r['closed'] and not P[r['iV']]['subalgebra'] and P[r['iW']]['subalgebra']),
        sum(1 for r in R if not r['closed'] and P[r['iV']]['subalgebra'] and not P[r['iW']]['subalgebra']),
        sum(1 for r in R if not r['closed'] and not P[r['iV']]['subalgebra'] and not P[r['iW']]['subalgebra']))

    form_ok = sum(1 for r in closed if r['index'] == P[r['iW']]['h'] * P[r['iV']]['h'] ** 2)
    log("formula [Lambda:S] = h(W) h(V)^2 holds on", form_ok, "of", len(closed), "closed pairs")
    log("2Lambda in S on", sum(r['twoLambda_in_S'] for r in closed), "of", len(closed), "closed pairs")
    log("index distribution over closed pairs (log2 index : count):",
        dict(sorted(Counter(r['log2'] for r in closed).items())))
    log("(2Lambda in S means Lambda/S = (Z/2)^k with 2^k the index, k = 24 - F_2-rank of S in Lambda/2Lambda)")
    log("index by family of (V,W):",
        {f: dict(sorted(Counter(r['log2'] for r in closed if P[r['iV']]['family'] == f).items()))
         for f in ('same', 'other')})
    log("(h(V), h(W)) combinations over closed pairs:",
        dict(sorted(Counter((P[r['iV']]['h'], P[r['iW']]['h']) for r in closed).items())))
    log("index by (V two-sided ideal?, W two-sided ideal?):",
        {k: dict(Counter(r['log2'] for r in closed if (P[r['iV']]['two_sided'], P[r['iW']]['two_sided']) == k))
         for k in [(a, b) for a in (True, False) for b in (True, False)]})
    log("index by (dim VV, dim WW):",
        {k: dict(Counter(r['log2'] for r in closed if (P[r['iV']]['dimUU'], P[r['iW']]['dimUU']) == k))
         for k in sorted(set((P[r['iV']]['dimUU'], P[r['iW']]['dimUU']) for r in closed))})
    w = next(r for r in R if (r['iV'], r['iW']) == (idxV0, idxW0))
    log("\nWilson's own pair (V0, W0): closed =", w['closed'], "failing products:", w['n_fail'], "of 576")

    with open('s3_lagrangians.json', 'w') as f:
        json.dump(dict(lagrangians=[sorted(U) for U in LAG], props=P, V0=idxV0, W0=idxW0,
                       pairs=[{k: v for k, v in r.items()} for r in R]), f)
    with open('s3_pairs.tsv', 'w') as f:
        f.write("iV\tiW\tclosed\tn_fail\tlog2_index\trank2\tV_sub\tW_sub\tV_2sided\tW_2sided\tV_family\th(V)\th(W)\n")
        for r in R:
            pv, pw = P[r['iV']], P[r['iW']]
            f.write("\t".join(map(str, [r['iV'], r['iW'], r['closed'], r['n_fail'], r.get('log2', ''),
                                         r.get('rank2', ''), pv['subalgebra'], pw['subalgebra'],
                                         pv['two_sided'], pw['two_sided'], pv['family'], pv['h'], pw['h']])) + "\n")
    log("total time %.0fs" % (time.time() - t0))
    with open('s3_lagrangians.out', 'w') as f:
        f.write("\n".join(lines) + "\n")
