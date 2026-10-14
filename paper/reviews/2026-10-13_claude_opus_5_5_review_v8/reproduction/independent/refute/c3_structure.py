#!/usr/bin/env python3
"""
c3_structure.py -- structural explanation of the span indices 2^16, 2^10 and 1,
and of which root reflections give which.

Run:   python3 c3_structure.py > c3_structure.out
Needs: leech_twist_lib.py.  Exact integer arithmetic throughout.
Conventions: main.tex v8, see the header of leech_twist_lib.py.

What is checked (the proofs are in the report; this script verifies every
ingredient on the twists in question):
  (L1) S_g contains 4*Lambda for EVERY twist g with g(L) = L that closes:
       eps = g^{-1}(e0) is the identity of ._g, 2*eps = g^{-1}(2e0) lies in L,
       so (4eps,0,0) lies in 2L^3, inside Lambda, and (4eps,0,0)*v = 4*(x',z',y'),
       4 times the block swap, which preserves Lambda.  Hence [Lambda:S_g] is a
       power of 2.
  (L2) 2-adic splitting along D = {(a,a,a)}, T = {sum 0}:
          [Lambda:S_g] = [B : B.B]_2 * [A : A.A]_2^2,
       A = g(L sbar), B = g(L s), X.X = Z-span of products (standard product).
  (L3) for sigma: S = Lambda_D + 2*Lambda, Lambda_D = {(a,a,a): a in L s}
       (so S/2Lambda is the image of the "diagonal" Lambda_D; the paper calls
       the identification of S/2Lambda open).
  (L4) mod-2 picture: V = A/2L, W = B/2L in L/2L; dims of V.V and W.W;
       ideal properties of A, B in L.
  (L5) which root reflections close / give which index, in terms of <r,s>.
"""
import sys
from fractions import Fraction as Fr
from collections import Counter, defaultdict

import leech_twist_lib as T


def say(*a):
    print(*a)
    sys.stdout.flush()


def v2(n):
    k = 0
    while n and n % 2 == 0:
        n //= 2
        k += 1
    return k


def span_index_in(rows, lat):
    return T.sub_index(rows, lat)


def square_index(rows):
    X = T.Lattice(rows)
    prods = [T.omul2(a, b) for a in rows for b in rows]
    closed = all(p in X for p in prods)
    if not closed:
        return None
    return T.sub_index(prods, X)


def ideal_sides(rows):
    """Is the sublattice X (rows) a left / right ideal of L (standard product)?"""
    X = T.Lattice(rows)
    left = all(T.omul2(l, a) in X for l in T.L_BASIS2 for a in rows)
    right = all(T.omul2(a, l) in X for l in T.L_BASIS2 for a in rows)
    return left, right


def mod2_square_dim(rows):
    """dim over F2 of span(V.V) inside L/2L, V = X/2L."""
    prods = [T.omul2(a, b) for a in rows for b in rows]
    C = [T.L_LAT.coords(p) for p in prods]
    return T.rank_mod_p(C, 2)


SIG = T.Twist(T.perm_matrix({1: 2, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}), "sigma=(1 2)")
R = T.Twist(T.reflection([1, 1, 1, 1, 1, 1, 1, -1]), "s_r, r=(e0+..+e6-e7)/2")
CYC = T.Twist(T.perm_matrix({1: 1, 2: 2, 3: 3, 4: 4, 5: 6, 6: 7, 7: 5}), "3-cycle (5 6 7)")
TWISTS = [SIG, R, CYC]

# ---------------------------------------------------------------- (L1)
say("=" * 78)
say("(L1) S_g contains 4*Lambda; [Lambda:S_g] is a power of 2")
say("=" * 78)
results = {}
for tw in TWISTS:
    o = T.closure_and_index(tw)
    results[tw.name] = o
    eps = tw.identity_element()
    two_eps = [int(4 * x) for x in eps]          # doubled coordinates of 2*eps
    four_eps_vec = [2 * t for t in two_eps] + [0] * 16    # doubled coords of (4eps,0,0)
    in_lam = T.in_Lambda2(four_eps_vec)
    # (4eps,0,0) * b = 4 * swap23(b) for every basis vector b
    swap_ok = True
    for b in T.LAMBDA_BASIS2:
        p = T.tprod2(tw, four_eps_vec, b)
        q = T.tprod2(tw, b, four_eps_vec)
        sw = [4 * t for t in b[0:8] + b[16:24] + b[8:16]]
        swap_ok &= (p == sw and q == sw)
    # 4*b_k in S and 2*b_k in S (HNF row reduction)
    H = o["hnf"]
    in_S = lambda vec: T.hnf_rows(H + [vec]) == H
    four_in = all(in_S([4 * int(i == j) for j in range(24)]) for i in range(24))
    two_in = all(in_S([2 * int(i == j) for j in range(24)]) for i in range(24))
    say(f"{tw.name:28s} eps={['%s' % x for x in eps]}")
    say(f"{'':28s} 2eps in L: {T.in_L2(two_eps)} | (4eps,0,0) in Lambda: {in_lam} | "
        f"(4eps,0,0)*v = v*(4eps,0,0) = 4*swap23(v) on the basis: {swap_ok}")
    say(f"{'':28s} 4*Lambda in S: {four_in} | 2*Lambda in S: {two_in} | [Lambda:S] = {o['index']}")

# ---------------------------------------------------------------- (L2)
say("\n" + "=" * 78)
say("(L2) 2-adic factorisation [Lambda:S] = [B:B.B] [A:A.A]^2")
say("=" * 78)
# the decomposition itself: Lambda_D = {(a,a,a): a in L s}, Lambda_T = {(x,y,-x-y): x,y in L sbar}
LD = [b + b + b for b in T.LS_LAT.rows]
LT = [x + [0] * 8 + [-t for t in x] for x in T.LSBAR_LAT.rows] + [[0] * 8 + y + [-t for t in y] for y in T.LSBAR_LAT.rows]
say("Lambda_D, Lambda_T inside Lambda:", all(T.in_Lambda2(v) for v in LD + LT))
idxDT = T.sub_index(LD + LT, T.LAMBDA_LAT)
say("[Lambda : Lambda_D + Lambda_T] =", idxDT, "= 3^%d" % (len(str(idxDT)) and next(k for k in range(40) if 3 ** k == idxDT)))
say("Lambda cap D = Lambda_D and Lambda cap T = Lambda_T: (x,x,x) in Lambda iff x in L s; "
    "(x,y,-x-y) in Lambda iff x,y in L sbar  [checked on the 2L-cosets below]")
# exhaustive check over representatives mod 2L (2L^3 lies in all lattices involved)
reps2L = []
for mask in range(256):
    v = [0] * 8
    for k in range(8):
        if (mask >> k) & 1:
            v = [a + b for a, b in zip(v, T.L_BASIS2[k])]
    reps2L.append(v)
okD = all(T.in_Lambda2(x + x + x) == (x in T.LS_LAT) for x in reps2L)
okT = all(T.in_Lambda2(x + y + [-a - b for a, b in zip(x, y)]) == (x in T.LSBAR_LAT and y in T.LSBAR_LAT)
          for x in reps2L for y in reps2L)
say("   verified over all classes mod 2L:", okD, okT)
for tw in TWISTS:
    A_rows = [tw.apply(b) for b in T.LSBAR_LAT.rows]
    B_rows = [tw.apply(b) for b in T.LS_LAT.rows]
    iA, iB = square_index(A_rows), square_index(B_rows)
    idx = results[tw.name]["index"]
    say(f"{tw.name:28s} [B:B.B] = {iB}, [A:A.A] = {iA}  ->  [B:B.B]*[A:A.A]^2 = {iB * iA * iA} ; [Lambda:S] = {idx}"
        f"  {'MATCH' if iB * iA * iA == idx else 'MISMATCH'}")

# ---------------------------------------------------------------- (L3)
say("\n" + "=" * 78)
say("(L3) sigma: S = Lambda_D + 2*Lambda ; s_r: S = Lambda")
say("=" * 78)
H = results[SIG.name]["hnf"]
gens = [T.LAMBDA_LAT.coords(v) for v in LD] + [[2 * int(i == j) for j in range(24)] for i in range(24)]
say("sigma: HNF(S) == HNF(Lambda_D + 2 Lambda):", T.hnf_rows(gens) == H)
Hc = results[CYC.name]["hnf"]
A_c = [CYC.apply(b) for b in T.LSBAR_LAT.rows]
AA_c = T.hnf_rows([T.L_LAT.coords(T.omul2(a, b)) for a in A_c for b in A_c])
# S for the 3-cycle: Lambda_D + {(x,y,-x-y): x,y in g^{-1}(A.A)} + 2 Lambda
AAc_rows = [CYC.apply_inv([sum(c[k] * T.L_BASIS2[k][i] for k in range(8)) for i in range(8)]) for c in AA_c]
gens_c = ([T.LAMBDA_LAT.coords(v) for v in LD]
          + [T.LAMBDA_LAT.coords(x + [0] * 8 + [-t for t in x]) for x in AAc_rows]
          + [T.LAMBDA_LAT.coords([0] * 8 + y + [-t for t in y]) for y in AAc_rows]
          + [[2 * int(i == j) for j in range(24)] for i in range(24)])
say("3-cycle: HNF(S) == HNF(Lambda_D + {(x,y,-x-y): x,y in pi^{-1}(A.A)} + 2 Lambda):", T.hnf_rows(gens_c) == Hc)

# ---------------------------------------------------------------- (L4)
say("\n" + "=" * 78)
say("(L4) mod-2 picture and ideal properties (standard product, A = g(L sbar), B = g(L s))")
say("=" * 78)
for tw in TWISTS:
    A_rows = [tw.apply(b) for b in T.LSBAR_LAT.rows]
    B_rows = [tw.apply(b) for b in T.LS_LAT.rows]
    lA, rA = ideal_sides(A_rows)
    lB, rB = ideal_sides(B_rows)
    AA = [T.omul2(a, b) for a in A_rows for b in A_rows]
    A_lat = T.Lattice(A_rows)
    A2A = all(all(t % 2 == 0 for t in A_lat.coords(p)) for p in AA)
    say(f"{tw.name:28s} A: left ideal {lA}, right ideal {rA}; dim span(V.V) = {mod2_square_dim(A_rows)} (dim V = 4); "
        f"A.A inside 2A: {A2A}")
    say(f"{'':28s} B: left ideal {lB}, right ideal {rB}; dim span(W.W) = {mod2_square_dim(B_rows)} (dim W = 4)")

# ---------------------------------------------------------------- (L5)
say("\n" + "=" * 78)
say("(L5) the 120 root reflections, classified by <r,e0>, <r,s>, <r,sbar> (r up to sign)")
say("=" * 78)
roots = []
for i in range(8):
    for j in range(i + 1, 8):
        for si in (2, -2):
            for sj in (2, -2):
                v = [0] * 8
                v[i], v[j] = si, sj
                roots.append(v)
for m in range(256):
    v = [1 if (m >> k) & 1 else -1 for k in range(8)]
    if T.in_L2(v):
        roots.append(v)
reps, seen = [], set()
for r in roots:
    if tuple(-t for t in r) in seen:
        continue
    seen.add(tuple(r))
    reps.append(r)
ip = lambda X, Y: Fr(sum(a * b for a, b in zip(X, Y)), 4)
table = defaultdict(Counter)
rl_identity = Counter()
for r in reps:
    tw = T.Twist(T.reflection(r))
    o = T.closure_and_index(tw)
    outcome = "no closure" if not o["closes"] else f"index 2^{v2(o['index'])}"
    A_rows = [tw.apply(b) for b in T.LSBAR_LAT.rows]
    B_rows = [tw.apply(b) for b in T.LS_LAT.rows]
    iA, iB = square_index(A_rows), square_index(B_rows)
    if o["closes"]:
        assert iA is not None and iB is not None and iB * iA * iA == o["index"], "factorisation fails"
    kind = "half-int" if all(t % 2 for t in r) else ("e0+-ej" if r[0] else "ei+-ej")
    key = (kind, abs(ip(r, [2, 0, 0, 0, 0, 0, 0, 0])), abs(ip(r, T.S2)), abs(ip(r, T.SBAR2)))
    table[key][(outcome, f"[A:A.A]={iA}", f"[B:B.B]={iB}")] += 1
    # identity s_r(L sbar) = r.L  (r the LEFT factor), for every root r
    rL = T.hnf_rows([T.L_LAT.coords(T.omul2(r, b)) for b in T.L_BASIS2])
    Ah = T.hnf_rows([T.L_LAT.coords(a) for a in A_rows])
    rl_identity[Ah == rL] += 1
say(f"{'type':9s} {'|<r,e0>|':>9s} {'|<r,s>|':>8s} {'|<r,sbar>|':>11s}   outcome counts (one per +-r pair)")
for key in sorted(table, key=lambda k: (k[0], k[1], k[2], k[3])):
    say(f"{key[0]:9s} {str(key[1]):>9s} {str(key[2]):>8s} {str(key[3]):>11s}   {dict(table[key])}")
say("(the factorisation [Lambda:S] = [B:B.B][A:A.A]^2 was asserted for every closing reflection;")
say(" [X:X.X] = None means X.X is not contained in X, i.e. X is not closed)")
say("identity s_r(L sbar) = r.L (as lattices, r the left factor) holds for", rl_identity[True],
    "of the 120 reflections, fails for", rl_identity[False])
