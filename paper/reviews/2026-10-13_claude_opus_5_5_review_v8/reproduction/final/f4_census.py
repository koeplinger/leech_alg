"""
f4_census.py -- exhaustive censuses of twists, computed directly on Wilson's Lambda.

Conventions: leech_final.py.  Permutations are tuples (pi(1),...,pi(7)) with
g_pi(e_a) = e_{pi(a)};  convention A twists by g_pi (the paper's
x ._pi y = pi^{-1}(pi(x).pi(y))), convention B by g_pi^{-1}.

For every twist g the 576 basis-pair products are formed, closure is decided by
integrality of their Lambda-coordinates (cross-checked by the mod-2 Wilson
conditions), and for closing twists the index [Lambda : S_g] is computed from
the Z/4 structure of the coordinate matrix.  That is exact because
4 Lambda <= S_g for every twist with <L ._g L> = L (proved in the review:
(2a,0,0) *_g (2b,0,0) = (4 a._g b, 0, 0) with (2L)^3 <= Lambda), and the
premise <L ._g L> = L is asserted for every twist.  All closing twists of
convention A, all closing reflections and a random sample of convention B are
re-done with the exact integer HNF (no theorem used), and the span lattices S
themselves are compared (one lattice for all index-2^16 twists).  For every closing twist the index formula
[Lambda:S_g] = [Ls:<Ls._g Ls>] * [Lsbar:<Lsbar._g Lsbar>]^2 is checked.

Parts:
 1. all 5,040 permutations, conventions A and B: closure per cycle type (vs the
    paper's Section 5.6 table), index by parity, the 3-cycle orientation rule;
 2. all 120 root reflections (one root per +-pair): closure and index by root type;
 3. all 645,120 signed permutations e_a -> eps_a e_{pi(a)} (e0 fixed): distinct
    twisted tables, closure and index per table, and their relation to (1).
"""
import sys
import random
from collections import Counter, defaultdict
from multiprocessing import Pool
from itertools import permutations, product
from fractions import Fraction
import numpy as np
from leech_final import (Setup, perm_map, perm_inverse, reflection_map, gamma_of, twisted_constants,
                         triple_products, lambda_coords, mod2_ok, f2_rank, f2_span_basis, z4_structure,
                         hnf, hnf_index, h_index, structure_constants, T_STD, C_STD, int_det, in_L,
                         in_Ls_def, in_Lsbar_def)

ST = Setup()
PAPER = {(2,): (21, 21), (3,): (70, 35), (2, 2): (105, 42), (4,): (210, 0), (3, 2): (420, 126),
         (5,): (504, 252), (2, 2, 2): (105, 21), (3, 3): (280, 112), (4, 2): (630, 294),
         (6,): (840, 210), (3, 2, 2): (210, 105), (4, 3): (420, 168), (5, 2): (504, 42),
         (7,): (720, 336)}


def analyze_C(Cg, G=None, exact=False):
    """Closure and index for twisted constants Cg (L-coords).  G (L-matrix of g)
    is used for the 8-dim factors when given; otherwise they are computed in the
    twisted frame from Cg directly (same numbers, see f2_witnesses)."""
    iLL, _ = hnf_index(Cg.reshape(64, 8).tolist(), 8)
    assert iLL == 1                        # premise of 4 Lambda <= S
    P = triple_products(ST.B, Cg)
    coords, ok = lambda_coords(P, ST.N)
    assert np.array_equal(ok, mod2_ok(P, ST.cond))
    res = {"closes": bool(ok.all()), "nfail": int((~ok).sum())}
    # closure criterion data (twisted frame): Lsbar and Ls closed under ._g ?
    res["V_sub"] = closed_under(ST.V0, Cg)
    res["W_sub"] = closed_under(ST.W0, Cg)
    if not res["closes"]:
        return res
    r2 = f2_rank(coords)
    a, b = z4_structure(coords)
    idx = 2 ** (48 - 2 * a - b)
    res.update(index=idx, r2=r2, z4=(a, b), elem_abelian=(a + b == 24),
               two_lambda_in_S=(idx == 2 ** (24 - r2)))
    hV, _, _ = h_index(ST.V0, Cg)
    hW, _, _ = h_index(ST.W0, Cg)
    res.update(hV=hV, hW=hW, formula_ok=(hW * hV * hV == idx))
    if exact:
        ex, _ = hnf_index(coords.tolist(), 24)
        res["exact_index"] = ex
    if exact == "span":
        res["S_hnf"] = tuple(tuple(r) for r in hnf(coords.tolist(), 24))
    return res


def closed_under(Ub, Cg):
    """Is the full preimage U~ (U = rows mod 2) closed under the product Cg?
    Equivalent to: U is a subalgebra of L/2L for the induced product."""
    U = np.array(Ub, dtype=np.int64)
    prods = np.einsum('pa,qb,abc->pqc', U, U, Cg).reshape(-1, 8) % 2
    return f2_rank(np.vstack([U, prods])) == f2_rank(U)


def cycle_type(pi):
    seen, ct = set(), []
    for a in range(1, 8):
        if a in seen:
            continue
        L, b = 0, a
        while b not in seen:
            seen.add(b)
            b = pi[b - 1]
            L += 1
        if L > 1:
            ct.append(L)
    return tuple(sorted(ct, reverse=True))


def sign_of(pi):
    return (-1) ** sum(L - 1 for L in cycle_type(pi))


def work_perm(args):
    pi, conv, exact = args
    g = pi if conv == "A" else perm_inverse(pi)
    G = gamma_of(perm_map(g))
    r = analyze_C(twisted_constants(G), G, exact)
    r.update(pi=pi, conv=conv, ct=cycle_type(pi), sign=sign_of(pi))
    return r


def work_refl(R):
    G = gamma_of(reflection_map(R))
    r = analyze_C(twisted_constants(G), G, exact="span")
    r.update(root=R, det=int_det(G.tolist()))
    # 2*1_g with 1_g = s_r(e0) = e0 - <e0,r> r ; doubled coords of 2*1_g = 4*1_g
    rr = [Fraction(t, 2) for t in R]
    one = [Fraction(int(i == 0)) - rr[0] * rr[i] for i in range(8)]
    X = [int(4 * t) for t in one]
    r["two1_in"] = ("Lsbar" if in_Lsbar_def(X) else "Ls" if in_Ls_def(X) else "neither") if in_L(X) else "notL"
    return r


def work_table(args):
    key, Tg = args
    C = structure_constants(Tg)
    r = analyze_C(C)
    r["key"] = key
    return r


out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


if __name__ == "__main__":
    perms = list(permutations(range(1, 8)))
    rng = random.Random(7)
    sample = set(rng.sample(range(5040), 60))
    jobs = [(p, c, ("span" if c == "A" else (i in sample))) for c in ("A", "B") for i, p in enumerate(perms)]
    with Pool(14) as pool:
        res = pool.map(work_perm, jobs, chunksize=40)
    # ---------------------------------------------------------- 1. permutations
    say("1. ALL 5,040 PERMUTATION TWISTS")
    for conv in ("A", "B"):
        R = [r for r in res if r["conv"] == conv]
        say("  convention %s (%s):" % (conv, "g = g_pi, the paper's" if conv == "A" else "g = g_pi^{-1}"))
        tab = defaultdict(lambda: [0, 0])
        for r in R:
            if r["ct"]:
                tab[r["ct"]][0] += 1
                tab[r["ct"]][1] += r["closes"]
        allok = True
        for ct in sorted(PAPER, key=lambda c: (sum(c), c)):
            n, c = tab[ct]
            ok = (n, c) == PAPER[ct]
            allok &= ok
            say("    cycle type %-10s %4d permutations, %4d close   paper: %4d / %4d   %s"
                % (str(ct), n, c, PAPER[ct][1], PAPER[ct][0], "ok" if ok else "MISMATCH"))
        idn = [r for r in R if not r["ct"]][0]
        say("    identity closes:", idn["closes"], "(failing basis products: %d)" % idn["nfail"])
        cl = [r for r in R if r["closes"]]
        say("    total closing: %d;  matches paper table in all 14 rows: %s" % (len(cl), allok))
        byp = Counter((r["sign"], r["index"]) for r in cl)
        say("    (parity, index) counts:", dict(byp))
        say("    all closing: Lambda/S elementary abelian:", all(r["elem_abelian"] for r in cl),
            "; 2Lambda <= S:", all(r["two_lambda_in_S"] for r in cl),
            "; index formula hW*hV^2 holds:", all(r["formula_ok"] for r in cl))
        say("    (hW, hV) values by parity:", dict(Counter((r["sign"], r["hW"], r["hV"]) for r in cl)))
        crit = all(r["closes"] == (r["V_sub"] and r["W_sub"]) for r in R)
        say("    closure <=> (Lsbar and Ls both closed under ._g): holds for all 5040:", crit)
        ex = [r for r in R if "exact_index" in r and r["closes"]]
        say("    exact integer HNF re-check on %d closing twists (%s): all agree with the Z/4 index: %s"
            % (len(ex), "all closing twists" if conv == "A" else "random sample", all(r["exact_index"] == r["index"] for r in ex)))
        # 3-cycles: one orientation per 3-subset
        three = defaultdict(list)
        for r in R:
            if r["ct"] == (3,):
                support = tuple(sorted(a for a in range(1, 8) if r["pi"][a - 1] != a))
                three[support].append(r["closes"])
        say("    3-cycles: %d subsets, exactly one orientation closes for all: %s"
            % (len(three), all(sorted(v) == [False, True] for v in three.values())))
    A = {r["pi"]: r for r in res if r["conv"] == "A"}
    B = {r["pi"]: r for r in res if r["conv"] == "B"}
    clA = [r for r in res if r["conv"] == "A" and r["closes"]]
    S16 = set(r["S_hnf"] for r in clA if r["index"] == 2 ** 16)
    S10 = set(r["S_hnf"] for r in clA if r["index"] == 2 ** 10)
    say("  distinct span lattices S (exact HNF, convention A): index 2^16 twists: %d;  index 2^10 twists: %d"
        % (len(S16), len(S10)))
    w = (1, 2, 3, 4, 6, 7, 5)
    say("  witness pi(5)=6, pi(6)=7, pi(7)=5: convention A closes=%s index=%s; convention B closes=%s"
        % (A[w]["closes"], A[w].get("index"), B[w]["closes"]))
    same = sum(A[p]["closes"] == B[p]["closes"] for p in perms)
    say("  permutations with the same closure verdict in A and B: %d / 5040" % same)

    # ---------------------------------------------------------- 2. reflections
    say("")
    say("2. ALL 120 ROOT REFLECTIONS s_r (one root per +-pair: first nonzero coordinate > 0)")
    roots = []
    for X in product(range(-2, 3), repeat=8):
        if sum(t * t for t in X) == 8 and in_L(X):
            first = next(t for t in X if t != 0)
            if first > 0:
                roots.append(X)
    assert len(roots) == 120
    with Pool(14) as pool:
        rr = pool.map(work_refl, roots)

    def rtype(R):
        if all(t % 2 == 0 for t in R):
            nz = [i for i in range(8) if R[i]]
            i, j = nz
            sgn = "+" if R[i] * R[j] > 0 else "-"
            return ("e0%se_j" % sgn) if i == 0 else ("e_i%se_j" % sgn)
        neg = sum(1 for t in R[1:] if t < 0)
        return "half, %d minus among e1..e7" % neg
    groups = defaultdict(list)
    for r in rr:
        groups[rtype(r["root"])].append(r)
    for k in sorted(groups):
        g_ = groups[k]
        ips = Counter(sum(a * b for a, b in zip(r["root"], (-1, 1, 1, 1, 1, 1, 1, 1))) // 4 for r in g_)
        say("  %-32s %3d roots; <r,s> values %s; close: %3d; indices %s; 2*1_g in %s"
            % (k, len(g_), dict(ips), sum(r["closes"] for r in g_),
               dict(Counter(r.get("index") for r in g_ if r["closes"])),
               dict(Counter(r["two1_in"] for r in g_))))
    cl = [r for r in rr if r["closes"]]
    say("  totals: %d do not close, %d close with index 2^16, %d close with index 1"
        % (sum(not r["closes"] for r in rr), sum(r.get("index") == 2 ** 16 for r in cl),
           sum(r.get("index") == 1 for r in cl)))
    say("  exact HNF agrees with the Z/4 index for all 56 closing reflections:",
        all(r["exact_index"] == r["index"] for r in cl))
    say("  index formula holds for all closing reflections:", all(r["formula_ok"] for r in cl),
        "; closure <=> criterion for all 120:", all(r["closes"] == (r["V_sub"] and r["W_sub"]) for r in rr))
    wit = [r for r in rr if r["root"] == (1, 1, 1, 1, 1, 1, 1, -1)][0]
    say("  witness r = (e0+...+e6-e7)/2: closes=%s index=%s hW=%s hV=%s"
        % (wit["closes"], wit["index"], wit["hW"], wit["hV"]))
    say("  sigma's root e1-e2: index", [r for r in rr if r["root"] == (0, 2, -2, 0, 0, 0, 0, 0)][0]["index"])
    S16r = set(r["S_hnf"] for r in cl if r["index"] == 2 ** 16)
    say("  the 28 index-2^16 reflections and the 588 odd closing permutations share ONE span lattice:",
        len(S16r | S16) == 1)

    # ---------------------------------------------------------- 3. signed permutations
    say("")
    say("3. ALL 645,120 SIGNED PERMUTATIONS g(e0)=e0, g(e_a) = eps_a e_{pi(a)}")
    tab_t = np.zeros((8, 8), dtype=np.int64)
    tab_c = np.zeros((8, 8), dtype=np.int64)
    for a in range(8):
        for b in range(8):
            c = int(np.nonzero(T_STD[a, b])[0][0])
            tab_c[a, b] = c
            tab_t[a, b] = T_STD[a, b, c]
    eps = np.array([[1] + [1 - 2 * ((m >> k) & 1) for k in range(7)] for m in range(128)], dtype=np.int64)
    nneg = np.array([bin(m).count("1") for m in range(128)])
    classes = {}
    for pi in perms:
        full = (0,) + pi                       # full[a] = pi(a), pi(0) = 0
        inv = [0] * 8
        for a in range(8):
            inv[full[a]] = a
        ca = tab_c[np.ix_(full, full)]         # e_{pi a} e_{pi b} = t e_c
        ta = tab_t[np.ix_(full, full)]
        d = np.vectorize(lambda c: inv[c])(ca)  # g^{-1}(e_c) = eps_d e_d with d = pi^{-1}(c)
        sgn = eps[:, :, None] * eps[:, None, :] * ta[None] * eps[:, d]
        keys = (sgn * (d[None] + 1)).reshape(128, 64).astype(np.int8)
        for m in range(128):
            k = keys[m].tobytes()
            if k not in classes:
                classes[k] = {"n": 0, "parities": set(), "perm_rep": None, "rep": (pi, m)}
            cl_ = classes[k]
            cl_["n"] += 1
            cl_["parities"].add(int(nneg[m]) % 2)
            if m == 0:
                cl_["perm_rep"] = pi
    say("  distinct twisted tables: %d; maps per table: %s" % (len(classes), dict(Counter(c["n"] for c in classes.values()))))
    say("  sign-count parity constant on each table:", all(len(c["parities"]) == 1 for c in classes.values()))
    even = [k for k, c in classes.items() if 0 in c["parities"]]
    odd = [k for k, c in classes.items() if 1 in c["parities"]]
    say("  tables from even sign count (these g preserve L): %d;  odd (g(L) != L): %d" % (len(even), len(odd)))
    say("  even-sign tables containing an unsigned permutation: %d" % sum(classes[k]["perm_rep"] is not None for k in even))

    def key_to_T(k):
        arr = np.frombuffer(k, dtype=np.int8).astype(np.int64).reshape(8, 8)
        T = np.zeros((8, 8, 8), dtype=np.int64)
        for a in range(8):
            for b in range(8):
                v = int(arr[a, b])
                T[a, b, abs(v) - 1] = 1 if v > 0 else -1
        return T

    def opp_key(k):
        arr = np.frombuffer(k, dtype=np.int8).reshape(8, 8)
        return np.ascontiguousarray(arr.T).tobytes()
    say("  every odd-sign table is the OPPOSITE table (T'[a,b] = T[b,a]) of an even-sign table:",
        sorted(opp_key(k) for k in odd) == sorted(even))
    with Pool(14) as pool:
        tr = pool.map(work_table, [(k, key_to_T(k)) for k in classes])
    trd = {r["key"]: r for r in tr}
    for lab, ks in (("even (L-preserving)", even), ("odd (opposite)", odd)):
        c = [trd[k] for k in ks if trd[k]["closes"]]
        say("  %-20s tables: %d close; indices %s" % (lab, len(c), dict(Counter(r["index"] for r in c))))
    allc = [r for r in tr if r["closes"]]
    say("  total closing tables: %d; indices %s" % (len(allc), dict(Counter(r["index"] for r in allc))))
    say("  opposite tables have identical closure and index:",
        all(trd[k]["closes"] == trd[opp_key(k)]["closes"] and trd[k].get("index") == trd[opp_key(k)].get("index")
            for k in odd))
    permA = {r["pi"]: r for r in res if r["conv"] == "A"}
    agree = all(trd[k]["closes"] == permA[classes[k]["perm_rep"]]["closes"] and
                trd[k].get("index") == permA[classes[k]["perm_rep"]].get("index") for k in even)
    say("  even-sign tables reproduce the permutation census (closure and index):", agree)
    say("  => C4's 168 closing tables = 84 permutation tables (1764/21) + their 84 opposites.")
    open("f4_census.out.txt", "w").write("\n".join(out) + "\n")
