"""
f5_polarizations.py -- the mod-2 picture: Lagrangians of L/2L, the closure
criterion on ALL complementary pairs, the index census over all twists in
W(E8) = Aut(L), and the rule that decides the index.

Conventions: leech_final.py.  Abar = L/2L = F2^8 in L-coordinates mod 2
(elements are 8-bit masks).  Product mu = the octonion product of L reduced
mod 2L (structure constants = Table LL mod 2).  q(x) = N(x)/2 mod 2.
For an ordered pair (V, W) of complementary Lagrangians,
    Lambda(V,W) = {(x,y,z) in L^3 : x+y, x+z, y+z in V~ ; x+y+z in W~}
(V~ = full preimage of V in L).  Wilson's lattice is Lambda(Lsbar/2L, Ls/2L),
and for a twist g in Aut(L), u -> g(u) blockwise is a ring isomorphism
(Lambda, *_g) -> (Lambda(gV0, gW0), *), so closure and index of *_g are those
of the UNTWISTED product on Lambda(V_g, W_g), V_g = g(Lsbar)/2L, W_g = g(Ls)/2L.

Computed here (all exhaustive, exact):
 1. 135 nonzero singular vectors, 270 Lagrangians (two families of 135),
    subalgebras, left/right/two-sided ideals.
 2. all 17,280 ordered complementary pairs: Lambda(V,W) built from scratch,
    576 untwisted products, closure, index (Z/4 method), compared with
    (i) the criterion "V and W both subalgebras", (ii) hW * hV^2.
 3. h(U) = [U~ : <U~ U~>] on the subalgebras; the rule
    h(U) = 1 if eps = [2 e0] not in U, else 2^5 or 2^8 by family; and the
    mod-2 invariant 2^(d+m) used by another agent, compared with h.
 4. W(E8) acts transitively on the 17,280 ordered pairs (BFS from Wilson's
    pair under the 120 reflections, acting on L/2L as transvections), and
    det g = +1 iff V_g lies in the family of Lsbar/2L (all perms and reflections).
"""
from collections import Counter, defaultdict, deque
from itertools import permutations, product
from multiprocessing import Pool
import numpy as np
from leech_final import (Setup, C_STD, from_Lcoords, Lcoords, norm2, f2_rref, f2_rank, annihilator,
                         f2_nullspace, triple_products, z4_structure, h_index, gamma_of, perm_map,
                         reflection_map, int_det, in_L)

ST = Setup()
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


def vec(m):
    return [(m >> k) & 1 for k in range(8)]


def mask(v):
    return sum(int(v[k]) % 2 << k for k in range(8))


# ------------------------------------------------------------------ mod-2 algebra
PM = [[mask(C_STD[i, j]) for j in range(8)] for i in range(8)]
MU = np.zeros((256, 256), dtype=np.int64)
for a in range(256):
    for b in range(256):
        r = 0
        for i in range(8):
            if (a >> i) & 1:
                for j in range(8):
                    if (b >> j) & 1:
                        r ^= PM[i][j]
        MU[a, b] = r
Q = [int(norm2(from_Lcoords(vec(m))) / 2) % 2 for m in range(256)]
B = lambda x, y: Q[x ^ y] ^ Q[x] ^ Q[y]
EPS = mask(Lcoords([4, 0, 0, 0, 0, 0, 0, 0]))           # class of 2 e0


def span(gens):
    """F2-span of a collection of 8-bit masks."""
    S = {0}
    for g in set(int(x) for x in gens):
        if g not in S:
            S |= {s ^ g for s in S}
    return frozenset(S)


V0 = span([mask(r) for r in ST.V0])
W0 = span([mask(r) for r in ST.W0])
singular = [m for m in range(1, 256) if Q[m] == 0]
say("1. L/2L: nonzero singular vectors:", len(singular), "; eps = [2e0] singular:", Q[EPS] == 0,
    "; eps in Lsbar/2L:", EPS in V0, "; eps in Ls/2L:", EPS in W0)
level = {span([v]) for v in singular}
sizes = [len(level)]
for _ in range(3):
    nxt = set()
    for U in level:
        for v in singular:
            if v not in U and all(B(v, u) == 0 for u in U):
                nxt.add(U | frozenset(u ^ v for u in U))
    level = nxt
    sizes.append(len(level))
LAG = sorted(level, key=lambda U: sorted(U))
say("   totally singular subspaces of dims 1..4:", sizes)
assert len(LAG) == 270 and V0 in level and W0 in level
LID = {U: i for i, U in enumerate(LAG)}
fam = [len(U & V0).bit_length() - 1 for U in LAG]        # log2 |U cap V0| = dim(U cap V0)
family = [d % 2 for d in fam]                            # 0 = family of Lsbar/2L
say("   Lagrangians: %d; family sizes %s (family 0 = that of Lsbar/2L; Ls/2L in family %d)"
    % (len(LAG), dict(Counter(family)), family[LID[W0]]))
ALL = range(256)


def is_sub(U):
    return all(MU[a, b] in U for a in U for b in U)


def is_left(U):
    return all(MU[x, a] in U for x in ALL for a in U)


def is_right(U):
    return all(MU[a, x] in U for x in ALL for a in U)


sub = [is_sub(U) for U in LAG]
left = [is_left(U) for U in LAG]
right = [is_right(U) for U in LAG]
say("   subalgebras: %d (by family %s); left ideals: %d; right ideals: %d; two-sided: %d"
    % (sum(sub), dict(Counter(family[i] for i in range(270) if sub[i])), sum(left), sum(right),
       sum(l and r for l, r in zip(left, right))))
say("   Lsbar/2L: subalgebra %s, left ideal %s, right ideal %s, square-zero %s; Ls/2L: subalgebra %s"
    % (sub[LID[V0]], left[LID[V0]], right[LID[V0]], all(MU[a, b] == 0 for a in V0 for b in V0), sub[LID[W0]]))

# ------------------------------------------------------------------ h-values and invariants
def basis_rows(U):
    R, piv = f2_rref(np.array([vec(m) for m in U], dtype=np.int64))
    return R


H = {}
for i, U in enumerate(LAG):
    if sub[i]:
        h, _, _ = h_index(basis_rows(U), C_STD)
        H[i] = h


def dmspan(U):
    d = 4 - (len(span([MU[a, b] for a in U for b in U])).bit_length() - 1)
    ext = span(list(U) + [MU[x, a] for x in ALL for a in U] + [MU[a, x] for x in ALL for a in U])
    m = 8 - (len(ext).bit_length() - 1)
    return d, m


rule_ok, dm_vs_h = True, Counter()
for i in H:
    U = LAG[i]
    pred = 1 if EPS not in U else (2 ** 5 if family[i] == 0 else 2 ** 8)
    rule_ok &= (pred == H[i])
    d, m = dmspan(U)
    dm_vs_h[(EPS in U, family[i], H[i], 2 ** (d + m))] += 1
say("3. h(U) = [U~ : <U~.U~>] over the %d subalgebra Lagrangians: %s" % (len(H), dict(Counter(H.values()))))
# eps lemma: eps in U  =>  U~ integral, Re(U~ . U~) even, while 2s in U~ has Re = -1, so h(U) >= 2
lem = Counter()
for i, U in enumerate(LAG):
    if EPS in U:
        Ub = basis_rows(U)
        R_, piv_ = f2_rref(np.array(Ub, dtype=np.int64))
        gens = [list(map(int, r)) for r in R_] + [[2 * int(j == c) for c in range(8)] for j in range(8) if j not in piv_]
        E = [from_Lcoords(g) for g in gens]                       # doubled e-coordinates
        integral = all(t % 2 == 0 for X in E for t in X)
        from leech_final import dmul as _dm
        re_even = all(_dm(X, Y)[0] % 4 == 0 for X in E for Y in E)  # Re(xy) even  <=> doubled Re = 0 mod 4
        lem[(integral, re_even)] += 1
say("   eps lemma on the %d Lagrangians containing eps: (U~ integral, Re(U~.U~) even) -> %s" % (sum(lem.values()), dict(lem)))
say("   rule  h = 1 if eps not in U;  2^5 if eps in U, family 0;  2^8 if eps in U, family 1:  holds for all:", rule_ok)
say("   (eps in U, family, h, 2^(d+m)) counts:", dict(dm_vs_h))
say("   h(Lsbar/2L) = %d, 2^(d+m)(Lsbar/2L) = %d" % (H[LID[V0]], 2 ** sum(dmspan(V0))))


# ------------------------------------------------------------------ pairs: direct computation
def wilson_basis(Vb, Wb):
    AV, AW = annihilator(Vb), annihilator(Wb)
    I8, Z8 = np.eye(8, dtype=np.int64), np.zeros((8, 8), dtype=np.int64)
    cond = np.vstack([np.hstack([AV @ P, AV @ Qm, AV @ R]) for P, Qm, R in ((I8, I8, Z8), (I8, Z8, I8), (Z8, I8, I8))]
                     + [np.hstack([AW, AW, AW])]) % 2
    K = f2_nullspace(cond)
    R, piv = f2_rref(K)
    nonpiv = [j for j in range(24) if j not in piv]
    rows = [list(map(int, r)) for r in R] + [[2 * int(j == c) for c in range(24)] for j in nonpiv]
    return np.array(rows, dtype=np.int64), R.astype(np.int64), piv, nonpiv


def coords_fast(P, R, piv, nonpiv):
    """Exact coordinates in the basis (lifted RREF rows, 2e_j for non-pivot j)."""
    C1 = P[:, piv]
    Rem = P - C1 @ R
    assert np.all(Rem[:, piv] == 0)
    ok = np.all(Rem[:, nonpiv] % 2 == 0, axis=1)
    return np.hstack([C1, Rem[:, nonpiv] // 2]), ok


def work_pair(ij):
    i, j = ij
    Bm, R, piv, nonpiv = wilson_basis(basis_rows(LAG[i]), basis_rows(LAG[j]))
    assert abs(int_det(Bm.tolist())) == 2 ** 12
    P = triple_products(Bm, C_STD)
    c, ok = coords_fast(P, R, piv, nonpiv)
    if not ok.all():
        return (i, j, False, int((~ok).sum()), None, None)
    a, b = z4_structure(c)
    return (i, j, True, 0, 2 ** (48 - 2 * a - b), (a + b == 24) and (2 ** (48 - 2 * a - b) == 2 ** (24 - f2_rank(c))))


pairs = [(i, j) for i in range(270) for j in range(270) if len(LAG[i] & LAG[j]) == 1]
say("2. ordered complementary pairs:", len(pairs), "; complements per Lagrangian:",
    dict(Counter(Counter(i for i, j in pairs).values())), "; same family always:",
    all(family[i] == family[j] for i, j in pairs))
both = [(i, j) for i, j in pairs if sub[i] and sub[j]]
say("   pairs with both members subalgebras:", len(both))
if __name__ == "__main__":
    with Pool(14) as pool:
        res = pool.map(work_pair, pairs, chunksize=50)
    closing = [r for r in res if r[2]]
    say("   direct computation on all %d pairs: %d close" % (len(res), len(closing)))
    say("   closure <=> both subalgebras, on every pair:", all(r[2] == (sub[r[0]] and sub[r[1]]) for r in res))
    say("   index census over the closing pairs:", dict(sorted(Counter(r[4] for r in closing).items())))
    say("   Lambda/S elementary abelian and 2Lambda <= S on every closing pair:", all(r[5] for r in closing))
    say("   index = h(W) * h(V)^2 on every closing pair:", all(r[4] == H[r[1]] * H[r[0]] ** 2 for r in closing))
    nV0 = sum(1 for i, j in both if i == LID[V0] or j == LID[V0])
    say("   closing pairs containing Lsbar/2L (as V or W):", nV0)
    wil = [r for r in res if r[0] == LID[V0] and r[1] == LID[W0]][0]
    say("   Wilson's own pair (Lsbar/2L, Ls/2L), untwisted: closes =", wil[2], "(failing basis products: %d)" % wil[3])

    # ------------------------------------------------------------------ W(E8) action
    roots = []
    for X in product(range(-2, 3), repeat=8):
        if sum(t * t for t in X) == 8 and in_L(X) and next(t for t in X if t) > 0:
            roots.append(mask(Lcoords(list(X))))
    assert len(set(roots)) == 120 and all(Q[r] == 1 for r in roots)
    trans = []
    for r in roots:
        img = [m ^ (r if B(m, r) else 0) for m in range(256)]
        trans.append([LID[frozenset(img[u] for u in U)] for U in LAG])
    start = (LID[V0], LID[W0])
    seen, dq = {start}, deque([start])
    while dq:
        i, j = dq.popleft()
        for t in trans:
            nb = (t[i], t[j])
            if nb not in seen:
                seen.add(nb)
                dq.append(nb)
    say("4. orbit of Wilson's pair under the 120 reflections (W(E8) acting on L/2L):", len(seen),
        "= all ordered complementary pairs:", seen == set(pairs))
    WE8 = 696729600
    say("   |W(E8)| / 17280 = %d elements of W(E8) realize each pair" % (WE8 // 17280))
    cnt = Counter(r[4] for r in closing)
    say("   twists g in W(E8) by index: " + ", ".join("%d: %d" % (k, v * (WE8 // 17280)) for k, v in sorted(cnt.items()))
        + "; closing: %d of %d" % (len(closing) * (WE8 // 17280), WE8))
    # det g versus family of V_g, for all permutations and reflections
    ok = True
    V0rows = np.array(ST.V0, dtype=np.int64)
    for pi in permutations(range(1, 8)):
        G = gamma_of(perm_map(pi))
        Vg = span([mask(r) for r in (V0rows @ G) % 2])
        ok &= ((family[LID[Vg]] == 0) == (int_det(G.tolist()) == 1))
    for X in product(range(-2, 3), repeat=8):
        if sum(t * t for t in X) == 8 and in_L(X) and next(t for t in X if t) > 0:
            G = gamma_of(reflection_map(X))
            Vg = span([mask(r) for r in (V0rows @ G) % 2])
            ok &= ((family[LID[Vg]] == 0) == (int_det(G.tolist()) == 1))
    say("   det g = +1  <=>  V_g in the family of Lsbar/2L: holds for all 5040 permutations and 120 reflections:", ok)
    open("f5_polarizations.out.txt", "w").write("\n".join(out) + "\n")
