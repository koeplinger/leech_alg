#!/usr/bin/env python3
"""
c3_refute_analysis.py -- adversarial, exact re-examination of claim C3
("the triple product twisted by the reflection in r = (e0+e1+...+e6-e7)/2
closes on Wilson's Leech lattice with span index 1").

Run:   python3 c3_refute_analysis.py > c3_refute_analysis.out
Needs: leech_twist_lib.py (same folder).  Pure Python ints/Fractions.

All conventions are those of main.tex v8 and are listed in the header of
leech_twist_lib.py; the essential ones:
  * e_a e_b = +e_c for the Fano triples (1,2,4),(2,3,5),(3,4,6),(4,5,7),
    (5,6,1),(6,7,2),(7,1,3) and their cyclic shifts; e_i^2 = -e0.
  * L = D8+ with the ODD-sum half-integer coset; s = (-e0+e1+..+e7)/2,
    sbar = (-e0-..-e7)/2; L sbar, L s are LEFT multiples {l.sbar}, {l.s}.
  * x ._g y = g^{-1}(g x . g y); for reflections g = g^{-1}.
  * Triple product P1 = x.x'+z.y'+y.z', P2 = y.y'+x.z'+z.x', P3 = z.z'+y.x'+x.y'.

Sections
  1. conventions pinned against numbers printed in the paper
  2. Lambda: Gram determinant, the three families of minimal vectors
  3. the reflection s_r: L -> L, octonion-algebra axioms for ._g, identity element
  4. closure and span index: untwisted, sigma, s_r (witness), C2's 3-cycle
  5. structural 2-adic factorisation  [Lambda:S]_2 = [B:B^2]_2 [A:A^2]_2^2
  6. census of all 120 root reflections (claim C3: 64 / 28 / 28)
  7. convention variant: right multiples sbar.l, s.l
"""
import random
import sys
import time
from fractions import Fraction as Fr
from collections import Counter

import leech_twist_lib as T

t0 = time.time()
random.seed(20261013)


def say(*a):
    print(*a)
    sys.stdout.flush()


def half(X):
    """doubled coords -> exact coordinates as strings"""
    return [str(Fr(t, 2)) for t in X]


# ==========================================================================
say("=" * 78)
say("1. CONVENTIONS PINNED AGAINST THE PAPER")
say("=" * 78)
G_L = T.L_LAT.gram()
say("L basis (Appendix A) in D8+:", all(T.in_L2(b) for b in T.L_BASIS2),
    "| Gram det:", T.det_bareiss([[int(x) for x in r] for r in G_L]))
say("s in L:", T.in_L2(T.S2), "| sbar in L:", T.in_L2(T.SBAR2))
say("[L : L sbar] =", T.LSBAR_LAT.index_in(T.L_LAT), "| [L : L s] =", T.LS_LAT.index_in(T.L_LAT))
both = T.hnf_rows([T.L_LAT.coords(r) for r in T.LSBAR_LAT.rows + T.LS_LAT.rows])
say("L sbar + L s = L:", len(both) == 8 and all(both[i][i] == 1 for i in range(8)))
twoL_in_both = all([2 * t for t in b] in T.LSBAR_LAT and [2 * t for t in b] in T.LS_LAT for b in T.L_BASIS2)
say("2L in L sbar and in L s:", twoL_in_both,
    "-> with 16*16/1 = 256 = [L:2L], L sbar cap L s = 2L")
# Remark 4.6 of the paper
u1 = [-2, 0, 0, 0, 2, 0, 2, 2]
u2 = [-2, 2, 0, 0, 0, 2, 0, 2]
p = T.omul2(u1, u2)
say("Remark 4.6: u1,u2 in L s:", u1 in T.LS_LAT, u2 in T.LS_LAT,
    "| u1.u2 =", half(p), "(paper (0,-2,2,1,-2,-1,-1,-1))",
    "| in L sbar:", p in T.LSBAR_LAT, "| in L s:", p in T.LS_LAT)
SIG = T.Twist(T.perm_matrix({1: 2, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}), "sigma=(1 2)")
M = [SIG.apply(T.omul2(b, T.S2)) for b in T.L_BASIS2]      # M_k = sigma(L_k . s)
N = [SIG.apply(T.omul2(b, T.SBAR2)) for b in T.L_BASIS2]   # N_k = sigma(L_k . sbar)
paper_M = [[-3, -1, -1, -1, -1, -1, -1, -1], [-2, 0, 0, 0, 2, 0, 2, 2], [0, 2, 2, 2, 0, 2, 0, 0],
           [-2, 0, 2, 0, 0, 2, 0, 2], [-2, 2, 2, 0, 0, 0, 2, 0], [-2, 2, 0, 2, 0, 0, 0, 2],
           [-2, 0, 2, 2, 2, 0, 0, 0], [-2, 2, 0, 0, 2, 2, 0, 0]]
paper_N = [[4, 0, 0, 0, 0, 0, 0, 0], [0, 0, -2, 0, -2, 0, -2, -2], [-2, -2, 0, -2, 0, -2, 0, 0],
           [0, -2, -2, 0, 0, -2, 0, -2], [0, -2, -2, -2, 0, 0, -2, 0], [0, -2, 0, -2, -2, 0, 0, -2],
           [0, 0, -2, -2, -2, -2, 0, 0], [0, -2, 0, 0, -2, -2, -2, 0]]
say("Appendix A: all 8 M_k = sigma(L_k.s) equal the paper's:", M == paper_M)
say("Appendix A: all 8 N_k = sigma(L_k.sbar) equal the paper's:", N == paper_N)
# a few entries of Table LL: L_i . L_j = sum c_k L_k
tblLL = {(1, 1): [-1, -1, -1, 0, 0, 0, 0, 0], (1, 2): [0, -2, -2, 1, 1, 0, 1, 0],
         (5, 2): [-2, -2, -4, 1, 2, 1, 1, 1], (8, 4): [-2, -3, -4, 2, 1, 1, 1, 2],
         (7, 6): [-2, -3, -4, 1, 1, 2, 2, 1]}
ok = all(T.L_LAT.coords(T.omul2(T.L_BASIS2[i - 1], T.L_BASIS2[j - 1])) == c for (i, j), c in tblLL.items())
say("Table LL entries (1,1),(1,2),(5,2),(8,4),(7,6) reproduced:", ok)

# ==========================================================================
say("\n" + "=" * 78)
say("2. WILSON'S LEECH LATTICE")
say("=" * 78)
G = T.LAMBDA_LAT.gram()
d = T.det_bareiss([[int(4 * x) for x in r] for r in G])
say("dim_F2 Lambda/2L^3 =", T._KERDIM, "| Gram det =", Fr(d, 4 ** 24))
say("all 24 basis vectors satisfy Wilson's conditions:", all(T.in_Lambda2(b) for b in T.LAMBDA_BASIS2))
# the three families of minimal vectors listed in the paper (Section 2.3)
roots = []
for i in range(8):
    for j in range(i + 1, 8):
        for si in (2, -2):
            for sj in (2, -2):
                v = [0] * 8
                v[i], v[j] = si, sj
                roots.append(v)
for signs in range(256):
    v = [1 if (signs >> k) & 1 else -1 for k in range(8)]
    if T.in_L2(v):
        roots.append(v)
assert len(roots) == 240 and all(T.nrm2(r) == 2 for r in roots)
J = [T.unit(t, sg * 2) for t in range(8) for sg in (1, -1)]   # +-e_t, doubled
def cyc(x, y, z):
    return [x + y + z, z + x + y, y + z + x]
minv = set()
Z8 = [0] * 8
for lam in roots:
    for v in cyc([2 * t for t in lam], Z8, Z8):            # (2 lam, 0, 0)
        minv.add(tuple(v))
    ls_bar = T.omul2(lam, T.SBAR2)
    ls = T.omul2(lam, T.S2)
    for j in J:
        for v in cyc(ls_bar, T.omul2(ls_bar, j), Z8):        # (lam sbar, (lam sbar) j, 0)
            minv.add(tuple(v))
        lsj = T.omul2(ls, j)
        lj = T.omul2(lam, j)
        for k in J:
            for v in cyc(lsj, T.omul2(lam, k), T.omul2(lj, k)):  # ((lam s) j, lam k, (lam j) k)
                minv.add(tuple(v))
minv = [list(v) for v in minv]
say("paper's three families: distinct vectors =", len(minv),
    "| all norm 8:", all(T.nrm2(v) == 8 for v in minv),
    "| all in Lambda:", all(T.in_Lambda2(v) for v in minv))
say("(minimum 8 and the count 196560 of norm-8 vectors are also checked by qfminim in c3_witness.gp)")

# ==========================================================================
say("\n" + "=" * 78)
say("3. THE REFLECTION s_r, r = (e0+e1+e2+e3+e4+e5+e6-e7)/2")
say("=" * 78)
r2 = [1, 1, 1, 1, 1, 1, 1, -1]
say("r in L:", T.in_L2(r2), "| N(r) =", T.nrm2(r2))
R = T.Twist(T.reflection(r2), "s_r")
Gm = R.G
say("s_r^2 = I:", T.matmul(Gm, Gm) == [[Fr(int(i == j)) for j in range(8)] for i in range(8)],
    "| s_r orthogonal:", T.matmul(T.transpose(Gm), Gm) == [[Fr(int(i == j)) for j in range(8)] for i in range(8)],
    "| det =", R.det)
say("s_r maps each L-basis vector into L, det +-1  =>  s_r(L) = L : True (asserted in Twist)")
def gfrac(tw, x):
    return [sum(a * b for a, b in zip(row, x)) for row in tw.G]
say("s_r(s) =", half(R.apply(T.S2)), "| s_r(sbar) =", [str(x) for x in gfrac(R, [Fr(-1, 2)] * 8)])
eps = R.identity_element()     # g^{-1}(e0) = s_r(e0)
say("identity element of ._g:  eps = s_r(e0) =", [str(x) for x in eps])
e = [[Fr(int(i == j)) for j in range(8)] for i in range(8)]
say("eps ._g e_k = e_k ._g eps = e_k for k = 0..7:",
    all(R.tmul_frac(eps, e[k]) == e[k] and R.tmul_frac(e[k], eps) == e[k] for k in range(8)))
say("e0 is NOT the identity of ._g:  e0 ._g e1 =", [str(x) for x in R.tmul_frac(e[0], e[1])])
twoeps_d = [int(4 * x) for x in eps]     # doubled coordinates of 2*eps (= 4*eps)
say("eps in L:", all((2 * x).denominator == 1 for x in eps) and T.in_L2([int(2 * x) for x in eps]),
    "| 2 eps in L:", T.in_L2(twoeps_d), "(so 4 eps in 2L and (4eps,0,0) in 2L^3 in Lambda)")
say("2 eps in L sbar:", twoeps_d in T.LSBAR_LAT, "| 2 eps in L s:", twoeps_d in T.LS_LAT)
Z8 = [0] * 8
cands = {"(2eps,0,0)": twoeps_d + Z8 + Z8, "(2eps,2eps,0)": twoeps_d + twoeps_d + Z8,
         "(2eps,-2eps,0)": twoeps_d + [-t for t in twoeps_d] + Z8, "(2eps,2eps,2eps)": twoeps_d * 3}
for k, v in cands.items():
    say(f"   {k:18s} in Lambda: {T.in_Lambda2(v)}")
# octonion-algebra axioms on random exact rational elements
def rnd():
    return [Fr(random.randint(-9, 9), random.choice([1, 2, 3, 4, 5])) for _ in range(8)]
def sub(a, b):
    return [p - q for p, q in zip(a, b)]
comp = altL = altR = flex = moufang = 0
NTEST = 300
for _ in range(NTEST):
    x, y, z = rnd(), rnd(), rnd()
    m = R.tmul_frac
    comp += T.nrm(m(x, y)) == T.nrm(x) * T.nrm(y)
    altL += m(m(x, x), y) == m(x, m(x, y))
    altR += m(m(y, x), x) == m(y, m(x, x))
    flex += m(m(x, y), x) == m(x, m(y, x))
    moufang += m(m(m(z, x), y), x) == m(z, m(m(x, y), x))
say(f"on {NTEST} random rational triples: composition N(x._g y)=N(x)N(y): {comp}/{NTEST}; "
    f"left alt: {altL}/{NTEST}; right alt: {altR}/{NTEST}; flexible: {flex}/{NTEST}; Moufang ((zx)y)x=z((xy)x): {moufang}/{NTEST}")
# non-associativity control (a genuine octonion algebra is not associative)
x, y, z = rnd(), rnd(), rnd()
say("associativity fails on a random triple (as it must):", R.tmul_frac(R.tmul_frac(x, y), z) != R.tmul_frac(x, R.tmul_frac(y, z)))
closedL = all(T.in_L2(R.tmul2(a, b)) for a in T.L_BASIS2 for b in T.L_BASIS2)
say("L ._g L in L (64 basis products):", closedL)

# ==========================================================================
say("\n" + "=" * 78)
say("4. CLOSURE AND SPAN INDEX (576 basis-pair products, HNF + SNF)")
say("=" * 78)
ID = T.Twist([[int(i == j) for j in range(8)] for i in range(8)], "identity")
C3CYC = T.Twist(T.perm_matrix({1: 1, 2: 2, 3: 3, 4: 4, 5: 6, 6: 7, 7: 5}),
                "3-cycle e5->e6->e7->e5 (C2 witness; pi(5)=6, pi(6)=7, pi(7)=5)")
C3CYCinv = T.Twist(T.perm_matrix({1: 1, 2: 2, 3: 3, 4: 4, 5: 7, 6: 5, 7: 6}),
                   "inverse 3-cycle e5->e7->e6->e5")
results = {}
for tw in (ID, SIG, R, C3CYC, C3CYCinv):
    o = T.closure_and_index(tw)
    results[tw.name] = o
    if o["closes"]:
        say(f"{tw.name:60s} closes; [Lambda:S] = {o['index']} ; Lambda/S = "
            + ("0" if not o["elementary_divisors"] else " + ".join(f"Z/{d}" for d in o["elementary_divisors"])))
        twoL = all(T.hnf_rows(o["hnf"] + [[2 * int(i == j) for j in range(24)]])  == o["hnf"] for i in range(24))
        say(f"{'':60s} 2*Lambda inside S: {twoL}")
    else:
        say(f"{tw.name:60s} does NOT close ({o['failing_pairs']} of 576 basis-pair products leave Lambda)")

# ==========================================================================
say("\n" + "=" * 78)
say("5. STRUCTURAL 2-ADIC FACTORISATION")
say("=" * 78)
say("""Facts used (proved in the report): R^24 = D + T, D = {(a,a,a)}, T = {sum 0},
D*T = T*D = 0 for ANY bilinear block product.  3*pi_D(Lambda) lies in Lambda, so
[Lambda : (Lambda cap D) + (Lambda cap T)] is a power of 3 and, 2-adically,
Lambda = Lambda_D + Lambda_T with Lambda_D = {(a,a,a): a in L s} and
Lambda_T = {(x,y,-x-y): x,y in L sbar}.  With A = g(L sbar), B = g(L s) and the
STANDARD product, (a,a,a)*(b,b,b) = 3(a.b)(1,1,1), and on Lambda_T the products
contain (a,-2a), (a,a), so (0,3a) for a in A.A; hence
      [Lambda : S]_2 = [B : span(B.B)]_2 * [A : span(A.A)]_2^2 .""")


def square_index(lat_rows):
    """[X : span(X.X)] for a sublattice X of L given by basis rows (doubled)."""
    X = T.Lattice(lat_rows)
    prods = [T.omul2(a, b) for a in lat_rows for b in lat_rows]
    assert all(p in X for p in prods), "X is not closed"
    return T.sub_index(prods, X), X


def two_part(n):
    k = 0
    while n % 2 == 0 and n:
        n //= 2
        k += 1
    return k, n


for tw in (SIG, R, C3CYC):
    A_rows = [tw.apply(b) for b in T.LSBAR_LAT.rows]
    B_rows = [tw.apply(b) for b in T.LS_LAT.rows]
    iA, _ = square_index(A_rows)
    iB, _ = square_index(B_rows)
    idx = results[tw.name]["index"]
    kA, oA = two_part(iA)
    kB, oB = two_part(iB)
    kI, oI = two_part(idx)
    say(f"{tw.name:60s} [B:B^2] = 2^{kB}*{oB}, [A:A^2] = 2^{kA}*{oA}  ->  predicted v2 = {kB}+2*{kA} = {kB + 2 * kA};"
        f" actual [Lambda:S] = 2^{kI}*{oI}")
# also the untwisted pair (closure fails there, but the sublattices are still defined)
iAu, _ = square_index(T.LSBAR_LAT.rows) if True else (None, None)
try:
    iBu, _ = square_index(T.LS_LAT.rows)
except AssertionError:
    iBu = "n/a (L s is not closed under the product: Remark 4.6)"
say(f"{'untwisted pair (A = L sbar, B = L s)':60s} [A:A^2] = {iAu}, [B:B^2] = {iBu}")

# ==========================================================================
say("\n" + "=" * 78)
say("6. CENSUS OF THE 120 ROOT REFLECTIONS (one per +-r pair)")
say("=" * 78)
reps = []
seen = set()
for r in roots:
    key = tuple(r)
    neg = tuple(-t for t in r)
    if neg in seen:
        continue
    seen.add(key)
    reps.append(r)
assert len(reps) == 120


def root_type(r):
    if all(t % 2 for t in r):
        return "half-integer (e0-coefficient +-1/2)"
    nz = [i for i in range(8) if r[i]]
    if 0 in nz:
        return "+-e0 +- e_j"
    return "+-e_i +- e_j (i,j>=1)"


census = Counter()
by_type = Counter()
index1 = []
for r in reps:
    tw = T.Twist(T.reflection(r), "")
    o = T.closure_and_index(tw)
    key = ("no closure",) if not o["closes"] else ("closes", o["index"])
    census[key] += 1
    by_type[(root_type(r), key)] += 1
    if o["closes"] and o["index"] == 1:
        index1.append(r)
for k, v in sorted(census.items(), key=lambda kv: str(kv[0])):
    say(f"  {str(k):30s} {v}")
say("  by root type:")
for (rt, k), v in sorted(by_type.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
    say(f"     {rt:40s} {str(k):22s} {v}")
say("  sigma's root e1-e2 is in the index-2^16 class:",
    T.closure_and_index(T.Twist(T.reflection([0, 2, -2, 0, 0, 0, 0, 0])))["index"] == 2 ** 16)
say("  witness root (e0+...+e6-e7)/2 is in the index-1 class:", any(r == r2 or r == [-t for t in r2] for r in index1))
say("  the 28 index-1 roots (doubled coordinates):")
for r in index1:
    say("     ", r)

# ==========================================================================
say("\n" + "=" * 78)
say("7. CONVENTION VARIANT: RIGHT MULTIPLES sbar.l, s.l INSTEAD OF LEFT MULTIPLES")
say("=" * 78)
LSBAR_R = T.Lattice(T.right_multiples(T.L_BASIS2, T.SBAR2))
LS_R = T.Lattice(T.right_multiples(T.L_BASIS2, T.S2))
say("sbar.L == L.sbar:", T.hnf_rows([T.L_LAT.coords(r) for r in LSBAR_R.rows]) ==
    T.hnf_rows([T.L_LAT.coords(r) for r in T.LSBAR_LAT.rows]),
    "| s.L == L.s:", T.hnf_rows([T.L_LAT.coords(r) for r in LS_R.rows]) ==
    T.hnf_rows([T.L_LAT.coords(r) for r in T.LS_LAT.rows]))
basisR, kd = T.leech_basis(LSBAR_R, LS_R)
latR = T.Lattice(basisR)
GR = latR.gram()
dR = T.det_bareiss([[int(4 * x) for x in r] for r in GR])
say("variant lattice: Gram det =", Fr(dR, 4 ** 24))


def closure_variant(tw):
    prods = T.basis_products(tw, basisR)
    bad = sum(not T.in_Lambda2(P, LSBAR_R, LS_R) for P in prods)
    if bad:
        return f"does not close ({bad} failing pairs)"
    H = T.hnf_rows([latR.coords(P) for P in prods])
    idx = 1
    for i in range(24):
        idx *= H[i][i]
    return f"closes, index {idx}"


for tw in (ID, SIG, R):
    say(f"  {tw.name:40s} {closure_variant(tw)}")
vcount = Counter()
for r in reps:
    vcount[closure_variant(T.Twist(T.reflection(r)))] += 1
say("  all 120 root reflections under the right-multiple convention:")
for k, v in sorted(vcount.items()):
    if k.startswith("does not close"):
        continue
    say(f"     {k:30s} {v}")
say(f"     {'does not close':30s} {sum(v for k, v in vcount.items() if k.startswith('does not close'))}")
say(f"\nruntime {time.time() - t0:.1f} s")
