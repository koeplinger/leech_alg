"""
f1_setup.py -- build O, L, Ls, Lsbar and Wilson's Lambda from the definitions
and pin every convention against numbers PRINTED in the paper (v8, read-only).

Conventions: see the header of leech_final.py (doubled coordinates, L = D8+
odd convention, LEFT multiples l.s and l.sbar, Appendix-A basis L_1..L_8).

Checks performed (all exact):
 1. The Fano table defines an octonion algebra: composition N(xy) = N(x)N(y),
    left/right alternativity, on 2,000 random integer pairs.
 2. L_1..L_8 lie in L and have covolume 1 (= covol D8+), so they are a Z-basis;
    all 240 roots have integer L-coordinates.  span(L.L) = L.
 3. The paper's Tables LL, LN, MM (Appendix A, 3 x 64 rows, parsed from
    main.tex) and the printed vectors M_k = sigma(L_k.s), N_k = sigma(L_k.sbar)
    are reproduced exactly; Remark 4.5's three membership claims hold.
 4. [L:Lsbar] = [L:Ls] = 16, Lsbar + Ls = L, Lsbar cap Ls = 2L; Lsbar/2L and
    Ls/2L are totally singular for q(x) = N(x)/2 mod 2.
 5. Lambda: [L^3 : Lambda] = 2^12; every basis vector satisfies Wilson's
    definition (definitional test, no basis); Gram det = 2^24 and Gram/2 is
    even unimodular.  The 196,560 vectors of the paper's three families of
    minimal vectors are distinct, of norm 8 and lie in Lambda.
    The Gram matrix is written to lambda_gram_final.gp for PARI's qfminim
    (f6_pari_independent.gp re-derives the minimum independently).
"""
import re
import random
import numpy as np
from fractions import Fraction
from leech_final import (T_STD, omul, dmul, conj, norm2, in_L, in_Lsbar_def, in_Ls_def,
                         in_Lambda_def, LB, Lcoords, from_Lcoords, C_STD, Setup,
                         hnf_index, int_det, f2_rank, perm_map, S2, SBAR2)

import os as _os
# the paper source (v8); default: the project root, three levels up from this folder.
# Override with the environment variable MAIN_TEX if the folder is moved.
MAIN_TEX = _os.environ.get("MAIN_TEX", _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "..", "main.tex"))
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


# ---------------------------------------------------------------- 1. octonions
random.seed(20261013)
for _ in range(2000):
    x = [random.randint(-5, 5) for _ in range(8)]
    y = [random.randint(-5, 5) for _ in range(8)]
    xy = omul(x, y)
    assert sum(t * t for t in xy) == sum(t * t for t in x) * sum(t * t for t in y)
    assert omul(omul(x, x), y) == omul(x, omul(x, y))
    assert omul(omul(x, y), y) == omul(x, omul(y, y))
say("[1] Fano table: composition and left/right alternativity hold on 2000 random pairs")

# ---------------------------------------------------------------- 2. L
assert all(in_L(v) for v in LB)
d = int_det(LB)
say("[2] det of doubled basis matrix =", d, "-> covolume of span(L_k) =", Fraction(d, 2 ** 8))
assert abs(d) == 2 ** 8
roots = []
for X in __import__('itertools').product(range(-2, 3), repeat=8):
    if sum(t * t for t in X) == 8 and in_L(X):
        roots.append(list(X))
assert len(roots) == 240
for r in roots:
    Lcoords(r)
say("    240 roots found in L (doubled norm 8 <=> N = 2); all have integer L-coordinates")
prods = [list(C_STD[i, j]) for i in range(8) for j in range(8)]
idx, _ = hnf_index(prods, 8)
say("    [L : span(L.L)] =", idx)
assert idx == 1

# ---------------------------------------------------------------- 3. paper tables
tex = open(MAIN_TEX).read()


def parse_table(label):
    end = tex.index("\\label{%s}" % label)
    start = tex.rindex("\\begin{table}", 0, end)
    body = tex[start:end]
    rows = {}
    for m in re.finditer(r"^\s*(\d+) & (\d+)((?: & -?\d+){8}) \\\\", body, re.M):
        i, j = int(m.group(1)), int(m.group(2))
        rows[(i, j)] = [int(t) for t in m.group(3).split("&")[1:]]
    assert len(rows) == 64, (label, len(rows))
    return rows


LL = parse_table("tbl:appendix-LL")
LN = parse_table("tbl:appendix-LN")
MM = parse_table("tbl:appendix-MM")
sig = perm_map((2, 1, 3, 4, 5, 6, 7))
M = [sig(dmul(LB[k], S2)) for k in range(8)]
Nv = []
for k in range(8):
    P = omul(LB[k], SBAR2)
    Nv.append(sig([p // 2 for p in P]))
printed_M = [[-3, -1, -1, -1, -1, -1, -1, -1], [-2, 0, 0, 0, 2, 0, 2, 2], [0, 2, 2, 2, 0, 2, 0, 0],
             [-2, 0, 2, 0, 0, 2, 0, 2], [-2, 2, 2, 0, 0, 0, 2, 0], [-2, 2, 0, 2, 0, 0, 0, 2],
             [-2, 0, 2, 2, 2, 0, 0, 0], [-2, 2, 0, 0, 2, 2, 0, 0]]          # doubled, from Appendix A
printed_N = [[4, 0, 0, 0, 0, 0, 0, 0], [0, 0, -2, 0, -2, 0, -2, -2], [-2, -2, 0, -2, 0, -2, 0, 0],
             [0, -2, -2, 0, 0, -2, 0, -2], [0, -2, -2, -2, 0, 0, -2, 0], [0, -2, 0, -2, -2, 0, 0, -2],
             [0, 0, -2, -2, -2, -2, 0, 0], [0, -2, 0, 0, -2, -2, -2, 0]]
assert M == printed_M, "M_k mismatch"
assert Nv == printed_N, "N_k mismatch"
say("[3] printed M_k = sigma(L_k.s) and N_k = sigma(L_k.sbar) reproduced (8 + 8 vectors)")


def coords_in(basis, X):
    Bi = __import__('leech_final').frac_inverse(basis)
    c = [sum(Fraction(X[i]) * Bi[i][k] for i in range(8)) for k in range(8)]
    assert all(t.denominator == 1 for t in c)
    return [int(t) for t in c]


bad = 0
for i in range(8):
    for j in range(8):
        bad += list(C_STD[i, j]) != LL[(i + 1, j + 1)]
        bad += coords_in(Nv, dmul(LB[i], Nv[j])) != LN[(i + 1, j + 1)]
        bad += coords_in(M, dmul(M[i], M[j])) != MM[(i + 1, j + 1)]
say("    Tables LL, LN, MM: 192 printed rows compared, mismatches =", bad)
assert bad == 0
u1 = [-2, 0, 0, 0, 2, 0, 2, 2]
u2 = [-2, 2, 0, 0, 0, 2, 0, 2]
u12 = dmul(u1, u2)
assert in_Ls_def(u1) and in_Ls_def(u2) and u12 == [0, -4, 4, 2, -4, -2, -2, -2]
assert in_Lsbar_def(u12) and not in_Ls_def(u12)
v = [-2, 0, 2, 0, 0, 2, 0, 2]
w = [0, 0, -2, 0, -2, 0, -2, -2]
assert not in_Ls_def(v) and in_Ls_def(sig(v))         # v in sigma(Ls) \ Ls
assert not in_Lsbar_def(w) and in_Lsbar_def(sig(w))   # w in sigma(Lsbar) \ Lsbar
say("    Remark 4.5 reproduced: u1,u2 in Ls, u1.u2 = (0,-2,2,1,-2,-1,-1,-1) in Lsbar\\Ls;",
    "v in sigma(Ls)\\Ls; w in sigma(Lsbar)\\Lsbar")

# ---------------------------------------------------------------- 4. Ls, Lsbar
st = Setup()
iLs, _ = hnf_index(st.Ls_L, 8)
iLsb, _ = hnf_index(st.Lsbar_L, 8)
isum, _ = hnf_index(st.Ls_L + st.Lsbar_L, 8)
say("[4] [L:Ls] =", iLs, " [L:Lsbar] =", iLsb, " [L : Ls+Lsbar] =", isum)
assert iLs == 16 and iLsb == 16 and isum == 1
# 2L inside both (each 2 L_k has integer coordinates in the bases), so by
# counting |L/2L| = 256 = 16*16/[L:Ls+Lsbar], the intersection is exactly 2L.
for k in range(8):
    e = [0] * 8
    e[k] = 2
    assert hnf_index(st.Ls_L + [e], 8)[0] == 16 and hnf_index(st.Lsbar_L + [e], 8)[0] == 16
say("    2L <= Ls and 2L <= Lsbar, hence Ls cap Lsbar = 2L (index count 16*16 = 256)")
# definitional test agrees with basis test on all 256 classes of L/2L
agree = 0
for m in range(256):
    c = [(m >> k) & 1 for k in range(8)]
    X = from_Lcoords(c)
    inV = f2_rank(np.vstack([st.V0, np.array(c)])) == 4
    inW = f2_rank(np.vstack([st.W0, np.array(c)])) == 4
    agree += (inV == in_Lsbar_def(X)) and (inW == in_Ls_def(X))
say("    definitional tests (t.s)/2 in L, (t.sbar)/2 in L agree with the bases on", agree, "/ 256 classes")
assert agree == 256


def q_of(c):
    X = from_Lcoords(c)
    return int(norm2(X) / 2) % 2


for U in (st.V0, st.W0):
    for m in range(16):
        c = np.zeros(8, dtype=np.int64)
        for k in range(4):
            if (m >> k) & 1:
                c = (c + U[k]) % 2
        assert q_of(list(c)) == 0
say("    Lsbar/2L and Ls/2L are totally singular for q = N/2 mod 2 (Lagrangians)")

# ---------------------------------------------------------------- 5. Lambda
dB = int_det(st.B.tolist())
say("[5] det of Lambda basis in L^3-coordinates =", dB, "= 2^12 -> [L^3:Lambda] = 4096")
assert dB == 2 ** 12
Be = st.basis_e()
assert all(in_Lambda_def(b) for b in Be)
say("    all 24 basis vectors satisfy Wilson's definition (definitional test)")
G = [[Fraction(sum(a * b for a, b in zip(u, v)), 4) for v in Be] for u in Be]
assert all(g.denominator == 1 for row in G for g in row)
Gi = [[int(g) for g in row] for row in G]
dG = int_det(Gi)
say("    Gram determinant =", dG, "= 2^24:", dG == 2 ** 24)
assert dG == 2 ** 24
assert all(Gi[i][j] % 2 == 0 for i in range(24) for j in range(24))
G2 = [[g // 2 for g in row] for row in Gi]
assert int_det(G2) == 1 and all(G2[i][i] % 2 == 0 for i in range(24))
say("    Gram/2 is integral, even, det 1 (Lambda/sqrt2 even unimodular)")
with open("lambda_gram_final.gp", "w") as f:
    f.write("\\\\ Gram matrix (norm = squared length) of the Lambda basis of f1_setup.py\n")
    f.write("G = [" + ";".join(",".join(str(x) for x in row) for row in Gi) + "];\n")

# minimal-vector families of the paper (Section 2.3) -- vectorized check
LBinvD = None
from leech_final import LB_INV
den = 1
for row in LB_INV:
    for x in row:
        den = den * x.denominator // __import__('math').gcd(den, x.denominator)
LBinvD = np.array([[int(x * den) for x in row] for row in LB_INV], dtype=np.int64)


def lc_np(X):
    Q = X @ LBinvD
    assert np.all(Q % den == 0)
    return Q // den


J = []
for t in range(8):
    for sgn in (1, -1):
        e = [0] * 8
        e[t] = 2 * sgn
        J.append(e)
fam = []
R = roots
for lam in R:
    fam.append([2 * a for a in lam] + [0] * 16)
for lam in R:
    ls = [p // 2 for p in omul(lam, SBAR2)]      # doubled lam.sbar
    for j in J:
        fam.append(ls + dmul(ls, j) + [0] * 8)
for lam in R:
    lsj_base = dmul(lam, S2)
    for j in J:
        a = dmul(lsj_base, j)
        lj = dmul(lam, j)
        for k in J:
            fam.append(a + dmul(lam, k) + dmul(lj, k))
cyc = []
for v in fam:
    x, y, z = v[0:8], v[8:16], v[16:24]
    cyc += [x + y + z, z + x + y, y + z + x]
A = np.array(cyc, dtype=np.int64)
assert len(set(map(tuple, A.tolist()))) == 196560
assert np.all((A * A).sum(axis=1) == 32)            # doubled norm 32 <=> N = 8
Lc = np.hstack([lc_np(A[:, 0:8]), lc_np(A[:, 8:16]), lc_np(A[:, 16:24])])
ok = np.all((st.cond @ (Lc.T % 2)) % 2 == 0, axis=0)
say("    paper's three families: 196560 distinct vectors, all of norm 8, in Lambda:", bool(ok.all()))
assert ok.all()

open("f1_setup.out.txt", "w").write("\n".join(out) + "\n")
