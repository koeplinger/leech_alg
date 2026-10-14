"""
f2_witnesses.py -- how each span index is obtained, twist by twist.

Twists examined (conventions in leech_final.py; x ._g y = g^{-1}(g x . g y)):
  identity                  pi = (1,2,3,4,5,6,7)            (untwisted, reference)
  sigma = (1 2)             pi = (2,1,3,4,5,6,7)            (the paper's twist)
  C2 witness  e5->e6->e7    pi = (1,2,3,4,6,7,5), g = g_pi  (pi(5)=6, pi(6)=7, pi(7)=5)
  its inverse e5->e7->e6    pi = (1,2,3,4,7,5,6)
  C3 witness  s_r,  r = (e0+e1+e2+e3+e4+e5+e6-e7)/2       (root reflection, involution)

For every twist:
  (a) the 576 products b_p *_g b_q of the Lambda basis (p = left factor) are
      formed in L^3-coordinates from the twisted structure constants; for
      permutations the constants are computed by TWO routes (L-basis matrix of g,
      and the twisted e-basis table) and compared;
  (b) closure is tested THREE ways, product by product: integrality of the
      Lambda-coordinates, the 12 mod-2 Wilson conditions, and the definitional
      test in e-coordinates ((t.s)/2 in L etc., no basis involved);
  (c) for a non-closing twist, the number of failing basis products (basis
      dependent!) and the F2-dimension of the obstruction subgroup
      (S + Lambda)/Lambda of L^3/Lambda (basis independent) are reported;
  (d) for a closing twist:
        r2   = F2-rank of the 576 x 24 coordinate matrix -> [Lambda : S+2Lambda] = 2^(24-r2)
        HNF  = exact Hermite normal form over Z (pure Python ints, no modular
               reduction, no theorem used) -> [Lambda : S] = product of pivots
        Z/4  = structure (Z/4)^a (Z/2)^b of (S + 4Lambda)/4Lambda; if
               2^(48-2a-b) equals the HNF index then 4Lambda <= S and
               Lambda/S = (Z/2)^b x (Z/4)^(24-a-b) exactly
        8-dim factors (index formula, proved in the review):
               hW = [Ls : <Ls ._g Ls>],  hV = [Lsbar : <Lsbar ._g Lsbar>],
               predicted index hW * hV^2
        the twisted identity 1_g = g^{-1}(e0) and whether 2*1_g lies in Lsbar / Ls
  (e) premise of the theorem 4 Lambda <= S: <L ._g L> = L.
For sigma the span is identified (open question of the paper's Sec. 5.5):
  S = Delta(Ls) + 2 Lambda = {(x,y,z) : x,y,z in Ls ; x-y, x-z in 2 Lsbar},
  Delta(w) = (w,w,w); both equalities are checked as equalities of exact HNFs.
The coordinate matrices are written to f2_matrices.gp for PARI (f6).
"""
import numpy as np
from fractions import Fraction
from leech_final import (Setup, perm_map, perm_inverse, reflection_map, gamma_of, gamma_inverse,
                         twisted_constants, twisted_table_e, structure_constants, triple_products,
                         lambda_coords, mod2_ok, f2_rank, f2_span_basis, z4_structure, hnf,
                         hnf_index, in_Lambda_def, from_Lcoords, in_L, in_Ls_def, in_Lsbar_def,
                         h_index, bil, int_det, C_STD)

st = Setup()
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


def to_e(P):
    """L^3-coordinates -> doubled e-coordinates (24 ints)."""
    return sum((from_Lcoords(P[8 * k:8 * k + 8]) for k in range(3)), [])


def factor_twos(n):
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k, n


TWISTS = [
    ("identity (untwisted)", "perm", (1, 2, 3, 4, 5, 6, 7)),
    ("sigma = (1 2)", "perm", (2, 1, 3, 4, 5, 6, 7)),
    ("C2 witness: 3-cycle e5->e6->e7->e5, g = g_pi", "perm", (1, 2, 3, 4, 6, 7, 5)),
    ("inverse 3-cycle e5->e7->e6->e5", "perm", (1, 2, 3, 4, 7, 5, 6)),
    ("C3 witness: reflection s_r, r = (e0+...+e6-e7)/2", "refl", (1, 1, 1, 1, 1, 1, 1, -1)),
]

gp_lines = ["\\\\ 576 x 24 integer Lambda-coordinate matrices of the basis-pair products (f2_witnesses.py)"]
summary = []
for name, kind, data in TWISTS:
    say("=" * 100)
    say("TWIST:", name, "  data =", data)
    if kind == "perm":
        f = perm_map(data)
        finv = perm_map(perm_inverse(data))
    else:
        f = reflection_map(data)
        finv = f
    G = gamma_of(f)
    detg = int_det(G.tolist())
    Cg = twisted_constants(G)
    if kind == "perm":
        Cg2 = structure_constants(twisted_table_e(f, finv))
        assert np.array_equal(Cg, Cg2), "two routes to the twisted constants disagree"
        say("  twisted structure constants: L-matrix route == e-table route (64 x 8 entries)")
    say("  det g =", detg)
    # premise of 4 Lambda <= S :  <L ._g L> = L
    iLL, _ = hnf_index(Cg.reshape(64, 8).tolist(), 8)
    say("  [L : <L ._g L>] =", iLL)
    assert iLL == 1
    # (a) products
    P = triple_products(st.B, Cg)
    coords, ok1 = lambda_coords(P, st.N)
    ok2 = mod2_ok(P, st.cond)
    ok3 = np.array([in_Lambda_def(to_e(list(map(int, p)))) for p in P])
    assert np.array_equal(ok1, ok2) and np.array_equal(ok1, ok3)
    nfail = int((~ok1).sum())
    obstruction = f2_rank((st.cond @ (P.T % 2)) % 2)
    say("  closure: failing basis products = %d / 576 (three tests agree product by product)" % nfail)
    say("  F2-dimension of the obstruction subgroup (S+Lambda)/Lambda in L^3/Lambda =", obstruction,
        "(basis independent)")
    if nfail:
        summary.append((name, "does not close", nfail, obstruction))
        continue
    # (d) index
    r2 = f2_rank(coords)
    say("  F2-rank of the 576x24 coordinate matrix r2 = %d  ->  [Lambda : S+2Lambda] = 2^%d" % (r2, 24 - r2))
    idx, H = hnf_index(coords.tolist(), 24)
    diag = [H[k][k] for k in range(24)]
    k2, odd = factor_twos(idx)
    say("  exact HNF (Python ints, no modular reduction): pivots =", diag)
    say("  [Lambda : S] = %d = 2^%d * %d" % (idx, k2, odd))
    a, b = z4_structure(coords)
    z4idx = 2 ** (48 - 2 * a - b)
    say("  Z/4 structure of (S+4Lambda)/4Lambda: (Z/4)^%d x (Z/2)^%d -> [Lambda : S+4Lambda] = %d" % (a, b, z4idx))
    assert z4idx == idx, "4 Lambda not inside S"
    grp = []
    if b:
        grp.append("(Z/2)^%d" % b)
    if 24 - a - b:
        grp.append("(Z/4)^%d" % (24 - a - b))
    say("  => 4Lambda <= S and Lambda/S =", " x ".join(grp) if grp else "0  (S = Lambda)")
    say("  2Lambda <= S:", idx == 2 ** (24 - r2))
    # 8-dimensional factors
    Vg = f2_span_basis((np.array(st.V0) @ G) % 2)     # g(Lsbar)/2L
    Wg = f2_span_basis((np.array(st.W0) @ G) % 2)     # g(Ls)/2L
    hV, _, _ = h_index(Vg, C_STD)
    hW, _, _ = h_index(Wg, C_STD)
    hV_tw, _, _ = h_index(st.V0, Cg)                  # same number in the twisted frame
    hW_tw, _, _ = h_index(st.W0, Cg)
    assert (hV, hW) == (hV_tw, hW_tw)
    say("  8-dim factors: hW = [Ls : <Ls ._g Ls>] = %d,  hV = [Lsbar : <Lsbar ._g Lsbar>] = %d" % (hW, hV))
    say("  index formula hW * hV^2 = %d   (direct HNF: %d)  agree: %s" % (hW * hV * hV, idx, hW * hV * hV == idx))
    assert hW * hV * hV == idx
    # twisted identity 1_g = g^{-1}(e0), exact rationals (e0 itself is not in L)
    if kind == "perm":
        one_g = [Fraction(1)] + [Fraction(0)] * 7
    else:
        r = [Fraction(t, 2) for t in data]
        ip = r[0]                                   # <e0, r>
        one_g = [Fraction(int(i == 0)) - ip * r[i] for i in range(8)]
    two1 = [2 * t for t in one_g]                   # 2*1_g
    assert all((2 * t).denominator == 1 for t in two1)
    two1_doubled = [int(2 * t) for t in two1]
    inL = in_L(two1_doubled)
    say("  twisted identity 1_g = g^{-1}(e0) =", "(" + ", ".join(str(t) for t in one_g) + ")")
    say("  2*1_g in L: %s,  in Lsbar: %s,  in Ls: %s" % (inL, inL and in_Lsbar_def(two1_doubled),
                                                         inL and in_Ls_def(two1_doubled)))
    summary.append((name, "closes", idx, diag.count(2), r2, (hW, hV)))
    tag = "t%d" % (len(gp_lines))
    gp_lines.append("\\\\ %s : %s" % (tag, name))
    gp_lines.append("%s = %s;" % (tag, "[" + ";".join(",".join(str(int(x)) for x in row) for row in coords) + "]"))
    if name.startswith("sigma"):
        # identify S for sigma: S = Delta(Ls) + 2 Lambda
        diagL = [list(w) * 3 for w in st.Ls_L]                    # (w,w,w) in L^3 coords
        twoLam = [[2 * int(x) for x in row] for row in st.B]
        gens = []
        for v in diagL + twoLam:
            c, ok = lambda_coords(np.array([v], dtype=np.int64), st.N)
            assert ok[0]
            gens.append([int(x) for x in c[0]])
        H2 = hnf(gens, 24)
        say("  sigma: HNF(S) == HNF(Delta(Ls) + 2 Lambda):", H2 == H)
        assert H2 == H
        # equivalently S = {(x,y,z) : x,y,z in Ls ; x-y, x-z in 2 Lsbar}
        gens = []
        for v in diagL + [[0] * 8 + [2 * t for t in x] + [0] * 8 for x in st.Lsbar_L] \
                       + [[0] * 16 + [2 * t for t in x] for x in st.Lsbar_L]:
            c, ok = lambda_coords(np.array([v], dtype=np.int64), st.N)
            assert ok[0]
            gens.append([int(x) for x in c[0]])
        H3 = hnf(gens, 24)
        say("  sigma: HNF(S) == HNF({x,y,z in Ls; x-y, x-z in 2Lsbar}):", H3 == H)
        assert H3 == H

say("=" * 100)
say("SUMMARY")
for s in summary:
    say("  ", s)
open("f2_witnesses.out.txt", "w").write("\n".join(out) + "\n")
open("f2_matrices.gp", "w").write("\n".join(gp_lines) + "\n")
