"""
f3_certificates.py -- explicit, independently checkable certificates for the
span indices, and re-verification of the certificates produced by the other
review agents, using only THIS folder's code (leech_final.py).

Conventions: leech_final.py.  The index-1 twist is g = s_r with
r = (e0+e1+e2+e3+e4+e5+e6-e7)/2 (an involution, so g = g^{-1}).

A. Own certificate for S = Lambda (index 1), built from basis-pair products:
   two (or more) sets of 24 products b_p *_g b_q whose 24x24 Lambda-coordinate
   determinants are odd and have gcd 1.  Since each set spans a sublattice of
   S of index |det|, [Lambda : S] divides every |det|, hence divides their gcd = 1.
   (A set with odd determinant is found by choosing products that are linearly
   independent modulo 2, which is possible because the F2-rank is 24.)
   Exported to cert_index1_final.gp, checked from scratch by f6_pari_independent.gp.
B. The refute agent's certificate: 24 single products u_k *_g v_k (general
   vectors u_k, v_k of Lambda) with determinant +-1, i.e. a Z-basis of Lambda
   consisting of products.  Re-verified here: factors in Lambda by Wilson's
   definition, products recomputed with this folder's octonion code, Lambda-
   coordinate determinant and Gram determinant.
C. The GAP agent's certificate (two coprime determinants 2^44 and odd).
D. The structure agent's certificates for sigma (2^16), the 3-cycle (2^10)
   and the reflection (1): gcd of the determinants of several 24-product sets.
"""
import re
import json
import random
from math import gcd
from fractions import Fraction
import numpy as np
from leech_final import (Setup, reflection_map, perm_map, perm_inverse, gamma_of, twisted_constants,
                         triple_products, lambda_coords, f2_rank, int_det, dmul, in_Lambda_def,
                         Lcoords, from_Lcoords)

import os as _os
# the four independent implementations whose certificates and bases are re-checked here
REPRO = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "independent") + _os.sep
st = Setup()
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


R_WIT = (1, 1, 1, 1, 1, 1, 1, -1)          # doubled r
refl = reflection_map(R_WIT)


def twisted_triple_e(u, v, f, finv):
    """u *_g v from the definition, doubled e-coordinates (blocks must lie in L)."""
    def tw(a, b):
        return finv(dmul(f(a), f(b)))
    x, y, z = u[0:8], u[8:16], u[16:24]
    xp, yp, zp = v[0:8], v[8:16], v[16:24]
    add = lambda *vs: [sum(t) for t in zip(*vs)]
    P1 = add(tw(x, xp), tw(z, yp), tw(y, zp))
    P2 = add(tw(y, yp), tw(x, zp), tw(z, xp))
    P3 = add(tw(z, zp), tw(y, xp), tw(x, yp))
    return P1 + P2 + P3


def lam_coords_e(U):
    """Lambda-coordinates (this folder's basis) of a doubled e-coordinate vector."""
    P = np.array([sum((Lcoords(U[8 * k:8 * k + 8]) for k in range(3)), [])], dtype=np.int64)
    c, ok = lambda_coords(P, st.N)
    assert ok[0]
    return [int(t) for t in c[0]]


def gram_det(vecs):
    G = [[Fraction(sum(a * b for a, b in zip(u, v)), 4) for v in vecs] for u in vecs]
    assert all(g.denominator == 1 for row in G for g in row)
    return int_det([[int(g) for g in row] for row in G])


def factorize(n):
    """Prime factorization as a readable string (sympy.factorint, exact)."""
    from sympy import factorint
    f = factorint(abs(int(n)))
    return " * ".join(("%d^%d" % (p, e)) if e > 1 else str(p) for p, e in sorted(f.items()))


# ------------------------------------------------------------------ A. own certificate
say("A. own certificate for [Lambda : S] = 1 under s_r, r = (e0+...+e6-e7)/2")
G = gamma_of(refl)
Cg = twisted_constants(G)
P = triple_products(st.B, Cg)
coords, ok = lambda_coords(P, st.N)
assert ok.all()
say("   all 576 basis-pair products lie in Lambda; F2-rank =", f2_rank(coords))
rng = random.Random(1)
sets, dets, g_all = [], [], 0
while g_all != 1:
    order = list(range(576))
    rng.shuffle(order)
    chosen, rows = [], []
    for i in order:
        trial = rows + [coords[i]]
        if f2_rank(np.array(trial)) == len(trial):
            chosen.append(i)
            rows = trial
            if len(rows) == 24:
                break
    d = int_det(np.array(rows).tolist())
    assert d % 2 == 1
    sets.append(chosen)
    dets.append(d)
    g_all = gcd(g_all, d)
    say("   set %d: 24 products independent mod 2, det = %d = %s;  running gcd = %d"
        % (len(sets), d, factorize(d), g_all))
say("   => [Lambda : S] divides gcd = 1, so S = Lambda (no HNF used).")
Be = st.basis_e()
for k, (chosen, d) in enumerate(zip(sets, dets)):
    prods_e = []
    for i in chosen:
        p, q = divmod(i, 24)
        pe = twisted_triple_e(Be[p], Be[q], refl, refl)
        assert pe == sum((from_Lcoords(P[i][8 * b:8 * b + 8]) for b in range(3)), [])
        assert in_Lambda_def(pe)
        prods_e.append(pe)
    assert gram_det(prods_e) == 2 ** 24 * d * d
say("   cross-check: products recomputed from the definition in e-coordinates match, lie in Lambda,")
say("   and Gram det of each set = 2^24 * det^2 (basis-free form of the index).")
with open("cert_index1_final.gp", "w") as fh:
    fh.write("\\\\ Certificate for S = Lambda under s_r, r = (e0+...+e6-e7)/2 (f3_certificates.py).\n")
    fh.write("\\\\ All vectors in DOUBLED e-coordinates (blocks x|y|z). Pairs are 1-based (left, right).\n")
    fh.write("R2 = [%s];\n" % ",".join(map(str, R_WIT)))
    fh.write("BASIS2 = [" + ";".join(",".join(map(str, b)) for b in Be) + "];\n")
    for k, (chosen, d) in enumerate(zip(sets, dets)):
        pairs = [[i // 24 + 1, i % 24 + 1] for i in chosen]
        fh.write("PAIRS%d = [%s];\n" % (k + 1, ";".join("%d,%d" % tuple(pq) for pq in pairs)))
        fh.write("DET%d = %d;\n" % (k + 1, d))
    fh.write("NSETS = %d;\n" % len(sets))

# ------------------------------------------------------------------ B. refute certificate
say("")
say("B. refute agent's 24-product Z-basis certificate (refute/c3_certificate.gp.out)")
txt = open(REPRO + "refute/c3_certificate.gp.out").read()
lefts = re.findall(r"left\s*=\s*\[([^\]]*)\]", txt)
rights = re.findall(r"right\s*=\s*\[([^\]]*)\]", txt)
assert len(lefts) == 24 and len(rights) == 24


def dbl(s):
    v = [Fraction(t.strip()) for t in s.split(",")]
    assert len(v) == 24 and all((2 * t).denominator == 1 for t in v)
    return [int(2 * t) for t in v]


prods, cs = [], []
for a, b in zip(lefts, rights):
    u, v = dbl(a), dbl(b)
    assert in_Lambda_def(u) and in_Lambda_def(v), "factor not in Lambda"
    p = twisted_triple_e(u, v, refl, refl)
    assert in_Lambda_def(p)
    prods.append(p)
    cs.append(lam_coords_e(p))
d = int_det(cs)
gd = gram_det(prods)
say("   48 factors and 24 products satisfy Wilson's definition; det(Lambda-coords of products) = %d" % d)
say("   Gram det of the 24 products = %d (= 2^24: %s)" % (gd, gd == 2 ** 24))
assert abs(d) == 1 and gd == 2 ** 24
say("   => the 24 single products form a Z-basis of Lambda: S = Lambda.  VERIFIED.")
nrm = lambda U: sum(t * t for t in U) // 4
say("   norms (left, right, product):", [(nrm(dbl(a_)), nrm(dbl(b_)), nrm(p_)) for a_, b_, p_ in zip(lefts, rights, prods)])

# ------------------------------------------------------------------ C. GAP certificate
say("")
say("C. GAP agent's certificate (gap/cert_index1_reflection.txt)")
t = open(REPRO + "gap/cert_index1_reflection.txt").read()


def grab_matrix(name):
    m = re.search(name + r"\s*=\s*\[(.*?)\n\s*\]", t, re.S)
    rows = re.findall(r"\[([^\[\]]*)\]", m.group(1))
    return [[int(x) for x in r.split(",")] for r in rows]


BASIS2 = grab_matrix("BASIS2")
assert len(BASIS2) == 24 and all(in_Lambda_def(b) for b in BASIS2)
dB = int_det([lam_coords_e(b) for b in BASIS2])
say("   their BASIS2: 24 vectors in Lambda, det in this folder's Lambda-coordinates = %d (a Z-basis)" % dB)
assert abs(dB) == 1
for tag in ("A", "B"):
    pairs = grab_matrix("PAIRS_" + tag)
    prods_listed = grab_matrix("PRODUCTS_%s2" % tag)
    detlisted = int(re.search(r"DET_%s\s*=\s*(\d+)" % tag, t).group(1))
    cs = []
    for (a, c), pl in zip(pairs, prods_listed):
        p = twisted_triple_e(BASIS2[a - 1], BASIS2[c - 1], refl, refl)
        assert p == pl, "listed product differs from recomputation"
        assert in_Lambda_def(p)
        cs.append(lam_coords_e(p))
    dd = int_det(cs)
    say("   set %s: 24 products recomputed == listed, in Lambda; |det| = %d (listed %d) = %s"
        % (tag, abs(dd), detlisted, factorize(dd)))
    assert abs(dd) == detlisted
    if tag == "A":
        dA = dd
    else:
        dBset = dd
say("   gcd = %d  => S = Lambda.  VERIFIED." % gcd(dA, dBset))
assert gcd(dA, dBset) == 1

# ------------------------------------------------------------------ D. structure certificates
say("")
say("D. structure agent's certificates (structure/cert_*_products_ecoords_doubled.txt)")
SB = [[int(x) for x in line.split()] for line in open(REPRO + "structure/lambda_basis_ecoords_doubled.txt")
      if line.strip() and not line.startswith("#")]
assert len(SB) == 24 and all(in_Lambda_def(b) for b in SB)
assert abs(int_det([lam_coords_e(b) for b in SB])) == 1
say("   their basis: 24 vectors in Lambda forming a Z-basis (det +-1 in this folder's coordinates)")
twists = {"sigma": (perm_map((2, 1, 3, 4, 5, 6, 7)),) * 2,
          "c3": (perm_map((1, 2, 3, 4, 6, 7, 5)), perm_map(perm_inverse((1, 2, 3, 4, 6, 7, 5)))),
          "refl": (refl, refl)}
for name in ("sigma", "c3", "refl"):
    f, finv = twists[name]
    lines = open(REPRO + "structure/cert_%s_products_ecoords_doubled.txt" % name).read().splitlines()
    sets_d, cur = [], None
    for ln in lines:
        m = re.match(r"# set \d+: det of coordinate matrix = (-?\d+)", ln)
        if m:
            cur = {"det": int(m.group(1)), "rows": []}
            sets_d.append(cur)
        elif ln.strip() and not ln.startswith("#"):
            ij, vec = ln.split(":")
            i, j = map(int, ij.split())
            cur["rows"].append((i, j, [int(x) for x in vec.split()]))
    g = 0
    dl = []
    for s_ in sets_d:
        cs = []
        for i, j, pl in s_["rows"]:
            p = twisted_triple_e(SB[i], SB[j], f, finv)
            assert p == pl and in_Lambda_def(p)
            cs.append(lam_coords_e(p))
        dd = int_det(cs)
        assert abs(dd) == abs(s_["det"])
        dl.append(abs(dd))
        g = gcd(g, dd)
    say("   %-5s: %d sets recomputed == listed; |det| = %s; gcd = %d" % (name, len(sets_d), dl, g))
say("   (the gcd is an UPPER bound for the index; the lower bounds 2^16, 2^10, 1 come from the")
say("    F2-ranks 8, 14, 24 in f2_witnesses.py; together: 2^16, 2^10, 1 exactly)")

open("f3_certificates.out.txt", "w").write("\n".join(out) + "\n")
