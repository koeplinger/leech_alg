#!/usr/bin/env python3
"""
c3_certificate_check.py -- independent Python check of the 24-product
certificate printed by c3_certificate.gp (PARI/GP).

Reads the 24 (left factor, right factor) pairs from c3_certificate.gp.out,
and with the octonion code of leech_twist_lib.py (no PARI involved):
  1. checks every factor against Wilson's three conditions (so u_k, v_k in Lambda),
  2. computes p_k = u_k *_g v_k for g = s_r, r = (e0+...+e6-e7)/2,
  3. checks every p_k against Wilson's conditions,
  4. computes the Gram determinant of p_1..p_24 exactly.
Since det Gram(Lambda) = 2^24 and the p_k lie in Lambda, Gram det = 2^24
means [Lambda : Z<p_1..p_24>]^2 = 1, i.e. the 24 products are a Z-basis of
Lambda, so S_g = Lambda.  This check does not use any basis of Lambda.

Run:   python3 c3_certificate_check.py [c3_certificate.gp.out] > c3_certificate_check.out
Conventions: main.tex v8, as in leech_twist_lib.py.
"""
import re
import sys
from fractions import Fraction as Fr

import leech_twist_lib as T

path = sys.argv[1] if len(sys.argv) > 1 else "c3_certificate.gp.out"
text = open(path).read()
vecs = re.findall(r"(left|right) += \[([^\]]*)\]", text)
assert len(vecs) == 48, f"expected 48 factor lines, found {len(vecs)}"


def doubled(s):
    return [int(Fr(t.strip()) * 2) for t in s.split(",")]


pairs = []
for k in range(24):
    (s1, a), (s2, b) = vecs[2 * k], vecs[2 * k + 1]
    assert s1 == "left" and s2 == "right"
    pairs.append((doubled(a), doubled(b)))

R = T.Twist(T.reflection([1, 1, 1, 1, 1, 1, 1, -1]), "s_r")
ok_f = all(T.in_Lambda2(u) and T.in_Lambda2(v) for u, v in pairs)
prods = [T.tprod2(R, u, v) for u, v in pairs]
ok_p = all(T.in_Lambda2(p) for p in prods)
gram4 = [[sum(a * b for a, b in zip(p, q)) for q in prods] for p in prods]   # 4 * Gram
gdet = Fr(T.det_bareiss(gram4), 4 ** 24)
print("factors read:", 2 * len(pairs))
print("all 48 factors satisfy Wilson's conditions:", ok_f)
print("all 24 products u_k *_g v_k satisfy Wilson's conditions:", ok_p)
print("Gram determinant of the 24 products:", gdet, "(det Lambda = 2^24 =", 2 ** 24, ")")
print("=> the 24 products form a Z-basis of Lambda:", ok_f and ok_p and gdet == 2 ** 24)
print("norms of the products:", sorted(int(T.nrm2(p)) for p in prods))
