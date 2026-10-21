"""The witnesses of paper v8, Section 5.1 (21 October 2026, prompt 318; item 1
of evidence_and_reasoning/2026-10-20_closed_form_plan.md): every standard
identity that fails on (Lambda, +, star) has a witness that can be checked by
hand.

Conventions are the paper's: the multiplication (2.1), sigma = (1 2),
x ._sigma y := sigma(sigma(x) . sigma(y)), and the triple product of
Definition 3.3,
    P1 = x.x' + z.y' + y.z',  P2 = y.y' + x.z' + z.x',  P3 = z.z' + y.x' + x.y',
with ._sigma in all nine products; eps_1 = (e0, 0, 0), eps_2 = (0, e0, 0),
eps_3 = (0, 0, e0) and omega = eps_1 + eps_2 + eps_3.  Exact arithmetic
throughout (fractions.Fraction; sympy for one rank computation).

(A) The routing rule on single blocks: for basis units e_i in block alpha and
    e_j in block beta, the product is e_i ._sigma e_j, in block alpha if
    alpha = beta and in the third block otherwise.  Checked on all 576 basis
    pairs, hence for all vectors with a single nonzero block (bilinearity).
(B) sigma fixes e0, e3, e4, e6; their span is closed under (2.1), and
    ._sigma equals the standard product on it (16 basis pairs).
(C) eps_a star eps_a = eps_a, and eps_a star eps_b = eps_c for {a, b, c} =
    {1, 2, 3}.
(D) The table of Section 5.1: for each identity, both sides at the witness,
    compared with the values the paper prints; they differ.  For the quartic
    identity all five parenthesizations of u star u star u star u are
    computed at u = eps_1 - eps_2 (the paper shows ((uu)u)u and (uu)(uu));
    u star u = eps_1 + eps_2 - 2 eps_3, (uu)u = 3u = u(uu), so the cube
    identity holds at this u.
(E) The scaled witnesses: four times each witness vector lies in Lambda
    (Wilson's three conditions, inLambda of verify_consecutive_twists_exact.py),
    and at four times the witness both sides of the identity are multiplied
    by 4 to the power of its degree, so the identity fails there too.
(F) Remark 5.1: u = (2 lam, 0, 0) with lam = e3 + e4 lies in Lambda, N(u) = 8,
    u star u = (-8 e0, 0, 0), of norm 64.
(G) Remark 5.2, cross-check: (R^24, star) has neither a left nor a right
    identity (the 576 linear equations e star b_k = b_k, resp. b_k star e =
    b_k, over the standard basis b_k are inconsistent).

Note (21 October 2026, prompt 320).  The paper's table in Section 5.1 now
lists only the identities and their witnesses, written as explicit vectors
(eps_1 = (e0, 0, 0), omega = (e0, e0, e0), eps_1 - eps_2 = (e0, -e0, 0));
it no longer prints the two sides, and the paper no longer states the
scaling into Lambda.  The witnesses are unchanged.  Item (D) checks the two
sides as given in the table of prompt 319
(evidence_and_reasoning/2026-10-21_plan.md, section 2.1), and item (E) the
scaling stated there; both checks stand as a record.

Usage: python3 verify_identity_witnesses.py   (a few seconds)
"""
from fractions import Fraction as F

import sympy as sp

from verify_ideal_decomposition import omul
from verify_consecutive_twists_exact import inLambda

ZERO = [F(0)] * 8


def unit(i, c=1):
    v = [F(0)] * 8
    v[i] = F(c)
    return v


def sig(x):
    return [x[0], x[2], x[1]] + list(x[3:])


def mul_s(x, y):                                   # x ._sigma y
    return sig(omul(sig(x), sig(y)))


def add(*vs):
    return [sum(t) for t in zip(*vs)]


def star(u, v):
    x, y, z = u
    xp, yp, zp = v
    m = mul_s
    return (add(m(x, xp), m(z, yp), m(y, zp)),
            add(m(y, yp), m(x, zp), m(z, xp)),
            add(m(z, zp), m(y, xp), m(x, yp)))


def vsum(*us):
    return tuple(add(*blocks) for blocks in zip(*us))


def smul(c, u):
    return tuple([c * t for t in b] for b in u)


def norm(u):
    return sum(t * t for b in u for t in b)


def inner(u, v):
    return sum(p * q for bu, bv in zip(u, v) for p, q in zip(bu, bv))


def blocks(a=None, b=None, c=None):
    return (a or ZERO, b or ZERO, c or ZERO)


EPS = [blocks(unit(0)), blocks(None, unit(0)), blocks(None, None, unit(0))]
OMEGA = vsum(*EPS)


def third(a, b):                                   # blocks indexed 0, 1, 2
    return a if a == b else 3 - a - b


def main():
    ok = True
    # ---------------------------------------------------------------- (A)
    a_ok = True
    for al in range(3):
        for be in range(3):
            for i in range(8):
                for j in range(8):
                    u = [ZERO, ZERO, ZERO]; u[al] = unit(i)
                    v = [ZERO, ZERO, ZERO]; v[be] = unit(j)
                    w = [ZERO, ZERO, ZERO]; w[third(al, be)] = mul_s(unit(i), unit(j))
                    a_ok &= star(tuple(u), tuple(v)) == tuple(w)
    print("(A) routing rule on single blocks, all 576 basis pairs:", a_ok)
    ok &= a_ok
    # ---------------------------------------------------------------- (B)
    H = (0, 3, 4, 6)
    fixed = all(sig(unit(i)) == unit(i) for i in H)
    span_closed = all(all(omul(unit(i), unit(j))[k] == 0 for k in range(8) if k not in H) for i in H for j in H)
    same = all(mul_s(unit(i), unit(j)) == omul(unit(i), unit(j)) for i in H for j in H)
    b_ok = fixed and span_closed and same
    print("(B) sigma fixes e0, e3, e4, e6:", fixed, "| their span is closed:", span_closed,
          "| ._sigma = standard product on it:", same)
    ok &= b_ok
    # ---------------------------------------------------------------- (C)
    c_ok = all(star(EPS[a], EPS[b]) == EPS[third(a, b)] for a in range(3) for b in range(3))
    print("(C) eps_a star eps_b = eps_(a o b):", c_ok)
    ok &= c_ok
    # ---------------------------------------------------------------- (D)
    e3, e4 = blocks(unit(3)), blocks(unit(4))
    x = vsum(EPS[0], smul(-1, EPS[1]))             # eps_1 - eps_2
    xx = star(x, x)
    rows = []
    # each row: (identity, both sides at the witness equal the values in the paper's table, and differ)
    l, r = star(star(EPS[0], EPS[1]), EPS[2]), star(EPS[0], star(EPS[1], EPS[2]))
    rows.append(("associativity", l == EPS[2] and r == EPS[0] and l != r))
    l, r = star(star(EPS[0], EPS[0]), EPS[1]), star(EPS[0], star(EPS[0], EPS[1]))
    rows.append(("left alternativity", l == EPS[2] and r == EPS[1] and l != r))
    l, r = star(star(EPS[1], EPS[0]), EPS[0]), star(EPS[1], star(EPS[0], EPS[0]))
    rows.append(("right alternativity", l == EPS[1] and r == EPS[2] and l != r))
    u = blocks(unit(3), unit(4))
    l, r = star(star(u, EPS[0]), u), star(u, star(EPS[0], u))
    rows.append(("flexibility", vsum(l, smul(-1, r)) == blocks(None, unit(6, -2), unit(6, 2))
                 and l == blocks(unit(0, -2), unit(6, -1), unit(6, 1))
                 and r == blocks(unit(0, -2), unit(6, 1), unit(6, -1))))
    u = blocks(unit(0), unit(3), unit(4))
    l, r = star(star(u, u), u), star(u, star(u, u))
    rows.append(("cube power-associativity", vsum(l, smul(-1, r)) == blocks(None, unit(6, -4), unit(6, 4))))
    p5 = [star(star(xx, x), x), star(star(x, xx), x), star(x, star(xx, x)), star(x, star(x, xx)), star(xx, xx)]
    quart = (xx == vsum(EPS[0], EPS[1], smul(-2, EPS[2])) and star(xx, x) == smul(3, x) == star(x, xx)
             and p5[0] == smul(3, xx) and p5[4] == smul(-3, xx) and len({str(t) for t in p5}) > 1)
    rows.append(("quartic power-associativity (cube holds at u)", quart))
    l, r = star(e3, e4), star(e4, e3)
    rows.append(("commutativity", l == blocks(unit(6)) and r == blocks(unit(6, -1))))
    rows.append(("norm multiplicativity", norm(star(OMEGA, x)) == 0 and norm(OMEGA) * norm(x) == 6))
    a = e3
    rows.append(("symmetric composition", inner(star(EPS[0], a), a) == 1 and inner(EPS[0], star(a, a)) == -1))
    d_ok = all(t for _, t in rows)
    for name, t in rows:
        print("(D)", name + ":", t)
    ok &= d_ok
    # ---------------------------------------------------------------- (E)
    four = lambda u: smul(4, u)
    vecs = [EPS[0], EPS[1], EPS[2], OMEGA, x, e3, e4, blocks(unit(3), unit(4)), blocks(unit(0), unit(3), unit(4))]
    in_lam = all(inLambda([list(b) for b in four(v)]) for v in vecs)
    # each row: (both sides of the identity as a function of its arguments, the witness, the degree)
    rows_e = [
        (lambda u, v, w: (star(star(u, v), w), star(u, star(v, w))), (EPS[0], EPS[1], EPS[2]), 3),
        (lambda u, v: (star(star(u, u), v), star(u, star(u, v))), (EPS[0], EPS[1]), 3),
        (lambda u, v: (star(star(u, v), v), star(u, star(v, v))), (EPS[1], EPS[0]), 3),
        (lambda u, v: (star(star(u, v), u), star(u, star(v, u))), (blocks(unit(3), unit(4)), EPS[0]), 3),
        (lambda u: (star(star(u, u), u), star(u, star(u, u))), (blocks(unit(0), unit(3), unit(4)),), 3),
        (lambda u: (star(star(star(u, u), u), u), star(star(u, u), star(u, u))), (x,), 4),
        (lambda u, v: (star(u, v), star(v, u)), (e3, e4), 2),
        (lambda u, v: (norm(star(u, v)), norm(u) * norm(v)), (OMEGA, x), 4),
        (lambda u, v, w: (inner(star(u, v), w), inner(u, star(v, w))), (EPS[0], a, a), 3),
    ]
    times = lambda c, t: smul(c, t) if isinstance(t, tuple) else c * t
    scaled = []
    for sides, wit, deg in rows_e:
        (l1, r1), (l4, r4) = sides(*wit), sides(*[four(t) for t in wit])
        scaled.append(l4 == times(4 ** deg, l1) and r4 == times(4 ** deg, r1) and l4 != r4)
    e_ok = in_lam and all(scaled)
    print("(E) four times each witness vector lies in Lambda:", in_lam,
          "| at four times the witness both sides scale by 4^degree and still differ:", all(scaled))
    ok &= e_ok
    # ---------------------------------------------------------------- (F)
    lam = add(unit(3), unit(4))
    u = blocks([2 * t for t in lam])
    f_ok = inLambda([list(b) for b in u]) and norm(u) == 8 and star(u, u) == blocks(unit(0, -8)) and norm(star(u, u)) == 64
    print("(F) Remark 5.1: u = (2 lam, 0, 0), lam = e3 + e4, in Lambda, N(u) = 8, u star u = (-8 e0, 0, 0), norm 64:", f_ok)
    ok &= f_ok
    # ---------------------------------------------------------------- (G)
    basis = []
    for bl in range(3):
        for i in range(8):
            v = [ZERO, ZERO, ZERO]; v[bl] = unit(i); basis.append(tuple(v))
    flat = lambda u: [t for b in u for t in b]
    def system(left):
        rows_, rhs = [], []
        for bk in basis:
            cols = [flat(star(bj, bk) if left else star(bk, bj)) for bj in basis]   # e = sum c_j b_j
            for comp in range(24):
                rows_.append([cols[j][comp] for j in range(24)])
                rhs.append(flat(bk)[comp])
        A = sp.Matrix(rows_)
        Ab = A.row_join(sp.Matrix(rhs))
        return A.rank() < Ab.rank()                 # inconsistent: no identity
    no_left, no_right = system(True), system(False)
    print("(G) Remark 5.2: no left identity:", no_left, "| no right identity:", no_right)
    ok &= no_left and no_right
    print("ALL PASS" if ok else "SOME CHECK FAILED")
    return ok


if __name__ == "__main__":
    main()
