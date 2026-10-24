"""The structural proof of Lemmas 4.2 to 4.4 (paper Section 4, label
sec:proof): the maximal order O := (1/2)(e1 + e2) L of Section 2.2, the two
lattice identities (4.1), the Moufang identity (4.2) for a purely imaginary
multiplier, and the three lemmas as their consequences.

Created: 2026-10-20.  Updated: 2026-10-24.

The paper cites this script in Section 2.2 at the maximal order O of (2.2)
(label eq:MO), with verify_sigma_sublattice_maximal_order.py; in Section 4.1
at the lattice identities (4.1) (label eq:multiples); and twice in Appendix C
(label sec:methodology, "Symbolic proof and numeric validation"): for the
finite content of the proof, and with symbolic_proof_checks.py and
gap_project/tests/test_lemmas.g for the lemmas on Z-bases.

Conventions are the paper's: the multiplication (2.1) (octonions.py), L = D_8^+
with s = (-e0 + e1 + ... + e7)/2, Wilson's sublattices L sbar = {l . sbar} and
L s = {l . s} (left multiples), sigma = (1 2).  Exact arithmetic throughout
(fractions.Fraction; sympy for the symbolic identity).

(A) O := (1/2) q L with q = e1 + e2, as the lattice spanned by (1/2) q b over a
    Z-basis b of L (left multiplication by q is injective and Z-linear).
    Checks: e0 in O; every basis unit e_i in O; q and p := e2 - e4 + e5 - e6
    in O, both purely imaginary; O closed under (2.1) on its 64 basis
    products; the trace form 2<,> integral on O with determinant 1 (O is
    unimodular, hence a maximal order); O = Z^8 + (1/2)(0 4)C_K = sigma(J_1),
    the descriptions of ledger claim 019 (verify_sigma_sublattice_maximal_order.py);
    and O = (1/2)(e0 + e4) B (e0 + e4) with B := (1/2) L sbar, Wilson's ring A_3
    (his A_t = (1/2)(1 + i_t) B (1 + i_t) with i_t = e_{t+1}; Wilson 2009, Sec. 2).
(B) L = q O: q O lies in L and has covolume 1 = covol(L).  (Immediate from
    the definition, since q^2 = -2 and L = -L; checked anyway.)
(C) sigma(L sbar) = 2 O and sigma(L s) = p O, each by eight products (the
    Z-basis of O) inside the target plus equal covolumes; and the covolume
    formula covol(x O) = N(x)^4 covol(O) = N(x)^4 / 16 for x = q, 2, p
    (left multiplication by x scales inner products by N(x)).
(D) The identity (x a)(x b) = x((a x) b) for purely imaginary x: symbolically
    in (2.1) (sympy, 23 symbols) and on 200 exact random rational triples;
    control: the defect (xa)(xb) - x((ax)b) equals 2 Re(x) [x, a, b], with
    [x, a, b] = (xa)b - x(ab), for every x (symbolically and on 200 random
    triples), so the identity fails for most x with a nonzero real part.  The
    laws the paper's proof uses (left and right alternative laws, flexibility,
    the left Moufang identity) are checked symbolically as well.
(E) The consequences, checked directly on basis pairs: (qO)(qO) in qO (Lemma
    4.2, L . L in L), L . 2O in 2O (Lemma 4.3), (pO)(pO) in pO (Lemma 4.4),
    and the three lemma inclusions on the paper's lattices L, sigma(L sbar),
    sigma(L s).
(F) The exclusion of Remark 4.5 (label rem:Ls-not-closed): sigma(L s) = p O
    is not an ideal of L on either side (some basis product of L . sigma(Ls)
    and of sigma(Ls) . L lies outside sigma(Ls)), while sigma(L sbar) = 2 O
    is a two-sided ideal of L.

Usage: python3 verify_structural_proof.py   (about 3 s)
"""
import random
from fractions import Fraction as F

import sympy as sp

from verify_ideal_decomposition import omul
from verify_consecutive_twists_exact import Lb, s, sbar, inL
import verify_sigma_sublattice_maximal_order as msb

E = [[F(1) if j == i else F(0) for j in range(8)] for i in range(8)]
ip = lambda a, b: sum(x * y for x, y in zip(a, b))
N = lambda a: ip(a, a)
sig = lambda x: [x[0], x[2], x[1]] + list(x[3:])
scale = lambda c, v: [c * x for x in v]
q = [F(0), F(1), F(1)] + [F(0)] * 5
p = [F(0)] * 8
p[2], p[4], p[5], p[6] = F(1), F(-1), F(1), F(-1)


def inside(vectors, basis):
    return all(msb.in_span(v, basis) for v in vectors)


def closed(basis):
    return all(msb.in_span(omul(a, b), basis) for a in basis for b in basis)


def main():
    ok = True
    # ---------------------------------------------------------------- (A)
    O = [scale(F(1, 2), omul(q, b)) for b in Lb]           # Z-basis of O = (1/2) q L
    units_in = all(msb.in_span(e, O) for e in E)
    q_in, p_in = msb.in_span(q, O), msb.in_span(p, O)
    imag = q[0] == 0 and p[0] == 0 and N(q) == 2 and N(p) == 4
    O_closed = closed(O)
    integral = all((2 * ip(a, b)).denominator == 1 for a in O for b in O) and all(N(a).denominator == 1 for a in O)
    unimod = integral and msb.covolume(O) ** 2 * 2 ** 8 == 1
    C_K = msb.span_code([[1] * 8] + [[1 if i in (0,) + tuple(l) else 0 for i in range(8)] for l in msb.FANO])
    swap = lambda c, n: tuple(c[n] if i == 0 else c[0] if i == n else c[i] for i in range(8))
    O4 = frozenset(swap(c, 4) for c in C_K)
    same_O4 = all(msb.in_ZC(v, O4) for v in O) and msb.covolume(O) == msb.covolume(msb.gens_ZC(O4))
    sC_K = frozenset(tuple(sig(list(c))) for c in C_K)        # sigma(J_1) = Z^8 + sigma(C_K)/2
    same_sJ1 = all(msb.in_ZC(v, sC_K) for v in O) and msb.covolume(O) == msb.covolume(msb.gens_ZC(sC_K))
    print("(A) O = (1/2)(e1+e2)L: e0 and all e_i in O:", units_in, "| q, p in O:", q_in and p_in,
          "| q, p purely imaginary of norms 2, 4:", imag)
    B = [scale(F(1, 2), omul(b, sbar)) for b in Lb]            # Wilson's B = (1/2) L sbar
    w = [a + b for a, b in zip(E[0], E[4])]                      # 1 + i_3 = e0 + e4
    A3 = [scale(F(1, 2), omul(omul(w, v), w)) for v in B]        # (w v) w = w (v w) by flexibility
    same_A3 = inside(A3, O) and msb.covolume(A3) == msb.covolume(O)
    print("    O closed on its 64 basis products:", O_closed, "| trace form integral with det 1 (maximal):", unimod,
          "| O = Z^8 + (0 4)C_K/2:", same_O4, "| O = sigma(J_1):", same_sJ1,
          "| O = Wilson's A_3 = (1/2)(e0+e4) B (e0+e4):", same_A3)
    ok &= units_in and q_in and p_in and imag and O_closed and unimod and same_O4 and same_sJ1 and same_A3
    # ---------------------------------------------------------------- (B)
    qO = [omul(q, b) for b in O]
    b_ok = all(inL(v) for v in qO) and msb.covolume(qO) == msb.covolume(Lb) == 1
    print("(B) L = q O: q O in L on the basis, covolume 1 = covol(L):", b_ok)
    ok &= b_ok
    # ---------------------------------------------------------------- (C)
    sLsbar = [sig(omul(b, sbar)) for b in Lb]
    sLs = [sig(omul(b, s)) for b in Lb]
    twoO = [scale(F(2), b) for b in O]
    pO = [omul(p, b) for b in O]
    c1 = inside(twoO, sLsbar) and msb.covolume(twoO) == msb.covolume(sLsbar) == 16
    c2 = inside(pO, sLs) and msb.covolume(pO) == msb.covolume(sLs) == 16
    formula = all(msb.covolume([omul(x, b) for b in O]) == N(x) ** 4 * msb.covolume(O) for x in (q, scale(F(2), E[0]), p)) \
        and msb.covolume(O) == F(1, 16)
    print("(C) sigma(L sbar) = 2 O (8 products + covolume 16):", c1, "| sigma(L s) = p O (8 products + covolume 16):", c2,
          "| covol(x O) = N(x)^4 / 16 for x = q, 2, p:", formula)
    ok &= c1 and c2 and formula
    # ---------------------------------------------------------------- (D)
    a = sp.symbols('a0:8'); b = sp.symbols('b0:8'); x = (sp.Integer(0),) + sp.symbols('x1:8')
    M = lambda u, v: [sp.expand(c) for c in omul(list(u), list(v))]
    sym_ok = all(sp.expand(l - r) == 0 for l, r in zip(M(M(x, a), M(x, b)), M(x, M(M(a, x), b))))
    rnd = random.Random(2026)
    rv = lambda imag_only: [F(0) if (imag_only and i == 0) else F(rnd.randint(-5, 5), rnd.randint(1, 3)) for i in range(8)]
    ex_ok = all(omul(omul(xx, aa), omul(xx, bb)) == omul(xx, omul(omul(aa, xx), bb))
                for xx, aa, bb in ((rv(True), rv(False), rv(False)) for _ in range(200)))
    xg = sp.symbols('x0:8')
    assoc_sym = [u - v for u, v in zip(M(M(xg, a), b), M(xg, M(a, b)))]
    defect_sym = all(sp.expand(l - r - 2 * xg[0] * c) == 0
                     for l, r, c in zip(M(M(xg, a), M(xg, b)), M(xg, M(M(a, xg), b)), assoc_sym))
    fails, defect_ok = 0, True
    for _ in range(200):
        xx, aa, bb = rv(False), rv(False), rv(False)
        lhs, rhs = omul(omul(xx, aa), omul(xx, bb)), omul(xx, omul(omul(aa, xx), bb))
        fails += lhs != rhs
        assoc = [u - v for u, v in zip(omul(omul(xx, aa), bb), omul(xx, omul(aa, bb)))]
        defect_ok &= [u - v for u, v in zip(lhs, rhs)] == scale(2 * xx[0], assoc)
    y = sp.symbols('y0:8'); z = sp.symbols('z0:8')
    laws = (all(sp.expand(l - r) == 0 for l, r in zip(M(M(xg, a), xg), M(xg, M(a, xg))))          # flexibility
            and all(sp.expand(l - r) == 0 for l, r in zip(M(M(a, xg), xg), M(a, M(xg, xg))))      # right alternative
            and all(sp.expand(l - r) == 0 for l, r in zip(M(xg, M(xg, b)), M(M(xg, xg), b)))      # left alternative
            and all(sp.expand(l - r) == 0 for l, r in zip(M(M(M(xg, y), xg), z), M(xg, M(y, M(xg, z))))))  # left Moufang
    print("(D) (xa)(xb) = x((ax)b) for purely imaginary x: symbolic:", sym_ok, "| 200 exact random triples:", ex_ok,
          "| defect = 2 Re(x)[x,a,b] for every x: symbolic:", defect_sym, "| on 200 random triples:", defect_ok,
          "(the identity fails for", fails, "of them)")
    print("    the laws used by the paper's proof hold symbolically (flexibility, both alternative laws, left Moufang):", laws)
    ok &= sym_ok and ex_ok and defect_sym and defect_ok and fails > 150 and laws
    # ---------------------------------------------------------------- (E)
    e1 = closed(qO) and all(msb.in_span(omul(u, v), twoO) for u in Lb for v in twoO) and closed(pO)
    e2 = all(inL(omul(u, v)) for u in Lb for v in Lb) \
        and all(msb.in_span(omul(u, v), sLsbar) for u in Lb for v in sLsbar) \
        and all(msb.in_span(omul(u, v), sLs) for u in sLs for v in sLs)
    print("(E) (qO)(qO) in qO, L.2O in 2O, (pO)(pO) in pO:", e1, "| Lemmas 4.2, 4.3, 4.4 on the paper's lattices:", e2)
    ok &= e1 and e2
    # ---------------------------------------------------------------- (F)
    left_out = sum(not msb.in_span(omul(u, v), sLs) for u in Lb for v in sLs)
    right_out = sum(not msb.in_span(omul(v, u), sLs) for u in Lb for v in sLs)
    ideal = all(msb.in_span(omul(u, v), sLsbar) and msb.in_span(omul(v, u), sLsbar) for u in Lb for v in sLsbar)
    print("(F) sigma(Ls) = pO not an ideal of L: products outside, left", left_out, "right", right_out,
          "(of 64 each; the counts depend on the bases) | sigma(L sbar) = 2O a two-sided ideal of L:", ideal)
    ok &= left_out > 0 and right_out > 0 and ideal
    print("ALL PASS" if ok else "FAILURE")


if __name__ == "__main__":
    main()
