"""Verify the Polarization paragraph of paper Appendix A.1: with
V = sigma(Ls_bar)/2L and W = sigma(Ls)/2L inside L_bar = L/2L,

  (1) V is a TWO-SIDED ideal of (L_bar, mu_bar): the right-ideal
      direction mu_bar(V, L_bar) <= V holds in addition to the
      left-ideal direction of Lemma 4.3;
  (2) both V and W are totally isotropic under the plus-type quadratic
      form q(v + 2L) = N(v)/2 mod 2, each of dimension 4 (the Witt
      index), and V (+) W = L_bar: a Witt decomposition of (L_bar, q)
      into a complementary pair of Lagrangians;
  (3) mu_bar is non-commutative on L_bar (the paper's footnote: L
      contains half-integer elements, for which the commutator need
      not lie in 2L).

Also certified along the way: q is the standard plus-type form on
F_2^8 (value distribution 136 zeros, 120 ones over the 256 cosets).

Machinery is reused unchanged from verify_mod2_quotient.py (octonion
product, L basis, sigma, mu_bar, F_2 spans); exact arithmetic
throughout.  Runtime: a few seconds.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fractions import Fraction as F

from verify_mod2_quotient import L, Lcoords, f2span, mubar, sLs, sLsbar


def f2_span_set(basis):
    out = {(0,) * 8}
    for b in basis:
        out |= {tuple(x ^ y for x, y in zip(t, b)) for t in out}
    return out


def q_form(coords_in_L):
    """q(v + 2L) = N(v)/2 mod 2, v reconstructed from L-basis F2-coords."""
    v = [F(0)] * 8
    for k in range(8):
        if coords_in_L[k]:
            for i in range(8):
                v[i] += L[k][i]
    n = sum(x * x for x in v)
    assert n.denominator == 1 and int(n) % 2 == 0, f"odd or fractional norm {n}"
    return (int(n) // 2) % 2


def main():
    print("=" * 70)
    print("Polarization of L/2L: V, W complementary Lagrangians (paper A.1)")
    print("=" * 70)

    V_basis = f2span([tuple(int(c) % 2 for c in Lcoords(v)) for v in sLsbar])
    W_basis = f2span([tuple(int(c) % 2 for c in Lcoords(v)) for v in sLs])
    V = f2_span_set(V_basis)
    W = f2_span_set(W_basis)
    e = [tuple(1 if k == i else 0 for k in range(8)) for i in range(8)]
    Lbar = f2_span_set(e)

    checks = []

    def check(name, ok):
        checks.append(ok)
        print(f"  {name}: {ok}")

    print(f"\ndim V = {len(V_basis)}, dim W = {len(W_basis)}, "
          f"|V| = {len(V)}, |W| = {len(W)}, |L/2L| = {len(Lbar)}")
    check("dim V = dim W = 4", len(V_basis) == 4 and len(W_basis) == 4)

    print("\n(1) V is a two-sided ideal of (L/2L, mu-bar):")
    check("left  ideal  mu-bar(L/2L, V) <= V",
          all(mubar(a, v) in V for a in Lbar for v in V))
    check("right ideal  mu-bar(V, L/2L) <= V",
          all(mubar(v, a) in V for v in V for a in Lbar))

    print("\n(2) Witt decomposition into complementary Lagrangians:")
    dist = {0: 0, 1: 0}
    for v in Lbar:
        dist[q_form(v)] += 1
    print(f"  q distribution on L/2L: q=0 -> {dist[0]}, q=1 -> {dist[1]}")
    check("q is plus-type (136 zeros, 120 ones)", dist == {0: 136, 1: 120})
    check("V totally isotropic (q = 0 on V)", all(q_form(v) == 0 for v in V))
    check("W totally isotropic (q = 0 on W)", all(q_form(w) == 0 for w in W))
    check("V and W complementary (V ^ W = 0, dims 4+4=8)",
          V & W == {(0,) * 8} and len(V_basis) + len(W_basis) == 8)

    print("\n(3) mu-bar is non-commutative on L/2L:")
    witness = None
    for a in Lbar:
        for b in Lbar:
            if mubar(a, b) != mubar(b, a):
                witness = (a, b)
                break
        if witness:
            break
    check("a commutator witness exists", witness is not None)
    if witness:
        a, b = witness
        print(f"  witness: a = {a}, b = {b}, "
              f"mu-bar(a,b) = {mubar(a, b)}, mu-bar(b,a) = {mubar(b, a)}")

    print()
    print("ALL PASS" if all(checks) else "FAIL -- INVESTIGATE")
    return 0 if all(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
