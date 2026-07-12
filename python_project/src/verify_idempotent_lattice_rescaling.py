"""Task A, part 4: how the eight idempotents of (R^24, +, star) sit relative
to the Leech lattice Lambda.

Uses the complete classification of verify_idempotent_classification.py
(exactly 8 idempotents: 0, eps_1, eps_2, eps_3, w/3, and eps_i - w/3) and the
exact Lambda membership test of verify_consecutive_twists_exact.py.

Facts established here (all exact):

 (A) e_0 is not in L (its coordinate sum is odd), so none of the seven nonzero
     idempotents lies in Lambda: every one of them has some block equal to a
     nonzero rational multiple of e_0 with a denominator or an odd numerator
     that L forbids.  Since the classification is COMPLETE, this PROVES
        Lambda contains no nonzero idempotent.
     (Not merely "the minimal shell contains none".)

 (B) RESCALING.  Each idempotent still spans a rational ray that MEETS
     Lambda.  The least positive integer n with n E in Lambda is
        eps_i        : n = 4    -> (4 e_0, 0, 0),          N = 16
        eps_i - w/3  : n = 6    -> (4 e_0, -2 e_0, -2 e_0), N = 24
        w/3          : n = 12   -> (4 e_0, 4 e_0, 4 e_0),   N = 48
     computed by scanning c over a rational grid and testing inLambda.

 (C) The scalar-square equation on the lattice.  For u in Lambda and a scalar
     c, u star u = c u with c =/= 0 holds IFF u = c E for one of the seven
     nonzero idempotents E (because then u/c is idempotent).  Combining with
     (B), the complete solution set in Lambda is
        u = n eps_i        (4 | n),   c = n;
        u = n (eps_i - w/3)(6 | n),   c = n;
        u = n (w/3)       (12 | n),   c = n.
     The smallest such u is (4 e_0, 0, 0) and its block permutations, of norm
     16, with c = 4 (and (-4 e_0, 0, 0) with c = -4).
     NO minimal vector of Lambda satisfies u star u = c u for ANY c: verified
     exhaustively over all 196,560 minimal vectors, and proved -- c = 0 would
     make u square-zero (impossible on the minimal shell, see
     verify_square_zero_classification.py: square-zero norms lie in 12Z),
     while c =/= 0 forces N(u) = c^2 N(E) in {c^2, 2c^2/3, c^2/3} = 8, i.e.
     c^2 in {8, 12, 24}, which has no rational root.

 (D) The orchestrator's witness, verified: for lambda a PURELY IMAGINARY root
     of L (Re(lambda) = 0; there are 84 of them),
        u = (2 lambda, 0, 0) in Min(Lambda)  and  u star u = (-8 e_0, 0, 0)
                                                           = -8 eps_1,
     and (-8 e_0, 0, 0) IS in Lambda.  For a root lambda with Re(lambda) =/= 0
     one gets instead u star u = (8 Re(lambda) lambda - 8 e_0, 0, 0), which is
     not a multiple of eps_1; so purely imaginary is exactly the right
     hypothesis.  Note -8 eps_1 is NOT the smallest lattice point on the ray:
     -4 eps_1 = (-4 e_0, 0, 0) is (item B).
"""
from fractions import Fraction as F
from itertools import combinations, product as iproduct

from verify_ideal_decomposition import (
    omul, odot, star, star_indep, add, zero8, E, eq, is_zero, Sigma,
)
from verify_consecutive_twists_exact import inLambda, inL, inLs, inLsbar, s, sbar
from verify_idempotent_classification import EIGHT, EPS1, EPS2, EPS3, OM3, P1, P2, P3


def N(v):   return sum(c*c for c in v)
def Nu(u):  return sum(N(b) for b in u)
def ip(u, v): return sum(a*b for a, b in zip(u, v))
def smul(c, u): return tuple([c*t for t in blk] for blk in u)


def main():
    h = F(1, 2)

    # ------------------------------------------------------------------ (A)
    print("=== (A) no nonzero idempotent lies in Lambda ===")
    e0 = E(0)
    print(f"  e_0 in L ? {inL(e0)}   (coordinate sum 1 is odd)")
    inl = []
    for name, u in EIGHT:
        inlam = inLambda(u)
        inl.append(inlam)
        print(f"  {name:<14} in Lambda : {inlam}")
    only_zero = (inl[0] is True) and not any(inl[1:])
    print(f"  => the ONLY idempotent in Lambda is 0 : {only_zero}")
    print( "  This is a PROOF (not a shell check): the 8 above are ALL the")
    print( "  idempotents of the algebra (verify_idempotent_classification.py).")
    print()

    # ------------------------------------------------------------------ (B)
    print("=== (B) least positive multiple of each idempotent lying in Lambda ===")
    least = {}
    for name, u in EIGHT[1:]:
        found = None
        for k in range(1, 241):                 # c = k/2, k = 1..240
            c = F(k, 2)
            if inLambda(smul(c, u)):
                found = c; break
        least[name] = found
        v = smul(found, u)
        blocks = ", ".join(
            ("0" if all(t == 0 for t in b)
             else f"{[str(t) for t in b if t != 0][0]}*e_0") for b in v)
        print(f"  {name:<14} least c > 0 with c*E in Lambda : c = {found}"
              f"   -> ({blocks}),  N = {Nu(v)}")
    okB = (least["eps_1"] == 4 and least["eps_2"] == 4 and least["eps_3"] == 4
           and least["eps_1 - w/3"] == 6 and least["eps_2 - w/3"] == 6
           and least["eps_3 - w/3"] == 6 and least["w/3"] == 12)
    print(f"  matches 4 / 6 / 12 : {okB}")
    # and the full set of c on the grid
    for name, u in EIGHT[1:]:
        good = [F(k, 2) for k in range(1, 121) if inLambda(smul(F(k, 2), u))]
        print(f"  {name:<14} all c in (0, 60] on the (1/2)Z grid with c*E in"
              f" Lambda : {[str(c) for c in good][:8]}...")
    print()

    # ------------------------------------------------------------------ (C)
    print("=== (C) u star u = c u  on Lambda ===")
    # the lattice points c*E really do satisfy u*u = c*u
    for name, u in EIGHT[1:]:
        c = least[name]
        v = smul(c, u)
        want = smul(c, v)
        good = eq(star(v, v), want) and eq(star_indep(v, v), want)
        print(f"  u = {c}*{name:<14}: u*u = {c}*u : {good}   (N(u) = {Nu(v)})")
    print()

    print("  exhaustive scan of the 196,560 minimal vectors for u*u = c*u ...")
    from verify_idempotents_min_shell import enumerate_min_shell
    from verify_section5_properties import star as star_i, eq as eq_i, Nsq
    vecs, counts = enumerate_min_shell()
    assert len(vecs) == 196560
    hits = []
    hist = {}
    for U in vecs:
        S2 = star_i(U, U)
        n = Nsq(S2) // 4
        hist[n] = hist.get(n, 0) + 1
        # is S2 a scalar multiple of U?  (doubled-int coords, same scale)
        flat_u = [t for b in U for t in b]
        flat_s = [t for b in S2 for t in b]
        k = next(i for i, t in enumerate(flat_u) if t != 0)
        cnum, cden = flat_s[k], flat_u[k]
        if all(a*cden == cnum*b for a, b in zip(flat_s, flat_u)):
            hits.append((U, F(cnum, cden)))
    print(f"  minimal vectors u with u*u = c*u for some scalar c : {len(hits)}")
    print(f"  histogram of N(u*u) over Min(Lambda): "
          f"{ {k: hist[k] for k in sorted(hist)} }")
    print(f"  min N(u*u) = {min(hist)} > 0  => no square-zero minimal vector")
    print( "  PROOF that the scan must come out empty:")
    print( "    c = 0 -> u square-zero -> N(u) in 12Z (square-zero script), but")
    print( "            N(u) = 8 for minimal vectors.  Impossible.")
    print( "    c =/= 0 -> u/c is an idempotent E =/= 0, so u = c E and")
    print( "            8 = N(u) = c^2 N(E) with N(E) in {1, 2/3, 1/3},")
    print( "            i.e. c^2 in {8, 12, 24} -- no rational root.  Impossible.")
    for nm, EE in [("eps_i", EPS1), ("eps_i - w/3", P1), ("w/3", OM3)]:
        print(f"    N({nm}) = {Nu(EE)}")
    print()

    # ------------------------------------------------------------------ (D)
    print("=== (D) the (2 lambda, 0, 0) witness ===")
    roots = []
    for i, j in combinations(range(8), 2):
        for si in (F(1), F(-1)):
            for sj in (F(1), F(-1)):
                v = [F(0)]*8; v[i] = si; v[j] = sj; roots.append(v)
    for sg in iproduct([h, -h], repeat=8):
        if sum(1 for t in sg if t < 0) % 2 == 1:
            roots.append(list(sg))
    assert len(roots) == 240
    im_roots = [lam for lam in roots if lam[0] == 0]     # Re(lambda) = 0
    print(f"  E_8 roots with Re = 0 (purely imaginary): {len(im_roots)}")

    target = smul(F(-8), EPS1)                            # (-8 e_0, 0, 0)
    ok_all = True
    for lam in im_roots:
        u = ([2*t for t in lam], zero8(), zero8())
        ok_all &= (Nu(u) == 8) and inLambda(u)
        ok_all &= eq(star(u, u), target) and eq(star_indep(u, u), target)
    print(f"  for ALL {len(im_roots)}: (2 lambda,0,0) in Min(Lambda) and")
    print(f"     (2 lambda,0,0)^2 = (-8 e_0, 0, 0) = -8 eps_1 : {ok_all}")
    print(f"  (-8 e_0, 0, 0) in Lambda : {inLambda(target)}   N = {Nu(target)}")
    print(f"  (-4 e_0, 0, 0) in Lambda : {inLambda(smul(F(-4), EPS1))}"
          f"   N = {Nu(smul(F(-4), EPS1))}   <- the smallest point on the ray")

    # sharpness: a root with Re =/= 0 does NOT give a multiple of eps_1
    prop = 0
    for lam in roots:
        u = ([2*t for t in lam], zero8(), zero8())
        S2 = star(u, u)
        blk = S2[0]
        if all(t == 0 for t in blk[1:]) and is_zero((S2[1], S2[2])):
            prop += 1
    print(f"  roots lambda for which (2 lambda,0,0)^2 is a multiple of eps_1 :"
          f" {prop}  (= the {len(im_roots)} purely imaginary ones: "
          f"{prop == len(im_roots)})")
    print()

    okall = (only_zero and okB and len(hits) == 0 and ok_all
             and prop == len(im_roots) and inLambda(target))
    print("ALL PASS" if okall else "FAILURE")


if __name__ == "__main__":
    main()
