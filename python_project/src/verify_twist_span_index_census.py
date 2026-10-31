"""Span index of twisted triple products: the permutation census and the
index-1 reflection.

Created: 2026-10-14.  Updated: 2026-10-31.

The paper cites this script in Appendix E (label sec:census, "Other twists:
closure and span of the image"): at the census table, with
verify_all_permutations_exact.py, for the span-index column, and at the
reflection rho, which is s_r of Part 1 below.

Written in this repository to check statements of the external review of v8
by Claude Opus 5.5 (paper/reviews/2026-10-13_claude_opus_5_5_review_v8/), from
the paper's definitions with the exact lattice machinery of
verify_consecutive_twists_exact.py (Z-basis B of Lambda, membership tests, L^3
coordinates) and the octonion product of verify_ideal_decomposition.py.  No
code of the review is used.

For an isometry g of R^8 with g(L) = L, the twisted octonion product is
x ._g y := g^{-1}(g(x) . g(y)); the triple product star_g is Definition 3.3
of the paper with every block product ._g; S_g := Z-span{u star_g v : u, v in
Lambda}.  The span index [Lambda : S_g] is |det| of the integer Hermite
normal form of the 576 basis-pair products in L^3 coordinates, divided by
[L^3 : Lambda] = 4096 (own HNF, Python integers).

Part 1, the reflection s_r in the root r = (e0 + e1 + ... + e6 - e7)/2 of L:
  s_r(L) = L; L ._r L <= L on the 64 basis products (hence L^3 is a
  star_r-closed overlattice of Lambda, of index 2^12); all 576 basis-pair
  products lie in Lambda; [Lambda : S_r] = 1.
Part 2, all 5,040 permutation twists pi (pi(e_a) = e_{pi(a)}, pi(e_0) = e_0,
  x ._pi y = pi^{-1}(pi(x) . pi(y)), the paper's Section 5.6 convention):
  closure, the span index of each closing twist, and two statements of the
  review:
    Theorem A  star_pi closes on Lambda  <=>  L sbar and L s are both closed
               under ._pi  (checked on the 64 basis products of each), on
               every permutation;
    Theorem B  [Lambda : S_pi] = h_W . h_V^2 with h_W = [Ls : <Ls ._pi Ls>],
               h_V = [L sbar : <L sbar ._pi L sbar>] (8-dimensional HNFs), on
               every closing permutation.
  Theorem A is the paper's Proposition 5.3 (label prop:criterion), stated and
  proved there for every twist, an invertible linear map g with g(L) = L.
  Theorem B is not in the paper.
  Result (750 s on 12 cores): 1,764 closing permutations; all 1,176 even ones
  have index 2^10, all 588 odd ones (sigma among them) index 2^16; Theorem A
  holds on all 5,040 permutations, Theorem B on all 1,764 closing ones.

Independent of this file: the review's own implementations (Python, GAP,
PARI/GP) obtained the same numbers.

Usage:  python3 verify_twist_span_index_census.py            # both parts
        python3 verify_twist_span_index_census.py --reflection-only
"""
import itertools, sys, time
from collections import Counter
from fractions import Fraction as F
from multiprocessing import Pool

from verify_consecutive_twists_exact import B, inLambda, inL, l3, detB, Lb, s, sbar
from verify_ideal_decomposition import omul


def hnf_det(rows, n):
    """|det| of the Z-span of the integer row vectors (length n), by integer
    row reduction; 0 if the span has rank < n."""
    mat = [row[:] for row in rows]
    det = 1
    for col in range(n):
        piv = None
        while True:
            nz = [i for i in range(len(mat)) if mat[i][col] != 0]
            if not nz:
                return 0
            i0 = min(nz, key=lambda i: abs(mat[i][col]))
            done = True
            for i in nz:
                if i == i0:
                    continue
                q = mat[i][col] // mat[i0][col]
                mat[i] = [a - q * b for a, b in zip(mat[i], mat[i0])]
                if mat[i][col] != 0:
                    done = False
            if done:
                piv = i0
                break
        det *= abs(mat[piv][col])
        mat.pop(piv)
    return det


def coords(v, basis):
    """Exact coordinates of v in an 8-dimensional basis (Gaussian elimination)."""
    n = 8
    M = [[basis[j][i] for j in range(n)] + [v[i]] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]


add = lambda *vs: [sum(c) for c in zip(*vs)]
LS_BASIS = [omul(b, s) for b in Lb]          # left multiples l.s   (Wilson's Ls)
LSBAR_BASIS = [omul(b, sbar) for b in Lb]    # left multiples l.sbar (Wilson's L sbar)


def sublattice_index(basis_sub, tw):
    """[U : <U ._g U>] for the sublattice U spanned by basis_sub, or None if U is not closed."""
    prods = [coords(tw(a, b), basis_sub) for a in basis_sub for b in basis_sub]
    if any(c.denominator != 1 for v in prods for c in v):
        return None
    return hnf_det([[int(c) for c in v] for v in prods], 8)


def triple(tw):
    def star(u, v):
        (x, y, z), (x2, y2, z2) = u, v
        return (add(tw(x, x2), tw(z, y2), tw(y, z2)),
                add(tw(y, y2), tw(x, z2), tw(z, x2)),
                add(tw(z, z2), tw(y, x2), tw(x, y2)))
    return star


def closure_and_index(tw):
    star = triple(tw)
    prods = []
    for bi in B:
        for bj in B:
            p = star(bi, bj)
            if not inLambda(p):
                return False, None
            prods.append([int(c) for c in l3(p)])
    return True, hnf_det(prods, 24) // detB


# ------------------------------------------------------------------ part 1
def reflection_part():
    r = [F(1, 2)] * 7 + [F(-1, 2)]
    ip = lambda a, b: sum(x * y for x, y in zip(a, b))
    refl = lambda x: [xi - ip(x, r) * ri for xi, ri in zip(x, r)]     # N(r) = 2
    tw = lambda a, b: refl(omul(refl(a), refl(b)))                     # s_r is an involution
    print("Part 1: the reflection s_r, r = (e0+...+e6-e7)/2")
    print("  r in L:", inL(r), "| s_r(L) = L on the basis of L:", all(inL(refl(b)) for b in Lb))
    print("  L ._r L <= L (64 basis products):", all(inL(tw(a, b)) for a in Lb for b in Lb),
          " => L^3 is star_r-closed, [L^3 : Lambda] =", detB)
    closed, idx = closure_and_index(tw)
    print("  all 576 basis-pair products in Lambda:", closed, "| [Lambda : S_r] =", idx)
    hW, hV = sublattice_index(LS_BASIS, tw), sublattice_index(LSBAR_BASIS, tw)
    print("  Theorem A/B factors: h_W =", hW, " h_V =", hV, " h_W*h_V^2 =", None if None in (hW, hV) else hW * hV ** 2)
    return closed and idx == 1


# ------------------------------------------------------------------ part 2
def run_perm(perm):
    inv = [0] * 8
    for a, b in enumerate(perm):
        inv[b] = a
    def P(x):
        w = [F(0)] * 8
        for a in range(8):
            w[perm[a]] = x[a]
        return w
    def Pi(x):
        w = [F(0)] * 8
        for a in range(8):
            w[inv[a]] = x[a]
        return w
    tw = lambda a, b: Pi(omul(P(a), P(b)))
    hW, hV = sublattice_index(LS_BASIS, tw), sublattice_index(LSBAR_BASIS, tw)
    closed, idx = closure_and_index(tw)
    odd = sum(1 for i in range(8) for j in range(i + 1, 8) if perm[i] > perm[j]) % 2
    return perm, closed, idx, hW, hV, odd


def census_part(workers=12):
    perms = [(0,) + p for p in itertools.permutations(range(1, 8))]
    t0 = time.time()
    with Pool(workers) as pool:
        out = pool.map(run_perm, perms, chunksize=8)
    closing = [r for r in out if r[1] and r[0] != tuple(range(8))]
    print(f"Part 2: all {len(perms)} permutation twists ({time.time() - t0:.0f} s)")
    print("  closing (non-identity):", len(closing), "of", len(perms) - 1)
    by_parity = Counter(("odd" if r[5] else "even", r[2]) for r in closing)
    for (par, idx), n in sorted(by_parity.items()):
        print(f"  {par:4s} permutations with index {idx}: {n}")
    thA = all(r[1] == (r[3] is not None and r[4] is not None) for r in out)
    thB = all(r[2] == r[3] * r[4] ** 2 for r in closing)
    print("  Theorem A (closes <=> L sbar and L s both ._pi-closed), all 5,040:", thA)
    print("  Theorem B ([Lambda:S] = h_W * h_V^2), all closing:", thB)
    return len(closing) == 1764 and by_parity == Counter({("even", 1024): 1176, ("odd", 65536): 588}) and thA and thB


if __name__ == "__main__":
    ok = reflection_part()
    if "--reflection-only" not in sys.argv:
        ok = census_part() and ok
    print("ALL PASS" if ok else "FAILURE")
