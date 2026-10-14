r"""
06_certificates.py -- HNF-free, numpy-free certificates of closure and of the
span index for
  * sigma = (1 2)                                            (claimed 2^16)
  * the C2 witness pi = (5 6 7): pi(5)=6, pi(6)=7, pi(7)=5    (claimed 2^10)
  * the inverse permutation (5 7 6)                          (expected: no closure)
  * the C3 witness, reflection in r = (e0+e1+...+e6-e7)/2     (claimed 1)
  * every other root reflection whose census index (03) is 1 (compact)

Everything is computed with fractions.Fraction directly from the definitions.
Shared with the census (03) are only the octonion table, the definitional
membership tests and the Lambda basis, which is read back from
lambda_basis.txt and re-verified here.

For each twist g, with A := Lambda / S_g (a finite abelian group):
 1. The 576 products b_i *_g b_j, with x ._g y = g^{-1}(g(x).g(y)) evaluated
    literally (g = coordinate permutation pi(e_a) = e_{pi(a)}, or the
    reflection x -> x - <x,r> r), in e-coordinates.
 2. CLOSURE: every product passes the definitional test (blocks in L by the
    coordinate rule; pairwise block sums a in L.sbar via (a.s)/2 in L; block
    sum a in L.s via (a.sbar)/2 in L).  Then the Lambda-coordinates of every
    product (exact rational Gauss-Jordan) are checked to be integers.
 3. |A/2A| = 2^(24 - r2), r2 = rank over F2 of the 576 coordinate rows.
 4. |A/4A| = 4^24 / |(S + 4 Lambda)/4 Lambda|, the order of the Z/4-span of
    the rows (leechlib.log2_order_mod4_span; no HNF).  If |A/4A| = |A/2A|,
    the 2-part of A is elementary abelian of order |A/2A|.
 5. ODD PART: pick 24 products T that are linearly independent (chosen
    F2-independent first, to keep det small in 2-adic terms); |A| divides
    D = |det T| because Z<T> <= S <= Lambda.  For every odd prime p | D,
    rank over F_p of the 576 rows = 24 means A has no p-torsion.  So the odd
    part of |A| is 1.
 => A = (Z/2)^(24 - r2),  [Lambda : S_g] = 2^(24 - r2)  (and 2 Lambda <= S_g).

A second, conceptual certificate (mod-2 algebra Abar = L/2L; proofs in the
report).  With V = g(L.sbar)/2L, W = g(L.s)/2L and, for X in {V, W},
   d(X) = 4 - dim span(XX),   m(X) = 8 - dim span(X + Abar X + X Abar),
 (T2) span(L.L) = L (checked here), hence (4L)^3 <= S_g: Lambda/S_g is a
      2-group (exponent dividing 4);
 (T3) [Lambda : S_g + (2L)^3] = 2^(d(W) + 2 d(V));
 (T4) [Lambda : S_g + 2 Lambda] <= 2^(d(W) + 2 d(V) + m(W) + 2 m(V)).
So d = m = 0 for V and W gives S_g + 2 Lambda = Lambda (odd index) and, with
(T2), S_g = Lambda.  These F2 invariants are recomputed here with Fractions.

Run:  python3 06_certificates.py > 06_certificates.out.txt
"""
import json
import random
from fractions import Fraction as F
from itertools import product as iproduct

import sympy
import leechlib as LL

P_MOD = 2 ** 61 - 1          # a prime; used only to choose Q-independent rows

# ---------------------------------------------------------------------------
# Lambda basis from lambda_basis.txt (doubled e-coordinates), re-verified
# ---------------------------------------------------------------------------
B = []
for line in open('lambda_basis.txt'):
    if line.startswith('#'):
        continue
    B.append([F(int(a), 2) for a in line.split('|L3|')[0].split()])
assert len(B) == 24 and all(LL.in_Lambda_def(b) for b in B)
GRAM = [[int(LL.ip(a, b)) for b in B] for a in B]
assert LL.det_bareiss(GRAM) == 2 ** 24
print("Lambda basis read from lambda_basis.txt: all 24 vectors pass the definitional test;"
      " Gram det = 2^24 = det(Lambda)  ->  Z-basis of Lambda.")

_unit = [[F(int(i == j)) for j in range(24)] for i in range(24)]
_Binv = LL.frac_solve_rows(B, _unit)     # row i: coordinates of the i-th unit vector


def lam_coords(v):
    return [sum(v[i] * _Binv[i][k] for i in range(24) if v[i] != 0) for k in range(24)]


def triple(mul, u, v):
    x, y, z = u[0:8], u[8:16], u[16:24]
    xp, yp, zp = v[0:8], v[8:16], v[16:24]
    P1 = LL.vadd(mul(x, xp), mul(z, yp), mul(y, zp))
    P2 = LL.vadd(mul(y, yp), mul(x, zp), mul(z, xp))
    P3 = LL.vadd(mul(z, zp), mul(y, xp), mul(x, yp))
    return P1 + P2 + P3


def f2_then_q_independent(rows, order):
    """24 row indices: greedily F2-independent first, then completed by rows
    independent mod the prime P_MOD (hence independent over Q)."""
    f2, chosen = {}, []
    for idx in order:
        m = 0
        for i, a in enumerate(rows[idx]):
            if a & 1:
                m |= 1 << i
        while m:
            hb = m.bit_length() - 1
            if hb in f2:
                m ^= f2[hb]
            else:
                f2[hb] = m
                chosen.append(idx)
                break
    basis = {}

    def add(v):
        v = [a % P_MOD for a in v]
        for c in sorted(basis):
            if v[c]:
                f = v[c]
                v = [(a - f * b) % P_MOD for a, b in zip(v, basis[c])]
        piv = next((c for c in range(24) if v[c]), None)
        if piv is None:
            return False
        inv = pow(v[piv], P_MOD - 2, P_MOD)
        basis[piv] = [(a * inv) % P_MOD for a in v]
        return True

    for idx in chosen:
        assert add(rows[idx])
    for idx in order:
        if len(chosen) == 24:
            break
        if idx not in chosen and add(rows[idx]):
            chosen.append(idx)
    return chosen


def mod2_algebra(gf, LC, BL):
    """F2 invariants (Fractions, from the definitions) of V_g = g(L.sbar)/2L and
    W_g = g(L.s)/2L inside Abar = L/2L: dimension, subalgebra?,
    d = 4 - dim span(XX), m = 8 - dim span(X + Abar X + X Abar)."""
    Vg = [LC.to_L(gf(LL.omul(b, LL.SBAR_VEC))) for b in BL]
    Wg = [LC.to_L(gf(LL.omul(b, LL.S_VEC))) for b in BL]
    out = {}
    for name, gens in (('V', Vg), ('W', Wg)):
        basis, _ = LL.f2_rowspace_basis(gens, 8)
        vecs = [LC.to_e(b) for b in basis]
        prods = [[c % 2 for c in LC.to_L(LL.omul(a, b))] for a in vecs for b in vecs]
        mods = [[c % 2 for c in LC.to_L(LL.omul(l, x))] for l in BL for x in vecs] + \
               [[c % 2 for c in LC.to_L(LL.omul(x, l))] for l in BL for x in vecs]
        out['dim_' + name] = len(basis)
        out['dim_span_' + name + name] = LL.f2_rank(prods)
        out[name + '_subalgebra'] = LL.f2_rank(basis + prods) == len(basis)
        out['d_' + name] = 4 - LL.f2_rank(prods)
        out['m_' + name] = 8 - LL.f2_rank(basis + mods)
    return out


def certify(label, gf, gi, rng, verbose=True, LC=None, BL=None):
    mul = (lambda x, y: gi(LL.omul(gf(x), gf(y))))
    prods = [triple(mul, B[i], B[j]) for i in range(24) for j in range(24)]
    n_in = sum(LL.in_Lambda_def(p) for p in prods)
    coords = [lam_coords(p) for p in prods]
    integral = all(c.denominator == 1 for row in coords for c in row)
    out = {'label': label, 'products_in_Lambda_by_definition': n_in, 'coords_integral': integral}
    if LC is not None:
        out.update(mod2_algebra(gf, LC, BL))
    if verbose:
        print(f"\n--- {label}")
        print(f"  1-2. products passing the definitional Lambda test: {n_in}/576;"
              f"  Lambda-coordinates integral: {integral}")
    if n_in != 576 or not integral:
        if verbose:
            print("  -> the product does NOT close on Lambda; no index is defined.")
        return out
    rows = [[int(c) for c in row] for row in coords]
    r2 = LL.f2_rank(rows)
    l2_A2 = 24 - r2
    l2_A4 = 48 - LL.log2_order_mod4_span(rows)
    order = list(range(576))
    rng.shuffle(order)
    T = f2_then_q_independent(rows, order)
    D = abs(LL.det_bareiss([rows[k] for k in T]))
    fac = sympy.factorint(D)
    odd = {p: e for p, e in fac.items() if p != 2}
    odd_ranks = {p: LL.rank_mod_p(rows, p) for p in odd}
    no_odd = all(v == 24 for v in odd_ranks.values())
    elem2 = (l2_A4 == l2_A2)
    index = 2 ** l2_A2 if (no_odd and elem2) else None
    out.update({'f2_rank': r2, 'log2_A_mod_2A': l2_A2, 'log2_A_mod_4A': l2_A4,
                'det_T': D, 'det_T_factorization': {str(p): e for p, e in fac.items()},
                'ranks_mod_odd_primes_of_det': {str(p): v for p, v in odd_ranks.items()},
                'T_pairs': [(k // 24 + 1, k % 24 + 1) for k in T],
                'index_certified': index})
    if verbose:
        print(f"  3. rank over F2 of the 576 coordinate rows: r2 = {r2}  ->  |A/2A| = 2^{l2_A2}")
        print(f"  4. order of the Z/4-span of the rows: 2^{48 - l2_A4}  ->  |A/4A| = 2^{l2_A4}"
              f"  ->  2-part of A elementary abelian: {elem2}")
        print(f"  5. 24 independent products T (pairs (i,j) of b_i*b_j): {out['T_pairs']}")
        print(f"     |det T| = {D} = " + ' * '.join(f'{p}^{e}' for p, e in sorted(fac.items())))
        print(f"     ranks of the 576 rows mod the odd primes of |det T|: {odd_ranks}"
              f"  ->  no odd torsion in A: {no_odd}")
        if index is not None:
            print(f"  => A = Lambda/S = (Z/2)^{l2_A2},  [Lambda : S] = {index}"
                  + ("  (S = Lambda)" if index == 1 else "  (and 2 Lambda <= S)"))
        if LC is not None:
            ub = out['d_W'] + 2 * out['d_V'] + out['m_W'] + 2 * out['m_V']
            print(f"  conceptual reading: V_g, W_g subalgebras: {out['V_subalgebra']}, {out['W_subalgebra']};"
                  f" d(V) = {out['d_V']}, m(V) = {out['m_V']}, d(W) = {out['d_W']}, m(W) = {out['m_W']}")
            print(f"     (T3) [Lambda : S + (2L)^3] = 2^{out['d_W'] + 2 * out['d_V']};"
                  f"  (T4) [Lambda : S + 2Lambda] <= 2^{ub}  (F2-rank gives exactly 2^{24 - r2})"
                  + ("  -> with (T2): S = Lambda" if ub == 0 else ""))
    return out


def perm_maps(pi):
    pinv = LL.perm_inverse(pi)
    return (lambda x: LL.perm_apply(pi, x)), (lambda x: LL.perm_apply(pinv, x))


def refl_maps(r2):
    r = [F(a, 2) for a in r2]
    f = (lambda x: LL.refl_apply(r, x))
    return f, f


BL = LL.l_basis()
LC = LL.LCoords(BL)
_LL = [LC.to_L(LL.omul(a, b)) for a in BL for b in BL]
print(f"(T2) [L : span(L.L)] = {LL.hnf_index(_LL, 8)[0]}  (64 products of the L-basis, Fractions)"
      "  ->  (4L)^3 <= S_g for every twist, so Lambda/S_g is a 2-group of exponent dividing 4")
rng = random.Random(20261013)
results = []
results.append(certify("sigma = (1 2): pi = (2,1,3,4,5,6,7)", *perm_maps((2, 1, 3, 4, 5, 6, 7)), rng, LC=LC, BL=BL))
results.append(certify("C2 witness pi = (5 6 7): pi(5)=6, pi(6)=7, pi(7)=5", *perm_maps((1, 2, 3, 4, 6, 7, 5)), rng, LC=LC, BL=BL))
results.append(certify("inverse of the C2 witness, (5 7 6)  [expected NOT to close]", *perm_maps((1, 2, 3, 4, 7, 5, 6)), rng, LC=LC, BL=BL))
results.append(certify("C3 witness: reflection in r = (1,1,1,1,1,1,1,-1)/2", *refl_maps((1, 1, 1, 1, 1, 1, 1, -1)), rng, LC=LC, BL=BL))

census = json.load(open('twist_results.json'))
idx1 = [r for r in census['reflections'] if r['closes'] and r['span_index'] == 1
        and tuple(r['root_doubled_e_coords']) != (1, 1, 1, 1, 1, 1, 1, -1)]
print(f"\n--- the other {len(idx1)} reflections with census index 1 (compact)")
n_cert = int(results[3].get('index_certified') == 1)
for rr in idx1:
    res = certify(f"reflection r = {rr['root']}", *refl_maps(tuple(rr['root_doubled_e_coords'])), rng,
                  verbose=False, LC=LC, BL=BL)
    ok = res.get('index_certified') == 1
    n_cert += ok
    print(f"  r = {rr['root']:44s} closes {res['products_in_Lambda_by_definition'] == 576 and res['coords_integral']}, "
          f"r2 = {res.get('f2_rank')}, |det T| = {res.get('det_T')} "
          f"(odd primes {list(res.get('ranks_mod_odd_primes_of_det', {}).keys())}, all ranks 24), "
          f"index 1 certified: {ok}; d,m of V,W = {res['d_V']},{res['m_V']},{res['d_W']},{res['m_W']}")
    results.append(res)
print(f"\nindex 1 certified without any HNF for {n_cert} of {len(idx1) + 1} reflections")
json.dump(results, open('06_certificates.json', 'w'), indent=1)
