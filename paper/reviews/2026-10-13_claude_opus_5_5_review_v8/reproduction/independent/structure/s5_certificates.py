"""
s5_certificates.py -- HNF-free certificates for the span index of the three
witness twists (conventions: leech_core.py):

  sigma   = transposition (1 2)                      (paper; index 2^16)
  c3      = 3-cycle pi(5)=6, pi(6)=7, pi(7)=5, g = g_pi  (convention A; 2^10)
  refl    = reflection s_r, r = (e0+e1+...+e6-e7)/2  (index 1)

For each twist g and the Z-basis u_1..u_24 of Lambda (lambda_basis_Lcoords.txt):
  closure : all 576 products u_i *_g u_j have integral coordinates;
  UPPER   : a few sets of 24 index pairs (i,j); for each set the products
            lie in S, so [Lambda : S] divides |det(coords)|; the sets are
            chosen so that gcd of the |det| is 2^k, hence [Lambda : S] | 2^k;
  LOWER   : k vectors f_1..f_k in F_2^24, linearly independent, with
            f . coords(u_i *_g u_j) = 0 mod 2 for all 576 products; then S lies
            in the index-2^k sublattice {x : f_t . x = 0 mod 2}, so
            [Lambda : S] >= 2^k.
Together: [Lambda : S] = 2^k exactly, verifiable with integer arithmetic only.
The sets are found by p-adic basis exchange: replacing row t of B by a
candidate c changes det(B) to (c . adj(B))_t exactly; if v_p(det B) exceeds
v_p([Lambda:S]) some product c is outside the Z_p-span of B, so some swap
lowers v_p(det).  One set is pushed to v_2(det) = k, and for every odd prime
p dividing its determinant a further set is pushed to v_p(det) = 0.
Outputs: s5_certificates.out, s5_certificates.json (pairs, functionals), and
cert_<name>_products_ecoords_doubled.txt (the 24 product vectors, e-coords x2).
"""
import json
import random
from fractions import Fraction as Fr
from math import gcd
import numpy as np
from sympy import factorint
import leech_core as C

lines = []


def log(*a):
    s_ = " ".join(str(t) for t in a)
    print(s_)
    lines.append(s_)


U = [[int(t) for t in l.split()] for l in open('lambda_basis_Lcoords.txt') if not l.startswith('#')]
U = np.array(U, dtype=np.int64)
assert U.shape == (24, 24)

TWISTS = {
    'sigma': C.perm_matrix_e({1: 2, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}),
    'c3': C.perm_matrix_e({1: 1, 2: 2, 3: 3, 4: 4, 5: 6, 6: 7, 7: 5}),
    'refl': C.reflection_matrix_e([Fr(1, 2)] * 7 + [Fr(-1, 2)]),
}


def int_det(rows):
    return int(C.mat_det([[int(t) for t in r] for r in rows]))


def adjugate(B):
    Bi = C.mat_inv([[Fr(int(t)) for t in r] for r in B])
    d = int(C.mat_det([[int(t) for t in r] for r in B]))
    return [[int(t * d) for t in row] for row in Bi], d


def vp(n, p):
    n = abs(n)
    if n == 0:
        return 10 ** 9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def independent_start(Cc, seed):
    """24 rows of Cc that are linearly independent over Q (random order)."""
    rng = random.Random(seed)
    rows = Cc.tolist()
    order = list(range(len(rows)))
    rng.shuffle(order)
    p = 2 ** 61 - 1
    chosen, basis = [], []
    for i in order:
        v = [x % p for x in rows[i]]
        for (piv, b) in basis:
            if v[piv]:
                f = v[piv] * pow(b[piv], -1, p) % p
                v = [(x - f * y) % p for x, y in zip(v, b)]
        nz = next((k for k in range(24) if v[k]), None)
        if nz is not None:
            basis.append((nz, v))
            chosen.append(i)
        if len(chosen) == 24:
            return chosen
    raise ValueError("rank < 24")


def padic_exchange(Cc, chosen, p, target):
    """Swap rows until v_p(det) <= target; returns (chosen, det)."""
    rows = Cc.tolist()
    chosen = list(chosen)
    while True:
        B = [rows[i] for i in chosen]
        adj, d = adjugate(B)
        if vp(d, p) <= target:
            return chosen, d
        best = None
        for c in range(len(rows)):
            if c in chosen:
                continue
            vals = [sum(rows[c][k] * adj[k][t] for k in range(24)) for t in range(24)]
            for t in range(24):
                if vals[t]:
                    v = vp(vals[t], p)
                    if v < vp(d, p) and (best is None or (v, abs(vals[t])) < best[0]):
                        best = ((v, abs(vals[t])), c, t)
        assert best is not None, "p-adic exchange stuck: v_p(det) is the true valuation"
        chosen[best[2]] = best[1]


def null_mod2(Cc):
    """Basis of {f in F_2^24 : Cc f = 0 mod 2} (vectors as int lists)."""
    M = [[int(t) & 1 for t in row] for row in Cc]
    # row-reduce M over F_2
    piv_cols, R = [], []
    for row in M:
        r = row[:]
        for (pc, pr) in zip(piv_cols, R):
            if r[pc]:
                r = [a ^ b for a, b in zip(r, pr)]
        nz = next((k for k in range(24) if r[k]), None)
        if nz is not None:
            # keep R in reduced form
            R = [[a ^ b for a, b in zip(q, r)] if q[nz] else q for q in R]
            R.append(r)
            piv_cols.append(nz)
    free = [k for k in range(24) if k not in piv_cols]
    out = []
    for fcol in free:
        f = [0] * 24
        f[fcol] = 1
        for (pc, pr) in zip(piv_cols, R):
            if pr[fcol]:
                f[pc] = 1
        out.append(f)
    return out


cert = {}
for name, G in TWISTS.items():
    T = C.tensor_int(C.structure_tensor_L(G))
    P = C.triple_products(U, T)
    Cc, ok = C.coords_in_basis(P, U)
    idx, _ = C.span_index(Cc)
    k = idx.bit_length() - 1
    log("%-5s closes: %s (576/576 integral: %d)   HNF index = %d = 2^%d" % (name, bool(ok.all()), int(ok.sum()), idx, k))
    # LOWER bound certificate
    fs = null_mod2(Cc)
    lower_ok = all(sum(f[t] * int(row[t]) for t in range(24)) % 2 == 0 for f in fs for row in Cc)
    # independence of the f's over F_2 (rank)
    rk = 0
    rows2 = [int("".join(map(str, f)), 2) for f in fs]
    while rows2:
        piv = max(rows2)
        top = piv.bit_length() - 1
        rows2 = [r ^ piv if (r >> top) & 1 else r for r in rows2 if r != piv]
        rows2 = [r for r in rows2 if r]
        rk += 1
    log("      LOWER: %d independent mod-2 functionals vanish on all 576 products: %s  => index >= 2^%d"
        % (rk, lower_ok, rk))
    # UPPER bound certificate: sets of 24 products, gcd of |det| = 2^k
    base = independent_start(Cc, 0)
    ch2, d2 = padic_exchange(Cc, base, 2, k)
    sets = [(ch2, d2)]
    for p in sorted(factorint(abs(d2))):
        if p == 2:
            continue
        chp, dp = padic_exchange(Cc, ch2, p, 0)
        sets.append((chp, dp))
    g = 0
    for chosen, d in sets:
        assert abs(int_det([Cc[i].tolist() for i in chosen])) == abs(d)
        g = gcd(g, abs(d))
    log("      UPPER: %d set(s) of 24 products; |det| = %s; gcd = %d = 2^%d  => index divides 2^%d"
        % (len(sets), [abs(d) for _, d in sets], g, g.bit_length() - 1, g.bit_length() - 1))
    assert g == idx
    pairs_sets = [[(i // 24, i % 24) for i in chosen] for chosen, _ in sets]
    cert[name] = dict(index=idx, k=k, functionals=fs,
                      sets=[dict(pairs=ps, det=int(d)) for ps, (_, d) in zip(pairs_sets, sets)])
    with open('cert_%s_products_ecoords_doubled.txt' % name, 'w') as f:
        f.write("# %s: certificate products u_i *_g u_j (i, j = 0-based rows of lambda_basis_*.txt),\n"
                "# rows = 2*(e-coords of P1|P2|P3); one block of 24 rows per set, with its determinant\n" % name)
        for sidx, ((chosen, d), ps) in enumerate(zip(sets, pairs_sets)):
            f.write("# set %d: det of coordinate matrix = %d\n" % (sidx, d))
            for (i, j), crow in zip(ps, chosen):
                vec = P[crow]
                e = []
                for blk in range(3):
                    e += [int(2 * t) for t in C.from_L(vec[8 * blk:8 * blk + 8].tolist())]
                f.write("%d %d : " % (i, j) + " ".join(map(str, e)) + "\n")
json.dump(cert, open('s5_certificates.json', 'w'), indent=1)
with open('s5_certificates.out', 'w') as f:
    f.write("\n".join(lines) + "\n")
