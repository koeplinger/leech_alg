"""
s5_verify_certificates.py -- STANDALONE verifier (does not import leech_core).

Reads   lambda_basis_ecoords_doubled.txt  (24 rows: 2 x e-coordinates of a
        claimed Z-basis u_1..u_24 of Wilson's Lambda, blocks x|y|z)
        s5_certificates.json             (functionals and 24-product sets)
and checks, with its own octonion product and Fractions only:
  (1) each u_i lies in Lambda by Wilson's definition: blocks in L = D8+ (odd
      half-integer coset); x+y, x+z, y+z in Lsbar; x+y+z in Ls.  Membership
      uses the inverse property of the octonions: t in Lsbar <=> t.(s/2) in L
      (sbar^{-1} = s/2), and t in Ls <=> t.(sbar/2) in L;
  (2) Gram determinant 2^24 = det(Lambda), so u_1..u_24 span Lambda;
  (3) for each twist g: all 576 products u_i *_g u_j have integral
      coordinates (closure), the listed mod-2 functionals vanish on all of them
      and are independent (index >= 2^k), and the listed 24-product sets have
      the listed determinants whose gcd is 2^k (index divides 2^k).
Conventions: e0 identity, Fano triples (1,2,4),(2,3,5),(3,4,6),(4,5,7),
(5,6,1),(6,7,2),(7,1,3) (cyclic +, anticyclic -), x ._g y = g^{-1}(g(x).g(y)),
P1 = x.x'+z.y'+y.z', P2 = y.y'+x.z'+z.x', P3 = z.z'+y.x'+x.y'.
Twists: sigma swaps e1,e2; c3: e5->e6->e7->e5 (g(e_a) = e_{pi(a)});
refl: x -> x - <x,r> r with r = (e0+e1+...+e6-e7)/2.
Output: s5_verify_certificates.out
"""
import json
from fractions import Fraction as Fr
from math import gcd

FANO = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]
TAB = {}
for a in range(8):
    TAB[(0, a)] = (1, a)
    TAB[(a, 0)] = (1, a)
for a in range(1, 8):
    TAB[(a, a)] = (-1, 0)
for (a, b, c) in FANO:
    for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
        TAB[(x, y)] = (1, z)
        TAB[(y, x)] = (-1, z)


def om(x, y):
    z = [Fr(0)] * 8
    for a in range(8):
        if x[a]:
            for b in range(8):
                if y[b]:
                    sg, c = TAB[(a, b)]
                    z[c] += sg * x[a] * y[b]
    return z


def inL(x):
    if all(t.denominator == 1 for t in x):
        return sum(x) % 2 == 0
    if all(t.denominator == 2 for t in x):
        return int(sum(x)) % 2 == 1
    return False


s = [Fr(-1, 2)] + [Fr(1, 2)] * 7
sbar = [Fr(-1, 2)] * 8
half_s = [t / 2 for t in s]
half_sbar = [t / 2 for t in sbar]
in_Lsbar = lambda t: inL(om(t, half_s))
in_Ls = lambda t: inL(om(t, half_sbar))
add = lambda a, b: [p + q for p, q in zip(a, b)]


def det(M):
    A = [[Fr(t) for t in r] for r in M]
    n, d = len(A), Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] / A[c][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return d


def inv(M):
    n = len(M)
    A = [[Fr(M[i][j]) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [t / pv for t in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [row[n:] for row in A]


lines = []


def log(*a):
    t = " ".join(str(x) for x in a)
    print(t)
    lines.append(t)


B = [[Fr(int(t), 2) for t in l.split()] for l in open('lambda_basis_ecoords_doubled.txt') if not l.startswith('#')]
assert len(B) == 24
ok1 = True
for u in B:
    x, y, z = u[:8], u[8:16], u[16:]
    ok1 &= inL(x) and inL(y) and inL(z)
    ok1 &= in_Lsbar(add(x, y)) and in_Lsbar(add(x, z)) and in_Lsbar(add(y, z)) and in_Ls(add(add(x, y), z))
log("(1) all 24 basis vectors satisfy Wilson's definition:", ok1)
gram = [[sum(a * b for a, b in zip(u, v)) for v in B] for u in B]
log("(2) Gram determinant:", det(gram), "(= 2^24 = det Lambda, so the vectors span Lambda)")
Binv = inv(B)


def twist(name):
    if name == 'sigma':
        perm = {0: 0, 1: 2, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}
        g = lambda x: [x[[k for k in perm if perm[k] == i][0]] for i in range(8)]   # (g x)_{pi(a)} = x_a
        return g, g
    if name == 'c3':
        perm = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 6, 6: 7, 7: 5}
        pinv = {v: k for k, v in perm.items()}
        g = lambda x: [x[pinv[i]] for i in range(8)]          # g(e_a) = e_{pi(a)}
        gi = lambda x: [x[perm[i]] for i in range(8)]
        return g, gi
    if name == 'refl':
        r = [Fr(1, 2)] * 7 + [Fr(-1, 2)]
        g = lambda x: [xi - sum(a * b for a, b in zip(x, r)) * ri for xi, ri in zip(x, r)]
        return g, g


def triple(u, v, g, gi):
    tw = lambda a, b: gi(om(g(a), g(b)))
    x, y, z = u[:8], u[8:16], u[16:]
    xp, yp, zp = v[:8], v[8:16], v[16:]
    P1 = add(add(tw(x, xp), tw(z, yp)), tw(y, zp))
    P2 = add(add(tw(y, yp), tw(x, zp)), tw(z, xp))
    P3 = add(add(tw(z, zp), tw(y, xp)), tw(x, yp))
    return P1 + P2 + P3


cert = json.load(open('s5_certificates.json'))
for name in ('sigma', 'c3', 'refl'):
    g, gi = twist(name)
    coords = {}
    integral = 0
    for i in range(24):
        for j in range(24):
            p = triple(B[i], B[j], g, gi)
            c = [sum(p[k] * Binv[k][m] for k in range(24)) for m in range(24)]
            integral += all(t.denominator == 1 for t in c)
            coords[(i, j)] = c
    cd = cert[name]
    k = cd['k']
    fs = cd['functionals']
    van = all(sum(f[m] * int(c[m]) for m in range(24)) % 2 == 0 for f in fs for c in coords.values())
    # rank of the functionals over F_2
    rows = [int("".join(map(str, f)), 2) for f in fs]
    rk = 0
    while rows:
        pv = max(rows)
        top = pv.bit_length() - 1
        rows = [r ^ pv if (r >> top) & 1 else r for r in rows if r != pv]
        rows = [r for r in rows if r]
        rk += 1
    dets = []
    for st in cd['sets']:
        d = det([coords[tuple(pq)] for pq in st['pairs']])
        assert d.denominator == 1 and abs(int(d)) == abs(st['det'])
        dets.append(abs(int(d)))
    gg = 0
    for d in dets:
        gg = gcd(gg, d)
    log("(3) %-5s: closure %d/576 integral | %d independent mod-2 functionals vanish: %s (index >= 2^%d) |"
        " gcd of %d certificate determinants = %d (index divides it)  ==> [Lambda:S] = %d"
        % (name, integral, rk, van, rk, len(dets), gg, gg if (van and gg == 2 ** rk) else -1))
with open('s5_verify_certificates.out', 'w') as f:
    f.write("\n".join(lines) + "\n")
