"""Core exact-arithmetic objects for the v8 mathematics review.

Everything here is built from the paper's definitions (Section 2):
 - octonion product from the Fano triples (1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)
 - L = D8+ : all-integer with even sum, or all-half-integer with odd sum
 - s = (-e0+e1+...+e7)/2, sbar = (-e0-e1-...-e7)/2
 - Wilson's Lambda: x,y,z in L; x+y,x+z,y+z in L sbar; x+y+z in L s
All arithmetic is integer arithmetic on doubled (or 4x) coordinates, or Fractions.
"""
from fractions import Fraction as Fr
import itertools
import numpy as np

FANO = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]


def octonion_tensor(fano=FANO):
    """T[a,b,c] with e_a e_b = sum_c T[a,b,c] e_c (integer)."""
    T = np.zeros((8, 8, 8), dtype=np.int64)
    for a in range(8):
        T[0, a, a] = 1
        T[a, 0, a] = 1
    for a in range(1, 8):
        T[a, a, 0] = -1
    for (a, b, c) in fano:
        for (p, q, r) in [(a, b, c), (b, c, a), (c, a, b)]:
            T[p, q, r] = 1
            T[q, p, r] = -1
    return T


T = octonion_tensor()


def twisted_tensor(perm, T=T):
    """perm: tuple of length 8 with perm[0]=0, a permutation of 0..7.
    x ._pi y := pi^{-1}(pi(x) . pi(y)), pi(e_a) = e_{pi(a)}.
    Then T_pi[a,b,c] = T[pi a, pi b, pi c]."""
    p = np.array(perm)
    return T[np.ix_(p, p, p)]


SIGMA = (0, 2, 1, 3, 4, 5, 6, 7)
TS = twisted_tensor(SIGMA)


def omul(x, y, TT=T):
    """product of coordinate vectors (any numeric dtype incl. object/Fraction)."""
    out = [0] * 8
    for a in range(8):
        if x[a] == 0:
            continue
        for b in range(8):
            if y[b] == 0:
                continue
            for c in range(8):
                t = TT[a, b, c]
                if t:
                    out[c] += t * x[a] * y[b]
    return out


def conj(x):
    return [x[0]] + [-v for v in x[1:]]


def norm(x):
    return sum(v * v for v in x)


def in_L(x):
    """x: sequence of Fractions/ints (real coords)."""
    xs = [Fr(v) for v in x]
    if all(v.denominator == 1 for v in xs):
        return sum(xs) % 2 == 0
    if all(v.denominator == 2 for v in xs):
        return (sum(xs)) % 2 == 1  # sum of 8 half-integers is an integer
    return False


half = Fr(1, 2)
s = [-half] + [half] * 7
sbar = [-half] + [-half] * 7

# paper's basis of L (Appendix A)
def e(i, coef=1):
    v = [Fr(0)] * 8
    v[i] = Fr(coef)
    return v

def add(*vs):
    return [sum(t) for t in zip(*vs)]

def scal(c, v):
    return [c * t for t in v]

LB = [s, add(e(0), e(1)), add(e(0), e(1, -1))] + [add(e(0), e(k)) for k in range(2, 7)]


def mat_inv_fr(M):
    """inverse of a square matrix of Fractions (Gauss-Jordan)."""
    n = len(M)
    A = [[Fr(v) for v in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [v / pv for v in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [a - f * b for a, b in zip(A[r], A[col])]
    return [row[n:] for row in A]


def det_fr(M):
    n = len(M)
    A = [[Fr(v) for v in row] for row in M]
    d = Fr(1)
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r][col] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
            d = -d
        pv = A[col][col]
        d *= pv
        for r in range(col + 1, n):
            if A[r][col] != 0:
                f = A[r][col] / pv
                A[r] = [a - f * b for a, b in zip(A[r], A[col])]
    return d


LBinv = mat_inv_fr(LB)  # rows: coords; coordinates a = x * LBinv  (row vector)


def L_coords(x):
    """coordinates of x (row vector) in basis LB; returns list of Fractions."""
    return [sum(Fr(x[i]) * LBinv[i][j] for i in range(8)) for j in range(8)]


# --- F2 linear algebra helpers ---------------------------------------------
def f2_rref(rows, ncols):
    rows = [list(r) for r in rows]
    piv_cols = []
    r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(rows)) if rows[i][c] % 2), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c] % 2:
                rows[i] = [(a + b) % 2 for a, b in zip(rows[i], rows[r])]
        piv_cols.append(c)
        r += 1
    return [[v % 2 for v in row] for row in rows[:r]], piv_cols


def f2_nullspace(rows, ncols):
    R, piv = f2_rref(rows, ncols)
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for f in free:
        v = [0] * ncols
        v[f] = 1
        for i, pc in enumerate(piv):
            v[pc] = R[i][f] % 2
        basis.append(v)
    return basis


def sublattice_mod2(gens):
    """gens: list of vectors in L; return F2 rref basis of their image in L/2L (in LB coords)."""
    rows = []
    for g in gens:
        a = L_coords(g)
        assert all(t.denominator == 1 for t in a)
        rows.append([int(t) % 2 for t in a])
    R, piv = f2_rref(rows, 8)
    return R


def int_rank_det(M):
    return det_fr(M)


# Lsbar, Ls generated by LB_k * sbar, LB_k * s
LsbarB = [omul(l, sbar) for l in LB]
LsB = [omul(l, s) for l in LB]

V0 = sublattice_mod2(LsbarB)   # F2 basis of Lsbar/2L in LB-coords
W0 = sublattice_mod2(LsB)      # F2 basis of Ls/2L
H_V = f2_nullspace(V0, 8)      # parity checks: a in V0  <=>  H_V a = 0 ; (nullspace of V0 as rows = dual)
H_W = f2_nullspace(W0, 8)


def f2_in(a2, H):
    return all(sum(h[i] * a2[i] for i in range(8)) % 2 == 0 for h in H)


def in_Lsbar(x):
    if not in_L(x):
        return False
    a = [int(t) % 2 for t in L_coords(x)]
    return f2_in(a, H_V)


def in_Ls(x):
    if not in_L(x):
        return False
    a = [int(t) % 2 for t in L_coords(x)]
    return f2_in(a, H_W)


def in_Lambda(v):
    x, y, z = v[0:8], v[8:16], v[16:24]
    if not (in_L(x) and in_L(y) and in_L(z)):
        return False
    if not (in_Lsbar(add(x, y)) and in_Lsbar(add(x, z)) and in_Lsbar(add(y, z))):
        return False
    return in_Ls(add(x, y, z))


# --- integer lattice index (row style) --------------------------------------
def xgcd(a, b):
    """returns (g, x, y) with x*a + y*b = g = gcd(a,b) >= 0"""
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return a, x0, y0


def bareiss_det(M):
    n = len(M)
    A = [list(map(int, r)) for r in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def independent_rows(rows, n, p=2147483647):
    """greedy choice of n rows independent mod p (hence over Q); None if rank < n."""
    basis = {}  # pivot col -> reduced row mod p
    chosen = []
    for idx, r in enumerate(rows):
        v = [int(x) % p for x in r]
        for c in range(n):
            if v[c] and c in basis:
                f = v[c]
                b = basis[c]
                v = [(x - f * y) % p for x, y in zip(v, b)]
        piv = next((c for c in range(n) if v[c]), None)
        if piv is None:
            continue
        inv = pow(v[piv], p - 2, p)
        v = [(x * inv) % p for x in v]
        # keep basis rows reduced on the new pivot as well
        for c, b in list(basis.items()):
            if b[piv]:
                f = b[piv]
                basis[c] = [(x - f * y) % p for x, y in zip(b, v)]
        basis[piv] = v
        chosen.append(idx)
        if len(chosen) == n:
            return chosen
    return None


def lattice_index(rows, n):
    """[Z^n : span_Z(rows)] exactly (0 if not full rank).
    D = |det| of n independent rows is a multiple of the index; then HNF mod D
    (column by column, with modulus R_k = D / prod of earlier pivots)."""
    rows = [list(map(int, r)) for r in rows]
    sel = independent_rows(rows, n)
    if sel is None:
        return 0
    D = abs(bareiss_det([rows[i] for i in sel]))
    assert D > 0
    R = D
    cur = [[x % R for x in r] for r in rows]
    index = 1
    for k in range(n):
        P = None
        nxt = []
        for r in cur:
            if r[k] % R == 0:
                nxt.append(r)
                continue
            if P is None:
                P = r
                continue
            g, a, b = xgcd(P[k], r[k])
            pk, rk = P[k] // g, r[k] // g
            newP = [(a * u + b * w) % R for u, w in zip(P, r)]
            rr = [(pk * w - rk * u) % R for u, w in zip(P, r)]
            P = newP
            nxt.append(rr)
        if P is None:
            d = R
        else:
            d = xgcd(P[k], R)[0]
        index *= d
        R //= d
        if R == 1:
            # remaining sub-lattice is everything
            return index
        cur = [[x % R for x in r] for r in nxt if any(x % R for x in r[k + 1:])]
    return index
