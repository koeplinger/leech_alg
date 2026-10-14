"""
leechlib.py -- exact-arithmetic library for the twisted triple-octonion product
on Wilson's Leech lattice (independent re-implementation for the review of
"An order on the Leech lattice from a Z3-symmetric triple-octonion product", v8).

Everything here is exact: Python ints (arbitrary precision) and
fractions.Fraction.  numpy is used only with dtype=int64 for bulk products of
small integers, and every such use is guarded by an explicit overflow bound
(see `_check_int64_bound`).

CONVENTIONS (stated once, used everywhere)
------------------------------------------
* Octonions O = R^8 with basis e0..e7, e0 the identity.  A vector x is the list
  [x0, ..., x7] of its coordinates in this basis.
* Multiplication from the seven Fano triples
      (1,2,4) (2,3,5) (3,4,6) (4,5,7) (5,6,1) (6,7,2) (7,1,3):
  for a triple (a,b,c): e_a e_b = +e_c, e_b e_c = +e_a, e_c e_a = +e_b, and
  the reversed orders carry a minus sign; e_i e_i = -e0 for i >= 1.
* Norm N(x) = sum x_i^2 (squared Euclidean length); <x,y> = sum x_i y_i.
* L = D8+ (odd convention): all coordinates integers with EVEN sum, or all
  coordinates half-integers with ODD sum.  Roots of L have norm 2.
* s = (-e0+e1+...+e7)/2 (in L), sbar = (-e0-e1-...-e7)/2 (not in L).
  L*sbar = { l . sbar : l in L },  L*s = { l . s : l in L }  (l on the LEFT).
* Wilson's Leech lattice
      Lambda = {(x,y,z) in L^3 : x+y, x+z, y+z in L*sbar ;  x+y+z in L*s }.
* Twisted product for a linear map g of R^8 with g(L) = L:
      x ._g y := g^{-1}( g(x) . g(y) ).
  For a permutation pi of the imaginary indices {1..7} (pi fixes 0) the map
  is pi(e_a) = e_{pi(a)}, i.e. (pi(x))_{pi(a)} = x_a, and the paper's
      x ._pi y := pi^{-1}( pi(x) . pi(y) )
  is exactly the g = pi case.  A permutation is stored as the tuple
  (pi(1), ..., pi(7)).  sigma = (1 2) is the tuple (2,1,3,4,5,6,7).
  For a root r of L the reflection is s_r(x) = x - <x,r> r (an involution).
* Triple product (Definition 3.3 of the paper), u=(x,y,z), v=(x',y',z'):
      P1 = x.x' + z.y' + y.z'
      P2 = y.y' + x.z' + z.x'
      P3 = z.z' + y.x' + x.y'
  with every block product the twisted product ._g.
* Closure: u *_g v in Lambda for all u, v in Lambda (checked on the 576
  products of a Z-basis, which suffices by bilinearity).
* Span index: S_g = Z-span{u *_g v : u, v in Lambda} = Z-span of the 576
  basis-pair products; [Lambda : S_g] from an exact Hermite normal form of
  the 576 x 24 matrix of Lambda-coordinates of those products.

COORDINATE SYSTEMS USED INTERNALLY
----------------------------------
* "e-coordinates": coordinates in e0..e7 (Fractions), or "doubled
  e-coordinates" 2x (ints; every element of L has 2x in Z^8).
* "L-coordinates": integer coefficients w.r.t. the Z-basis of L built in
  `l_basis()` (rows of BL).  L^3 elements are 24-vectors of L-coordinates
  (x block, then y block, then z block).
* "Lambda-coordinates": integer coefficients w.r.t. the Z-basis of Lambda
  built in `lambda_basis()` (rows of BLam, which are given in L^3 coordinates).
"""

from fractions import Fraction
from itertools import product as iproduct
from math import gcd
import numpy as np

# --------------------------------------------------------------------------
# 1. Octonion multiplication
# --------------------------------------------------------------------------

FANO = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]


def _build_mul_table():
    """MUL[a][b] = (sign, c) meaning e_a e_b = sign * e_c."""
    mul = [[None] * 8 for _ in range(8)]
    for a in range(8):
        mul[0][a] = (1, a)
        mul[a][0] = (1, a)
    for a in range(1, 8):
        mul[a][a] = (-1, 0)
    for (a, b, c) in FANO:
        for (p, q, r) in ((a, b, c), (b, c, a), (c, a, b)):   # cyclic: +
            mul[p][q] = (1, r)
            mul[q][p] = (-1, r)                                 # anticyclic: -
    assert all(mul[a][b] is not None for a in range(8) for b in range(8))
    return mul


MUL = _build_mul_table()

# Structure tensor of the standard octonion product in e-coordinates:
# (x.y)_c = sum_{a,b} x_a y_b OCT[a,b,c]
OCT = np.zeros((8, 8, 8), dtype=np.int64)
for _a in range(8):
    for _b in range(8):
        _s, _c = MUL[_a][_b]
        OCT[_a, _b, _c] = _s


def omul(x, y):
    """Standard octonion product of two coordinate lists (ints or Fractions)."""
    z = [0] * 8
    for a in range(8):
        xa = x[a]
        if xa == 0:
            continue
        for b in range(8):
            yb = y[b]
            if yb == 0:
                continue
            sgn, c = MUL[a][b]
            z[c] += sgn * xa * yb
    return z


def conj(x):
    return [x[0]] + [-t for t in x[1:]]


def ip(x, y):
    return sum(a * b for a, b in zip(x, y))


def norm(x):
    return ip(x, x)


def vadd(*vs):
    return [sum(t) for t in zip(*vs)]


def vsub(x, y):
    return [a - b for a, b in zip(x, y)]


def vscale(k, x):
    return [k * a for a in x]


F = Fraction
HALF = Fraction(1, 2)
S_VEC = [F(-1, 2)] + [F(1, 2)] * 7            # s
SBAR_VEC = [F(-1, 2)] * 8                     # sbar


# --------------------------------------------------------------------------
# 2. The lattice L = D8+ (odd convention), definitional membership tests
# --------------------------------------------------------------------------

def in_L(x):
    """Definitional test of x in L = D8+ (integers/even sum or halves/odd sum)."""
    x = [F(t) for t in x]
    if all(t.denominator == 1 for t in x):
        return sum(x) % 2 == 0
    if all(t.denominator == 2 for t in x):
        return sum(x) % 2 == 1
    return False


def in_L_sbar(x):
    """x in L.sbar  <=>  (x . s)/2 in L.

    Reason: sbar^{-1} = conj(sbar)/N(sbar) = s/2, and in the alternative
    algebra O, (l . sbar) . sbar^{-1} = l (Artin's theorem), so
    x = l . sbar  <=>  l = x . (s/2)."""
    return in_L(vscale(HALF, omul(x, S_VEC)))


def in_L_s(x):
    """x in L.s  <=>  (x . sbar)/2 in L   (s^{-1} = conj(s)/N(s) = sbar/2)."""
    return in_L(vscale(HALF, omul(x, SBAR_VEC)))


def in_Lambda_def(u):
    """Definitional membership test for Wilson's Lambda (u = 24 Fractions)."""
    x, y, z = u[0:8], u[8:16], u[16:24]
    return (in_L(x) and in_L(y) and in_L(z)
            and in_L_sbar(vadd(x, y)) and in_L_sbar(vadd(x, z))
            and in_L_sbar(vadd(y, z)) and in_L_s(vadd(x, y, z)))


def roots_of_L():
    """All 240 roots of L (as Fraction lists): +-e_i +- e_j (112) and the 128
    half-integer vectors (+-1/2,...,+-1/2) with odd coordinate sum."""
    roots = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [F(0)] * 8
                    v[i] = F(si)
                    v[j] = F(sj)
                    roots.append(v)
    for signs in iproduct((1, -1), repeat=8):
        v = [F(sg, 2) for sg in signs]
        if sum(v) % 2 == 1:
            roots.append(v)
    assert len(roots) == 240 and all(norm(r) == 2 and in_L(r) for r in roots)
    return roots


# --------------------------------------------------------------------------
# 3. Exact integer lattice tools: HNF, SNF, determinant, F2 rank
# --------------------------------------------------------------------------

def xgcd(a, b):
    """Return (g, p, q) with p*a + q*b = g = gcd(a, b) >= 0."""
    old_r, r = a, b
    old_p, p = 1, 0
    old_q, q = 0, 1
    while r != 0:
        t = old_r // r
        old_r, r = r, old_r - t * r
        old_p, p = p, old_p - t * p
        old_q, q = q, old_q - t * q
    if old_r < 0:
        old_r, old_p, old_q = -old_r, -old_p, -old_q
    return old_r, old_p, old_q


def hnf(rows, n):
    """Exact row-style Hermite normal form of the Z-span of `rows`.

    rows : iterable of length-n integer sequences (Python ints).
    Returns the list of HNF rows (upper "echelon" form, positive pivots, each
    entry above a pivot reduced into [0, pivot)), one row per pivot column.
    Pure Python integers: no overflow, no modular shortcuts.  Every operation
    is a unimodular row operation (or the addition of a generator), so the
    Z-span is preserved exactly.
    """
    piv = {}                       # pivot column -> row (list of ints)

    def reduce_row_against_later(v, start):
        # make entries of v in pivot columns > start lie in [0, pivot)
        for k in range(start + 1, n):
            if k in piv and v[k] != 0:
                hk = piv[k]
                q = v[k] // hk[k]
                if q:
                    v = [a - q * b for a, b in zip(v, hk)]
        return v

    for r in rows:
        v = [int(t) for t in r]
        j = 0
        while True:
            while j < n and v[j] == 0:
                j += 1
            if j == n:
                break
            if j not in piv:
                if v[j] < 0:
                    v = [-a for a in v]
                piv[j] = reduce_row_against_later(v, j)
                break
            h = piv[j]
            if v[j] % h[j] == 0:
                q = v[j] // h[j]
                v = [a - q * b for a, b in zip(v, h)]
            else:
                g, p, q = xgcd(h[j], v[j])
                hj, vj = h[j] // g, v[j] // g
                newh = [p * a + q * b for a, b in zip(h, v)]      # pivot g
                newv = [vj * a - hj * b for a, b in zip(h, v)]    # zero at j
                # [[p, q], [vj, -hj]] has determinant -p*hj - q*vj = -1: unimodular
                piv[j] = reduce_row_against_later(newh, j)
                v = newv
    # final full Hermite reduction (entries above each pivot into [0, pivot)).
    # For each row, reduce at the later pivot columns in increasing order:
    # subtracting piv[k] only changes entries in columns >= k.
    cols = sorted(piv)
    for i in cols:
        piv[i] = reduce_row_against_later(piv[i], i)
    return [piv[c] for c in cols]


def hnf_index(rows, n):
    """Return (index, H) with index = [Z^n : span(rows)] (None if rank < n)."""
    H = hnf(rows, n)
    if len(H) < n:
        return None, H
    idx = 1
    for i, h in enumerate(H):
        assert h[i] > 0 and all(h[k] == 0 for k in range(i))
        idx *= h[i]
    return idx, H


def snf_diagonal(M):
    """Elementary divisors of a square nonsingular integer matrix (Python ints),
    by explicit unimodular row and column operations.  Returns the sorted list
    d_1 | d_2 | ... | d_n (all > 0)."""
    A = [list(map(int, r)) for r in M]
    n = len(A)
    m = len(A[0])
    diag = []
    t = 0
    while t < min(n, m):
        # find a nonzero entry of minimal absolute value in the submatrix
        best = None
        for i in range(t, n):
            for j in range(t, m):
                if A[i][j] != 0 and (best is None or abs(A[i][j]) < abs(A[best[0]][best[1]])):
                    best = (i, j)
        if best is None:
            break
        i, j = best
        A[t], A[i] = A[i], A[t]
        for row in A:
            row[t], row[j] = row[j], row[t]
        while True:
            done = True
            p = A[t][t]
            # clear column t
            for i in range(t + 1, n):
                if A[i][t] != 0:
                    q = A[i][t] // p
                    A[i] = [a - q * b for a, b in zip(A[i], A[t])]
                    if A[i][t] != 0:
                        done = False
            # clear row t
            for j in range(t + 1, m):
                if A[t][j] != 0:
                    q = A[t][j] // p
                    for row in A:
                        row[j] -= q * row[t]
                    if A[t][j] != 0:
                        done = False
            if done:
                # divisibility condition: p must divide every remaining entry
                bad = None
                for i in range(t + 1, n):
                    for j in range(t + 1, m):
                        if A[i][j] % p != 0:
                            bad = i
                            break
                    if bad is not None:
                        break
                if bad is None:
                    break
                A[t] = [a + b for a, b in zip(A[t], A[bad])]
                done = False
            if not done:
                # move the smallest nonzero entry of row t / column t to (t,t)
                best = None
                for i in range(t, n):
                    if A[i][t] != 0 and (best is None or abs(A[i][t]) < abs(A[best[0]][best[1]])):
                        best = (i, t)
                for j in range(t, m):
                    if A[t][j] != 0 and (best is None or abs(A[t][j]) < abs(A[best[0]][best[1]])):
                        best = (t, j)
                i, j = best
                A[t], A[i] = A[i], A[t]
                for row in A:
                    row[t], row[j] = row[j], row[t]
        diag.append(abs(A[t][t]))
        t += 1
    return sorted(diag)


def det_bareiss(M):
    """Exact determinant of a square integer matrix (fraction-free Bareiss)."""
    A = [list(map(int, r)) for r in M]
    n = len(A)
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


def f2_rank(rows):
    """Rank over GF(2) of integer row vectors (reduced mod 2)."""
    basis = {}                     # leading bit -> bitmask
    for r in rows:
        m = 0
        for i, a in enumerate(r):
            if int(a) & 1:
                m |= 1 << i
        while m:
            hb = m.bit_length() - 1
            if hb in basis:
                m ^= basis[hb]
            else:
                basis[hb] = m
                break
    return len(basis)


def frac_solve_rows(B, targets):
    """Solve c . B = t exactly over Q for each target row t (B square, rows are
    basis vectors).  B, targets: lists of Fraction/int lists.  Returns list of
    coefficient rows (Fractions).  Gauss-Jordan on the transposed system."""
    n = len(B)
    # we need c with sum_k c_k B[k][j] = t_j  ->  (B^T) c^T = t^T
    BT = [[F(B[k][j]) for k in range(n)] for j in range(n)]
    # invert BT once (exact)
    A = [row[:] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(BT)]
    for col in range(n):
        p = next(i for i in range(col, n) if A[i][col] != 0)
        A[col], A[p] = A[p], A[col]
        pv = A[col][col]
        A[col] = [a / pv for a in A[col]]
        for i in range(n):
            if i != col and A[i][col] != 0:
                f = A[i][col]
                A[i] = [a - f * b for a, b in zip(A[i], A[col])]
    inv = [row[n:] for row in A]          # inv = (B^T)^{-1}
    out = []
    for t in targets:
        tt = [F(a) for a in t]
        out.append([sum(inv[i][j] * tt[j] for j in range(n)) for i in range(n)])
    return out


# --------------------------------------------------------------------------
# 4. Z-bases of L, L.sbar, L.s and Lambda
# --------------------------------------------------------------------------

def l_basis():
    """Z-basis of L built from scratch: HNF of the doubled coordinates of all
    240 roots.  Returns BL as a list of 8 Fraction rows (e-coordinates).
    (The roots of the unimodular lattice E8 generate it.)"""
    roots = roots_of_L()
    H = hnf([[int(2 * t) for t in r] for r in roots], 8)
    assert len(H) == 8
    BL = [[F(a, 2) for a in h] for h in H]
    assert all(in_L(b) for b in BL)
    return BL


def int_inverse_scaled(M):
    """For a square matrix M with Fraction entries, return (Minv, den) where
    Minv is an integer matrix and M^{-1} = Minv / den (den > 0 minimal)."""
    n = len(M)
    A = [[F(a) for a in row] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for col in range(n):
        p = next(i for i in range(col, n) if A[i][col] != 0)
        A[col], A[p] = A[p], A[col]
        pv = A[col][col]
        A[col] = [a / pv for a in A[col]]
        for i in range(n):
            if i != col and A[i][col] != 0:
                f = A[i][col]
                A[i] = [a - f * b for a, b in zip(A[i], A[col])]
    inv = [row[n:] for row in A]
    den = 1
    for row in inv:
        for a in row:
            den = den * a.denominator // gcd(den, a.denominator)
    return [[int(a * den) for a in row] for row in inv], den


class LCoords:
    """Conversion between e-coordinates and L-coordinates (exact)."""

    def __init__(self, BL):
        self.BL = BL
        self.inv, self.den = int_inverse_scaled(BL)     # BL^{-1} = inv/den

    def to_L(self, x):
        """L-coordinates (ints) of x in e-coordinates; raises if x not in L."""
        c = [sum(F(x[a]) * self.inv[a][k] for a in range(8)) / self.den for k in range(8)]
        assert all(t.denominator == 1 for t in c), "vector not in L"
        return [int(t) for t in c]

    def to_e(self, c):
        return [sum(c[k] * self.BL[k][a] for k in range(8)) for a in range(8)]


def lambda_basis(BL, LC):
    """Build a Z-basis of Lambda from scratch.

    Returns dict with
      'BLam'   : 24x24 int matrix, rows = basis of Lambda in L^3-coordinates
      'Vrows'  : L-coordinates of the basis {l_k . sbar} of L.sbar
      'Wrows'  : L-coordinates of the basis {l_k . s} of L.s
      'Q'      : 12 x 24 GF(2) condition matrix (rows of 0/1) with
                 u in Lambda <=> u in Z^24 and Q u = 0 mod 2
    Construction: Lambda/(2L)^3 is the 12-dimensional F2-space of residues
    (w+v1+v2, w+v2, w+v1), w in W = L.s/2L, v1, v2 in V = L.sbar/2L.  We take the
    generators (w,w,w), (v,0,v), (v,v,0) for w, v running over the bases of
    L.s and L.sbar, plus 2*(L^3), and reduce to an HNF basis.  The result is
    verified independently afterwards (index in L^3 from the rank of Q, and
    the definitional membership test on every basis vector)."""
    Vrows = [LC.to_L(omul(b, SBAR_VEC)) for b in BL]     # L.sbar in L-coords
    Wrows = [LC.to_L(omul(b, S_VEC)) for b in BL]        # L.s   in L-coords
    gens = []
    Z8 = [0] * 8
    for w in Wrows:
        gens.append(w + w + w)
    for v in Vrows:
        gens.append(v + Z8 + v)
        gens.append(v + v + Z8)
    for k in range(24):
        e = [0] * 24
        e[k] = 2
        gens.append(e)
    BLam = hnf(gens, 24)
    assert len(BLam) == 24
    # GF(2) condition matrix.  For a sublattice A of L with 2L <= A, let
    # its dual-type functionals be the rows of an F2 matrix whose kernel is A/2L.
    Qv = f2_kernel_functionals(Vrows)       # 4 x 8 over F2, kernel = V
    Qw = f2_kernel_functionals(Wrows)       # 4 x 8 over F2, kernel = W
    Q = []
    for f in Qv:                             # x+y in V, x+z in V (y+z follows)
        Q.append(f + f + [0] * 8)
        Q.append(f + [0] * 8 + f)
        Q.append([0] * 8 + f + f)
    for f in Qw:                             # x+y+z in W
        Q.append(f + f + f)
    return {'BLam': BLam, 'Vrows': Vrows, 'Wrows': Wrows, 'Q': Q}


def f2_rowspace_basis(rows, n):
    """Reduced row echelon basis (list of 0/1 lists) of the F2 row space."""
    M = [[int(a) & 1 for a in r] for r in rows]
    basis = []
    pivcols = []
    for r in M:
        r = r[:]
        for b, pc in zip(basis, pivcols):
            if r[pc]:
                r = [x ^ y for x, y in zip(r, b)]
        if any(r):
            pc = r.index(1)
            # eliminate pc from existing basis rows
            for i in range(len(basis)):
                if basis[i][pc]:
                    basis[i] = [x ^ y for x, y in zip(basis[i], r)]
            basis.append(r)
            pivcols.append(pc)
    return basis, pivcols


def f2_kernel_functionals(rows):
    """Given integer rows spanning (mod 2) a subspace U of F2^8, return a list of
    functionals f (0/1 lists of length 8) whose common kernel is exactly U."""
    n = len(rows[0])
    basis, pivcols = f2_rowspace_basis(rows, n)
    free = [c for c in range(n) if c not in pivcols]
    # annihilator: vectors f with <f, b> = 0 for all basis rows b
    funcs = []
    for fc in free:
        f = [0] * n
        f[fc] = 1
        for b, pc in zip(basis, pivcols):
            if b[fc]:
                f[pc] = 1
        funcs.append(f)
    for f in funcs:
        for b in basis:
            assert sum(x * y for x, y in zip(f, b)) % 2 == 0
    return funcs


# --------------------------------------------------------------------------
# 5. Twists g (as integer matrices in L-coordinates) and twisted structure
#    constants
# --------------------------------------------------------------------------

def perm_apply(pi, x):
    """pi = (pi(1),...,pi(7)); returns pi(x) with pi(e_a) = e_{pi(a)}, pi(0)=0."""
    full = (0,) + tuple(pi)
    y = [None] * 8
    for a in range(8):
        y[full[a]] = x[a]
    return y


def perm_inverse(pi):
    full = (0,) + tuple(pi)
    inv = [0] * 8
    for a in range(8):
        inv[full[a]] = a
    return tuple(inv[1:])


def refl_apply(r, x):
    """Reflection in the root r: s_r(x) = x - <x,r> r   (N(r) = 2)."""
    t = ip(x, r)
    return [a - t * b for a, b in zip(x, r)]


def twist_matrix(gfun, BL, LC):
    """Integer matrix G (8x8, rows) with g(l_i) = sum_k G[i][k] l_k.
    Asserts g(L) <= L (integrality) and |det G| = 1 (so g(L) = L)."""
    G = [LC.to_L(gfun(b)) for b in BL]
    d = det_bareiss(G)
    assert abs(d) == 1, "g does not map L onto L"
    return G


def std_structure_constants(BL, LC):
    """C[i][j][k] (ints): l_i . l_j = sum_k C[i][j][k] l_k for the standard
    product.  LC.to_L asserts integrality, which proves L.L <= L."""
    C = np.zeros((8, 8, 8), dtype=np.int64)
    for i in range(8):
        for j in range(8):
            C[i, j, :] = LC.to_L(omul(BL[i], BL[j]))
    return C


def _check_int64_bound(*arrays_and_factors):
    """Guard: sum over terms of |products| must stay far below 2^63."""
    total = 1
    for f in arrays_and_factors:
        total *= f
    assert total < 2 ** 62, "int64 overflow risk"


def twisted_structure_constants(C, G):
    """C_g[i][j][:] = L-coords of g^{-1}( g(l_i) . g(l_j) ).

    In L-coordinates, g(l_i) = G[i] (row), the standard product of
    sum_a G[i,a] l_a and sum_b G[j,b] l_b is sum_{a,b} G[i,a] G[j,b] C[a,b,:],
    and applying g^{-1} multiplies a coordinate row by G^{-1} on the right."""
    Gm = np.array(G, dtype=np.int64)
    Ginv_int, Gden = int_inverse_scaled([[F(a) for a in r] for r in G])
    assert Gden == 1          # G unimodular -> inverse is integral
    Ginv = np.array(Ginv_int, dtype=np.int64)
    mg = int(np.abs(Gm).max()); mc = int(np.abs(C).max()); mi = int(np.abs(Ginv).max())
    _check_int64_bound(64 * mg * mg * mc, 8 * mi)
    T = np.einsum('ia,jb,abk->ijk', Gm, Gm, C)
    Cg = np.einsum('ijk,kl->ijl', T, Ginv)
    # sanity: integer by construction; verify g(l_i ._g l_j) = g(l_i).g(l_j)
    _check_int64_bound(int(np.abs(Cg).max()) + 1, 8 * mg)
    chk = np.einsum('ijk,kl->ijl', Cg, Gm)
    assert np.array_equal(chk, T)
    return Cg


# --------------------------------------------------------------------------
# 6. Triple product on L^3 coordinates, closure and span index
# --------------------------------------------------------------------------

def block_products(A, B, Cg):
    """For row-sets A (m x 8) and B (m x 8) of L-coordinates, return the
    m x m x 8 array of twisted products A_i ._g B_j (L-coordinates)."""
    ma = int(np.abs(A).max()) if A.size else 0
    mb = int(np.abs(B).max()) if B.size else 0
    mc = int(np.abs(Cg).max())
    _check_int64_bound(64 * ma * mb * mc)
    K = np.einsum('ia,abk->ibk', A, Cg)
    return np.einsum('ibk,jb->ijk', K, B)


def all_basis_products(BLam_np, Cg):
    """576 products b_i *_g b_j of the Lambda basis, in L^3-coordinates.
    Returns an array of shape (24, 24, 24): [i, j, :] = b_i * b_j."""
    X, Y, Z = BLam_np[:, 0:8], BLam_np[:, 8:16], BLam_np[:, 16:24]
    P1 = block_products(X, X, Cg) + block_products(Z, Y, Cg) + block_products(Y, Z, Cg)
    P2 = block_products(Y, Y, Cg) + block_products(X, Z, Cg) + block_products(Z, X, Cg)
    P3 = block_products(Z, Z, Cg) + block_products(Y, X, Cg) + block_products(X, Y, Cg)
    return np.concatenate([P1, P2, P3], axis=2)


def lambda_coordinates(P, Ninv2):
    """Given L^3-coordinate rows P (k x 24, int64) and N = 2*BLam^{-1} (int),
    return (ok_mask, coords2) where coords2 = P @ N (= 2 * Lambda-coordinates)
    and ok_mask[i] says whether row i lies in Lambda (coords2 all even)."""
    mp = int(np.abs(P).max()); mn = int(np.abs(Ninv2).max())
    _check_int64_bound(24 * mp * mn)
    coords2 = P @ Ninv2
    ok = np.all(coords2 % 2 == 0, axis=1)
    return ok, coords2


def cycle_type(pi):
    """Cycle type of a permutation of {1..7} given as (pi(1),...,pi(7)),
    as a tuple of cycle lengths >= 2 in decreasing order."""
    seen = set()
    lens = []
    for a in range(1, 8):
        if a in seen:
            continue
        L = 0
        b = a
        while b not in seen:
            seen.add(b)
            b = pi[b - 1]
            L += 1
        if L >= 2:
            lens.append(L)
    return tuple(sorted(lens, reverse=True))


def parity(pi):
    """0 for even, 1 for odd permutations."""
    return sum(L - 1 for L in cycle_type(pi)) % 2


CYCLE_TYPE_NAMES = {
    (2,): '2-cycle', (3,): '3-cycle', (2, 2): '(2,2)', (4,): '4-cycle',
    (3, 2): '(3,2)', (5,): '5-cycle', (2, 2, 2): '(2,2,2)', (3, 3): '(3,3)',
    (4, 2): '(4,2)', (6,): '6-cycle', (3, 2, 2): '(3,2,2)', (4, 3): '(4,3)',
    (5, 2): '(5,2)', (7,): '7-cycle', (): 'identity',
}


def fmt_vec(x):
    """Pretty-print a Fraction vector."""
    return '(' + ', '.join(str(t) for t in x) + ')'


# --------------------------------------------------------------------------
# 7. One-call setup used by the driver scripts
# --------------------------------------------------------------------------

def build_all():
    """Build and minimally verify all fixed data.  Returns a dict:
      BL      : Z-basis of L (8 Fraction rows, e-coordinates)
      LC      : LCoords converter
      C       : standard structure constants in L-coordinates (8x8x8 int64)
      Vrows, Wrows : L-coordinates of bases of L.sbar and L.s
      BLam    : 24x24 Python-int matrix, Lambda basis in L^3-coordinates
      BLam_np : same as int64 numpy array
      N2      : int64 numpy array, N2 = 2 * BLam^{-1} (integral)
      Q       : F2 condition matrix (12 independent rows among 16)
      BLam_e  : Lambda basis in e-coordinates (24 Fraction rows of length 24)
    """
    BL = l_basis()
    LC = LCoords(BL)
    C = std_structure_constants(BL, LC)
    lam = lambda_basis(BL, LC)
    BLam = lam['BLam']
    inv, den = int_inverse_scaled([[F(a) for a in r] for r in BLam])
    # 2 * BLam^{-1} must be integral because (2L)^3 <= Lambda
    assert all((2 * a) % den == 0 for row in inv for a in row)
    N2 = np.array([[2 * a // den for a in row] for row in inv], dtype=np.int64)
    BLam_np = np.array(BLam, dtype=np.int64)
    assert np.array_equal(BLam_np @ N2, 2 * np.eye(24, dtype=np.int64))
    BLam_e = [LC.to_e(r[0:8]) + LC.to_e(r[8:16]) + LC.to_e(r[16:24]) for r in BLam]
    return {'BL': BL, 'LC': LC, 'C': C, 'Vrows': lam['Vrows'], 'Wrows': lam['Wrows'],
            'BLam': BLam, 'BLam_np': BLam_np, 'N2': N2, 'Q': lam['Q'], 'BLam_e': BLam_e}


# --------------------------------------------------------------------------
# 8. HNF-free invariants of a finite-index sublattice S = span(rows) of Z^n
#    (used by 06_certificates.py).  With A = Z^n / S:
#      |A/pA| = p^(n - rank_p(rows))                       (any prime p)
#      |A/4A| = 4^n / |M|,  M = (S + 4Z^n)/4Z^n  (order computed below)
# --------------------------------------------------------------------------

def rank_mod_p(rows, p):
    """Rank over GF(p) of integer rows (Gaussian elimination mod p)."""
    basis = {}                                   # pivot column -> normalized row
    for r in rows:
        v = [int(a) % p for a in r]
        for c in sorted(basis):
            if v[c]:
                f = v[c]
                v = [(a - f * b) % p for a, b in zip(v, basis[c])]
        piv = next((c for c in range(len(v)) if v[c]), None)
        if piv is None:
            continue
        inv = pow(v[piv], p - 2, p)
        basis[piv] = [(a * inv) % p for a in v]
    return len(basis)


def log2_order_mod4_span(rows):
    """log2 of |M|, M = Z/4-submodule of (Z/4)^n generated by the integer rows.

    M -> F2^n (reduction mod 2) has image R = span(rows mod 2), dim r1, and
    kernel M /\\ 2(Z/4)^n = 2 * (R + Y), where Y is spanned by the vectors
    y(k) = (sum_{j in k} rows_j)/2 mod 2 for k running over a basis of the
    F2-relations k among the rows mod 2 (0/1 lift of k).  [Changing the lift
    of k changes y(k) by an element of R, so R + Y is well defined.]
    Hence log2|M| = r1 + dim(R + Y)."""
    n = len(rows[0])
    piv = {}                     # leading bit -> (mask of row, combination bitmask)
    relations = []
    for idx, r in enumerate(rows):
        m = 0
        for i, a in enumerate(r):
            if int(a) & 1:
                m |= 1 << i
        comb = 1 << idx
        while m:
            hb = m.bit_length() - 1
            if hb in piv:
                m ^= piv[hb][0]
                comb ^= piv[hb][1]
            else:
                piv[hb] = (m, comb)
                break
        if m == 0:
            relations.append(comb)
    r1 = len(piv)
    Rrows = [[(pm >> i) & 1 for i in range(n)] for pm, _ in piv.values()]
    Yrows = []
    for comb in relations:
        acc = [0] * n
        j = 0
        c = comb
        while c:
            if c & 1:
                acc = [a + int(b) for a, b in zip(acc, rows[j])]
            c >>= 1
            j += 1
        assert all(a % 2 == 0 for a in acc)
        Yrows.append([(a // 2) % 2 for a in acc])
    return r1 + f2_rank(Rrows + Yrows)
