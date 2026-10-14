"""
leech_twist_lib.py -- exact-arithmetic library for twisted Z3 triple products
on Wilson's Leech lattice.  Written from the definitions in main.tex (v8) of
"An order on the Leech lattice from a Z3-symmetric triple-octonion product";
it shares no code with any other implementation.

CONVENTIONS (each stated in main.tex v8):
  * Octonions O, basis e0..e7, e0 the identity.  A vector of R^8 is a list of
    8 numbers, index k holding the coefficient of e_k.
  * Fano triples (1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3): for a
    triple (a,b,c), e_a e_b = e_c, e_b e_c = e_a, e_c e_a = e_b, reversed
    products negative, e_i^2 = -e0 (i >= 1).
  * L = D8+ : all coordinates integers with even sum, or all half-odd-integers
    with ODD sum.  Roots have norm 2.
  * s = (-e0+e1+...+e7)/2 in L, sbar = (-e0-...-e7)/2 not in L.
  * L sbar = {l.sbar : l in L}, L s = {l.s : l in L}: LEFT multiples (lattice
    element is the left factor).
  * Lambda = {(x,y,z) in L^3 : x+y, x+z, y+z in L sbar ; x+y+z in L s}.
  * Twist: x ._g y := g^{-1}(g(x) . g(y)) for a linear g with g(L) = L.
    For a coordinate permutation pi of the imaginary units, pi(e_a) = e_{pi(a)}
    and x ._pi y := pi^{-1}(pi(x) . pi(y)).  Reflections s_r(x) = x - <x,r> r
    (r a root) are involutions, so for them g = g^{-1}.
  * Triple product (Definition 3.3): u=(x,y,z), v=(x',y',z'),
        P1 = x.x' + z.y' + y.z',  P2 = y.y' + x.z' + z.x',  P3 = z.z' + y.x' + x.y'
    with every block product equal to ._g.
  * S_g = Z-span{u *_g v : u,v in Lambda}; [Lambda : S_g] from the 576
    basis-pair products (bilinearity), by Hermite/Smith normal form.

REPRESENTATION.  Every lattice vector handled here (elements of L, L sbar,
L s, Lambda and their images under the twists used) has coordinates in
(1/2)Z.  Such vectors are stored as Python-int lists of DOUBLED coordinates
("X = 2x").  A twist g is stored as the integer matrix 4g (all elements of
W(E8) in these coordinates have entries in (1/4)Z), and g(x) is recovered as
(4g . X)/4 with an exactness assertion.  Every division in this file is
asserted to be exact, so no rounding can happen silently.  Python ints are
unbounded, so there is no overflow anywhere.
"""
from fractions import Fraction as Fr
import itertools

# --------------------------------------------------------------------------
# octonion multiplication table
# --------------------------------------------------------------------------
FANO = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]


def _build_table():
    """MT[a][b] = (sign, c) meaning e_a e_b = sign * e_c."""
    MT = [[None] * 8 for _ in range(8)]
    for a in range(8):
        MT[0][a] = (1, a)
        MT[a][0] = (1, a)
    for a in range(1, 8):
        MT[a][a] = (-1, 0)
    for (a, b, c) in FANO:
        for (p, q, r) in ((a, b, c), (b, c, a), (c, a, b)):
            assert MT[p][q] is None and MT[q][p] is None, "Fano lines overlap"
            MT[p][q] = (1, r)
            MT[q][p] = (-1, r)
    assert all(MT[a][b] is not None for a in range(8) for b in range(8))
    return MT


MT = _build_table()
# flat list of (a, b, sign, c) for fast loops
_TERMS = [(a, b, MT[a][b][0], MT[a][b][1]) for a in range(8) for b in range(8)]


def omul(x, y):
    """Octonion product of two coefficient lists (any exact number type)."""
    z = [0] * 8
    for (a, b, sg, c) in _TERMS:
        xa = x[a]
        if xa:
            yb = y[b]
            if yb:
                z[c] += sg * xa * yb
    return z


def omul2(X, Y):
    """Product in doubled coordinates: X = 2x, Y = 2y  ->  2(x.y).
    Requires x.y in (1/2)Z^8 (true for x, y in L); asserted."""
    Z = omul(X, Y)  # = 4 x.y
    assert all(t % 2 == 0 for t in Z), "product leaves (1/2)Z^8"
    return [t // 2 for t in Z]


def nrm(x):
    """Squared Euclidean norm of a coefficient list."""
    return sum(t * t for t in x)


def nrm2(X):
    """Norm of x from doubled coordinates X = 2x (returned exactly)."""
    return Fr(sum(t * t for t in X), 4)


def unit(k, scale=1):
    v = [0] * 8
    v[k] = scale
    return v


# --------------------------------------------------------------------------
# small exact integer linear algebra
# --------------------------------------------------------------------------
def det_bareiss(M):
    """Exact determinant of a square integer matrix (fraction-free Bareiss)."""
    A = [list(r) for r in M]
    n = len(A)
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    A[k], A[i] = A[i], A[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def hnf_rows(rows):
    """Row-style Hermite normal form of an integer matrix given as a list of
    rows; returns the nonzero rows (a Z-basis of the row lattice), upper
    triangular with positive pivots and reduced entries above each pivot."""
    A = [list(r) for r in rows if any(r)]
    if not A:
        return []
    m, n = len(A), len(A[0])
    pr = 0
    for col in range(n):
        if pr == m:
            break
        while True:
            nz = [i for i in range(pr, m) if A[i][col] != 0]
            if not nz:
                break
            imin = min(nz, key=lambda i: abs(A[i][col]))
            A[pr], A[imin] = A[imin], A[pr]
            p = A[pr][col]
            clean = True
            for i in range(pr + 1, m):
                if A[i][col]:
                    q = A[i][col] // p
                    if q:
                        Ai, Ap = A[i], A[pr]
                        A[i] = [a - q * b for a, b in zip(Ai, Ap)]
                    if A[i][col]:
                        clean = False
            if clean:
                break
        if any(A[i][col] for i in range(pr, m)):
            if A[pr][col] < 0:
                A[pr] = [-a for a in A[pr]]
            p = A[pr][col]
            for i in range(pr):
                q = A[i][col] // p
                if q:
                    A[i] = [a - q * b for a, b in zip(A[i], A[pr])]
            pr += 1
    return [r for r in A[:pr]]


def smith_diag(M):
    """Elementary divisors (diagonal of the Smith normal form) of an integer
    matrix, sorted increasingly, zeros last."""
    A = [list(r) for r in M]
    m, n = len(A), len(A[0])
    diag = []
    t = 0
    while t < min(m, n):
        # find a nonzero entry of minimal absolute value in the submatrix
        best = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] and (best is None or abs(A[i][j]) < abs(A[best[0]][best[1]])):
                    best = (i, j)
        if best is None:
            break
        i, j = best
        A[t], A[i] = A[i], A[t]
        for r in A:
            r[t], r[j] = r[j], r[t]
        while True:
            p = A[t][t]
            changed = False
            for i in range(t + 1, m):
                if A[i][t]:
                    q = A[i][t] // p
                    A[i] = [a - q * b for a, b in zip(A[i], A[t])]
                    if A[i][t]:
                        changed = True
            for j in range(t + 1, n):
                if A[t][j]:
                    q = A[t][j] // p
                    for r in A:
                        r[j] -= q * r[t]
                    if A[t][j]:
                        changed = True
            if changed:
                # move the smallest nonzero entry of row/col t to the pivot
                cand = [(abs(A[i][t]), i, t) for i in range(t, m) if A[i][t]] + \
                       [(abs(A[t][j]), t, j) for j in range(t, n) if A[t][j]]
                _, i, j = min(cand)
                A[t], A[i] = A[i], A[t]
                for r in A:
                    r[t], r[j] = r[j], r[t]
                continue
            # divisibility condition: pivot divides the rest of the submatrix
            bad = None
            for i in range(t + 1, m):
                for j in range(t + 1, n):
                    if A[i][j] % p:
                        bad = i
                        break
                if bad is not None:
                    break
            if bad is None:
                break
            A[t] = [a + b for a, b in zip(A[t], A[bad])]
        diag.append(abs(A[t][t]))
        t += 1
    diag += [0] * (min(m, n) - len(diag))
    nz = sorted(d for d in diag if d)
    return nz + [d for d in diag if d == 0]


def inv_frac(M):
    """Exact inverse of a square rational matrix (Gauss-Jordan on Fractions)."""
    n = len(M)
    A = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(i for i in range(c, n) if A[i][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    return [row[n:] for row in A]


def matmul(A, B):
    return [[sum(a * b for a, b in zip(row, col)) for col in zip(*B)] for row in A]


def transpose(A):
    return [list(r) for r in zip(*A)]


def rank_mod_p(rows, p):
    """Rank over F_p of an integer matrix given as rows."""
    A = [[x % p for x in r] for r in rows]
    if not A:
        return 0
    n = len(A[0])
    rk = 0
    for c in range(n):
        piv = next((i for i in range(rk, len(A)) if A[i][c]), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], -1, p)
        A[rk] = [(x * inv) % p for x in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c]:
                f = A[i][c]
                A[i] = [(a - f * b) % p for a, b in zip(A[i], A[rk])]
        rk += 1
    return rk


def kernel_mod2(rows):
    """Basis (0/1 lists) of the right kernel over F_2 of the matrix `rows`."""
    A = [[x % 2 for x in r] for r in rows]
    n = len(A[0])
    pivcols = []
    rk = 0
    for c in range(n):
        piv = next((i for i in range(rk, len(A)) if A[i][c]), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        for i in range(len(A)):
            if i != rk and A[i][c]:
                A[i] = [(a + b) % 2 for a, b in zip(A[i], A[rk])]
        pivcols.append(c)
        rk += 1
    free = [c for c in range(n) if c not in pivcols]
    basis = []
    for f in free:
        v = [0] * n
        v[f] = 1
        for i, pc in enumerate(pivcols):
            v[pc] = A[i][f] % 2
        basis.append(v)
    return basis


# --------------------------------------------------------------------------
# lattices given by a basis in doubled coordinates
# --------------------------------------------------------------------------
class Lattice:
    """Full-rank lattice in R^n given by basis rows in doubled coordinates."""

    def __init__(self, rows):
        self.rows = [list(r) for r in rows]
        self.n = len(rows[0])
        assert len(self.rows) == self.n
        inv = inv_frac(transpose(self.rows))  # coords c = inv . X
        den = 1
        for r in inv:
            for x in r:
                den = den * x.denominator // _gcd(den, x.denominator)
        self.den = den
        self.K = [[int(x * den) for x in r] for r in inv]  # integer matrix

    def coords(self, X):
        """Integer coordinates of X (doubled coords) w.r.t. the basis, or None."""
        out = []
        for r in self.K:
            t = sum(a * b for a, b in zip(r, X))
            if t % self.den:
                return None
            out.append(t // self.den)
        return out

    def __contains__(self, X):
        return self.coords(X) is not None

    def gram(self):
        """Gram matrix of the basis (exact Fractions; doubled coords /4)."""
        return [[Fr(sum(a * b for a, b in zip(u, v)), 4) for v in self.rows] for u in self.rows]

    def index_in(self, other):
        """[other : self] for self a sublattice of other."""
        C = []
        for r in self.rows:
            c = other.coords(r)
            assert c is not None, "not a sublattice"
            C.append(c)
        return abs(det_bareiss(C))


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return abs(a)


def in_L2(X):
    """D8+ membership from doubled coordinates X = 2x."""
    if all(t % 2 == 0 for t in X):      # integer coordinates
        return sum(X) % 4 == 0          # even coordinate sum
    if all(t % 2 == 1 for t in X):      # half-odd-integer coordinates (Python: -1 % 2 == 1)
        return sum(X) % 4 == 2          # odd coordinate sum
    return False


# --------------------------------------------------------------------------
# the fixed objects of the paper
# --------------------------------------------------------------------------
S2 = [-1, 1, 1, 1, 1, 1, 1, 1]          # 2s
SBAR2 = [-1] * 8                         # 2 sbar
# basis of L from Appendix A: L1 = s, L2 = e0+e1, L3 = e0-e1, L_{k+2} = e0+e_k (k=2..6)
L_BASIS2 = [S2,
            [2, 2, 0, 0, 0, 0, 0, 0], [2, -2, 0, 0, 0, 0, 0, 0],
            [2, 0, 2, 0, 0, 0, 0, 0], [2, 0, 0, 2, 0, 0, 0, 0],
            [2, 0, 0, 0, 2, 0, 0, 0], [2, 0, 0, 0, 0, 2, 0, 0],
            [2, 0, 0, 0, 0, 0, 2, 0]]


def left_multiples(basis2, a2):
    """Basis of {l . a : l in lattice}: products L_k . a (lattice element LEFT)."""
    return [omul2(b, a2) for b in basis2]


def right_multiples(basis2, a2):
    """Basis of {a . l} (used only for the convention-variant check)."""
    return [omul2(a2, b) for b in basis2]


L_LAT = Lattice(L_BASIS2)
LSBAR_LAT = Lattice(left_multiples(L_BASIS2, SBAR2))
LS_LAT = Lattice(left_multiples(L_BASIS2, S2))


def leech_basis(A_lat=None, B_lat=None):
    """Z-basis (doubled coordinates, 24-lists) of
       {(x,y,z) in L^3 : x+y, x+z, y+z in A ; x+y+z in B}
    with defaults A = L sbar, B = L s (Wilson's Lambda).  Uses 2L^3 inside
    the lattice: Lambda/2L^3 is the F2-kernel of the conditions, written in
    L-coordinates; the result is the HNF basis in L^3-coordinates mapped to R^24."""
    A_lat = A_lat or LSBAR_LAT
    B_lat = B_lat or LS_LAT
    # 2 * (coordinates of L-basis vectors in A, B): integral because 2L in A, B
    def two_inv(lat):
        M = []
        for k in range(8):
            c = lat.coords([2 * t for t in L_BASIS2[k]])
            assert c is not None, "2L not inside the sublattice"
            M.append(c)
        return transpose(M)   # column k = coords of 2 L_k
    QA, QB = two_inv(A_lat), two_inv(B_lat)
    Z = [[0] * 8 for _ in range(8)]
    def hcat(*Ms):
        return [sum((M[i] for M in Ms), []) for i in range(8)]
    K = hcat(QA, QA, Z) + hcat(QA, Z, QA) + hcat(QB, QB, QB)   # 24 x 24, conditions mod 2
    ker = kernel_mod2(K)
    gens = ker + [[2 * int(i == j) for j in range(24)] for i in range(24)]
    H = hnf_rows(gens)
    assert len(H) == 24
    basis = []
    for c in H:
        v = []
        for blk in range(3):
            cc = c[8 * blk: 8 * blk + 8]
            v += [sum(cc[k] * L_BASIS2[k][i] for k in range(8)) for i in range(8)]
        basis.append(v)
    return basis, len(ker)


LAMBDA_BASIS2, _KERDIM = leech_basis()
LAMBDA_LAT = Lattice(LAMBDA_BASIS2)


def in_Lambda2(U, A_lat=None, B_lat=None):
    """Wilson's three conditions, tested directly (doubled coordinates)."""
    A_lat = A_lat or LSBAR_LAT
    B_lat = B_lat or LS_LAT
    x, y, z = U[0:8], U[8:16], U[16:24]
    if not (in_L2(x) and in_L2(y) and in_L2(z)):
        return False
    add = lambda a, b: [p + q for p, q in zip(a, b)]
    return (add(x, y) in A_lat and add(x, z) in A_lat and add(y, z) in A_lat
            and add(add(x, y), z) in B_lat)


# --------------------------------------------------------------------------
# twists
# --------------------------------------------------------------------------
class Twist:
    """A linear map g of R^8 with g(L) = L, stored exactly (Fractions) and as
    the integer matrices 4g, 4g^{-1}.  Twisted product x ._g y = g^{-1}(g x . g y)."""

    def __init__(self, G, name="", require_L=True):
        """require_L=True: assert g(L) = L (all twists of the paper and of
        claim C3).  require_L=False (used only for claim C4's signed
        permutations, half of which map L onto the other D8+ coset L'):
        assert instead that the twisted product maps L x L into L."""
        self.name = name
        self.G = [[Fr(x) for x in r] for r in G]
        self.Gi = inv_frac(self.G)
        self.G4 = [[int(4 * x) for x in r] for r in self.G]
        self.Gi4 = [[int(4 * x) for x in r] for r in self.Gi]
        assert all(4 * x == int(4 * x) for r in self.G + self.Gi for x in r), "entries not in (1/4)Z"
        d = det_frac(self.G)
        assert abs(d) == 1
        self.det = d
        self.maps_L_to_L = all(in_L2(self.apply(b)) for b in L_BASIS2)
        if require_L:
            # g(L) = L: images of the L basis lie in L and det g = +-1
            assert self.maps_L_to_L, "g does not map L into L"
        else:
            assert all(in_L2(self.tmul2(a, b)) for a in L_BASIS2 for b in L_BASIS2), \
                "twisted product does not preserve L"

    @staticmethod
    def _app(M4, X):
        out = []
        for r in M4:
            t = sum(a * b for a, b in zip(r, X))
            assert t % 4 == 0, "image leaves (1/2)Z^8"
            out.append(t // 4)
        return out

    def apply(self, X):
        return self._app(self.G4, X)

    def apply_inv(self, X):
        return self._app(self.Gi4, X)

    def tmul2(self, X, Y):
        """Twisted product in doubled coordinates."""
        return self.apply_inv(omul2(self.apply(X), self.apply(Y)))

    def tmul_frac(self, x, y):
        """Twisted product on exact Fraction vectors (any element of O)."""
        gx = [sum(a * b for a, b in zip(r, x)) for r in self.G]
        gy = [sum(a * b for a, b in zip(r, y)) for r in self.G]
        p = omul(gx, gy)
        return [sum(a * b for a, b in zip(r, p)) for r in self.Gi]

    def identity_element(self):
        """g^{-1}(e0): the two-sided identity of ._g."""
        return [r[0] for r in self.Gi]


def det_frac(M):
    n = len(M)
    A = [[Fr(x) for x in r] for r in M]
    d = Fr(1)
    for c in range(n):
        p = next((i for i in range(c, n) if A[i][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for i in range(c + 1, n):
            f = A[i][c] / A[c][c]
            if f:
                A[i] = [a - f * b for a, b in zip(A[i], A[c])]
    return d


def reflection(r2):
    """Matrix of s_r(x) = x - <x,r> r for a root r given in doubled coords."""
    r = [Fr(t, 2) for t in r2]
    assert sum(t * t for t in r) == 2, "not a root"
    return [[Fr(int(i == j)) - r[i] * r[j] for j in range(8)] for i in range(8)]


def perm_matrix(pi):
    """pi: dict on {1..7} (pi(e_a) = e_{pi(a)}), e0 fixed.  Matrix P with P e_a = e_{pi(a)}."""
    P = [[0] * 8 for _ in range(8)]
    P[0][0] = 1
    for a in range(1, 8):
        P[pi[a]][a] = 1
    return P


def signed_perm_matrix(pi, eps):
    """e_a -> eps[a] e_{pi(a)} for a = 1..7, e0 fixed."""
    P = [[0] * 8 for _ in range(8)]
    P[0][0] = 1
    for a in range(1, 8):
        P[pi[a]][a] = eps[a]
    return P


# --------------------------------------------------------------------------
# triple product, closure, span index
# --------------------------------------------------------------------------
def tprod2(tw, U, V):
    """u *_g v in doubled coordinates (Definition 3.3 routing)."""
    x, y, z = U[0:8], U[8:16], U[16:24]
    xp, yp, zp = V[0:8], V[8:16], V[16:24]
    m = tw.tmul2
    add3 = lambda a, b, c: [p + q + r for p, q, r in zip(a, b, c)]
    P1 = add3(m(x, xp), m(z, yp), m(y, zp))
    P2 = add3(m(y, yp), m(x, zp), m(z, xp))
    P3 = add3(m(z, zp), m(y, xp), m(x, yp))
    return P1 + P2 + P3


def basis_products(tw, basis=None):
    """All 576 products b_i *_g b_j (doubled coordinates), row-major in (i, j)."""
    basis = basis or LAMBDA_BASIS2
    return [tprod2(tw, bi, bj) for bi in basis for bj in basis]


def closure_and_index(tw, basis=None, lat=None):
    """Return dict: closes, number of failing pairs, index, elementary divisors,
    coordinate rows of the products (if closing)."""
    basis = basis or LAMBDA_BASIS2
    lat = lat or LAMBDA_LAT
    prods = basis_products(tw, basis)
    coords, bad = [], 0
    for P in prods:
        if not in_Lambda2(P):
            bad += 1
            coords.append(None)
            continue
        c = lat.coords(P)
        assert c is not None, "Wilson test and basis coordinates disagree"
        coords.append(c)
    out = {"closes": bad == 0, "failing_pairs": bad}
    if bad == 0:
        H = hnf_rows(coords)
        out["rank"] = len(H)
        index = 1
        for i in range(len(H)):
            index *= H[i][i] if len(H) == 24 else 0
        out["index"] = abs(index) if len(H) == 24 else 0
        out["elementary_divisors"] = [d for d in smith_diag(H) if d != 1]
        out["coords"] = coords
        out["hnf"] = H
    return out


def sub_index(rows_doubled, ambient_lat):
    """[ambient : Z-span(rows)] (0 if not of full rank)."""
    C = [ambient_lat.coords(r) for r in rows_doubled]
    assert all(c is not None for c in C)
    H = hnf_rows(C)
    if len(H) < ambient_lat.n:
        return 0
    d = 1
    for i in range(len(H)):
        d *= H[i][i]
    return abs(d)
