"""
leech_core.py -- shared exact-arithmetic definitions for the structural review
of the twisted Z3 triple product on Wilson's Leech lattice
(paper: "An order on the Leech lattice from a Z3-symmetric triple-octonion
product", v8).  Written from the paper's definitions only.

CONVENTIONS (stated once; every script imports them from here)
---------------------------------------------------------------
1. Octonions.  Basis e0..e7, e0 the identity.  Fano triples
       (1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3).
   For a triple (a,b,c):  e_a e_b = +e_c, e_b e_c = +e_a, e_c e_a = +e_b
   (cyclic = positive), e_b e_a = -e_c etc. (anticyclic = negative);
   e_i e_i = -e0 for i >= 1.  An element x = sum_i x_i e_i is a length-8
   list of ints/Fractions ("e-coordinates").

2. L = D8+ : all-integer vectors with EVEN coordinate sum, together with
   all-half-integer vectors with ODD coordinate sum.  N(x) = sum x_i^2
   (roots have N = 2).  Z-basis (the paper's Appendix A basis):
       L1 = s = (-1,1,1,1,1,1,1,1)/2,  L2 = e0+e1,  L3 = e0-e1,
       L_{k+2} = e0+e_k  (k = 2..6).
   "L-coordinates" of x in L: the integer row vector c with x = sum_k c_k L_k.

3. s = (-e0+e1+...+e7)/2 (in L), sbar = (-e0-e1-...-e7)/2 (not in L).
   Ls = {l.s : l in L}, Lsbar = {l.sbar : l in L}  (LEFT multiples, standard
   product).

4. Wilson's Leech lattice
       Lambda = {(x,y,z) in L^3 : x+y, x+z, y+z in Lsbar;  x+y+z in Ls}.
   Elements of L^3 are stored as 24 integers: L-coords of x, then y, then z.

5. Twists.  A linear map g of R^8 is an 8x8 rational matrix G acting on
   COLUMN vectors of e-coordinates: g(x) = G x.
       twisted product:   x ._g y := g^{-1}( g(x) . g(y) ).
   * permutation pi of {1..7} (pi(0) = 0):  g_pi(e_a) = e_{pi(a)},
     i.e. G[pi(a), a] = 1.  The paper's census writes x ._pi y =
     pi^{-1}(pi(x).pi(y)), i.e. g = g_pi ("convention A").  Convention B
     (g = g_pi^{-1}) is also reported where it matters.
   * root reflection: s_r(x) = x - <x,r> r, G = I - r r^T (<r,r> = 2);
     s_r is an involution, so g^{-1} = g.
   * signed permutation: g(e_a) = eps_a e_{pi(a)}, eps_0 = +1.

6. Triple product with Z3 routing (Definition 3.3), every block product ._g:
       P1 = x.x' + z.y' + y.z',  P2 = y.y' + x.z' + z.x',
       P3 = z.z' + y.x' + x.y'.

7. Mod-2 quotient A = L/2L = F_2^8 through L-coordinates mod 2.  An element
   of A is an int 0..255 whose bit k is c_k mod 2 (bit 0 <-> L1, ..., bit 7
   <-> L8).  A subspace is stored as a frozenset of its elements.

Arithmetic is exact: Python ints / fractions.Fraction everywhere, numpy int64
only where an explicit bound is asserted (see check_int64_bound).
"""
from fractions import Fraction as Fr
import itertools
import numpy as np

# ---------------------------------------------------------------------------
# 1. Octonion multiplication
# ---------------------------------------------------------------------------
FANO = [(1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 7), (5, 6, 1), (6, 7, 2), (7, 1, 3)]


def _build_mult():
    M = [[None] * 8 for _ in range(8)]
    for a in range(8):
        M[0][a] = (1, a)
        M[a][0] = (1, a)
    for a in range(1, 8):
        M[a][a] = (-1, 0)
    for (a, b, c) in FANO:
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            M[x][y] = (1, z)
            M[y][x] = (-1, z)
    assert all(M[a][b] is not None for a in range(8) for b in range(8))
    return M


MULT = _build_mult()                      # MULT[a][b] = (sign, c): e_a e_b = sign e_c
OCT = np.zeros((8, 8, 8), dtype=np.int64)  # OCT[a,b,c] = coefficient of e_c in e_a e_b
for _a in range(8):
    for _b in range(8):
        _s, _c = MULT[_a][_b]
        OCT[_a, _b, _c] = _s


def omul(x, y):
    """Standard octonion product in e-coordinates (exact for ints/Fractions)."""
    z = [0] * 8
    for a in range(8):
        if x[a] == 0:
            continue
        for b in range(8):
            if y[b] == 0:
                continue
            sg, c = MULT[a][b]
            z[c] += sg * x[a] * y[b]
    return z


def conj(x):
    return [x[0]] + [-t for t in x[1:]]


def norm(x):
    return sum(t * t for t in x)


def inner(x, y):
    return sum(a * b for a, b in zip(x, y))


# ---------------------------------------------------------------------------
# small exact linear algebra over Q (Fractions)
# ---------------------------------------------------------------------------
def mat_inv(M):
    """Inverse of a square rational matrix (list of lists), Gauss-Jordan."""
    n = len(M)
    A = [[Fr(M[i][j]) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [t / pv for t in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [a - f * b for a, b in zip(A[r], A[col])]
    return [row[n:] for row in A]


def mat_det(M):
    """Determinant of a square rational matrix (Fractions, exact)."""
    n = len(M)
    A = [[Fr(t) for t in row] for row in M]
    det = Fr(1)
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r][col] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
            det = -det
        det *= A[col][col]
        for r in range(col + 1, n):
            if A[r][col] != 0:
                f = A[r][col] / A[col][col]
                A[r] = [a - f * b for a, b in zip(A[r], A[col])]
    return det


def matvec(G, x):
    """G x for a matrix G (list of rows) and a column vector x (list)."""
    return [sum(G[i][j] * x[j] for j in range(len(x))) for i in range(len(G))]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]


def vecmat(x, M):
    """Row vector x times matrix M."""
    return [sum(x[k] * M[k][j] for k in range(len(x))) for j in range(len(M[0]))]


def is_integral(v):
    return all(Fr(t).denominator == 1 for t in v)


# ---------------------------------------------------------------------------
# exact integer HNF (Python ints; used for small constructions)
# ---------------------------------------------------------------------------
def hnf(rows):
    """Row-style Hermite normal form of an integer matrix (list of lists).
    Returns the list of nonzero rows: upper echelon, positive pivots, entries
    above each pivot reduced into [0, pivot).  Exact (Python ints)."""
    A = [list(map(int, r)) for r in rows if any(r)]
    if not A:
        return []
    m = len(A[0])
    out = []
    r0 = 0
    for c in range(m):
        # bring a gcd into position r0 by Euclid on the column
        while True:
            nz = [i for i in range(r0, len(A)) if A[i][c] != 0]
            if not nz:
                break
            p = min(nz, key=lambda i: abs(A[i][c]))
            A[r0], A[p] = A[p], A[r0]
            done = True
            for i in range(r0 + 1, len(A)):
                if A[i][c] != 0:
                    q = A[i][c] // A[r0][c]
                    A[i] = [a - q * b for a, b in zip(A[i], A[r0])]
                    if A[i][c] != 0:
                        done = False
            if done:
                break
        if r0 < len(A) and A[r0][c] != 0:
            if A[r0][c] < 0:
                A[r0] = [-a for a in A[r0]]
            for i in range(r0):
                q = A[i][c] // A[r0][c]
                if q:
                    A[i] = [a - q * b for a, b in zip(A[i], A[r0])]
            r0 += 1
        A = A[:r0] + [r for r in A[r0:] if any(r)]
    return A[:r0]


def hnf_det(rows, m):
    """Index in Z^m of the lattice spanned by integer rows (0 if not full rank)."""
    H = hnf(rows)
    if len(H) < m:
        return 0
    d = 1
    for i in range(m):
        d *= H[i][i]
    return d


# ---------------------------------------------------------------------------
# fast exact index computation (numpy int64 with an explicit overflow guard,
# falling back to the Python-int HNF if the guard trips)
# ---------------------------------------------------------------------------
INT64_GUARD = 2 ** 31   # entries kept below 2^31, so q*row products stay below 2^62


def check_int64_bound(arr, bound=INT64_GUARD):
    mx = int(np.abs(arr).max()) if arr.size else 0
    assert mx < bound, f"int64 guard tripped: max |entry| = {mx}"
    return mx


def span_index(M):
    """[Z^m : span(rows of M)] for an integer numpy array M (n x m).
    Returns (index, pivots) with pivots the diagonal of the row-HNF, or
    (0, None) if the rows do not span a full-rank lattice.  Column-by-column
    Euclid, vectorized over rows; exact as long as the int64 guard holds,
    otherwise recomputed with the Python-int hnf()."""
    R = np.array(M, dtype=np.int64)
    R = R[np.any(R != 0, axis=1)]
    n, m = R.shape if R.size else (0, M.shape[1])
    pivots = []
    try:
        for c in range(m):
            while True:
                nz = np.nonzero(R[:, c])[0]
                if nz.size == 0:
                    return 0, None
                if nz.size == 1:
                    break
                vals = R[nz, c]
                p = nz[np.argmin(np.abs(vals))]
                q = R[:, c] // R[p, c]
                q[p] = 0
                R = R - q[:, None] * R[p][None, :]
                check_int64_bound(R)
            p = nz[0]
            pivots.append(abs(int(R[p, c])))
            R = np.delete(R, p, axis=0)
            R = R[np.any(R != 0, axis=1)] if R.size else R
        idx = 1
        for d in pivots:
            idx *= d
        return idx, pivots
    except AssertionError:
        H = hnf(M.tolist())
        if len(H) < m:
            return 0, None
        piv = [H[i][i] for i in range(m)]
        idx = 1
        for d in piv:
            idx *= d
        return idx, piv


# ---------------------------------------------------------------------------
# 2. The lattice L and coordinates
# ---------------------------------------------------------------------------
H2 = Fr(1, 2)
L_BASIS = [
    [-H2, H2, H2, H2, H2, H2, H2, H2],      # L1 = s
    [1, 1, 0, 0, 0, 0, 0, 0],               # L2 = e0+e1
    [1, -1, 0, 0, 0, 0, 0, 0],              # L3 = e0-e1
    [1, 0, 1, 0, 0, 0, 0, 0],               # L4 = e0+e2
    [1, 0, 0, 1, 0, 0, 0, 0],               # L5 = e0+e3
    [1, 0, 0, 0, 1, 0, 0, 0],               # L6 = e0+e4
    [1, 0, 0, 0, 0, 1, 0, 0],               # L7 = e0+e5
    [1, 0, 0, 0, 0, 0, 1, 0],               # L8 = e0+e6
]
L_BASIS = [[Fr(t) for t in row] for row in L_BASIS]
L_BASIS_INV = mat_inv(L_BASIS)
GRAM_L = [[int(inner(a, b)) for b in L_BASIS] for a in L_BASIS]


def in_L_e(x):
    """Membership in D8+ (odd half-integer coset) straight from the definition."""
    x = [Fr(t) for t in x]
    if all(t.denominator == 1 for t in x):
        return sum(x) % 2 == 0
    if all((2 * t).denominator == 1 and t.denominator == 2 for t in x):
        return sum(x).denominator == 1 and int(sum(x)) % 2 == 1
    return False


def to_L(x):
    """e-coordinates -> L-coordinates (Fractions; integral iff x in L)."""
    return vecmat([Fr(t) for t in x], L_BASIS_INV)


def to_L_int(x):
    c = to_L(x)
    assert is_integral(c), f"not in L: {x}"
    return [int(t) for t in c]


def from_L(c):
    return vecmat(list(c), L_BASIS)


S_VEC = [Fr(-1, 2)] + [Fr(1, 2)] * 7          # s
SBAR_VEC = [Fr(-1, 2)] * 8                     # sbar = conj(s)


def roots_of_L():
    """The 240 roots of L (e-coordinates, Fractions)."""
    out = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [Fr(0)] * 8
                v[i], v[j] = Fr(si), Fr(sj)
                out.append(v)
    for signs in itertools.product((1, -1), repeat=8):
        if sum(1 for t in signs if t < 0) % 2 == 1:      # odd sum <=> odd # of minus signs
            out.append([Fr(t, 2) for t in signs])
    assert len(out) == 240 and all(norm(v) == 2 and in_L_e(v) for v in out)
    return out


# ---------------------------------------------------------------------------
# structure tensors in L-coordinates
# ---------------------------------------------------------------------------
def structure_tensor_L(G=None):
    """T[i,j,:] = L-coordinates of L_i ._g L_j  (g = identity if G is None).
    Returned as an (8,8,8) object array of Fractions."""
    if G is None:
        Gm = Gi = None
    else:
        Gm = [[Fr(t) for t in row] for row in G]
        Gi = mat_inv(Gm)
    T = np.empty((8, 8, 8), dtype=object)
    imgs = [L_BASIS[i] if Gm is None else matvec(Gm, L_BASIS[i]) for i in range(8)]
    for i in range(8):
        for j in range(8):
            p = omul(imgs[i], imgs[j])
            if Gi is not None:
                p = matvec(Gi, p)
            T[i, j, :] = to_L(p)
    return T


def tensor_int(T):
    """Convert an object tensor of Fractions to int64, asserting integrality."""
    flat = [Fr(t) for t in T.flat]
    assert all(t.denominator == 1 for t in flat), "tensor is not integral"
    return np.array([int(t) for t in flat], dtype=np.int64).reshape(T.shape)


T_STD = tensor_int(structure_tensor_L())        # standard product on L (L.L in L asserted)


def lmul(c1, c2, T=T_STD):
    """Product of two elements of L given in L-coordinates (int64 vectors)."""
    return np.einsum('i,j,ijk->k', np.asarray(c1, dtype=np.int64),
                     np.asarray(c2, dtype=np.int64), T)


# ---------------------------------------------------------------------------
# twists as L-coordinate matrices
# ---------------------------------------------------------------------------
def perm_matrix_e(pi):
    """pi: dict or list with pi[a] for a in 1..7 (pi[0] = 0 implied).
    Returns G with g(e_a) = e_{pi(a)}, i.e. G[pi(a)][a] = 1."""
    G = [[0] * 8 for _ in range(8)]
    for a in range(8):
        b = 0 if a == 0 else pi[a]
        G[b][a] = 1
    return G


def signed_perm_matrix_e(pi, eps):
    """g(e_a) = eps[a] e_{pi(a)}; eps[0] must be +1."""
    G = [[0] * 8 for _ in range(8)]
    for a in range(8):
        b = 0 if a == 0 else pi[a]
        G[b][a] = eps[a]
    return G


def reflection_matrix_e(r):
    """s_r = I - r r^T for a root r (<r,r> = 2)."""
    assert norm(r) == 2
    return [[Fr(int(i == j)) - r[i] * r[j] for j in range(8)] for i in range(8)]


def L_matrix(G):
    """Row-convention L-coordinate matrix of g: L-coords(g(x)) = c(x) @ GL.
    Asserts integrality (g(L) in L) and returns a list of int rows."""
    rows = []
    for i in range(8):
        img = matvec(G, L_BASIS[i])
        rows.append(to_L_int(img))
    return rows


# ---------------------------------------------------------------------------
# 7. mod-2 quotient A = L/2L
# ---------------------------------------------------------------------------
def to_bits(c):
    v = 0
    for k, t in enumerate(c):
        if int(t) % 2:
            v |= 1 << k
    return v


def from_bits(v):
    return [(v >> k) & 1 for k in range(8)]


def span_bits(gens):
    """F_2-span of a list of bitmasks (as a frozenset)."""
    S = {0}
    for g in gens:
        if g not in S:
            S |= {x ^ g for x in S}
    return frozenset(S)


def basis_bits(S):
    """A basis (list) of an F_2-subspace given as a set of bitmasks."""
    basis, cur = [], {0}
    for v in sorted(S):
        if v not in cur:
            basis.append(v)
            cur |= {x ^ v for x in cur}
    return basis


def dim_bits(S):
    return len(S).bit_length() - 1


def mul_table_mod2(T):
    """256x256 table of the product induced on A = L/2L by an integral
    structure tensor T (bilinear over F_2)."""
    T2 = [[to_bits(T[i, j]) for j in range(8)] for i in range(8)]
    # rowprod[i][v] = L_i * v  (mod 2), by linearity in v
    rowprod = [[0] * 256 for _ in range(8)]
    for i in range(8):
        for v in range(256):
            acc = 0
            for j in range(8):
                if (v >> j) & 1:
                    acc ^= T2[i][j]
            rowprod[i][v] = acc
    tab = np.zeros((256, 256), dtype=np.int64)
    for u in range(256):
        for v in range(256):
            acc = 0
            for i in range(8):
                if (u >> i) & 1:
                    acc ^= rowprod[i][v]
            tab[u, v] = acc
    return tab


def q_form(v):
    """Quadratic form q(v) = N(v~)/2 mod 2 on A, v~ the 0/1 lift in L-coords."""
    c = from_bits(v)
    N = sum(c[i] * GRAM_L[i][j] * c[j] for i in range(8) for j in range(8))
    assert N % 2 == 0
    return (N // 2) % 2


def b_form(u, v):
    """Polar form b(u,v) = q(u+v) - q(u) - q(v) = <u~,v~> mod 2."""
    cu, cv = from_bits(u), from_bits(v)
    return sum(cu[i] * GRAM_L[i][j] * cv[j] for i in range(8) for j in range(8)) % 2


# ---------------------------------------------------------------------------
# sublattices given by subspaces of A, and the Leech lattices Lambda_{V,W}
# ---------------------------------------------------------------------------
def preimage_basis(U):
    """Z-basis (HNF rows, L-coords) of the full preimage in L of a subspace
    U of A = L/2L (it contains 2L)."""
    gens = [from_bits(u) for u in basis_bits(U)] + [[2 * int(i == j) for j in range(8)] for i in range(8)]
    H = hnf(gens)
    assert len(H) == 8
    return H


def leech_basis_from_pair(V, W):
    """HNF Z-basis (24 x 24, rows in L^3-coords) of
       Lambda_{V,W} = {(x,y,z) in L^3 : x+y, x+z, y+z in V~ ; x+y+z in W~}
    (V~, W~ = preimages of V, W in L), generated by 2L^3, (a,0,a), (0,a,a)
    for a in V and (w,w,w) for w in W."""
    z8 = [0] * 8
    gens = []
    for blk in range(3):
        for k in range(8):
            v = [0] * 24
            v[8 * blk + k] = 2
            gens.append(v)
    for a in basis_bits(V):
        ca = from_bits(a)
        gens.append(ca + z8 + ca)
        gens.append(z8 + ca + ca)
    for w in basis_bits(W):
        cw = from_bits(w)
        gens.append(cw + cw + cw)
    H = hnf(gens)
    assert len(H) == 24
    return H


def in_pair_lattice_mod2(P, Vset, Wset):
    """Vectorized membership test in Lambda_{V,W} for integer rows P (n x 24):
    x+y, x+z in V (mod 2) and x+y+z in W (mod 2) [y+z then follows]."""
    P = np.asarray(P, dtype=np.int64) & 1
    w8 = (1 << np.arange(8)).astype(np.int64)
    bx = (P[:, 0:8] * w8).sum(1)
    by = (P[:, 8:16] * w8).sum(1)
    bz = (P[:, 16:24] * w8).sum(1)
    inV = np.zeros(256, dtype=bool)
    inV[list(Vset)] = True
    inW = np.zeros(256, dtype=bool)
    inW[list(Wset)] = True
    return inV[bx ^ by] & inV[bx ^ bz] & inW[bx ^ by ^ bz]


def triple_products(U, T):
    """All 576 products U[p] * U[q] (row p*24+q) of the rows of U (24 x 24,
    L^3-coords) under the Z3-routed triple product with block structure
    tensor T (8x8x8 int64).  Returns a (576, 24) int64 array."""
    U = np.asarray(U, dtype=np.int64)
    X, Y, Z = U[:, 0:8], U[:, 8:16], U[:, 16:24]

    def bp(A, B):   # bp[p,q,:] = A[p] . B[q]
        return np.einsum('pi,qj,ijk->pqk', A, B, T, optimize=True)
    P1 = bp(X, X) + bp(Z, Y) + bp(Y, Z)
    P2 = bp(Y, Y) + bp(X, Z) + bp(Z, X)
    P3 = bp(Z, Z) + bp(Y, X) + bp(X, Y)
    out = np.concatenate([P1, P2, P3], axis=2).reshape(576, 24)
    check_int64_bound(out)
    return out


def coords_in_basis(P, U):
    """Coordinates of the rows of P (n x 24 integer) in the basis U (rows),
    exact.  Returns (C, ok) with C an int64 array and ok a boolean array
    marking the rows whose coordinates are integral (row in the lattice)."""
    Ui = mat_inv([[Fr(int(t)) for t in row] for row in U])
    den = 1
    for row in Ui:
        for t in row:
            den = den * t.denominator // np.gcd(den, t.denominator)
    adj = np.array([[int(t * den) for t in row] for row in Ui], dtype=np.int64)
    P = np.asarray(P, dtype=np.int64)
    # explicit overflow bound for P @ adj
    bound = int(np.abs(P).max()) * int(np.abs(adj).max()) * P.shape[1]
    assert bound < 2 ** 62, "int64 bound exceeded in coords_in_basis"
    num = P @ adj
    ok = np.all(num % den == 0, axis=1)
    C = num // den
    return C, ok

