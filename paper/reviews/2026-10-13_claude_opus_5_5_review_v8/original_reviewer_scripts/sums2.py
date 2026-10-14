import sys, json, itertools
sys.path.insert(0, '.')
from vec import *
from core import lattice_index
from multiprocessing import Pool
from collections import Counter

Cj = np.diag([1, -1, -1, -1, -1, -1, -1, -1]).astype(np.int64)


def phi_arr(U, V, TT=T):
    """Baez-Egan projected Jordan product as in the repo script (standard product):
    phi(u,v) = 2(ybar zpbar + ypbar zbar, zbar xpbar + zpbar xbar, xbar ypbar + xpbar ybar)."""
    x, y, z = U[:, 0:8] @ Cj, U[:, 8:16] @ Cj, U[:, 16:24] @ Cj
    xp, yp, zp = V[:, 0:8] @ Cj, V[:, 8:16] @ Cj, V[:, 16:24] @ Cj
    m = lambda a, b: np.einsum('na,nb,abc->nc', a, b, TT)
    return 2 * np.concatenate([m(y, zp) + m(yp, z), m(z, xp) + m(zp, x), m(x, yp) + m(xp, y)], 1)


n = 24
I, J = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
Ub, Vb = B2[I.ravel()], B2[J.ravel()]


def span_index_of(P, scale=4):
    ok, C = lam_coords_scaled(P, scale)
    if not ok.all():
        return None
    return lattice_index(C.tolist(), 24)


Pstar = basis_products(B2, TS)


def sum_with_star(p):
    return tuple(p), span_index_of(Pstar + basis_products(B2, twisted_tensor(p)))


def f2_rank_rows(C):
    basis = {}
    for row in C:
        v = 0
        for b in row:
            v = (v << 1) | (int(b) & 1)
        while v:
            h = v.bit_length() - 1
            if h in basis:
                v ^= basis[h]
            else:
                basis[h] = v
                break
    return len(basis)


if __name__ == '__main__':
    Pphi = phi_arr(Ub, Vb)
    print("phi closed on Lambda:", in_Lambda_scaled(Pphi, 4).all(), "; commutative:", (phi_arr(Vb, Ub) == Pphi).all(), flush=True)
    phi_idx = span_index_of(Pphi)
    print("phi span index:", phi_idx, "= 2^%d" % (phi_idx.bit_length() - 1), flush=True)
    okp, Cp = lam_coords_scaled(Pphi, 4)
    two = [[2 if i == j else 0 for i in range(24)] for j in range(24)]
    print("  index of S_phi + 2Lambda:", lattice_index(Cp.tolist() + two, 24), flush=True)
    print("span index of star + phi:", span_index_of(Pstar + Pphi), flush=True)
    res = json.load(open('census_results.json'))
    closing = [tuple(p) for p, c, i in res if c]
    rows = []
    for p in closing:
        ok, C = lam_coords_scaled(basis_products(B2, twisted_tensor(p)), 4)
        rows.append(C % 2)
    allrows = np.concatenate(rows, 0)
    r = f2_rank_rows(allrows.tolist())
    print("F2-rank of all closing products mod 2Lambda:", r,
          " => index of the sum of all S_pi (contains S_sigma, hence 2Lambda) = 2^%d" % (24 - r), flush=True)
    with Pool(15) as pool:
        out = pool.map(sum_with_star, closing, chunksize=10)
    cnt = Counter(k for _, k in out)
    print("span index of (star + star_pi) over all 1764 closing pi:", sorted(cnt.items()), flush=True)
    best = min(out, key=lambda t: t[1])
    print("smallest:", best, flush=True)
    json.dump([(list(p), k) for p, k in out], open('sum_results.json', 'w'))
