"""Random Weyl-group twists g in W(E8) = O(L): x ._g y := g^{-1}(g(x) . g(y)), blockwise, Z3 routing.
Tests (i) criterion: closure on Lambda <=> g(Lsbar)/2L and g(Ls)/2L both subalgebras of L/2L,
(ii) span-index distribution over closing g.  Exact integer arithmetic on scaled coordinates."""
import sys, itertools, random
sys.path.insert(0, '.')
from vec import *
from core import lattice_index
from multiprocessing import Pool
from collections import Counter

# roots of L at scale 2
roots2 = []
for i, j in itertools.combinations(range(8), 2):
    for si in (2, -2):
        for sj in (2, -2):
            v = [0] * 8; v[i] = si; v[j] = sj; roots2.append(v)
for sg in itertools.product((1, -1), repeat=8):
    if (sum(sg) // 2) % 2 == 1 or (sum(sg) // 2) % 2 == -1:
        roots2.append(list(sg))
roots2 = np.array(roots2, dtype=np.int64)
assert len(roots2) == 240 and in_L_scaled(roots2, 2).all()
pos = roots2[:120] if False else roots2  # use all (reflection s_r = s_{-r})


def refl4(r2):
    """4 * (I - r r^T) for root r = r2/2 (norm 2): 4I - r2 r2^T."""
    return 4 * np.eye(8, dtype=np.int64) - np.outer(r2, r2)


def random_g(rng, k=40):
    """returns G (exact Fractions not needed): integer matrix G4 = 4g, with g a product of k reflections."""
    G = np.eye(8, dtype=object)
    for _ in range(k):
        r = roots2[rng.randrange(240)]
        R = refl4(r).astype(object)
        G = (R.dot(G))
        G = np.vectorize(lambda t: Fr(t))(G) if False else G
        # renormalize: R is 4*s_r, so divide by 4 exactly
        G = np.vectorize(lambda t: Fr(t) / 4)(G)
    return G   # Fractions, entries in (1/4)Z


MU = None


def setup_mu():
    global MU
    LB2 = np.array([[int(2 * v) for v in b] for b in LB], dtype=np.int64)
    I, J = np.meshgrid(range(8), range(8), indexing='ij')
    prod = omul_arr(LB2[I.ravel()], LB2[J.ravel()])
    ok, A = Lcoords_scaled(prod, 4)
    MU = (A % 2).reshape(8, 8, 8)


setup_mu()


def span(gs):
    S = {(0,) * 8}
    for g in gs:
        S |= {tuple((x + y) % 2 for x, y in zip(s_, g)) for s_ in S}
    return frozenset(S)


def is_subalg(S):
    arr = np.array(sorted(S))
    P = np.einsum('ni,mj,ijk->nmk', arr, arr, MU) % 2
    return all(tuple(int(x) for x in v) in S for v in P.reshape(-1, 8))


Lsb4 = np.array([[int(4 * v) for v in b] for b in LsbarB], dtype=np.int64)
Ls4 = np.array([[int(4 * v) for v in b] for b in LsB], dtype=np.int64)


def work(seed):
    rng = random.Random(seed)
    G = random_g(rng)
    G4 = np.array([[int(4 * t) for t in row] for row in G], dtype=np.int64)
    # g(L) = L check (rows of LB): g acts on column vectors; images of basis vectors
    LB2 = np.array([[int(2 * v) for v in b] for b in LB], dtype=np.int64)
    imgs = (G4 @ LB2.T).T          # scale 8
    assert in_L_scaled(imgs, 8).all()
    # pair (V, W) = (g(Lsbar)/2L, g(Ls)/2L)
    _, Av = Lcoords_scaled((G4 @ Lsb4.T).T, 16)
    _, Aw = Lcoords_scaled((G4 @ Ls4.T).T, 16)
    V = span([tuple(int(x) % 2 for x in r) for r in Av.tolist()])
    W = span([tuple(int(x) % 2 for x in r) for r in Aw.tolist()])
    crit = is_subalg(V) and is_subalg(W)
    # twisted tensor: T_g[a,b,c] = sum g[i,a] g[j,b] g[k,c] T[i,j,k]  (g orthogonal); scaled by 64
    Tg64 = np.einsum('ia,jb,kc,ijk->abc', G4, G4, G4, T)
    P = basis_products(B2, Tg64)       # scale 4*64 = 256
    ok, C = lam_coords_scaled(P, 256)
    closed = bool(ok.all())
    idx = lattice_index(C.tolist(), 24) if closed else None
    return seed, crit, closed, idx, (V, W)


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    with Pool(15) as pool:
        out = pool.map(work, range(N), chunksize=4)
    agree = sum(1 for o in out if o[1] == o[2])
    print(f"random g in W(E8): {N}; criterion==closure: {agree}/{N}; closing: {sum(o[2] for o in out)}", flush=True)
    cnt = Counter(o[3] for o in out if o[2])
    print("span index distribution over closing g:", sorted((("2^%d" % (k.bit_length() - 1)) if k & (k - 1) == 0 else str(k), v) for k, v in cnt.items()))
    print("distinct (V,W) pairs hit:", len(set(o[4] for o in out)), "; closing pairs:", len(set(o[4] for o in out if o[2])))
    best = min((o for o in out if o[2]), key=lambda o: o[3])
    print("smallest index:", best[3], "at seed", best[0])
