import sys, json, random
sys.path.insert(0,'.')
from vec import *
from core import lattice_index
from sums2 import f2_rank_rows, span_index_of
from multiprocessing import Pool
from collections import Counter
res = json.load(open('census_results.json'))
closing = [(tuple(p), i) for p,c,i in res if c]
def check2L(item):
    p, idx = item
    ok, C = lam_coords_scaled(basis_products(B2, twisted_tensor(p)),4)
    r = f2_rank_rows((C % 2).tolist())
    return idx == 2**(24-r)
def pair(ab):
    a,b = ab
    return ab, span_index_of(basis_products(B2,twisted_tensor(a)) + basis_products(B2,twisted_tensor(b)))
if __name__ == '__main__':
    with Pool(15) as pool:
        flags = pool.map(check2L, closing, chunksize=20)
    print("2Lambda contained in S_pi for all closing pi (index == 2^(24 - F2rank)):", all(flags), sum(flags), "/", len(flags), flush=True)
    even = [p for p,i in closing if i == 1024]
    odd = [p for p,i in closing if i == 65536]
    rng = random.Random(7)
    sample = [tuple(rng.sample(even,2)) for _ in range(300)] + [ (rng.choice(even), rng.choice(odd)) for _ in range(150)]
    with Pool(15) as pool:
        out = pool.map(pair, sample, chunksize=5)
    ce = Counter(k for (ab,k) in out[:300]); co = Counter(k for (ab,k) in out[300:])
    print("star_a + star_b, a,b even-closing (300 random pairs):", sorted(ce.items()), flush=True)
    print("star_a + star_b, a even, b odd (150 random pairs):", sorted(co.items()), flush=True)
    best = min(out, key=lambda t: t[1]); print("smallest:", best)
