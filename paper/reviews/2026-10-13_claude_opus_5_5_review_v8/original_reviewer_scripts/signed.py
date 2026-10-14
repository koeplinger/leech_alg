import sys, json, itertools
sys.path.insert(0,'.')
from vec import *
from core import lattice_index
from multiprocessing import Pool
from collections import Counter, defaultdict
signs = np.array(list(itertools.product((1,-1),repeat=7)),dtype=np.int64)
signs = np.concatenate([np.ones((128,1),dtype=np.int64), signs],1)   # s_0 = +1
tables = {}
for q in itertools.permutations(range(1,8)):
    p = (0,)+q
    Tp = twisted_tensor(p)
    S = signs[:,:,None,None]*signs[:,None,:,None]*signs[:,None,None,:]   # (128,8,8,8)
    Tg = (S * Tp[None]).astype(np.int8)
    for k in range(128):
        key = Tg[k].tobytes()
        if key not in tables:
            tables[key] = (p, tuple(signs[k].tolist()))
print("distinct multiplication tables from signed permutations:", len(tables), flush=True)
def work(item):
    key,(p,sg) = item
    TT = np.frombuffer(key,dtype=np.int8).reshape(8,8,8).astype(np.int64)
    P = basis_products(B2, TT)
    ok, C = lam_coords_scaled(P,4)
    if not ok.all(): return p, sg, False, None
    return p, sg, True, lattice_index(C.tolist(),24)
if __name__ == '__main__':
    with Pool(15) as pool:
        out = pool.map(work, list(tables.items()), chunksize=4)
    cl = [o for o in out if o[2]]
    print("tables closing on Lambda:", len(cl), "of", len(out))
    cnt = Counter(o[3] for o in cl)
    print("span index distribution over closing tables:", sorted(("2^%d"%(k.bit_length()-1), v) for k,v in cnt.items()))
    # unsigned-permutation tables only (all signs +1)
    un = [o for o in out if all(x==1 for x in o[1])]
    print("tables from unsigned permutations:", len(un), " closing:", sum(o[2] for o in un),
          " index dist:", sorted(Counter(o[3] for o in un if o[2]).items()))
    best = min(cl, key=lambda o:o[3]); print("smallest-index closing table: perm", best[0], "signs", best[1], "index", best[3])
