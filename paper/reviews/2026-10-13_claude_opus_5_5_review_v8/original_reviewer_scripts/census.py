import sys, time, itertools, json
sys.path.insert(0,'.')
from vec import *
from core import lattice_index
from collections import Counter, defaultdict
from multiprocessing import Pool

def cycle_type(p):   # p: tuple length 8, p[0]=0; cycle type on 1..7, non-trivial cycles
    seen=set(); ct=[]
    for i in range(1,8):
        if i in seen: continue
        L=0; j=i
        while j not in seen:
            seen.add(j); j=p[j]; L+=1
        if L>1: ct.append(L)
    return tuple(sorted(ct, reverse=True))

def work(p):
    TT = twisted_tensor(p)
    P = basis_products(B2, TT)          # scale 4
    ok, C = lam_coords_scaled(P, 4)
    closed = bool(ok.all())
    # cross-check closure by Wilson membership
    mem = bool(in_Lambda_scaled(P,4).all())
    assert closed == mem
    idx = lattice_index(C.tolist(), 24) if closed else None
    return p, closed, idx

if __name__ == '__main__':
    perms = [(0,)+q for q in itertools.permutations(range(1,8))]
    t=time.time(); r = work(SIGMA); print("sigma:", r[1], r[2], "time %.2fs" % (time.time()-t))
    if len(sys.argv) > 1 and sys.argv[1] == 'time': sys.exit()
    t=time.time()
    with Pool(15) as pool:
        res = pool.map(work, perms, chunksize=20)
    print("total time %.1fs" % (time.time()-t))
    by_ct = defaultdict(lambda: [0,0]); idx_by_ct = defaultdict(Counter); allidx = Counter()
    for p, closed, idx in res:
        ct = cycle_type(p)
        by_ct[ct][0]+=1
        if closed:
            by_ct[ct][1]+=1; idx_by_ct[ct][idx]+=1; allidx[idx]+=1
    json.dump([(list(p),c,i) for p,c,i in res], open('census_results.json','w'))
    names = {():'identity'}
    for ct in sorted(by_ct, key=lambda c:(sum(c),c)):
        print(f"{str(ct):12s} perms {by_ct[ct][0]:4d} close {by_ct[ct][1]:4d}   indices: " +
              ", ".join(f"2^{i.bit_length()-1}x{n}" if i & (i-1)==0 else f"{i}x{n}" for i,n in sorted(idx_by_ct[ct].items())))
    tot = sum(v[0] for k,v in by_ct.items() if k!=()); cl = sum(v[1] for k,v in by_ct.items() if k!=())
    print("non-identity total:", tot, "closing:", cl)
    print("overall index distribution:", sorted((("2^%d"%(i.bit_length()-1)) if i&(i-1)==0 else i, n) for i,n in allidx.items()))
    print("minimum index over closing products:", min(allidx), " (2^16 =", 2**16, ")")
