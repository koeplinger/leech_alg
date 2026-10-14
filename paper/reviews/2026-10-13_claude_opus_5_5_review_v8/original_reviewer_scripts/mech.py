import sys, json, itertools
sys.path.insert(0,'.')
from vec import *
from collections import Counter
res = json.load(open('census_results.json'))
LB2 = np.array([[int(2*v) for v in b] for b in LB],dtype=np.int64)
Lsb4 = np.array([[int(4*v) for v in b] for b in LsbarB],dtype=np.int64)
Ls4 = np.array([[int(4*v) for v in b] for b in LsB],dtype=np.int64)
def permute(X, p):
    Y = np.zeros_like(X)
    for a in range(8): Y[:,p[a]] = X[:,a]
    return Y
def span_set(rows):
    S={(0,)*8}
    for r in rows:
        r = tuple(int(x)%2 for x in r)
        S |= {tuple((a+b)%2 for a,b in zip(s_,r)) for s_ in S}
    return S
def inside(prod, scale, S):
    ok, A = Lcoords_scaled(prod, scale)
    return ok.all() and all(tuple(int(x)%2 for x in r) in S for r in A.tolist())
tab = Counter()
for p,closed,idx in res:
    p=tuple(p)
    Vp = permute(Lsb4,p); Wp = permute(Ls4,p)
    _,Av = Lcoords_scaled(Vp,4); _,Aw = Lcoords_scaled(Wp,4)
    Vs, Ws = span_set(Av.tolist()), span_set(Aw.tolist())
    left = inside(omul_arr(np.repeat(LB2,8,0), np.tile(Vp,(8,1))), 8, Vs)
    right = inside(omul_arr(np.tile(Vp,(8,1)), np.repeat(LB2,8,0)), 8, Vs)
    sub = inside(omul_arr(np.repeat(Wp,8,0), np.tile(Wp,(8,1))), 16, Ws)
    tab[(closed, idx, left, right, sub)] += 1
print("(closed, span index, V_pi left ideal, V_pi right ideal, W_pi subalgebra) : #perms")
for k,v in sorted(tab.items(), key=lambda kv: (not kv[0][0], str(kv[0]))):
    print("  ", k, ":", v)
