import sys, json, itertools
sys.path.insert(0,'.')
from vec import *
from collections import Counter
res = json.load(open('census_results.json'))
# mod-2 product table of L-bar in LB coordinates: mu(e_i, e_j) for basis classes
LB2 = np.array([[int(2*v) for v in b] for b in LB],dtype=np.int64)
I,J = np.meshgrid(range(8),range(8),indexing='ij')
prod = omul_arr(LB2[I.ravel()], LB2[J.ravel()])   # scale 4
ok, A = Lcoords_scaled(prod, 4); assert ok.all()
MU = (A % 2).reshape(8,8,8)        # MU[i,j,:] = mu(L_i, L_j) mod 2
def mu(a, b):
    a=np.array(a); b=np.array(b)
    return tuple(int(x) for x in (np.einsum('i,j,ijk->k', a, b, MU) % 2))
allv = list(itertools.product((0,1),repeat=8))
def vec_of(a): return [sum(Fr(a[k])*LB[k][i] for k in range(8)) for i in range(8)]
q = {a: int(norm(vec_of(a))/2) % 2 for a in allv}
def span(gs):
    S={(0,)*8}
    for g in gs: S |= {tuple((x+y)%2 for x,y in zip(s_,g)) for s_ in S}
    return frozenset(S)
def is_subalg(S): return all(mu(a,b) in S for a in S for b in S)
def is_lideal(S): return all(mu(a,b) in S for a in allv for b in S)
def is_rideal(S): return all(mu(b,a) in S for a in allv for b in S)
# (a) permutations
def permute(X,p):
    Y=np.zeros_like(X)
    for a in range(8): Y[:,p[a]]=X[:,a]
    return Y
Lsb4 = np.array([[int(4*v) for v in b] for b in LsbarB],dtype=np.int64)
Ls4 = np.array([[int(4*v) for v in b] for b in LsB],dtype=np.int64)
tab=Counter()
for p,closed,idx in res:
    p=tuple(p)
    _,Av = Lcoords_scaled(permute(Lsb4,p),4); _,Aw = Lcoords_scaled(permute(Ls4,p),4)
    V = span([tuple(int(x)%2 for x in r) for r in Av.tolist()]); W = span([tuple(int(x)%2 for x in r) for r in Aw.tolist()])
    tab[(closed, is_subalg(V), is_subalg(W))]+=1
print("(closed, V_pi subalgebra, W_pi subalgebra): #perms")
for k,v in sorted(tab.items(), key=str): print("  ",k,":",v)
# (b) all 270 Lagrangians
sing=[a for a in allv if any(a) and q[a]==0]
def bform(a,b):
    ab=tuple((x+y)%2 for x,y in zip(a,b)); return (q[ab]-q[a]-q[b])%2
lag=set()
def extend(S,start):
    if len(S)==16: lag.add(S); return
    for i in range(start,len(sing)):
        v=sing[i]
        if v in S: continue
        if all(bform(v,w)==0 for w in S):
            T_=frozenset(S|{tuple((x+y)%2 for x,y in zip(v,w)) for w in S})
            if all(q[t]==0 for t in T_): extend(T_,i+1)
extend(frozenset({(0,)*8}),0)
lag=list(lag)
sub=[X for X in lag if is_subalg(X)]
li=[X for X in lag if is_lideal(X)]; ri=[X for X in lag if is_rideal(X)]
print("Lagrangians:", len(lag), "; subalgebras of (Lbar, mu):", len(sub), "; left ideals:", len(li), "; right ideals:", len(ri),
      "; two-sided:", len([X for X in li if X in ri]))
pairs_sub = [(V,W) for V in sub for W in sub if len(V & W)==1]
print("ordered complementary pairs (V,W) with both subalgebras:", len(pairs_sub), " of", 270*64)
pairs_sigma_type = [(V,W) for V in li if V in ri for W in sub if len(V & W)==1]
print("ordered complementary pairs (V two-sided ideal, W subalgebra) [sigma-type]:", len(pairs_sigma_type))
