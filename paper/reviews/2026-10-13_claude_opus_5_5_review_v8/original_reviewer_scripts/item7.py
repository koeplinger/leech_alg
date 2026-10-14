import sys; sys.path.insert(0,'.')
from vec import *
import itertools
sig = lambda v: [v[0], v[2], v[1]] + list(v[3:])
Nb = [sig(omul(l, sbar)) for l in LB]     # sigma(L sbar) basis (paper's N_k)
Mb = [sig(omul(l, s)) for l in LB]        # sigma(L s) basis (paper's M_k)
# compare with the paper's printed N_k, M_k
paperN = [(2,0,0,0,0,0,0,0),(0,0,-1,0,-1,0,-1,-1),(-1,-1,0,-1,0,-1,0,0),(0,-1,-1,0,0,-1,0,-1),(0,-1,-1,-1,0,0,-1,0),(0,-1,0,-1,-1,0,0,-1),(0,0,-1,-1,-1,-1,0,0),(0,-1,0,0,-1,-1,-1,0)]
h=Fr(1,2)
paperM = [(-3*h,-h,-h,-h,-h,-h,-h,-h),(-1,0,0,0,1,0,1,1),(0,1,1,1,0,1,0,0),(-1,0,1,0,0,1,0,1),(-1,1,1,0,0,0,1,0),(-1,1,0,1,0,0,0,1),(-1,0,1,1,1,0,0,0),(-1,1,0,0,1,1,0,0)]
print("N_k match paper:", [tuple(x) for x in Nb]==[tuple(Fr(v) for v in p) for p in paperN],
      " M_k match paper:", [tuple(x) for x in Mb]==[tuple(Fr(v) for v in p) for p in paperM])
def gram(B): return [[sum(a*b for a,b in zip(u,v)) for v in B] for u in B]
for name,B in [("sigma(L sbar)",Nb),("sigma(L s)",Mb)]:
    G=gram(B)
    print(f"{name}: det Gram = {det_fr(G)}, all norms = 0 mod 4: {all(G[i][i]%4==0 for i in range(8))}, all inner products even: {all(g%2==0 for r in G for g in r)}"
          f"  => (1/sqrt2)*lattice even unimodular rank 8 = E8, i.e. isometric to sqrt(2) L")
V = sublattice_mod2(Nb); W = sublattice_mod2(Mb)
R,_ = f2_rref(V+W,8)
print("dim V =",len(V)," dim W =",len(W)," dim(V+W) =",len(R), "-> V cap W = 0 (meet only in 2L)")
# quadratic form q on L/2L: q(a) = N(sum a_k L_k)/2 mod 2
def vec_of(a): return [sum(Fr(a[k])*LB[k][i] for k in range(8)) for i in range(8)]
allv = list(itertools.product((0,1),repeat=8))
q = {a: int(norm(vec_of(a))/2) % 2 for a in allv}
def bform(a,b):
    ab = tuple((x+y)%2 for x,y in zip(a,b)); return (q[ab]-q[a]-q[b])%2
sing = [a for a in allv if any(a) and q[a]==0]
print("nonzero singular vectors:", len(sing))
def spanf2(gs):
    S={(0,)*8}
    for g in gs: S |= {tuple((x+y)%2 for x,y in zip(s_,g)) for s_ in S}
    return frozenset(S)
Vs = spanf2([tuple(v) for v in V]); Ws = spanf2([tuple(v) for v in W])
print("q vanishes on V:", all(q[a]==0 for a in Vs), "; on W:", all(q[a]==0 for a in Ws))
# enumerate all totally singular 4-spaces
lag=set()
singset=set(sing)
def extend(S, start):
    if len(S)==16:
        lag.add(S); return
    for i in range(start,len(sing)):
        v=sing[i]
        if v in S: continue
        if all(bform(v,w)==0 for w in S):
            T_=frozenset(S | {tuple((x+y)%2 for x,y in zip(v,w)) for w in S})
            if all(q[t]==0 for t in T_):
                extend(T_, i+1)
extend(frozenset({(0,)*8}),0)
print("totally singular 4-dim subspaces (Lagrangians):", len(lag))
comp = [X for X in lag if len(X & Vs)==1]
print("Lagrangians meeting V trivially:", len(comp), "; W among them:", Ws in comp, "; V is a Lagrangian:", Vs in lag)
# ideal / subalgebra under mu-bar: integer products of bases, then mod 2
def coords2(x):
    a=L_coords(x); assert all(t.denominator==1 for t in a); return tuple(int(t)%2 for t in a)
Lg = LB
left_ideal = all(coords2(omul(l,n)) in Vs for l in Lg for n in Nb)
right_ideal = all(coords2(omul(n,l)) in Vs for l in Lg for n in Nb)
subalg_W = all(coords2(omul(m1,m2)) in Ws for m1 in Mb for m2 in Mb)
W_left = all(coords2(omul(l,m)) in Ws for l in Lg for m in Mb)
W_right = all(coords2(omul(m,l)) in Ws for l in Lg for m in Mb)
print("V left ideal:", left_ideal, " V right ideal:", right_ideal, " W subalgebra:", subalg_W, " W left ideal:", W_left, " W right ideal:", W_right)
# integer-level: sigma(L sbar) two-sided ideal of L, sigma(Ls) closed
def inlat(x, B):
    Binv_ = mat_inv_fr(B); c=[sum(Fr(x[i])*Binv_[i][j] for i in range(8)) for j in range(8)]; return all(t.denominator==1 for t in c)
print("integer level: L.N in N:", all(inlat(omul(l,n),Nb) for l in LB for n in Nb), " N.L in N:", all(inlat(omul(n,l),Nb) for l in LB for n in Nb),
      " M.M in M:", all(inlat(omul(a,b),Mb) for a in Mb for b in Mb))
# Also: Lsbar/2L and Ls/2L (untwisted) - Lagrangians too? complementary?
Vu = spanf2([tuple(v) for v in V0]); Wu = spanf2([tuple(v) for v in W0])
print("untwisted Lsbar/2L, Ls/2L Lagrangian:", Vu in lag, Wu in lag, " complementary:", len(Vu & Wu)==1)
