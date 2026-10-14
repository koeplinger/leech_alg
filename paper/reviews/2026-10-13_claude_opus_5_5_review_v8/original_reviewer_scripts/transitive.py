import sys, itertools
sys.path.insert(0,'.')
from vec import *
from weyl import roots2
# work in LB coordinates mod 2, vectors as 8-bit ints
def vec_of(a): return [sum(Fr(a[k])*LB[k][i] for k in range(8)) for i in range(8)]
allv = list(itertools.product((0,1),repeat=8))
toint = lambda a: sum(b<<i for i,b in enumerate(a))
Q = {}; 
for a in allv: Q[toint(a)] = int(norm(vec_of(a))/2) % 2
def B(x,y): return (Q[x^y]^Q[x]^Q[y])
# root classes mod 2
_, Ar = Lcoords_scaled(roots2, 2)
rcls = sorted(set(toint(tuple(int(t)%2 for t in row)) for row in Ar.tolist()))
print("root classes mod 2L:", len(rcls))
# reflection s_r mod 2: v -> v + <v,r> r ; <v,r> mod 2 = polar form B(v, r)
maps = [[(v ^ r) if B(v,r) else v for v in range(256)] for r in rcls]
def sub_of(rows):
    S={0}
    for r in rows: S |= {s^r for s in S}
    return tuple(sorted(S))
_,Av = Lcoords_scaled(np.array([[int(4*v) for v in b] for b in LsbarB]),4)
_,Aw = Lcoords_scaled(np.array([[int(4*v) for v in b] for b in LsB]),4)
V0 = sub_of([toint(tuple(int(t)%2 for t in row)) for row in Av.tolist()])
W0 = sub_of([toint(tuple(int(t)%2 for t in row)) for row in Aw.tolist()])
start=(V0,W0); seen={start}; frontier=[start]
while frontier:
    nxt=[]
    for (V,W) in frontier:
        for m in maps:
            st=(tuple(sorted(m[v] for v in V)), tuple(sorted(m[w] for w in W)))
            if st not in seen: seen.add(st); nxt.append(st)
    frontier=nxt
print("ordered pairs reachable from (Lsbar/2L, Ls/2L) under W(E8):", len(seen), " (270*64 =", 270*64, ")")
