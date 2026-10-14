import sys; sys.path.insert(0,'.')
from vec import *
import itertools
from collections import Counter
# all doubly-even self-dual [8,4,4] codes (extended Hamming codes) as frozensets of 8-bit tuples
w4 = [tuple(1 if i in S else 0 for i in range(8)) for S in itertools.combinations(range(8),4)]
one = (1,)*8; zero=(0,)*8
def span(gs):
    S={zero}
    for g in gs:
        S |= {tuple((a+b)%2 for a,b in zip(x,g)) for x in S}
    return frozenset(S)
codes=set()
for a,b,c in itertools.combinations(w4,3):
    C = span([one,a,b,c])
    if len(C)==16 and all(sum(x) in (0,4,8) for x in C):
        codes.add(C)
codes=list(codes)
print("number of [8,4,4] doubly-even codes:", len(codes))

def lattice_closed(C, TT):
    """O_C = Z^8 + (1/2)C (unit scale). Generators: e_i and c/2 (c in C). Closure: all generator products in O_C.
    Work at scale 2: generator g2 integer; product at scale 4; membership: p4/4 in (1/2)Z^8 and (p4/2 mod 2) in C."""
    gens2 = [tuple(2 if i==j else 0 for i in range(8)) for j in range(8)] + [c for c in C if c!=zero]
    G = np.array(gens2,dtype=np.int64)
    n=len(G)
    I,J=np.meshgrid(range(n),range(n),indexing='ij')
    P = np.einsum('na,nb,abc->nc', G[I.ravel()], G[J.ravel()], TT)   # scale 4
    if (P % 2 != 0).any(): return False
    Q = (P//2) % 2
    return all(tuple(q) in C for q in Q.tolist())

def left_unit_stable(C, TT):
    gens2 = [c for c in C if c!=zero]
    for t in range(8):
        for v in gens2:
            et = np.zeros(8,dtype=np.int64); et[t]=1
            p = np.einsum('a,b,abc->c', et, np.array(v), TT)   # scale 2
            if tuple(int(x)%2 for x in p) not in C: return False
    return True
def right_unit_stable(C, TT):
    gens2 = [c for c in C if c!=zero]
    for t in range(8):
        for v in gens2:
            et = np.zeros(8,dtype=np.int64); et[t]=1
            p = np.einsum('a,b,abc->c', np.array(v), et, TT)
            if tuple(int(x)%2 for x in p) not in C: return False
    return True

Top = np.transpose(T,(1,0,2)).copy()   # opposite multiplication x.op y = y.x
for name,TT in [("paper's Fano triples", T), ("opposite multiplication", Top), ("sigma-twisted (1 2)", TS)]:
    closed = [C for C in codes if lattice_closed(C,TT)]
    lstab = [C for C in codes if left_unit_stable(C,TT)]
    rstab = [C for C in codes if right_unit_stable(C,TT)]
    print(f"{name:26s}: closed E8 overlattices of Z^8: {len(closed)};  left-unit-stable: {len(lstab)}; right-unit-stable: {len(rstab)};"
          f" closed subset of left-stable: {all(C in lstab for C in closed)}")
cl_std = {C for C in codes if lattice_closed(C,T)}
cl_op = {C for C in codes if lattice_closed(C,Top)}
cl_tw = {C for C in codes if lattice_closed(C,TS)}
print("closed for both orientations:", len(cl_std & cl_op), "; std vs sigma-twisted overlap:", len(cl_std & cl_tw))
# the left-stable-but-not-closed one(s) (Kirmse's J_1 analogue)
lst = [C for C in codes if left_unit_stable(C,T)]
bad = [C for C in lst if C not in cl_std]
print("left-unit-stable but not closed (paper's orientation):", len(bad))
for C in bad:
    print("   weight-4 supports:", sorted(tuple(i for i in range(8) if x[i]) for x in C if sum(x)==4))
# each closed lattice: is it an E8 (unimodular for 2<x,y>, min norm 1, 240 units)?
for C in list(cl_std)[:1]:
    units = [v for v in itertools.product(*([range(-1,2)]*8)) if sum(x*x for x in v)==1]
    halfs = [v for v in itertools.product((-1,1),repeat=8)]
    n_units = 16 + sum(1 for c in C if sum(c)==4)*16
    print("units in a closed one:", n_units)
# Do the seven closed ones contain the Gravesian integers? (by construction yes) -- and are they the only closed orders strictly
# between Z^8 and (1/2)Z^8 that are maximal?  Check: every closed lattice Z^8 + (1/2)S with S a linear code of dim<=4,
# S subset of some closed C, is contained in a closed E8; list closed codes of all dimensions that are 'even' (all weights = 0 mod 4)
