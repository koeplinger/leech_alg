import sys; sys.path.insert(0,'.')
from vec import *
zero=(0,)*8
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

