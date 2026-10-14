import sys, itertools, random, pickle
sys.path.insert(0, '.')
from core import *

# --- octonion sanity: alternative, composition, e_i^2=-1 -------------------
import numpy as np
rng = random.Random(0)
for trial in range(200):
    x = [Fr(rng.randint(-3,3),2) for _ in range(8)]
    y = [Fr(rng.randint(-3,3),2) for _ in range(8)]
    xy = omul(x,y)
    assert norm(xy) == norm(x)*norm(y)
    assert omul(omul(x,x),y) == omul(x,omul(x,y))
    assert omul(omul(y,x),x) == omul(y,omul(x,x))
    # same for the sigma-twisted product
    assert omul(omul(x,x,TS),y,TS) == omul(x,omul(x,y,TS),TS)
    assert norm(omul(x,y,TS)) == norm(x)*norm(y)
print("octonion + twisted: composition and alternativity OK on 200 random pairs")

# --- L basis --------------------------------------------------------------
assert all(in_L(b) for b in LB)
print("det LB =", det_fr(LB))
# L contains all 240 roots? enumerate roots of D8+ (odd half-int coset) and check integrality of coords
roots = []
for i,j in itertools.combinations(range(8),2):
    for si in (1,-1):
        for sj in (1,-1):
            v=[Fr(0)]*8; v[i]=Fr(si); v[j]=Fr(sj); roots.append(v)
for signs in itertools.product((1,-1),repeat=8):
    v=[Fr(sg,2) for sg in signs]
    if in_L(v): roots.append(v)
print("roots in L:", len(roots), " all norm 2:", all(norm(r)==2 for r in roots))
assert all(all(c.denominator==1 for c in L_coords(r)) for r in roots)
print("every root has integer LB-coordinates: LB spans L (det 1, so equality)")
print("s in L:", in_L(s), " sbar in L:", in_L(sbar), " e0 in L:", in_L(e(0)))

# --- Lsbar, Ls ------------------------------------------------------------
print("det(Lsbar basis) =", det_fr(LsbarB), " det(Ls basis) =", det_fr(LsB))
print("dim Lsbar/2L =", len(V0), " dim Ls/2L =", len(W0))
# all of Lsbar, Ls basis in L
assert all(in_L(b) for b in LsbarB) and all(in_L(b) for b in LsB)
# intersection mod 2 trivial and sum = L/2L
R,_ = f2_rref(V0+W0, 8)
print("dim (V0 + W0) =", len(R), " => Lsbar + Ls = L and Lsbar cap Ls = 2L (since both contain 2L, index 16 each)")
# 2L subset Lsbar: check 2*LB_k in Lsbar
assert all(in_Lsbar(scal(2,b)) and in_Ls(scal(2,b)) for b in LB)

# --- Lambda basis via F2 kernel --------------------------------------------
# variables: a_x, a_y, a_z in F2^8 (LB coords mod 2)
conds = []
def blockvec(cx, cy, cz):
    return cx + cy + cz
for h in H_V:
    conds.append(blockvec(h, h, [0]*8))   # x+y
    conds.append(blockvec(h, [0]*8, h))   # x+z
    conds.append(blockvec([0]*8, h, h))   # y+z
for h in H_W:
    conds.append(blockvec(h, h, h))
K = f2_nullspace(conds, 24)
print("dim K = Lambda/2L^3 =", len(K))
# generators of Lambda in LB^3 coordinates: K lifts + 2*unit
gens = [list(k) for k in K] + [[2 if i==j else 0 for i in range(24)] for j in range(24)]
# HNF to get a basis (exact, small): use sympy for the basis (we only need some basis; index checked separately)
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
Hm = hermite_normal_form(Matrix(gens).T)   # column-style HNF of the column lattice
Bc = [list(Hm.col(j)) for j in range(Hm.shape[1])]
assert len(Bc) == 24
# convert to R^24
def to_R24(c):
    out=[]
    for blk in range(3):
        a = c[8*blk:8*blk+8]
        out += [sum(Fr(int(a[k]))*LB[k][i] for k in range(8)) for i in range(8)]
    return out
LamB = [to_R24(c) for c in Bc]
assert all(in_Lambda(v) for v in LamB)
G = [[sum(a*b for a,b in zip(u,v)) for v in LamB] for u in LamB]
dG = det_fr(G)
print("Gram det of Lambda basis =", dG, "= 2^", dG.numerator.bit_length()-1)
print("all Gram entries even integers:", all(g.denominator==1 and g%2==0 for row in G for g in row))
print("all norms divisible by 4:", all(G[i][i] % 4 == 0 for i in range(24)))
pickle.dump(LamB, open('LamB.pkl','wb'))
