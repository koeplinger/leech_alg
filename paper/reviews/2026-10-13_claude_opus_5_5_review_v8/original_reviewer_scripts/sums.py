import sys, json, itertools
sys.path.insert(0,'.')
from vec import *
from core import lattice_index
Cj = np.diag([1,-1,-1,-1,-1,-1,-1,-1]).astype(np.int64)
def phi_arr(U, V, TT=T):
    x,y,z = U[:,0:8]@Cj, U[:,8:16]@Cj, U[:,16:24]@Cj        # conjugates
    xp,yp,zp = V[:,0:8]@Cj, V[:,8:16]@Cj, V[:,16:24]@Cj
    m = lambda a,b: np.einsum('na,nb,abc->nc', a, b, TT)
    return 2*np.concatenate([m(y,zp)+m(yp,z), m(z,xp)+m(zp,x), m(x,yp)+m(xp,y)],1)
n=24
I,J = np.meshgrid(np.arange(n),np.arange(n),indexing='ij')
Ub, Vb = B2[I.ravel()], B2[J.ravel()]
def span_index_of(P, scale=4):
    ok, C = lam_coords_scaled(P, scale)
    if not ok.all(): return None
    return lattice_index(C.tolist(), 24)
Pphi = phi_arr(Ub, Vb)
print("phi closed on Lambda:", in_Lambda_scaled(Pphi,4).all(), "; commutative:", (phi_arr(Vb,Ub)==Pphi).all())
ip = phi_idx = span_index_of(Pphi)
print("phi span index:", phi_idx, "= 2^%d" % (phi_idx.bit_length()-1))
okp, Cp = lam_coords_scaled(Pphi,4)
two = [[2 if i==j else 0 for i in range(24)] for j in range(24)]
print("  index of S_phi + 2Lambda:", lattice_index(Cp.tolist()+two,24))

res = json.load(open('census_results.json'))
closing = [tuple(p) for p,c,i in res if c]
idx_of = {tuple(p):i for p,c,i in res if c}
Pstar = basis_products(B2, TS)
# 1) union of all spans
allC = []
for p in closing:
    ok, C = lam_coords_scaled(basis_products(B2, twisted_tensor(p)),4); allC += C.tolist()
print("index of sum of ALL closing spans S_pi:", lattice_index(allC, 24))
# 2) sums of two closed products: star + star_pi for all closing pi
from collections import Counter
cnt = Counter(); best=None
for p in closing:
    P = Pstar + basis_products(B2, twisted_tensor(p))
    k = span_index_of(P)
    cnt[k]+=1
    if best is None or k < best[0]: best=(k,p)
print("span index of (star + star_pi) over closing pi:", sorted(cnt.items()))
print("smallest:", best)
# 3) star + phi
print("span index of star + phi:", span_index_of(Pstar + Pphi))
# 4) pairs of 3-cycle products
three = [p for p in closing if idx_of[p]==1024][:60]
cnt2=Counter(); best2=None
for a,b in itertools.combinations(three,2):
    k = span_index_of(basis_products(B2,twisted_tensor(a)) + basis_products(B2,twisted_tensor(b)))
    cnt2[k]+=1
    if best2 is None or k<best2[0]: best2=(k,a,b)
print("span index of star_a + star_b over pairs among 60 index-2^10 perms:", sorted(cnt2.items()), " smallest:", best2)
