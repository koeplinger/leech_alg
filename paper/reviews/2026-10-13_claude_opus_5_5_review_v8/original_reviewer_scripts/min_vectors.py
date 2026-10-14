import sys; sys.path.insert(0,'.')
from vec import *
import itertools
# roots of L, scaled by 2
roots2 = []
for i,j in itertools.combinations(range(8),2):
    for si in (2,-2):
        for sj in (2,-2):
            v=[0]*8; v[i]=si; v[j]=sj; roots2.append(v)
for signs in itertools.product((1,-1),repeat=8):
    if sum(signs) % 4 == 2 or sum(signs) % 4 == -2:  # half-int sum odd  <=> sum(signs)/2 odd
        roots2.append(list(signs))
roots2 = np.array(roots2, dtype=np.int64)
assert len(roots2) == 240 and in_L_scaled(roots2, 2).all()
s2 = np.array([-1]+[1]*7, dtype=np.int64); sb2 = np.array([-1]*8, dtype=np.int64)
J = []
for t in range(8):
    for sg in (1,-1):
        v = [0]*8; v[t] = sg; J.append(v)
J = np.array(J, dtype=np.int64)
def cyc(x,y,z):
    return [np.concatenate([x,y,z],1), np.concatenate([z,x,y],1), np.concatenate([y,z,x],1)]
vecs = []
Z0 = np.zeros((240,8),dtype=np.int64)
# Type 1: (2 lam,0,0), scale 4: 2*lam = roots2 (scale 2) -> *2 for scale 4

vecs = cyc(4*roots2, Z0, Z0)   # 2*lam at scale 4 = 8 lam = 4*roots2
# Type 2: (lam sbar, (lam sbar) j, 0), scale 4
lsb = omul_arr(roots2, np.broadcast_to(sb2,(240,8)))      # scale 4
A = np.repeat(lsb, 16, axis=0); Jr = np.tile(J, (240,1))
B = omul_arr(A, Jr)                                        # scale 4
vecs += cyc(A, B, np.zeros_like(A))
# Type 3: ((lam s) j, lam k, (lam j) k), scale 4
ls = omul_arr(roots2, np.broadcast_to(s2,(240,8)))         # scale 4
idx = np.array(list(itertools.product(range(240), range(16), range(16))))
L_ = ls[idx[:,0]]; jj = J[idx[:,1]]; kk = J[idx[:,2]]; lam = roots2[idx[:,0]]
X = omul_arr(L_, jj)                       # scale 4
Y = 2*omul_arr(lam, kk)                    # lam k at scale 2 -> scale 4
Zt = 2*omul_arr(omul_arr(lam, jj), kk)     # scale 2 -> 4
vecs += cyc(X, Y, Zt)
V = np.concatenate(vecs, 0)
print("generated", V.shape[0], "vectors")
norms = (V*V).sum(1)
print("norms (x16):", np.unique(norms), "-> real norm", np.unique(norms)/16)
U = np.unique(V, axis=0)
print("distinct:", U.shape[0])
print("all in Lambda:", in_Lambda_scaled(V, 4).all())
ok, C = lam_coords_scaled(V, 4)
print("all have integer Lambda-basis coordinates:", ok.all())
np.save('minvecs4.npy', U)
