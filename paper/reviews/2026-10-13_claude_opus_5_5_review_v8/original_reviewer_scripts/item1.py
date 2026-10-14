import sys; sys.path.insert(0,'.')
from vec import *
from core import lattice_index
print("--- Theorem sanity: Lambda * Lambda in Lambda (576 basis products) ---")
P = basis_products(B2)                       # scale 4
ok, C = lam_coords_scaled(P, 4)
print("all 576 basis products have integer Lambda-coordinates:", ok.all(), "; membership test:", in_Lambda_scaled(P,4).all())
idx = lattice_index(C.tolist(), 24)
print("span index [Lambda : S] =", idx, "= 2^%d" % (idx.bit_length()-1) if idx & (idx-1) == 0 else "")
# 2 Lambda subset S ?
Cs = C.tolist()
twoL = [[2 if i==j else 0 for i in range(24)] for j in range(24)]
print("index of S + 2Lambda:", lattice_index(Cs + twoL, 24), " (equals index of S iff 2Lambda in S)")

print("--- (1) L^3 * L^3 in L^3 ---")
# basis of L^3 (scale 2)
LB2 = np.array([[int(2*v) for v in b] for b in LB], dtype=np.int64)
L3 = np.zeros((24,24), dtype=np.int64)
for blk in range(3):
    L3[8*blk:8*blk+8, 8*blk:8*blk+8] = LB2
P3 = basis_products(L3)   # scale 4
okL = in_L_scaled(P3[:,0:8],4) & in_L_scaled(P3[:,8:16],4) & in_L_scaled(P3[:,16:24],4)
print("all 576 L^3-basis products have blocks in L:", okL.all())
# even integral: L^3 Gram
G3 = (L3 @ L3.T)
print("L^3 Gram (x4) entries divisible by 4 (integral):", (G3 % 4 == 0).all(), "; diag divisible by 8 (even):", (np.diag(G3) % 8 == 0).all())
# index [L^3 : Lambda]: Lambda-basis in L^3 coords
okc, A = Lcoords_scaled(B2.reshape(24*3,8), 2)
A = A.reshape(24,24)
print("Lambda basis in L^3: integral coords:", okc.all(), " index [L^3:Lambda] =", lattice_index(A.tolist(),24))
print("2 L^3 in Lambda:", in_Lambda_scaled(2*L3, 2).all())

print("--- (1) epsilon_1 ---")
eps1 = np.zeros(24, dtype=np.int64); eps1[0] = 2      # scale 2
rng = np.random.default_rng(5)
# random lattice vectors of Lambda: integer combos of basis
coef = rng.integers(-3,4,size=(2000,24))
Vv = coef @ B2                                         # scale 2
E = np.broadcast_to(eps1,(2000,24))
left = star_arr(E, Vv)       # scale 4
right = star_arr(Vv, E)      # scale 4
swap = np.concatenate([Vv[:,0:8], Vv[:,16:24], Vv[:,8:16]],1)*2   # scale 4
print("eps1*v == swap23(v) on 2000 random v:", (left==swap).all(), "; v*eps1 == swap23(v):", (right==swap).all())
# also on arbitrary R^24 vectors (not only lattice): random integer vectors
W = rng.integers(-5,6,size=(2000,24))
print("eps1*w == swap23(w) for 2000 random integer w in R^24:", (star_arr(E,W)==2*np.concatenate([W[:,0:8],W[:,16:24],W[:,8:16]],1)).all(),
      "; w*eps1:", (star_arr(W,E)==2*np.concatenate([W[:,0:8],W[:,16:24],W[:,8:16]],1)).all())
ee = star_arr(eps1[None,:], eps1[None,:])[0]           # scale 4
print("eps1*eps1 == eps1:", (ee == 2*eps1).all())
# swap preserves Lambda: apply to basis
sw = np.concatenate([B2[:,0:8], B2[:,16:24], B2[:,8:16]],1)
print("swap23 maps Lambda basis into Lambda:", in_Lambda_scaled(sw,2).all())
# multiples of eps1 in Lambda: t*eps1 for t in (1/2)Z, 0<t<=4 (coords of t eps1 must be in (1/2)Z)
for t2 in range(1,9):    # t = t2/2
    v = np.zeros((1,24),dtype=np.int64); v[0,0] = t2      # scale 2 => t = t2/2
    print("  t=%s: t*eps1 in Lambda? %s" % (Fr(t2,2), in_Lambda_scaled(v,2)[0]))
# index of Lambda + Z eps1 over Lambda: express in a common Z^24 frame (scale-2 coords as integer rows)
rows = B2.tolist() + [eps1.tolist()]
idx_big = lattice_index(rows, 24); idx_lam = lattice_index(B2.tolist(), 24)
print("[Z^24 : 2Lambda-frame] =", idx_lam, " [Z^24 : Lambda+Z eps1] =", idx_big, " ratio =", Fr(idx_lam, idx_big))
# closure of Lambda + Z eps1 under star: products of generators
G = np.vstack([B2, eps1[None,:]])
n = G.shape[0]
I,Jj = np.meshgrid(np.arange(n),np.arange(n),indexing='ij')
PP = star_arr(G[I.ravel()], G[Jj.ravel()])            # scale 4
# membership in Lambda + Z eps1: v - k*eps1 in Lambda for some k in {0,1,2,3}
def in_ext(Vs, scale):
    res = np.zeros(len(Vs), bool)
    for k in range(4):
        sh = Vs.copy(); sh[:,0] -= k*scale
        res |= in_Lambda_scaled(sh, scale)
    return res
print("all", len(PP), "generator products lie in Lambda + Z eps1:", in_ext(PP,4).all())
# not in L^3:
print("eps1 in L^3:", in_L_scaled(np.array([[2,0,0,0,0,0,0,0]]),2)[0], "(e0 in L)")
# inner products of eps1 with Lambda vectors: x_0 coordinate
ip = B2[:,0]   # scale 2 => <eps1, b> = B2[:,0]/2
print("inner products <eps1, basis vector> (x2):", sorted(set(ip.tolist())), "-> half-integers present:", any(v % 2 for v in ip))
