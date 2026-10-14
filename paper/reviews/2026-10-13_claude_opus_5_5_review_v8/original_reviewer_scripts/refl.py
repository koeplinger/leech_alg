import sys
sys.path.insert(0,'.')
from vec import *
from weyl import roots2, MU, span, is_subalg, Lsb4, Ls4
from core import lattice_index
from collections import Counter
print("max |Binv_int| =", np.abs(Binv_int).max(), " dB =", dB)
# one reflection per +- pair
seen=set(); refls=[]
for r in roots2:
    key = tuple(r) if tuple(r) > tuple(-r) else tuple(-r)
    if key in seen: continue
    seen.add(key); refls.append(np.array(key,dtype=np.int64))
print("reflections:", len(refls))
res=[]
for r in refls:
    G4 = 4*np.eye(8,dtype=np.int64) - np.outer(r,r)          # 4*s_r, r = r2/2 (norm 2)
    assert (G4 @ G4 == 16*np.eye(8,dtype=np.int64)).all()     # s_r^2 = 1
    Tg64 = np.einsum('ia,jb,kc,ijk->abc', G4, G4, G4, T)
    P = basis_products(B2, Tg64)
    assert np.abs(P).max() < 2**40
    ok, C = lam_coords_scaled(P, 256)
    closed = bool(ok.all())
    assert closed == bool(in_Lambda_scaled(P,256).all())
    idx = lattice_index(C.tolist(),24) if closed else None
    _,Av = Lcoords_scaled((G4 @ Lsb4.T).T,16); _,Aw = Lcoords_scaled((G4 @ Ls4.T).T,16)
    V = span([tuple(int(x)%2 for x in row) for row in Av.tolist()]); W = span([tuple(int(x)%2 for x in row) for row in Aw.tolist()])
    crit = is_subalg(V) and is_subalg(W)
    res.append((tuple(r), closed, idx, crit))
print("criterion==closure:", sum(c==cl for _,cl,_,c in res), "/", len(res))
cnt = Counter((cl, idx) for _,cl,idx,_ in res)
print("single reflections (closed, index): ", sorted(cnt.items(), key=str))
def rstr(r):   # root r = r2/2 in e-basis
    terms=[]
    for i,c in enumerate(r):
        if c: terms.append(("+" if c>0 else "-") + ("" if abs(c)==2 else "1/2 ") + f"e{i}")
    return " ".join(terms)
for r,cl,idx,_ in res:
    if cl and idx is not None and idx < 2**16:
        print("  witness reflection root r =", rstr(r), " index", idx)
print("sigma = reflection in e1 - e2:", [ (cl,idx) for r,cl,idx,_ in res if r==(0,2,-2,0,0,0,0,0) or r==(0,-2,2,0,0,0,0,0)])
