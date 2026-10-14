import sys, pickle, subprocess, itertools
sys.path.insert(0,'.')
from core import *          # Fraction-based only
h = Fr(1,2)
r = [h,h,h,h,h,h,h,-h]      # r = (e0+e1+...+e6-e7)/2
assert in_L(r) and norm(r) == 2
def ip(a,b): return sum(x*y for x,y in zip(a,b))
def g(x): c = ip(x,r); return [xi - c*ri for xi,ri in zip(x,r)]      # reflection s_r (involution, orthogonal)
def mul_g(x,y): return g(omul(g(x), g(y)))
print("s_r(e0) =", g(e(0)), "(identity of the twisted algebra; not e0)")
assert all(in_L(g(b)) for b in LB)   # g(L) = L
def star_g(u,v):
    x,y,z = u[0:8],u[8:16],u[16:24]; xp,yp,zp = v[0:8],v[8:16],v[16:24]
    return add(mul_g(x,xp),mul_g(z,yp),mul_g(y,zp)) + add(mul_g(y,yp),mul_g(x,zp),mul_g(z,xp)) + add(mul_g(z,zp),mul_g(y,xp),mul_g(x,yp))
LamB = pickle.load(open('LamB.pkl','rb'))
Binv = mat_inv_fr(LamB)
prods = [star_g(u,v) for u in LamB for v in LamB]
print("all 576 basis-pair products in Lambda (Wilson's conditions, Fractions):", all(in_Lambda(w) for w in prods))
C = [[sum(w[i]*Binv[i][j] for i in range(24)) for j in range(24)] for w in prods]
assert all(c.denominator==1 for row in C for c in row)
Ci = [[int(c) for c in row] for row in C]
M = "[" + ";".join(",".join(map(str,row)) for row in Ci) + "]"
out = subprocess.run(['gp','-q','-f'], input=f"A=({M})~; H=mathnf(A); print(matsize(H)); print(abs(matdet(H))); print(matsnf(H));\nquit\n", capture_output=True, text=True).stdout
print("PARI/GP (size, index, SNF):", out.split("\n")[:3])
# random-pair check beyond the basis (bilinearity sanity)
import random
rng = random.Random(3)
ok=True
for _ in range(40):
    cu=[rng.randint(-2,2) for _ in range(24)]; cv=[rng.randint(-2,2) for _ in range(24)]
    u=[sum(c*b[i] for c,b in zip(cu,LamB)) for i in range(24)]; v=[sum(c*b[i] for c,b in zip(cv,LamB)) for i in range(24)]
    ok &= in_Lambda(star_g(u,v))
print("40 random lattice pairs: products in Lambda:", ok)
# ideal / subalgebra properties of V_r = s_r(Lsbar)/2L, W_r = s_r(Ls)/2L
Vb=[g(b) for b in LsbarB]; Wb=[g(b) for b in LsB]
def inlat(x,B):
    Bi=mat_inv_fr(B); c=[sum(Fr(x[i])*Bi[i][j] for i in range(8)) for j in range(8)]; return all(t.denominator==1 for t in c)
print("L.V in V:", all(inlat(omul(l,v),Vb) for l in LB for v in Vb), " V.L in V:", all(inlat(omul(v,l),Vb) for l in LB for v in Vb),
      " V.V in V:", all(inlat(omul(a,b),Vb) for a in Vb for b in Vb), " W.W in W:", all(inlat(omul(a,b),Wb) for a in Wb for b in Wb),
      " L.W in W:", all(inlat(omul(l,w),Wb) for l in LB for w in Wb))
