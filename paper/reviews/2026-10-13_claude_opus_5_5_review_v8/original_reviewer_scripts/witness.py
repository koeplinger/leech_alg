import sys, json, itertools, subprocess, pickle
sys.path.insert(0,'.')
from core import *     # Fraction-based: omul, in_Lambda, LB, ...
res = json.load(open('census_results.json'))
def cyc_str(p):
    seen=set(); out=[]
    for i in range(1,8):
        if i in seen or p[i]==i: continue
        c=[]; j=i
        while j not in seen:
            seen.add(j); c.append(j); j=p[j]
        out.append("("+" ".join(map(str,c))+")")
    return "".join(out) or "id"
three = [(p,i) for p,c,i in res if c and sorted(cyc_str(p).count(' ') for _ in [0])==[2] and len(cyc_str(p))==7]
print("closing 3-cycles (pi(e_a)=e_{pi(a)} convention):", [cyc_str(p) for p,i in three][:40])
p, idx = three[0]
print("witness:", cyc_str(p), "census index", idx)
# literal definition with permutation matrices and Fractions
Pm = [[Fr(1 if p[a]==b else 0) for a in range(8)] for b in range(8)]   # column a -> e_{p(a)}
def apply(M, x): return [sum(M[i][j]*x[j] for j in range(8)) for i in range(8)]
Pinv = [[Pm[j][i] for j in range(8)] for i in range(8)]                  # orthogonal
def mul_pi(x,y): return apply(Pinv, omul(apply(Pm,x), apply(Pm,y)))
def star_pi(u,v):
    x,y,z = u[0:8],u[8:16],u[16:24]; xp,yp,zp = v[0:8],v[8:16],v[16:24]
    P1 = add(mul_pi(x,xp), mul_pi(z,yp), mul_pi(y,zp))
    P2 = add(mul_pi(y,yp), mul_pi(x,zp), mul_pi(z,xp))
    P3 = add(mul_pi(z,zp), mul_pi(y,xp), mul_pi(x,yp))
    return P1+P2+P3
LamB = pickle.load(open('LamB.pkl','rb'))
Binv = mat_inv_fr(LamB)
prods = [star_pi(u,v) for u in LamB for v in LamB]
print("all 576 products in Lambda (Wilson membership, Fractions):", all(in_Lambda(w) for w in prods))
C = [[sum(w[i]*Binv[i][j] for i in range(24)) for j in range(24)] for w in prods]
print("all coordinates integral:", all(c.denominator==1 for row in C for c in row))
Ci = [[int(c) for c in row] for row in C]
json.dump({'perm':p,'C':Ci}, open('witness_C.json','w'))
M = "[" + ";".join(",".join(map(str,r)) for r in Ci) + "]"
gp = f"A=({M})~; H=mathnf(A); print(matsize(H)); print(abs(matdet(H))); print(matsnf(H));\nquit\n"
out = subprocess.run(['gp','-q','-f'], input=gp, capture_output=True, text=True).stdout
print("PARI/GP: size, index, elementary divisors:\n", out)
# also sigma for comparison
prods_s = [None]
