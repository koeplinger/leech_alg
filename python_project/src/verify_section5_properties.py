"""Item (H): Section-5 algebraic-properties of the order (Lambda, star),
re-run on large samples, with power-associativity split into cube and quartic
and the Elduque symmetric-composition property added.

All vectors are stored in *doubled-integer* coordinates (= 2 x the actual
1/2-integer coordinate), so the whole hot loop is exact integer arithmetic.

Properties tested on minimal vectors u, v, w of Lambda (norm 8):
  commutativity        u*v = v*u
  norm multiplicativity N(u*v) = N(u)N(v) = 64
  left alternativity   (u*u)*v = u*(u*v)
  right alternativity  (u*v)*v = u*(v*v)
  flexibility          (u*v)*u = u*(v*u)
  cube power-assoc.    (u*u)*u = u*(u*u)
  quartic power-assoc. all 5 degree-4 parenthesisations of u coincide
  symmetric (Elduque)  <u*v, w> = <u, v*w>

The product star is the paper's sigma-twisted triple-octonion product,
sigma = (1 2).
"""
import time, random

# --- octonion multiplication, paper's standard Fano triples ---
FANO = [(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)]
MUL = {}
for a in range(8):
    MUL[(0,a)] = (1,a); MUL[(a,0)] = (1,a)
for a in range(1,8):
    MUL[(a,a)] = (-1,0)
for (a,b,c) in FANO:
    for (x,y,z) in [(a,b,c),(b,c,a),(c,a,b)]:
        MUL[(x,y)] = (1,z); MUL[(y,x)] = (-1,z)
# flatten to a list of (i, j, sign, c)
TERMS = [(i,j,s,c) for (i,j),(s,c) in MUL.items()]

def omul2(U,V):
    """doubled-int octonion product: returns 2*omul(U/2, V/2)."""
    r = [0]*8
    for (i,j,s,c) in TERMS:
        ui = U[i]
        if ui:
            vj = V[j]
            if vj:
                r[c] += s*ui*vj
    return [x//2 for x in r]          # exact: the sum is even for L-elements

SIG = (0,2,1,3,4,5,6,7)               # sigma = (1 2)
def perm(P,v):
    return [v[P[k]] for k in range(8)]   # (P involutive here)
def dots(a,b):                        # sigma-twisted octonion product, doubled
    return perm(SIG, omul2(perm(SIG,a), perm(SIG,b)))

def star(u,v):
    x,y,z = u; xp,yp,zp = v
    def s3(p,q,r): return [p[k]+q[k]+r[k] for k in range(8)]
    P1 = s3(dots(x,xp), dots(z,yp), dots(y,zp))
    P2 = s3(dots(y,yp), dots(x,zp), dots(z,xp))
    P3 = s3(dots(z,zp), dots(y,xp), dots(x,yp))
    return (P1,P2,P3)
def eq(u,v): return u[0]==v[0] and u[1]==v[1] and u[2]==v[2]
def Nsq(u):  return sum(c*c for blk in u for c in blk)        # = 4*N(actual)
def ip(u,v): return sum(a*b for bu,bv in zip(u,v) for a,b in zip(bu,bv))  # 4*<.,.>

# --- octonion units, roots, s, sbar in doubled-int coordinates ---
def e(t):  v=[0]*8; v[t]=2; return v          # 2 * e_t
ds    = [-1,1,1,1,1,1,1,1]                    # 2 * s
dsbar = [-1,-1,-1,-1,-1,-1,-1,-1]             # 2 * sbar
droots = []
from itertools import combinations, product as iproduct
for i,j in combinations(range(8),2):
    for si in (2,-2):
        for sj in (2,-2):
            v=[0]*8; v[i]=si; v[j]=sj; droots.append(v)      # 2*lambda, integer root
for sg in iproduct((1,-1),repeat=8):
    if sum(1 for x in sg if x<0) % 2 == 1:
        droots.append(list(sg))               # 2*lambda, half-integer root
assert len(droots)==240
dunits = [e(t) for t in range(8)] + [[-c for c in e(t)] for t in range(8)]  # 2*(+-e_t)

Z = [0]*8
def cyc(trip, sh): return tuple(trip[(k-sh)%3] for k in range(3))

def rand_min():
    """a uniformly-ish random minimal vector of Lambda, doubled-int coords."""
    r = random.random()
    lam = random.choice(droots)
    sh  = random.randint(0,2)
    if r < 720/196560:                                   # Type 1: (2 lambda,0,0)
        blk = [2*c for c in lam]
        return cyc([blk, Z, Z], sh)
    elif r < 12240/196560:                               # Type 2
        j = random.choice(dunits)
        lsb = omul2(lam, dsbar)
        return cyc([lsb, omul2(lsb,j), Z], sh)
    else:                                                # Type 3
        j = random.choice(dunits); k = random.choice(dunits)
        return cyc([omul2(omul2(lam,ds),j), omul2(lam,k),
                    omul2(omul2(lam,j),k)], sh)

if __name__ == "__main__":
    random.seed(20260522)
    # sanity: sampled vectors have norm 8  (Nsq = 4*8 = 32)
    assert all(Nsq(rand_min())==32 for _ in range(500)), "minimal-vector norm check"

    budget = 540.0
    t0 = time.time()
    keys = ["commutativity","norm multiplicativity","left alternativity",
            "right alternativity","flexibility","cube power-assoc",
            "quartic power-assoc","symmetric (Elduque)"]
    cnt = {k:0 for k in keys}
    N = 0
    while time.time()-t0 < budget:
        for _ in range(2000):                  # check the clock every 2000 iters
            u = rand_min(); v = rand_min(); w = rand_min()
            uv = star(u,v); vu = star(v,u)
            uu = star(u,u); vv = star(v,v)
            if eq(uv,vu):                       cnt["commutativity"] += 1
            if Nsq(uv)==256:                    cnt["norm multiplicativity"] += 1
            if eq(star(uu,v), star(u,uv)):      cnt["left alternativity"] += 1
            if eq(star(uv,v), star(u,vv)):      cnt["right alternativity"] += 1
            if eq(star(uv,u), star(u,vu)):      cnt["flexibility"] += 1
            A = star(uu,u); B = star(u,uu)      # the two cubes
            if eq(A,B):                         cnt["cube power-assoc"] += 1
            p1=star(A,u); p2=star(B,u); p3=star(uu,uu)
            p4=star(u,A); p5=star(u,B)
            if eq(p1,p2) and eq(p1,p3) and eq(p1,p4) and eq(p1,p5):
                cnt["quartic power-assoc"] += 1
            if ip(uv,w)==ip(u,star(v,w)):       cnt["symmetric (Elduque)"] += 1
            N += 1
    el = time.time()-t0
    print(f"samples per property: {N}     runtime: {el:.0f} s")
    print(f"{'property':<26}{'holds':>12}{'of':>10}{'percent':>11}")
    for k in keys:
        print(f"{k:<26}{cnt[k]:>12}{N:>10}{100.0*cnt[k]/N:>10.3f}%")
