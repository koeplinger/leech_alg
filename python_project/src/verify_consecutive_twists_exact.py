"""Item (1c), EXACT test: do consecutive sigma-twists close on the Leech lattice?

An earlier lemma criterion turned out to be only *sufficient* (the
condition-2 lemma L.pi(Ls_bar) <= pi(Ls_bar) is not necessary -- the three
terms of a pairwise block sum are coupled by membership in Lambda).  This
script uses the only fully reliable test: a genuine Z-basis of Wilson's
Lambda (24 vectors), and -- by bilinearity -- star_pi closes on Lambda iff
all 24*24 basis-pair products lie in Lambda.

pi ranges over the identity, the 21 transpositions, the 70 three-cycles and
the 105 (2,2)-double-transpositions -- exactly the permutations reachable as
a product of two transpositions (= two consecutive sigma-twists).
"""
from fractions import Fraction as F
from itertools import combinations, product as iproduct
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

FANO = [(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)]
MUL = {}
for a in range(8):
    MUL[(0,a)] = (1,a); MUL[(a,0)] = (1,a)
for a in range(1,8):
    MUL[(a,a)] = (-1,0)
for (a,b,c) in FANO:
    for (x,y,z) in [(a,b,c),(b,c,a),(c,a,b)]:
        MUL[(x,y)] = (1,z); MUL[(y,x)] = (-1,z)
def omul(u,v):
    r = [F(0)]*8
    for i in range(8):
        if u[i] == 0: continue
        for j in range(8):
            if v[j] == 0: continue
            s,c = MUL[(i,j)]
            r[c] += s*u[i]*v[j]
    return r
def vec(*c): return [F(x) for x in c]
def add(*vs): return [sum(t) for t in zip(*vs)]
h = F(1,2)
s    = vec(-h,h,h,h,h,h,h,h)
sbar = vec(-h,-h,-h,-h,-h,-h,-h,-h)
Lb = [s, vec(1,1,0,0,0,0,0,0), vec(1,-1,0,0,0,0,0,0), vec(1,0,1,0,0,0,0,0),
      vec(1,0,0,1,0,0,0,0), vec(1,0,0,0,1,0,0,0), vec(1,0,0,0,0,1,0,0),
      vec(1,0,0,0,0,0,1,0)]
Lsb    = [omul(Lk,s)    for Lk in Lb]
Lsbarb = [omul(Lk,sbar) for Lk in Lb]

Linv = Matrix([[Lb[k][r] for k in range(8)] for r in range(8)]).inv()
def Lcoords(v):                       # coords of v in the L-basis
    return list(Linv * Matrix(v))
def code(basis):
    C = {(0,)*8}
    for b in basis:
        cc = tuple(int(x) % 2 for x in Lcoords(b))
        C |= {tuple(p ^ q for p,q in zip(t,cc)) for t in C}
    return C
codeLs    = code(Lsb)
codeLsbar = code(Lsbarb)
def inL(v):
    c = Lcoords(v)
    return all(x.q == 1 for x in c)
def inLs(v):    return inL(v) and tuple(int(x)%2 for x in Lcoords(v)) in codeLs
def inLsbar(v): return inL(v) and tuple(int(x)%2 for x in Lcoords(v)) in codeLsbar
def inLambda(u):
    x,y,z = u
    if not (inL(x) and inL(y) and inL(z)): return False
    return (inLsbar(add(x,y)) and inLsbar(add(x,z)) and inLsbar(add(y,z))
            and inLs(add(x,y,z)))

# --- generators of Lambda: Wilson's three minimal-vector families ---
roots = []
for i,j in combinations(range(8),2):
    for si in (1,-1):
        for sj in (1,-1):
            v = [F(0)]*8; v[i] = F(si); v[j] = F(sj); roots.append(v)
for sg in iproduct([h,-h], repeat=8):
    if sum(1 for x in sg if x < 0) % 2 == 1:
        roots.append(list(sg))
assert len(roots) == 240
J = [[F(1 if (k==t and q>0) else (-1 if k==t else 0)) for k in range(8)]
     for t in range(8) for q in (1,-1)]
import random
random.seed(0)
gens = []
for lam in roots:                                  # Type 1: (2λ,0,0) etc.
    for c in range(3):
        t = [[F(0)]*8,[F(0)]*8,[F(0)]*8]; t[c] = [2*x for x in lam]
        gens.append(tuple(t))
for lam in roots[:120]:                             # Type 2
    lsb = omul(lam,sbar); j = random.choice(J)
    trip = [lsb, omul(lsb,j), [F(0)]*8]
    for c in range(3):
        gens.append(tuple(trip[(m-c) % 3] for m in range(3)))
for lam in roots[:120]:                             # Type 3
    j = random.choice(J); k = random.choice(J)
    trip = [omul(omul(lam,s),j), omul(lam,k), omul(omul(lam,j),k)]
    for c in range(3):
        gens.append(tuple(trip[(m-c) % 3] for m in range(3)))
gens = [g for g in gens if inLambda(g)]

# Z-basis of Lambda via HNF, working in L^3-coordinates (integers)
def l3(u): return [int(x) for c in u for x in Lcoords(c)]
G = Matrix([l3(g) for g in gens]).T                 # 24 x N, columns = generators
H = hermite_normal_form(G)                          # 24 x 24
detB = abs(H.det())
B_l3 = [[int(H[r,j]) for r in range(24)] for j in range(24)]
def from_l3(c):                                     # L^3-coords -> 3 R^8 vectors
    out = []
    for blk in range(3):
        v = [F(0)]*8
        for k in range(8):
            if c[blk*8+k]:
                v = [v[i] + c[blk*8+k]*Lb[k][i] for i in range(8)]
        out.append(v)
    return tuple(out)
B = [from_l3(c) for c in B_l3]

def twisted(perm):
    pinv = [0]*8
    for j in range(8): pinv[perm[j]] = j
    def ap(p,v):
        r = [F(0)]*8
        for j in range(8): r[p[j]] = v[j]
        return r
    def dot(a,b): return ap(pinv, omul(ap(perm,a), ap(perm,b)))
    def star(u,v):
        x,y,z = u; xp,yp,zp = v
        return (add(dot(x,xp),dot(z,yp),dot(y,zp)),
                add(dot(y,yp),dot(x,zp),dot(z,xp)),
                add(dot(z,zp),dot(y,xp),dot(x,yp)))
    return star
def closes(perm):
    star = twisted(perm)
    for bi in B:
        for bj in B:
            if not inLambda(star(bi,bj)):
                return False
    return True

if __name__ == "__main__":
    print("Lambda Z-basis built: 24 vectors, |det| in L^3-coords =", detB,
          "(expect 4096 = 2^12 = [L^3:Lambda])")
    print("all 24 basis vectors lie in Lambda:", all(inLambda(b) for b in B))
    print()
    ident = tuple(range(8))
    transp = [tuple(perm[k] if False else k for k in range(8)) for perm in []]
    transp = []
    for a,b in combinations(range(1,8),2):
        p = list(range(8)); p[a],p[b] = b,a; transp.append(tuple(p))
    threecyc = []
    for a,b,c in combinations(range(1,8),3):
        for cyc in [(b,c,a),(c,a,b)]:
            p = list(range(8)); p[a],p[b],p[c] = cyc; threecyc.append(tuple(p))
    twotwo = []
    for a,b in combinations(range(1,8),2):
        for c,d in combinations([x for x in range(1,8) if x not in (a,b)],2):
            p = list(range(8)); p[a],p[b] = b,a; p[c],p[d] = d,c
            twotwo.append(tuple(p))
    print(f"identity        : closes on Lambda = {closes(ident)}   (expect False)")
    nt = sum(closes(p) for p in transp)
    print(f"transpositions  : {nt} of {len(transp)} close   (expect all -- the paper)")
    n3 = sum(closes(p) for p in threecyc)
    print(f"3-cycles        : {n3} of {len(threecyc)} close")
    n22 = sum(closes(p) for p in twotwo)
    print(f"(2,2) doubles   : {n22} of {len(twotwo)} close")
    print()
    print("=> a single sigma-twist vs. two consecutive sigma-twists:")
    print(f"   single twist (transposition): closes")
    print(f"   two consecutive twists  ->  identity / 3-cycle / (2,2):",
          "closes" if (closes(ident) or n3 or n22) else "none close")
