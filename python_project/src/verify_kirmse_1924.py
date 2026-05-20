"""Verify Kirmse 1924 from scratch: his multiplication table (1) and his order J_1.
Sources, both from Kirmse's paper (pp.63-82, scans seq_325-344):
  - table (1), p.64 (seq_326): 7 cyclic Fano triples below.
  - J_1, p.70 (seq_332): J_1 = [J_0, a1,a2,a3,a4].
"""
from fractions import Fraction as F
from itertools import product as iproduct

# --- Kirmse's octonion multiplication, table (1), p.64 ---
# cyclic triples (a,b,c): i_a i_b = +i_c (and cyclic); anti-cyclic = -.
TRIPLES=[(1,2,3),(1,5,4),(1,6,7),(2,4,7),(2,6,5),(3,4,6),(3,5,7)]
mul={}
for a in range(8):
    mul[(0,a)]=(1,a); mul[(a,0)]=(1,a)
for a in range(1,8):
    mul[(a,a)]=(-1,0)
for (a,b,c) in TRIPLES:
    for (x,y,z) in [(a,b,c),(b,c,a),(c,a,b)]:
        mul[(x,y)]=(1,z); mul[(y,x)]=(-1,z)
assert len(mul)==64, len(mul)

def kmul(u,v):
    r=[F(0)]*8
    for i in range(8):
        if u[i]==0: continue
        for j in range(8):
            if v[j]==0: continue
            s,c=mul[(i,j)]
            r[c]+=s*u[i]*v[j]
    return r

def N(u): return sum(x*x for x in u)
# sanity: composition law N(uv)=N(u)N(v)
import random; random.seed(0)
ok=all(N(kmul([F(random.randint(-3,3)) for _ in range(8)],
               [F(random.randint(-3,3)) for _ in range(8)]))
       == N(u)*N(v) for u,v in [([F(random.randint(-2,2)) for _ in range(8)],
                                 [F(random.randint(-2,2)) for _ in range(8)]) for _ in range(1)])
# proper composition check
good=True
for _ in range(300):
    u=[F(random.randint(-3,3)) for _ in range(8)]
    v=[F(random.randint(-3,3)) for _ in range(8)]
    if N(kmul(u,v))!=N(u)*N(v): good=False
print("Kirmse table (1) is a composition algebra (N(uv)=N(u)N(v)):", good)

def half(*idx):
    v=[F(0)]*8
    for k in idx: v[k]=F(1,2)
    return v
# --- Kirmse's J_1, p.70:  J_1 = [J_0, a1,a2,a3,a4] ---
a1=half(0,1,2,3)            # 1/2(1+i1+i2+i3)
a2=half(4,5,6,7)            # 1/2(i4+i5+i6+i7)
a3=half(0,1,4,5)            # 1/2(1+i1+i4+i5)
a4=half(0,2,4,7)            # 1/2(1+i2+i4+i7)
units=[[F(1) if k==j else F(0) for k in range(8)] for j in range(8)]
gens=units+[a1,a2,a3,a4]    # 12 Z-module generators of J_1

# code C = F2-span of the half-integer support patterns
def supp(v): return tuple(int(2*x)%2 for x in v)
basisC=[supp(a1),supp(a2),supp(a3),supp(a4)]
C=set()
for coeffs in iproduct([0,1],repeat=4):
    w=[0]*8
    for c,b in zip(coeffs,basisC):
        if c:
            w=[(x^y) for x,y in zip(w,b)]
    C.add(tuple(w))
print("dim of code C =", len(C).bit_length()-1, " |C| =",len(C),"  => [J_1:J_0] =",len(C))

def inJ1(v):
    if any((2*x).denominator!=1 for x in v): return False   # must be in (1/2)Z^8
    return supp(v) in C

# --- closure test: all 144 generator products ---
fails=[]
for gi in gens:
    for gj in gens:
        p=kmul(gi,gj)
        if not inJ1(p): fails.append((gi,gj,p))
print()
print("J_1 . J_1  subseteq J_1  ?  ", "YES (closed)" if not fails else "NO -- NOT closed")
print("number of generator-pair products that leave J_1:", len(fails))
if fails:
    gi,gj,p=fails[0]
    def s(v): return "("+",".join(str(x) for x in v)+")"
    print("  example:  x =",s(fails[0][0]))
    print("            y =",s(fails[0][1]))
    print("          x*y =",s(p),"  -- leaves J_1")

print()
print("="*64)
print("Enumerating maximal orders for Kirmse's table (1)")
print("="*64)
from itertools import permutations

# E8-type candidates: doubly-even self-dual [8,4] codes = extended Hamming
# variants.  Start from RM(1,3) and apply all coordinate permutations.
H0=[(1,1,1,1,1,1,1,1),(0,1,0,1,0,1,0,1),(0,0,1,1,0,0,1,1),(0,0,0,0,1,1,1,1)]
def span(gens):
    S=set()
    for co in iproduct([0,1],repeat=len(gens)):
        w=[0]*8
        for c,g in zip(co,gens):
            if c: w=[x^y for x,y in zip(w,g)]
        S.add(tuple(w))
    return frozenset(S)
codes=set()
for p in permutations(range(8)):
    codes.add(span([tuple(g[p[k]] for k in range(8)) for g in H0]))
print("distinct doubly-even self-dual [8,4] codes (E8 candidates):",len(codes))

def Lclosed(code):
    """is L_code = J_0 + {1/2 c : c in code} closed under Kirmse mult?"""
    # basis of the code
    cb=[]
    for c in code:
        if c==tuple([0]*8): continue
        # gaussian over F2
        v=list(c)
        red=v[:]
        for b in cb:
            p=next(k for k in range(8) if b[k])
            if red[p]: red=[x^y for x,y in zip(red,b)]
        if any(red): cb.append(red)
    cb=cb[:4]
    hs=[[F(x,2) for x in c] for c in cb]              # 1/2-basis vectors
    gens=units+hs
    for gi in gens:
        for gj in gens:
            p=kmul(gi,gj)
            if any((2*x).denominator!=1 for x in p): return False
            if tuple(int(2*x)%2 for x in p) not in code: return False
    return True

closed=[c for c in codes if Lclosed(c)]
print("of these, closed under Kirmse's table (1)  =>  maximal orders:", len(closed))

# where does Kirmse's named J_1 sit?
J1code=frozenset(C)
print()
print("Kirmse's J_1 code is among the 30 E8 candidates:", J1code in codes)
print("Kirmse's J_1 is closed (a genuine maximal order):", J1code in closed)
print()
print("VERDICT:  Kirmse claimed 8 maximal orders; the true number is",
      len(closed),".")
print("          His named, used order J_1 is one of the non-closed candidates.")
