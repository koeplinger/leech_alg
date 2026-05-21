"""Forensic reconstruction of Kirmse's count of EIGHT maximal orders (1924).

Uses only Kirmse's own multiplication table (1), p.64, and his order J_1, p.70.
No later author enters.

Question (prompt 103): Kirmse's paper asserts "Man findet acht solche"
(eight maximal orders) but prints no computation.  Is there a natural
necessary-but-insufficient closure test that yields exactly eight, and that
would explain both his count and his choice of the defective representative
J_1?

Result: yes.  Of the 30 E8-type candidate lattices over J_0,
  - exactly  7  are closed under multiplication           -> genuine orders;
  - exactly  8  are invariant under multiplication by the units;
  - the 8 unit-invariant ones are precisely the 7 genuine orders together
    with Kirmse's J_1;
  - J_1's failure of closure is located entirely in the products of two
    half-integer elements (10 of 16), i.e. exactly the products a
    unit-invariance test does not examine.
"""
from fractions import Fraction as F
from itertools import product as iproduct, permutations

# --- Kirmse's table (1), p.64: seven cyclic triples ---
TRIPLES = [(1,2,3),(1,5,4),(1,6,7),(2,4,7),(2,6,5),(3,4,6),(3,5,7)]
mul = {}
for a in range(8):
    mul[(0,a)] = (1,a); mul[(a,0)] = (1,a)
for a in range(1,8):
    mul[(a,a)] = (-1,0)
for (a,b,c) in TRIPLES:
    for (x,y,z) in [(a,b,c),(b,c,a),(c,a,b)]:
        mul[(x,y)] = (1,z); mul[(y,x)] = (-1,z)

def kmul(u,v):
    r = [F(0)]*8
    for i in range(8):
        if u[i]==0: continue
        for j in range(8):
            if v[j]==0: continue
            s,c = mul[(i,j)]
            r[c] += s*u[i]*v[j]
    return r

# --- the 30 E8-type candidate codes (doubly-even self-dual [8,4]) ---
H0 = [(1,1,1,1,1,1,1,1),(0,1,0,1,0,1,0,1),(0,0,1,1,0,0,1,1),(0,0,0,0,1,1,1,1)]
def span(gens):
    S=set()
    for co in iproduct([0,1],repeat=len(gens)):
        w=[0]*8
        for c,g in zip(co,gens):
            if c: w=[x^y for x,y in zip(w,g)]
        S.add(tuple(w))
    return frozenset(S)
codes = sorted({span([tuple(g[p[k]] for k in range(8)) for g in H0])
                for p in permutations(range(8))}, key=lambda c: sorted(c))

units = [[F(1) if k==j else F(0) for k in range(8)] for j in range(8)]
def half_basis(code):
    cb=[]
    for c in sorted(code):
        v=list(c); red=v[:]
        for b in cb:
            p=next(k for k in range(8) if b[k])
            if red[p]: red=[x^y for x,y in zip(red,b)]
        if any(red): cb.append(red)
    return [[F(x,2) for x in c] for c in cb[:4]]

def closed(code):
    hs = half_basis(code)
    for gi in units+hs:
        for gj in units+hs:
            p = kmul(gi,gj)
            if any((2*x).denominator!=1 for x in p): return False
            if tuple(int(2*x)%2 for x in p) not in code: return False
    return True

# unit support-permutations  (sigma: left, tau: right)
sigma = [[mul[(k,j)][1] for j in range(8)] for k in range(8)]
tau   = [[mul[(j,k)][1] for j in range(8)] for k in range(8)]
def apply_perm(perm,w):
    r=[0]*8
    for j in range(8): r[perm[j]] = w[j]
    return tuple(r)
def unit_invariant(code):
    for perm in [sigma[k] for k in range(8)]+[tau[k] for k in range(8)]:
        if frozenset(apply_perm(perm,w) for w in code) != code: return False
    return True

genuine        = [c for c in codes if closed(c)]
unitinv        = [c for c in codes if unit_invariant(c)]

# Kirmse's J_1, p.70
J1 = span([(1,1,1,1,0,0,0,0),(0,0,0,0,1,1,1,1),
           (1,1,0,0,1,1,0,0),(1,0,1,0,1,0,0,1)])

if __name__ == "__main__":
    print("E8-type candidate lattices over J_0 :", len(codes))
    print("  closed under multiplication (genuine maximal orders):", len(genuine))
    print("  invariant under multiplication by the units         :", len(unitinv))
    print()
    print("the 8 unit-invariant lattices == the 7 genuine orders + J_1 :",
          set(unitinv) == set(genuine) | {J1})
    print("J_1 is unit-invariant            :", J1 in unitinv)
    print("J_1 is a genuine order           :", J1 in genuine)
    bad = [c for c in unitinv if c not in genuine]
    print("unique unit-invariant non-order  == J_1 :", bad == [J1])
    print()
    # where do J_1's products fail?  -- using Kirmse's own generators a_1..a_4 (p.70)
    def half(*idx):
        v=[F(0)]*8
        for k in idx: v[k]=F(1,2)
        return v
    A = [half(0,1,2,3), half(4,5,6,7), half(0,1,4,5), half(0,2,4,7)]
    gens = [('e',u) for u in units] + [('h',a) for a in A]
    blk = {'ee':[0,0],'eh':[0,0],'he':[0,0],'hh':[0,0]}
    for ti,gi in gens:
        for tj,gj in gens:
            k=ti+tj; blk[k][1]+=1
            p=kmul(gi,gj)
            ok = (not any((2*x).denominator!=1 for x in p)) and \
                 tuple(int(2*x)%2 for x in p) in J1
            if not ok: blk[k][0]+=1
    print("J_1 generator-product failures by block:")
    for k in ('ee','eh','he','hh'):
        print("  %s : %2d fail / %3d" % (k, blk[k][0], blk[k][1]))
