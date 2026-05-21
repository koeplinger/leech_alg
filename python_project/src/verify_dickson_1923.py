"""Verification of Dickson 1923, Section 19, "The Integers of Cayley's Algebra".

Source: L. E. Dickson, "A new simple theory of hypercomplex integers,"
J. Math. Pures Appl. (9) 2 (1923), 281-326; Section 19, pp. 319-324.

From Dickson's own data only -- his Cayley-Dickson product (his eq. 24) and his
order (30) with basis (29) -- this script checks:
  - the product is a composition algebra; basal units close to +/- units;
  - Dickson's order O_D (system (30)) IS closed -- a genuine maximal order;
  - there are exactly 7 maximal orders for his product (not 3);
  - Dickson's THEOREM XV ("the only maximal systems having C,U,N are (30)
    and its two (i j k)-cyclic images") undercounts: 3 stated, 7 actual;
  - the 3 of Theorem XV are exactly the maximal orders containing the Hurwitz
    quaternion p = 1/2(1+i+j+k) -- the very restriction Dickson states in his
    derivation ("We shall assume that the latter occur in our set in all
    cases," p. 321), but which holds for only 3 of the 7 maximal orders.

Dickson's axioms (his Section 6): C = closure, U = contains the unit 1,
N = norm is a rational integer.  All 7 maximal orders satisfy C, U, N.
"""
from fractions import Fraction as F
from itertools import product as iproduct, permutations
import random

# --- quaternions 1,i,j,k ---
def qmul(a, b):
    a0,a1,a2,a3 = a; b0,b1,b2,b3 = b
    return (a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
            a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0)
def qconj(a): return (a[0], -a[1], -a[2], -a[3])

# --- Dickson's octonion product, eq. (24): (q+Qe)(r+Re) = t + Te,
#     t = qr - R'Q,  T = Rq + Q r'   (prime = quaternion conjugate) ---
def omul(x, y):
    q, Q = tuple(x[:4]), tuple(x[4:])
    r, R = tuple(y[:4]), tuple(y[4:])
    t = [u-v for u, v in zip(qmul(q, r), qmul(qconj(R), Q))]
    T = [u+v for u, v in zip(qmul(R, q), qmul(Q, qconj(r)))]
    return t + T
def N(x): return sum(c*c for c in x)

units = [[F(1) if k == t else F(0) for k in range(8)] for t in range(8)]

# --- Dickson's basis (29): i,j,k,p,e,W,Z,V  in coords (1,i,j,k,e,ie,je,ke) ---
h = F(1, 2)
B = [[F(c) for c in v] for v in (
    [0,1,0,0,0,0,0,0],            # i
    [0,0,1,0,0,0,0,0],            # j
    [0,0,0,1,0,0,0,0],            # k
    [h,h,h,h,0,0,0,0],            # p = 1/2(1+i+j+k)
    [0,0,0,0,1,0,0,0],            # e
    [h,0,h,0,h,h,0,0],            # W = 1/2(1+j) + 1/2(1+i)e
    [h,h,0,0,h,0,h,0],            # Z = 1/2(1+i) + 1/2(1+j)e
    [h,0,0,h,h,0,0,h])]           # V = 1/2(1+k) + 1/2(1+k)e

def supp(v): return tuple(int(2*c) % 2 for c in v)

def code_of(halfvecs):
    gens = [supp(v) for v in halfvecs]
    C = set()
    for co in iproduct([0, 1], repeat=len(gens)):
        w = [0]*8
        for c, g in zip(co, gens):
            if c: w = [a ^ b for a, b in zip(w, g)]
        C.add(tuple(w))
    return frozenset(C)

CD = code_of([v for v in B if any(c.denominator == 2 for c in v)])

def in_lattice(v, code):
    if any((2*c).denominator != 1 for c in v): return False
    return supp(v) in code

# --- the 30 E8-type candidate codes (doubly-even self-dual [8,4]) ---
H0 = [(1,1,1,1,1,1,1,1),(0,1,0,1,0,1,0,1),(0,0,1,1,0,0,1,1),(0,0,0,0,1,1,1,1)]
def span(gs):
    S = set()
    for co in iproduct([0,1], repeat=len(gs)):
        w = [0]*8
        for c, g in zip(co, gs):
            if c: w = [a ^ b for a, b in zip(w, g)]
        S.add(tuple(w))
    return frozenset(S)
codes = sorted({span([tuple(g[p[k]] for k in range(8)) for g in H0])
                for p in permutations(range(8))}, key=lambda c: sorted(c))

def half_basis(code):
    cb = []
    for c in sorted(code):
        v = list(c); red = v[:]
        for b in cb:
            pp = next(k for k in range(8) if b[k])
            if red[pp]: red = [a ^ b for a, b in zip(red, b)]
        if any(red): cb.append(red)
    return [[F(x, 2) for x in c] for c in cb[:4]]

def closed(code):
    hs = half_basis(code)
    for gi in units + hs:
        for gj in units + hs:
            pr = omul(gi, gj)
            if not in_lattice(pr, code): return False
    return True

# (i j k) cyclic permutation of coordinates (1,i,j,k,e,ie,je,ke):
# fixes 1 and e; i->j->k->i; ie->je->ke->ie
PERM = [0, 2, 3, 1, 4, 6, 7, 5]
def cyc(code):
    out = set()
    for w in code:
        r = [0]*8
        for t in range(8): r[PERM[t]] = w[t]
        out.add(tuple(r))
    return frozenset(out)

if __name__ == "__main__":
    random.seed(0)
    comp = all(N(omul([F(random.randint(-3,3)) for _ in range(8)],
                      [F(random.randint(-3,3)) for _ in range(8)]))
               == N(x)*N(y) for x, y in
               [([F(random.randint(-3,3)) for _ in range(8)],
                 [F(random.randint(-3,3)) for _ in range(8)]) for _ in range(300)])
    # (recompute cleanly)
    comp = True
    for _ in range(300):
        x = [F(random.randint(-3,3)) for _ in range(8)]
        y = [F(random.randint(-3,3)) for _ in range(8)]
        if N(omul(x, y)) != N(x)*N(y): comp = False
    print("Dickson's product is a composition algebra:", comp)
    print("basal units close to +/- basal units:",
          all(sorted(abs(c) for c in omul(units[a], units[b])) == [0]*7+[1]
              for a in range(8) for b in range(8)))

    odclosed = all(in_lattice(omul(a, b), CD) for a in B for b in B)
    print("Dickson's order O_D (system (30)) is closed [a genuine order]:", odclosed)

    genuine = [c for c in codes if closed(c)]
    print("candidate E8-type lattices:", len(codes),
          "   maximal orders for Dickson's product:", len(genuine))
    # all 7 satisfy C, U, N:
    cun = all(closed(c) and (0,)*8 in c        # 0 in code <=> J_0 (hence 1) inside
              and all(sum(b) % 4 == 0 for b in c)   # doubly-even => integral norm
              for c in genuine)
    print("every maximal order satisfies Dickson's C, U, N:", cun)

    p_supp = (1,1,1,1,0,0,0,0)
    withp = [c for c in genuine if p_supp in c]
    orbit = {CD, cyc(CD), cyc(cyc(CD))}
    print("Dickson's O_D is among the maximal orders:", CD in genuine)
    print("maximal orders containing p=1/2(1+i+j+k):", len(withp),
          " == O_D and its two (ijk)-cyclic images:", set(withp) == orbit)
    print()
    print("THEOREM XV states 3 ;  verified count of maximal orders = ", len(genuine))
    print("Dickson's three = exactly the maximal orders containing the Hurwitz")
    print("quaternion p; the four he missed are the maximal orders without p.")
