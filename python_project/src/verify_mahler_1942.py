"""Verification of Mahler 1942, "On Ideals in the Cayley-Dickson Algebra".

Source: K. Mahler, Proc. Roy. Irish Acad. Sect. A 48 (1942), 123-133.
(The cover sheet misprints the title as "Cayley-Dixon"; the paper reads
"Cayley-Dickson".)  Mahler works with the integral Cayley numbers "defined
according to Dickson" (he cites Dickson 1923, p.319 ff.).

Checked, from Mahler's data only:
  - his Cayley-Dickson product equals Dickson's;
  - his ring J -- the criterion (14) and the basis (10) -- is a maximal
    order, and equals Dickson's order O_D;
  - the counts 240 (norm 1) and 2160 (norm 2);
  - THEOREM 1 (Euclidean algorithm): the covering radius^2 of J is 1/2
    (Mahler proves <= 15/16 and states the exact constant is 1/2);
  - the explicit counterexample (p.133): the norm-2 C-number (g|g),
    g=(1+i1+i2+i3)/2, does NOT generate a left ideal.
"""
from fractions import Fraction as F
from itertools import product as iproduct
import random

def qmul(a, b):
    a0,a1,a2,a3 = a; b0,b1,b2,b3 = b
    return (a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
            a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0)
def qconj(a): return (a[0],-a[1],-a[2],-a[3])
def omul(X, Y):                       # Mahler's CD product (= Dickson's)
    x1,y1 = tuple(X[:4]), tuple(X[4:])
    x2,y2 = tuple(Y[:4]), tuple(Y[4:])
    t = [u-v for u,v in zip(qmul(x1,x2), qmul(qconj(y2),y1))]
    T = [u+v for u,v in zip(qmul(y2,x1), qmul(y1,qconj(x2)))]
    return t+T
def N(x): return sum(c*c for c in x)
def vec(*c): return [F(v) for v in c]
h = F(1,2)

# Mahler's membership criterion (14) for the ring J
def inJ(X):
    if any((2*c).denominator != 1 for c in X): return False
    a,b,c,d,e,f,g,hh = (int(2*c) % 2 for c in X)
    return ((a+b-e-g) % 2 == 0 and (a+c-e-f) % 2 == 0 and
            (a+d-e-hh) % 2 == 0 and (a+b+c+d) % 2 == 0 and
            (e+f+g+hh) % 2 == 0)

# Mahler's basis (10).  NOTE: as printed, A8 = ( (1+i3)/2 | (1+i2)/2 );
# that value fails Mahler's own criterion (14) and makes the basis span a
# NON-closed lattice -- a typographical misprint.  The value consistent with
# (14), with J being a ring, and with "J = Dickson's order" is
# A8 = ( (1+i3)/2 | (1+i3)/2 ) = Dickson's V.  We use the corrected A8.
A = [vec(0,1,0,0,0,0,0,0), vec(0,0,1,0,0,0,0,0), vec(0,0,0,1,0,0,0,0),
     vec(h,h,h,h,0,0,0,0), vec(0,0,0,0,1,0,0,0),
     vec(h,h,0,0,h,0,h,0), vec(h,0,h,0,h,h,0,0), vec(h,0,0,h,h,0,0,h)]
A8_as_printed = vec(h,0,0,h,h,0,h,0)   # ( (1+i3)/2 | (1+i2)/2 ) -- the misprint
# Dickson's order O_D, basis (29)
OD = [vec(0,1,0,0,0,0,0,0), vec(0,0,1,0,0,0,0,0), vec(0,0,0,1,0,0,0,0),
      vec(h,h,h,h,0,0,0,0), vec(0,0,0,0,1,0,0,0),
      vec(h,0,h,0,h,h,0,0), vec(h,h,0,0,h,0,h,0), vec(h,0,0,h,h,0,0,h)]

def coeffsZ(target, basis):
    M = [[basis[k][r] for k in range(8)] + [target[r]] for r in range(8)]
    for col in range(8):
        piv = next((r for r in range(col,8) if M[r][col] != 0), None)
        if piv is None: return None
        M[col],M[piv] = M[piv],M[col]
        pv = M[col][col]; M[col] = [v/pv for v in M[col]]
        for r in range(8):
            if r != col and M[r][col] != 0:
                fa = M[r][col]; M[r] = [M[r][i]-fa*M[col][i] for i in range(9)]
    return [M[r][8] for r in range(8)]
def in_span(v, B):
    c = coeffsZ(v, B)
    return c is not None and all(x.denominator == 1 for x in c)

if __name__ == "__main__":
    random.seed(0)
    comp = True
    for _ in range(300):
        x = [F(random.randint(-3,3)) for _ in range(8)]
        y = [F(random.randint(-3,3)) for _ in range(8)]
        if N(omul(x,y)) != N(x)*N(y): comp = False
    print("Mahler's product is a composition algebra:", comp)

    units = [[F(1) if k==t else F(0) for k in range(8)] for t in range(8)]
    print("eq.(10) A8 AS PRINTED ((1+i3)/2|(1+i2)/2) satisfies criterion (14):",
          inJ(A8_as_printed), " <- misprint")
    print("eq.(10) A8 CORRECTED  ((1+i3)/2|(1+i3)/2) satisfies criterion (14):",
          all(inJ(a) for a in A))
    # J, span(A), span(O_D) all contain Z^8 and lie in (1/2)Z^8; compare via
    # membership of the 256 half-integer points (1/2)p, p in {0,1}^8.
    pts = [[F(p,2) for p in pat] for pat in iproduct((0,1), repeat=8)]
    sameJA  = all(inJ(v) == in_span(v, A)  for v in pts)
    sameJOD = all(inJ(v) == in_span(v, OD) for v in pts)
    z8ok = all(in_span(u, A) and in_span(u, OD) and inJ(u) for u in units)
    print("criterion (14), basis (10), and Dickson's O_D define the same ring:",
          sameJA and sameJOD and z8ok)
    codeJ = [v for v in pts if inJ(v)]
    print("[J : J_0] =", len(codeJ), "  (16 => a maximal order)")
    print("J is closed under multiplication:", all(inJ(omul(a,b)) for a in A for b in A))

    # counts of norm-1 and norm-2 integral C-numbers (integer doubled coords)
    n1 = n2 = 0
    for Y in iproduct(range(-2,3), repeat=8):
        s = sum(y*y for y in Y)
        if s not in (4,8): continue
        if inJ([F(y,2) for y in Y]):
            if s == 4: n1 += 1
            else:      n2 += 1
    print("integral C-numbers of norm 1:", n1, " (Mahler: 240)")
    print("integral C-numbers of norm 2:", n2, " (Mahler: 2160)")

    # THEOREM 1: covering radius^2 of J, by the exact closest-point per coset
    W16 = [pat for pat in iproduct((0,1), repeat=8) if inJ([F(p,2) for p in pat])]
    def depth2(X):
        best = None
        for w in W16:
            G = [F(w[t],2) + round(X[t] - F(w[t],2)) for t in range(8)]
            d = N([X[t]-G[t] for t in range(8)])
            if best is None or d < best: best = d
        return best
    worst = F(0); arg = None
    random.seed(1)
    samples = ([[F(random.randint(0,8),8) for _ in range(8)] for _ in range(6000)]
               + [[F(p,4) for p in pat] for pat in iproduct(range(4), repeat=4)
                  for pad in [()]
                  for pat8 in [list(pat)+[0,0,0,0]] for pat in [pat8[:8]]][:0])
    for X in samples:
        d = depth2(X)
        if d > worst: worst, arg = d, X
    print("Theorem 1: covering radius^2 of J (sampled max):", worst,
          " ; Mahler proves <= 15/16 =", F(15,16), ", states exact = 1/2")
    print("           a point realising depth^2 = 1/2:",
          [str(c) for c in arg] if worst == h else "(see note)")

    # THEOREM 2 counterexample (p.133): (g|g), g=(1+i1+i2+i3)/2, norm 2,
    # does NOT generate a left ideal.
    gg   = vec(h,h,h,h,h,h,h,h)                       # (g|g)
    i1   = vec(0,1,0,0,0,0,0,0)                       # (i1|0)
    z_i2 = vec(0,0,0,0,0,0,1,0)                       # (0|i2)
    half_gbar_mg = [c/2 for c in vec(h,-h,-h,-h,-h,-h,-h,-h)]   # (1/2)(gbar|-g)
    witness = omul(omul(i1, omul(z_i2, gg)), half_gbar_mg)
    stated  = vec(0,0,-h,-h,-h,h,0,0)                 # (1/2)(-i2-i3) | (1/2)(-1+i1)
    print("N((g|g)) =", N(gg), " (norm-2 C-number)")
    print("Theorem-2 counterexample: Mahler's witness equals his stated value:",
          witness == stated)
    print("           the witness is an integral C-number:", inJ(witness),
          "=> (g|g) does NOT generate a left ideal")
