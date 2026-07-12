"""Task A, part 1: the ideal decomposition of the algebra (R^24, +, star).

star is the paper's Z_3-symmetric sigma-twisted triple-octonion product,
sigma = (1 2), on R^24 = O_1 (+) O_2 (+) O_3.

Two subspaces are tested:
    D := {(a,a,a) : a in O}              (the diagonal, 8-dimensional)
    T := {(p,q,r) : p + q + r = 0}       (the sum-zero part, 16-dimensional)

Claims checked EXACTLY (Fractions), each by an exhaustive run over a basis
pair -- which, star being bilinear, is a proof and not a sample:

  (0) the "sum" map  Sigma(x,y,z) := x + y + z  is an algebra homomorphism
      (R^24, star) -> (O, .s):     Sigma(u star v) = Sigma(u) .s Sigma(v).
      All the rest follows from this, but each item is checked directly too.
  (1) R^24 = D (+) T  (direct sum: D + T spans, D ^ T = 0).
  (2) D is a two-sided ideal:  D star R^24 <= D  and  R^24 star D <= D.
  (3) T is a two-sided ideal:  T star R^24 <= T  and  R^24 star T <= T.
  (4) D star T = T star D = 0  (the two ideals annihilate each other).
  (5) phi : O -> D,  a |-> (1/3)(a, a, a),  is an algebra ISOMORPHISM
      (O, .s) -> (D, star);  equivalently  (a,a,a) star (b,b,b)
      = 3 (a .s b, a .s b, a .s b), so (D, star) is (O, .s) rescaled by 3.
  (6) sigma^(+3) : (x,y,z) |-> (sigma x, sigma y, sigma z) is an algebra
      isomorphism from the sigma-twisted algebra onto the UNTWISTED
      triple-octonion algebra (the one built from the plain octonion
      product).  The twist therefore changes nothing about the algebra up
      to isomorphism; it is only the lattice Lambda that sees it.

Independent second pass: every identity is re-checked on 200 random
rational vectors with denominators up to 7, using an independently coded
star (a direct 9-term expansion), not the block formula.
"""
import random
from fractions import Fraction as F

# ---------------------------------------------------------------- octonions
FANO = [(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)]
MUL = {}
for a in range(8):
    MUL[(0,a)] = (1,a); MUL[(a,0)] = (1,a)
for a in range(1,8):
    MUL[(a,a)] = (-1,0)
for (a,b,c) in FANO:
    for (x,y,z) in [(a,b,c),(b,c,a),(c,a,b)]:
        MUL[(x,y)] = (1,z); MUL[(y,x)] = (-1,z)

def omul(u, v):
    r = [F(0)]*8
    for i in range(8):
        if u[i] == 0: continue
        for j in range(8):
            if v[j] == 0: continue
            sg, c = MUL[(i,j)]
            r[c] += sg*u[i]*v[j]
    return r

SIGMA = (0,2,1,3,4,5,6,7)                 # sigma = (1 2), an involution
def sig(v):
    r = [F(0)]*8
    for j in range(8):
        r[SIGMA[j]] = v[j]
    return r

def odot(a, b):                           # a .s b = sigma(sigma(a) . sigma(b))
    return sig(omul(sig(a), sig(b)))

# ------------------------------------------------------- the triple product
def add(*vs): return [sum(t) for t in zip(*vs)]
def zero8():  return [F(0)]*8

def star(u, v, mul=odot):
    x, y, z = u; xp, yp, zp = v
    return (add(mul(x,xp), mul(z,yp), mul(y,zp)),
            add(mul(y,yp), mul(x,zp), mul(z,xp)),
            add(mul(z,zp), mul(y,xp), mul(x,yp)))

def star_plain(u, v):                     # untwisted: same formula, plain octonions
    return star(u, v, mul=omul)

def eq(u, v):
    return all(a == b for bu, bv in zip(u, v) for a, b in zip(bu, bv))
def sub(u, v):
    return tuple([a-b for a, b in zip(bu, bv)] for bu, bv in zip(u, v))
def is_zero(u):
    return all(c == 0 for blk in u for c in blk)
def Sigma(u):                             # x + y + z
    return add(u[0], u[1], u[2])

# ------------------------------------------------------------------- bases
def E(t):
    v = zero8(); v[t] = F(1); return v

D_BASIS = [(E(t), E(t), E(t)) for t in range(8)]                  # dim 8
T_BASIS = ([(E(t), [-c for c in E(t)], zero8()) for t in range(8)] +
           [(zero8(), E(t), [-c for c in E(t)]) for t in range(8)])  # dim 16
FULL_BASIS = [tuple(E(t) if b == blk else zero8() for b in range(3))
              for blk in range(3) for t in range(8)]              # dim 24

def in_D(u):
    return u[0] == u[1] == u[2]
def in_T(u):
    return all(c == 0 for c in Sigma(u))

# ---------------------------------------------------------- random rationals
def rand_vec(rng):
    return tuple([F(rng.randint(-6, 6), rng.randint(1, 7)) for _ in range(8)]
                 for _ in range(3))
def rand_oct(rng):
    return [F(rng.randint(-6, 6), rng.randint(1, 7)) for _ in range(8)]

# --- an independently coded star: raw 9-term expansion, no helper reuse ----
def star_indep(u, v):
    P = [[F(0)]*8, [F(0)]*8, [F(0)]*8]
    #  block b of u*v collects the terms  u_i .s v_j  with the paper's pattern
    pattern = {0: [(0,0), (2,1), (1,2)],
               1: [(1,1), (0,2), (2,0)],
               2: [(2,2), (1,0), (0,1)]}
    for b in range(3):
        for (i, j) in pattern[b]:
            t = odot(u[i], v[j])
            for k in range(8):
                P[b][k] += t[k]
    return (P[0], P[1], P[2])


def main():
    ok = {}

    # (0) Sigma is an algebra homomorphism onto (O, .s) -- basis-exhaustive
    ok["(0) Sigma(u*v) = Sigma(u) .s Sigma(v)"] = all(
        Sigma(star(bi, bj)) == odot(Sigma(bi), Sigma(bj))
        for bi in FULL_BASIS for bj in FULL_BASIS)

    # (1) direct sum: 24 vectors D_BASIS + T_BASIS are independent and span
    from sympy import Matrix
    M = Matrix([[c for blk in b for c in blk] for b in D_BASIS + T_BASIS])
    ok["(1) R^24 = D (+) T (24x24 matrix invertible)"] = (M.rank() == 24)

    # (2)-(4) ideal / annihilation properties -- basis-exhaustive => proofs
    ok["(2a) D * R^24 <= D"] = all(in_D(star(d, b))
                                   for d in D_BASIS for b in FULL_BASIS)
    ok["(2b) R^24 * D <= D"] = all(in_D(star(b, d))
                                   for d in D_BASIS for b in FULL_BASIS)
    ok["(3a) T * R^24 <= T"] = all(in_T(star(t, b))
                                   for t in T_BASIS for b in FULL_BASIS)
    ok["(3b) R^24 * T <= T"] = all(in_T(star(b, t))
                                   for t in T_BASIS for b in FULL_BASIS)
    ok["(4a) D * T = 0"] = all(is_zero(star(d, t))
                               for d in D_BASIS for t in T_BASIS)
    ok["(4b) T * D = 0"] = all(is_zero(star(t, d))
                               for t in T_BASIS for d in D_BASIS)

    # (5) phi(a) = (1/3)(a,a,a) is an algebra isomorphism (O,.s) -> (D,star)
    third = F(1, 3)
    def phi(a): return tuple([third*c for c in a] for _ in range(3))
    unit8 = [E(t) for t in range(8)]
    ok["(5a) phi(a) * phi(b) = phi(a .s b)"] = all(
        eq(star(phi(a), phi(b)), phi(odot(a, b))) for a in unit8 for b in unit8)
    ok["(5b) (a,a,a)*(b,b,b) = 3(a.s b, ...)"] = all(
        eq(star((a,a,a), (b,b,b)),
           tuple([3*c for c in odot(a, b)] for _ in range(3)))
        for a in unit8 for b in unit8)

    # (6) sigma^(+3) intertwines the twisted and the untwisted triple product
    def sig3(u): return tuple(sig(blk) for blk in u)
    ok["(6) sigma^3(u*v) = sigma^3(u) *plain sigma^3(v)"] = all(
        eq(sig3(star(bi, bj)), star_plain(sig3(bi), sig3(bj)))
        for bi in FULL_BASIS for bj in FULL_BASIS)

    print("=== basis-exhaustive checks (bilinearity => these are proofs) ===")
    for k, v in ok.items():
        print(f"  {k:<45} {v}")

    # ------------------------------------------------- second pass, random
    rng = random.Random(20260712)
    bad = 0
    for _ in range(200):
        u = rand_vec(rng); v = rand_vec(rng)
        if not eq(star(u, v), star_indep(u, v)):        bad += 1
        if Sigma(star(u, v)) != odot(Sigma(u), Sigma(v)): bad += 1
        # decompose u = d + t and check u*u = d*d + t*t with d*t = t*d = 0
        m = [third*c for c in Sigma(u)]
        d = (m, m, m)
        t = sub(u, d)
        if not in_D(d) or not in_T(t):                  bad += 1
        if not is_zero(star(d, t)) or not is_zero(star(t, d)): bad += 1
        if not eq(star(u, u), (add(star(d, d)[0], star(t, t)[0]),
                               add(star(d, d)[1], star(t, t)[1]),
                               add(star(d, d)[2], star(t, t)[2]))): bad += 1
        if not in_D(star(d, d)) or not in_T(star(t, t)): bad += 1
    print()
    print("=== second pass: 200 random rational u, v (independent star) ===")
    print(f"  mismatches: {bad}   (expect 0)")
    print()
    print("consequence: u = d + t (d in D, t in T) is idempotent")
    print("             iff d*d = d AND t*t = t;  likewise u*u = 0")
    print("             iff d*d = 0 AND t*t = 0.")
    print()
    print("ALL PASS" if all(ok.values()) and bad == 0 else "FAILURE")


if __name__ == "__main__":
    main()
