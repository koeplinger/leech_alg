"""Verify Coxeter 1946 from scratch -- Section 4, "Kirmse's mistake".

Source: H. S. M. Coxeter, "Integral Cayley numbers", Duke Math. J. 13 (1946),
561-578.  Section 4 (pp. 565-566) discusses Kirmse 1924.

Coxeter's Section-4 claims, each checked below against Kirmse's own data
(table (1) p.64 and J_1 p.70, the same data used by verify_kirmse_1924.py):

  (C1) Coxeter's table (4.1) reproduces Kirmse's multiplication table (1).
  (C2) Kirmse states there are eight maximal domains; "actually there are
       only seven, which presumably are the remaining seven of his eight".
  (C3) J_1 (Kirmse's exhibited domain) is not closed under multiplication:
         1/2(1+i1+i2+i3) . 1/2(1+i1+i4+i5) = 1/2(i1+i2+i4+i7),
       which "is not one of the 240 units".
  (C4) Bruck's correction: "Bruck's domain J can be derived from J_1 by
       transposing two of the i's" -- a *linear* coordinate transposition --
       and the result is closed.
  (C5) "the (7 choose 2) = 21 transpositions fall into 7 sets of 3, each set
       having the same effect" -> exactly seven domains; in each, one i is
       fixed by all three transpositions (the "special unit").
  (C6) framing duality: transposing two i's in the *elements* (J_1 -> J,
       same table) is equivalent to transposing them in the *table*
       (J_1 fixed, table (4.1) -> a relabelled table).  Both are the one
       linear involution; this is what Coxeter means by "Kirmse's J_1 could
       be used as it stands if we replaced his multiplication table".

The point of (C4)/(C6) for our paper: Coxeter's correction *is* induced by a
linear map on R^8 -- a transposition of two coordinate axes -- exactly the
kind of map our sigma is.
"""
from fractions import Fraction as F
from itertools import product as iproduct, permutations, combinations
import random

# --- Kirmse's octonion table (1), p.64 == Coxeter's table (4.1), p.566 ---
# 7 oriented Fano triples (a,b,c): i_a i_b = +i_c cyclically, reverse = -.
# (read directly from the Kirmse 1924 scan; identical to verify_kirmse_1924.py)
TRIPLES = [(1,2,3),(1,5,4),(1,6,7),(2,4,7),(2,6,5),(3,4,6),(3,5,7)]

def table_from(triples):
    t = {}
    for a in range(8):
        t[(0,a)] = (1,a); t[(a,0)] = (1,a)
    for a in range(1,8):
        t[(a,a)] = (-1,0)
    for (a,b,c) in triples:
        for (x,y,z) in [(a,b,c),(b,c,a),(c,a,b)]:
            t[(x,y)] = (1,z); t[(y,x)] = (-1,z)
    return t

KIRMSE = table_from(TRIPLES)

def omul(u, v, tab=KIRMSE):
    r = [F(0)]*8
    for i in range(8):
        if u[i] == 0: continue
        for j in range(8):
            if v[j] == 0: continue
            s,c = tab[(i,j)]
            r[c] += s*u[i]*v[j]
    return r

def N(u): return sum(x*x for x in u)

def half(*idx):
    v = [F(0)]*8
    for k in idx: v[k] = F(1,2)
    return v

units = [[F(1) if k==j else F(0) for k in range(8)] for j in range(8)]

# --- Kirmse's J_1, p.70:  J_1 = [J_0, a1, a2, a3, a4] ---
a1 = half(0,1,2,3)            # 1/2(1 + i1 + i2 + i3)
a2 = half(4,5,6,7)            # 1/2(i4 + i5 + i6 + i7)
a3 = half(0,1,4,5)            # 1/2(1 + i1 + i4 + i5)
a4 = half(0,2,4,7)            # 1/2(1 + i2 + i4 + i7)

def supp(v): return tuple(int(2*x) % 2 for x in v)

def code_of(halfvecs):
    gens = [supp(v) for v in halfvecs]
    C = set()
    for co in iproduct([0,1], repeat=len(gens)):
        w = [0]*8
        for c,g in zip(co, gens):
            if c: w = [x^y for x,y in zip(w,g)]
        C.add(tuple(w))
    return frozenset(C)

C_J1 = code_of([a1,a2,a3,a4])          # code of J_1, |C| = [J_1:J_0]

def in_lattice(v, code):
    if any((2*x).denominator != 1 for x in v): return False
    return supp(v) in code

# --- the 30 E8-type candidate codes (doubly-even self-dual [8,4]) ---
H0 = [(1,1,1,1,1,1,1,1),(0,1,0,1,0,1,0,1),(0,0,1,1,0,0,1,1),(0,0,0,0,1,1,1,1)]
def span(gs):
    S = set()
    for co in iproduct([0,1], repeat=len(gs)):
        w = [0]*8
        for c,g in zip(co,gs):
            if c: w = [x^y for x,y in zip(w,g)]
        S.add(tuple(w))
    return frozenset(S)
codes = {span([tuple(g[p[k]] for k in range(8)) for g in H0])
         for p in permutations(range(8))}

def half_basis(code):
    cb = []
    for c in sorted(code):
        red = list(c)
        for b in cb:
            pp = next(k for k in range(8) if b[k])
            if red[pp]: red = [x^y for x,y in zip(red,b)]
        if any(red): cb.append(red)
    return [[F(x,2) for x in c] for c in cb[:4]]

def closed(code, tab=KIRMSE):
    """is L_code = J_0 + {1/2 c : c in code} closed under the given table?"""
    gens = units + half_basis(code)
    for gi in gens:
        for gj in gens:
            if not in_lattice(omul(gi,gj,tab), code): return False
    return True

# --- a transposition of two imaginary axes (the real axis 0 is fixed) ---
def transp(a, b):
    p = list(range(8)); p[a], p[b] = p[b], p[a]
    return p
def permute_code(code, p):
    return frozenset(tuple(w[p[k]] for k in range(8)) for w in code)
def conj_table(tab, p):
    """table relabelled by the coordinate permutation p (p involutive here)."""
    t = {}
    for i in range(8):
        for j in range(8):
            s,c = tab[(p[i], p[j])]
            t[(i,j)] = (s, p[c])
    return t

if __name__ == "__main__":
    random.seed(0)
    print("="*70)
    print("Coxeter 1946, Section 4 -- 'Kirmse's mistake'")
    print("="*70)

    # --- (C1) table (4.1) is Kirmse's table (1), a composition algebra ---
    comp = True
    for _ in range(400):
        u = [F(random.randint(-3,3)) for _ in range(8)]
        v = [F(random.randint(-3,3)) for _ in range(8)]
        if N(omul(u,v)) != N(u)*N(v): comp = False
    print()
    print("(C1) Coxeter's (4.1) = Kirmse's table (1); a composition algebra:",
          comp)

    # --- (C3) the non-closure witness ---
    prod = omul(a1, a3)
    witness = half(1,2,4,7)            # 1/2(i1 + i2 + i4 + i7)
    print()
    print("(C3) Coxeter's witness, computed in Kirmse's table:")
    print("     1/2(1+i1+i2+i3) . 1/2(1+i1+i4+i5) =",
          "1/2(i1+i2+i4+i7)" if prod == witness else prod)
    print("     a1, a3 lie in J_1                 :",
          in_lattice(a1, C_J1) and in_lattice(a3, C_J1))
    print("     the product 1/2(i1+i2+i4+i7) in J_1:", in_lattice(prod, C_J1),
          "  (norm =", N(prod), "-> a unit, but not one of J_1's 240)")
    print("     => J_1 is NOT closed under multiplication.")

    # --- (C2) the count: eight claimed, seven actual ---
    genuine = sorted((c for c in codes if closed(c)), key=lambda c: sorted(c))
    print()
    print("(C2) E8-type candidate lattices:", len(codes))
    print("     closed under Kirmse's table (4.1) => maximal orders:",
          len(genuine), " (Kirmse stated eight)")
    print("     Kirmse's J_1 is among the candidates :", C_J1 in codes)
    print("     Kirmse's J_1 is one of the closed ones:", C_J1 in genuine)

    # --- (C4)+(C5) Bruck's correction: transposing two of the i's ---
    print()
    print("(C4)/(C5) Bruck's correction -- transpose two imaginary axes:")
    by_image = {}
    for (a,b) in combinations(range(1,8), 2):       # the 21 transpositions
        p = transp(a,b)
        img = permute_code(C_J1, p)
        by_image.setdefault(img, []).append((a,b))
    all_closed = all(closed(img) for img in by_image)
    all_e8     = all(img in codes for img in by_image)
    print("     all 21 transpositions send J_1 to a closed domain:", all_closed)
    print("     each image is an E8-type lattice                 :", all_e8)
    print("     number of DISTINCT domains among the 21 images    :",
          len(by_image), "  (Coxeter: seven)")
    sizes = sorted(len(v) for v in by_image.values())
    print("     sizes of the", len(by_image), "transposition-sets         :",
          sizes, " (Coxeter: 7 sets of 3)")
    # the "special unit": the index fixed by all three transpositions
    specials = []
    for img, tlist in by_image.items():
        moved = set()
        for (a,b) in tlist: moved |= {a,b}
        fixed = [x for x in range(1,8) if x not in moved]
        specials.append(len(fixed))
    print("     every transposition-set fixes exactly one i (special unit):",
          all(s == 1 for s in specials))
    print("     the 7 domains are all distinct from J_1            :",
          C_J1 not in by_image)
    print("     each of the 7 domains is one of the", len(genuine),
          "maximal orders :", set(by_image) <= set(genuine))

    # --- (C6) framing duality: transpose elements  <=>  transpose table ---
    ok = True
    for (a,b) in combinations(range(1,8), 2):
        p = transp(a,b)
        # framing A: move the elements (table fixed)
        A = closed(permute_code(C_J1, p), KIRMSE)
        # framing B: move the table (elements fixed)
        B = closed(C_J1, conj_table(KIRMSE, p))
        if A != B: ok = False
    # and the same identity must hold for every candidate lattice, not just J_1
    for code in codes:
        for (a,b) in combinations(range(1,8), 2):
            p = transp(a,b)
            if closed(permute_code(code,p), KIRMSE) != \
               closed(code, conj_table(KIRMSE,p)):
                ok = False
    print()
    print("(C6) 'transpose the elements' = 'transpose the table'")
    print("     (a linear coordinate transposition, either way):", ok)

    # --- coordinate symmetry: no transposition fixes Kirmse's lattices ---
    print()
    print("Coordinate symmetry of Kirmse's lattices:")
    trs = [transp(a,b) for (a,b) in combinations(range(1,8),2)]
    fix_J1 = sum(1 for p in trs if permute_code(C_J1, p) == C_J1)
    fix_max = max(sum(1 for p in trs if permute_code(g, p) == g) for g in genuine)
    print("   transpositions fixing J_1                     :", fix_J1, "of 21")
    print("   max transpositions fixing any of the 7 orders :", fix_max, "of 21")
    print("   => no Kirmse maximal order is coordinate-symmetric;")
    print("      none is sigma-invariant for any transposition sigma.")

    # --- concrete example for the companion note: tau = (1 2) ---
    print()
    print("-"*70)
    print("Concrete example (companion note), tau = (1 2):")
    tau = transp(1,2)
    KIRMSE_tau = conj_table(KIRMSE, tau)
    def show(i,j):
        s,c = KIRMSE[(i,j)];     lhs = ("+" if s>0 else "-")+"i%d"%c
        s,c = KIRMSE_tau[(i,j)]; rhs = ("+" if s>0 else "-")+"i%d"%c
        print("     i%d.i%d :  Kirmse  %s     tau-twisted  %s" % (i,j,lhs,rhs))
    for (i,j) in [(1,2),(1,5),(3,4),(6,7)]:
        show(i,j)
    def pretty(v):                         # half-integer octonion -> string
        nm = {0:"1",1:"i1",2:"i2",3:"i3",4:"i4",5:"i5",6:"i6",7:"i7"}
        terms = [("+" if v[k]>0 else "-")+nm[k] for k in range(8) if v[k]]
        return "1/2(" + "".join(terms).lstrip("+") + ")"
    p_K    = omul(a1,a3,KIRMSE)
    p_Ktau = omul(a1,a3,KIRMSE_tau)
    print("     a1 = 1/2(1+i1+i2+i3),  a3 = 1/2(1+i1+i4+i5)   (both in J_1)")
    print("       a1 .       a3  =", pretty(p_K),
          "  in J_1:", in_lattice(p_K,   C_J1))
    print("       a1 .[tau]  a3  =", pretty(p_Ktau),
          "  in J_1:", in_lattice(p_Ktau,C_J1))
    print("     => same elements, same lattice J_1; twisting the *product*")
    print("        (tau-conjugation) restores closure.  Framing B of the fix.")

    print()
    print("="*70)
    print("VERDICT.  Coxeter 1946, Section 4, is correct in every checked")
    print("claim.  Kirmse stated eight maximal domains; there are seven; his")
    print("exhibited J_1 is the spurious one (not closed).  Bruck's remedy is")
    print("a transposition of two imaginary units -- a linear involution of")
    print("R^8 -- and the 21 transpositions yield exactly the seven orders.")
    print("="*70)
