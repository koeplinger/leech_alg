"""Item (2): the mod-2-quotient nature of appendix Tables A.2 and A.3.

Tables A.2 and A.3 are the explicit integer certificates for Lemmas 4.3 and
4.4 (L . sigma(Ls_bar) <= sigma(Ls_bar)  and  sigma(Ls).sigma(Ls) <= sigma(Ls)).
This script shows what those two lemmas *are*, structurally: statements about
the 8-dimensional F2-bilinear algebra obtained from L by the mod-2 quotient
L -> L/2L.

  - L_bar s ∩ Ls = 2L (Wilson), so 2L is contained in sigma(Ls) and
    sigma(Ls_bar);
  - L.L <= L and Z-bilinearity make the octonion product descend to an
    F2-bilinear product mu-bar on L_bar := L/2L = F2^8;
  - sigma(Ls_bar) and sigma(Ls) are full preimages of 4-dim subspaces
    V, W of L_bar, with L_bar = V (+) W;
  - Lemma 4.3  <=>  mu-bar(L_bar, V) <= V   (V a LEFT IDEAL of L_bar);
  - Lemma 4.4  <=>  mu-bar(W, W)     <= W   (W a SUBALGEBRA of L_bar).

The paper's Section A.1 records two further structural facts (verified
empirically in verify_discovery2_V_isotropic.py): V is in fact a TWO-SIDED
ideal of L_bar; under the natural plus-type quadratic form
q(v + 2L) := N(v)/2 mod 2, both V and W are LAGRANGIANS (maximal totally
isotropic subspaces). The decomposition L_bar = V (+) W is a Witt
decomposition into a complementary pair of Lagrangians -- a polarization
of the orthogonal F2-space (L_bar, q).

This is a mod-2 *quotient* (L -> L/2L) -- distinct from a 2-adic scaling or
saturation.
"""
from fractions import Fraction as F
from itertools import product as iproduct

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
        if u[i]==0: continue
        for j in range(8):
            if v[j]==0: continue
            s,c = MUL[(i,j)]
            r[c] += s*u[i]*v[j]
    return r
def vec(*c): return [F(x) for x in c]
h = F(1,2)
s    = vec(-h,h,h,h,h,h,h,h)
sbar = vec(-h,-h,-h,-h,-h,-h,-h,-h)
L = [s, vec(1,1,0,0,0,0,0,0), vec(1,-1,0,0,0,0,0,0), vec(1,0,1,0,0,0,0,0),
     vec(1,0,0,1,0,0,0,0), vec(1,0,0,0,1,0,0,0), vec(1,0,0,0,0,1,0,0),
     vec(1,0,0,0,0,0,1,0)]
Ls    = [omul(Lk,s)    for Lk in L]
Lsbar = [omul(Lk,sbar) for Lk in L]
def sigma(v): r=list(v); r[1],r[2]=r[2],r[1]; return r
sLs    = [sigma(v) for v in Ls]
sLsbar = [sigma(v) for v in Lsbar]

def coeffs(t,B):
    M = [[B[k][r] for k in range(8)] + [t[r]] for r in range(8)]
    for col in range(8):
        piv = next((r for r in range(col,8) if M[r][col]!=0), None)
        if piv is None: return None
        M[col],M[piv] = M[piv],M[col]
        pv = M[col][col]; M[col] = [x/pv for x in M[col]]
        for r in range(8):
            if r!=col and M[r][col]!=0:
                fa = M[r][col]; M[r] = [M[r][i]-fa*M[col][i] for i in range(9)]
    return [M[r][8] for r in range(8)]
def Lcoords(v):
    c = coeffs(v, L); assert all(x.denominator==1 for x in c); return [int(x) for x in c]
def integer_subset(A,B):
    for a in A:
        for b in B:
            c = coeffs(omul(a,b), B)
            if c is None or any(x.denominator!=1 for x in c): return False
    return True

# --- F2 linear algebra in L_bar = L/2L ---
def f2span(gens):
    basis = []
    for g in gens:
        v = list(g)
        for b in basis:
            p = next(k for k in range(8) if b[k])
            if v[p]: v = [x^y for x,y in zip(v,b)]
        if any(v):
            basis.append(v); basis.sort(key=lambda r: next(k for k in range(8) if r[k]))
    return basis
def in_f2span(v, basis):
    v = list(v)
    for b in basis:
        p = next(k for k in range(8) if b[k])
        if v[p]: v = [x^y for x,y in zip(v,b)]
    return not any(v)

# mu-bar structure constants: [L_i].[L_j] = Lcoords(L_i.L_j) mod 2
MUBAR = [[tuple(c%2 for c in Lcoords(omul(L[i],L[j]))) for j in range(8)]
         for i in range(8)]
def mubar(a,b):
    r = [0]*8
    for i in range(8):
        if not a[i]: continue
        for j in range(8):
            if not b[j]: continue
            r = [x^y for x,y in zip(r, MUBAR[i][j])]
    return tuple(r)

if __name__ == "__main__":
    # 2L is contained in sigma(Ls) and sigma(Ls_bar)
    twoL = [[2*c for c in Lk] for Lk in L]
    in_sLs    = all(coeffs(t,sLs)    and all(x.denominator==1 for x in coeffs(t,sLs))
                    for t in twoL)
    in_sLsbar = all(coeffs(t,sLsbar) and all(x.denominator==1 for x in coeffs(t,sLsbar))
                    for t in twoL)
    print("2L is contained in sigma(Ls):", in_sLs,
          "   and in sigma(Ls_bar):", in_sLsbar)

    # V, W as F2-subspaces of L_bar
    V = f2span([tuple(c%2 for c in Lcoords(v)) for v in sLsbar])   # sigma(Ls_bar)/2L
    W = f2span([tuple(c%2 for c in Lcoords(v)) for v in sLs])      # sigma(Ls)/2L
    spanV = {(0,)*8}
    for b in V: spanV |= {tuple(x^y for x,y in zip(t,b)) for t in spanV}
    spanW = {(0,)*8}
    for b in W: spanW |= {tuple(x^y for x,y in zip(t,b)) for t in spanW}
    directsum = (len(V)==4 and len(W)==4 and (spanV & spanW)=={(0,)*8})
    print(f"dim V = {len(V)}, dim W = {len(W)};  V (+) W = L/2L:", directsum)

    e = [tuple(1 if k==i else 0 for k in range(8)) for i in range(8)]
    # Lemma 4.3
    lem43_int = integer_subset(L, sLsbar)
    lem43_f2  = all(in_f2span(mubar(ei,v), V) for ei in e for v in spanV)
    print()
    print("Lemma 4.3  L . sigma(Ls_bar) <= sigma(Ls_bar)  [Table A.2]:")
    print("   integer certificate holds :", lem43_int)
    print("   mod-2 quotient: V is a LEFT IDEAL of (L/2L, mu-bar):", lem43_f2)
    # Lemma 4.4
    lem44_int = integer_subset(sLs, sLs)
    lem44_f2  = all(in_f2span(mubar(w1,w2), W) for w1 in spanW for w2 in spanW)
    print()
    print("Lemma 4.4  sigma(Ls).sigma(Ls) <= sigma(Ls)  [Table A.3]:")
    print("   integer certificate holds :", lem44_int)
    print("   mod-2 quotient: W is a SUBALGEBRA of (L/2L, mu-bar):", lem44_f2)
    print()
    print("=> Tables A.2 and A.3 certify, respectively, that V is a left")
    print("   ideal and W is a subalgebra of the octonion F2-algebra L/2L.")
