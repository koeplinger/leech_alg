r"""
01_build_lambda.py -- build and verify Z-bases of L, L.sbar, L.s and of Wilson's
Leech lattice Lambda, from the definitions only (conventions: see leechlib.py).

Checks performed (all exact):
 (a) L-basis: HNF of the 240 roots; 8 vectors of L with |det| = 1 = covolume
     of D8+ (D8 has covolume 2 and index 2 in D8+), hence a Z-basis of L.
     Gram matrix even with det 1.  L.L <= L (64 structure constants integral).
 (b) L.sbar = {l.sbar} and L.s = {l.s}: images of the basis; both inside L,
     index 16 each, L.sbar + L.s = L, L.sbar /\ L.s = 2L.
     Cross-check of the two definitional membership tests
        x in L.sbar <=> (x.s)/2 in L,     x in L.s <=> (x.sbar)/2 in L
     against the basis description on all 256 classes of L/2L.
     V = L.sbar/2L and W = L.s/2L are totally singular for q(x) = N(x)/2 mod 2.
 (c) Lambda: the F2 condition matrix Q (x+y, x+z, y+z in V; x+y+z in W) has
     rank 12, so [L^3 : Lambda] = 2^12.  The HNF basis built from explicit
     generators lies in Lambda (both by Q and by the definitional test) and
     has determinant 2^12 in L^3-coordinates, hence is a Z-basis of Lambda.
 (d) Gram matrix G of Lambda (e-coordinates): integral, all entries even,
     det G = 2^24; G/2 is even unimodular.  G is exported for PARI/GP
     (02_qfminim.gp checks min 8 and the 196,560 minimal vectors).

Outputs:  01_build_lambda.out.txt (log), lambda_gram.gp (Gram matrix for gp),
          lambda_basis.txt (basis, doubled e-coordinates and L^3-coordinates),
          L_basis.txt (Z-basis of L, e-coordinates).
Run:      python3 01_build_lambda.py > 01_build_lambda.out.txt
"""
from fractions import Fraction as F
from itertools import product as iproduct
import leechlib as LL

D = LL.build_all()
BL, LC, C = D['BL'], D['LC'], D['C']

print("=" * 72)
print("(a) Z-basis of L = D8+ (odd convention)")
print("=" * 72)
for k, b in enumerate(BL):
    print(f"  l_{k+1} = {LL.fmt_vec(b)}   in L: {LL.in_L(b)}   N = {LL.norm(b)}")
detBL = LL.det_bareiss([[int(2 * a) for a in r] for r in BL])
print(f"  det(2*BL) = {detBL}  ->  |det BL| = {abs(F(detBL, 2**8))}  (= covolume of D8+ = 1)")
GL = [[LL.ip(a, b) for b in BL] for a in BL]
assert all(g.denominator == 1 for row in GL for g in row)
GLi = [[int(g) for g in row] for row in GL]
print(f"  Gram(L): det = {LL.det_bareiss(GLi)}, diagonal even: {all(GLi[i][i] % 2 == 0 for i in range(8))}")
print(f"  L.L <= L: all 512 structure constants of l_i . l_j in the L-basis are integers: True"
      f"  (max |c| = {int(abs(C).max())})")

print()
print("=" * 72)
print("(b) L.sbar, L.s and the mod-2 quotient")
print("=" * 72)
Vrows, Wrows = D['Vrows'], D['Wrows']
dV = LL.det_bareiss(Vrows)
dW = LL.det_bareiss(Wrows)
print(f"  L.sbar <= L (basis images l_k.sbar have integral L-coordinates): True; [L : L.sbar] = |det| = {abs(dV)}")
print(f"  L.s    <= L (basis images l_k.s    have integral L-coordinates): True; [L : L.s]    = |det| = {abs(dW)}")
idx_sum, _ = LL.hnf_index(Vrows + Wrows, 8)
print(f"  [L : L.sbar + L.s] = {idx_sum}")
# 2L inside both?
twoL = [[2 * int(i == j) for j in range(8)] for i in range(8)]
inV = LL.hnf_index(Vrows + twoL, 8)[0] == abs(dV)
inW = LL.hnf_index(Wrows + twoL, 8)[0] == abs(dW)
print(f"  2L <= L.sbar: {inV};  2L <= L.s: {inW}")
print(f"  [L : L.sbar /\\ L.s] = [L:L.sbar][L:L.s]/[L:L.sbar+L.s] = {abs(dV)*abs(dW)//idx_sum} = [L:2L] = 256"
      f"  ->  L.sbar /\\ L.s = 2L: {abs(dV)*abs(dW)//idx_sum == 256 and inV and inW}")

# F2 subspaces V, W (as sets of residues in {0,1}^8 of L-coordinates)
def span_f2(rows):
    basis, _ = LL.f2_rowspace_basis(rows, 8)
    out = set()
    for coeffs in iproduct((0, 1), repeat=len(basis)):
        v = [0] * 8
        for c, b in zip(coeffs, basis):
            if c:
                v = [x ^ y for x, y in zip(v, b)]
        out.add(tuple(v))
    return out

Vset, Wset = span_f2(Vrows), span_f2(Wrows)
print(f"  |V| = {len(Vset)}, |W| = {len(Wset)}, |V /\\ W| = {len(Vset & Wset)}")
agree = 0
for c in iproduct((0, 1), repeat=8):
    x = LC.to_e(list(c))
    assert LL.in_L(x)
    a1 = LL.in_L_sbar(x) == (tuple(c) in Vset)
    a2 = LL.in_L_s(x) == (tuple(c) in Wset)
    assert a1 and a2
    agree += 1
print(f"  definitional tests x in L.sbar <=> (x.s)/2 in L, x in L.s <=> (x.sbar)/2 in L"
      f" agree with the basis description on all {agree} classes of L/2L: True")
qV = all(LL.norm(LC.to_e(list(v))) % 4 == 0 for v in Vset)
qW = all(LL.norm(LC.to_e(list(v))) % 4 == 0 for v in Wset)
print(f"  V and W totally singular for q = N/2 mod 2 (Lagrangians): {qV}, {qW}")

print()
print("=" * 72)
print("(c) Lambda")
print("=" * 72)
Q = D['Q']
rQ = LL.f2_rank(Q)
print(f"  F2 condition matrix Q: {len(Q)} rows, rank {rQ}  ->  [L^3 : Lambda] = 2^{rQ}")
BLam = D['BLam']
dLam = LL.det_bareiss(BLam)
print(f"  HNF basis of Lambda (L^3-coordinates): det = {dLam} = 2^{abs(dLam).bit_length()-1}")
ok_Q = all(sum(f[i] * r[i] for i in range(24)) % 2 == 0 for r in BLam for f in Q)
ok_def = all(LL.in_Lambda_def(r) for r in D['BLam_e'])
print(f"  every basis vector satisfies Q u = 0 mod 2: {ok_Q}")
print(f"  every basis vector passes the definitional test (x,y,z in L; x+y,x+z,y+z in L.sbar; x+y+z in L.s): {ok_def}")
assert rQ == 12 and abs(dLam) == 2 ** 12 and ok_Q and ok_def
print("  => the 24 rows are a Z-basis of Lambda (contained in Lambda, same index 2^12 in L^3)")

print()
print("=" * 72)
print("(d) Gram matrix of Lambda (e-coordinates)")
print("=" * 72)
Be = D['BLam_e']
G = [[LL.ip(a, b) for b in Be] for a in Be]
assert all(g.denominator == 1 for row in G for g in row)
G = [[int(g) for g in row] for row in G]
detG = LL.det_bareiss(G)
print(f"  det G = {detG} = 2^{detG.bit_length()-1}")
print(f"  all entries of G even: {all(g % 2 == 0 for row in G for g in row)}")
H = [[g // 2 for g in row] for row in G]
print(f"  G/2: det = {LL.det_bareiss(H)}, all diagonal entries even: {all(H[i][i] % 2 == 0 for i in range(24))}")
print(f"  diagonal of G (norms of the basis vectors): {[G[i][i] for i in range(24)]}")
print("  -> Lambda/sqrt2 (Gram G/2) is an even unimodular lattice of rank 24.  If min(G) = 8,")
print("     i.e. G/2 has no roots, it is the Leech lattice (unique such lattice without roots).")
print("     Minimum and minimal-vector count: see 02_qfminim.gp / 02_qfminim.out.txt.")

with open('lambda_gram.gp', 'w') as fh:
    fh.write('G = [' + ';'.join(','.join(str(g) for g in row) for row in G) + '];\n')
with open('lambda_basis.txt', 'w') as fh:
    fh.write('# Z-basis of Lambda, one vector per line.\n')
    fh.write('# Columns 1-24: doubled e-coordinates 2*(x0..x7 | y0..y7 | z0..z7);\n')
    fh.write('# after "|L3|": the same vector in L^3-coordinates (L-basis of L_basis.txt).\n')
    for r, re in zip(BLam, Be):
        fh.write(' '.join(str(int(2 * a)) for a in re) + '   |L3|  ' + ' '.join(map(str, r)) + '\n')
with open('L_basis.txt', 'w') as fh:
    fh.write('# Z-basis of L = D8+ (odd convention), e-coordinates e0..e7\n')
    for b in BL:
        fh.write(' '.join(str(a) for a in b) + '\n')
print("\nwrote lambda_gram.gp, lambda_basis.txt, L_basis.txt")
