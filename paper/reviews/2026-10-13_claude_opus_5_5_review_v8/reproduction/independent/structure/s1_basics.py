"""
s1_basics.py -- basic objects and the facts the closure criterion needs.

Conventions: see leech_core.py (octonion table from the paper's Fano triples,
L = D8+ with the ODD half-integer coset, s, sbar, LEFT multiples Ls, Lsbar,
Wilson's Lambda in L^3, L-coordinates w.r.t. the paper's Appendix-A basis).

Checks (all exact):
  (a) the multiplication table is an octonion algebra (alternative laws and
      N(xy) = N(x)N(y) on random integer vectors; sbar = conj(s), s^2+s+2 = 0);
  (b) the 8 basis vectors lie in D8+ and have determinant +-1, so they are a
      Z-basis of L; the Gram matrix is even unimodular (E8);
  (c) L.L in L: all 64 basis products are integral in L-coordinates, and all
      240 x 240 root products lie in D8+ by the coordinate definition;
  (d) M := Z-span(L.L) and its index in L;
  (e) Lsbar, Ls in L, index 16 each, 2L in both, Lsbar + Ls = L,
      Lsbar cap Ls = 2L;
  (f) V0 = Lsbar/2L, W0 = Ls/2L: complementary 4-dim, totally singular for
      q(v) = N(v)/2 mod 2;
  (g) Wilson's Lambda: Z-basis, every basis vector satisfies Wilson's
      conditions by exact lattice membership, the 12 mod-2 conditions are
      independent (so the basis spans all of Lambda), Gram det 2^24, Gram/2
      even unimodular;
  (h) the mod-2 algebra A = L/2L: commutative? unital?
Outputs: s1_basics.out (this log), lambda_basis_Lcoords.txt,
lambda_basis_ecoords_doubled.txt, lambda_gram.gp (for the PARI check).
"""
import random
import sys
from fractions import Fraction as Fr
import numpy as np
import leech_core as C

out_lines = []


def log(*a):
    s = " ".join(str(t) for t in a)
    print(s)
    out_lines.append(s)


# (a) octonion sanity -------------------------------------------------------
random.seed(1)
ok_alt = ok_norm = True
for _ in range(300):
    x = [random.randint(-3, 3) for _ in range(8)]
    y = [random.randint(-3, 3) for _ in range(8)]
    if C.omul(C.omul(x, x), y) != C.omul(x, C.omul(x, y)):
        ok_alt = False
    if C.omul(C.omul(y, x), x) != C.omul(y, C.omul(x, x)):
        ok_alt = False
    if C.norm(C.omul(x, y)) != C.norm(x) * C.norm(y):
        ok_norm = False
log("(a) left+right alternative on 300 random pairs:", ok_alt, "| N(xy)=N(x)N(y):", ok_norm)
s, sb = C.S_VEC, C.SBAR_VEC
log("    sbar == conj(s):", sb == C.conj(s),
    "| s^2 + s + 2 == 0:", [a + b + (2 if i == 0 else 0) for i, (a, b) in enumerate(zip(C.omul(s, s), s))] == [0] * 8)
log("    s in L:", C.in_L_e(s), "| sbar in L:", C.in_L_e(sb), "| e0 in L:", C.in_L_e([1, 0, 0, 0, 0, 0, 0, 0]))

# (b) the basis -------------------------------------------------------------
detB = C.mat_det(C.L_BASIS)
log("(b) basis vectors in D8+:", all(C.in_L_e(v) for v in C.L_BASIS), "| det(basis) =", detB)
G = np.array(C.GRAM_L, dtype=np.int64)
log("    Gram(L) det =", round(np.linalg.det(G)), "(exact:", C.mat_det(C.GRAM_L), ") even diagonal:", all(G[i, i] % 2 == 0 for i in range(8)))

# (c) L.L in L ----------------------------------------------------------------
roots = C.roots_of_L()
bad = 0
for r1 in roots:
    for r2 in roots:
        if not C.in_L_e(C.omul(r1, r2)):
            bad += 1
log("(c) L.L in L: 64 basis products integral in L-coords: True (asserted when building T_STD);",
    "root x root products outside D8+:", bad, "of", 240 * 240)

# (d) span of L.L -------------------------------------------------------------
prods = [list(C.T_STD[i, j]) for i in range(8) for j in range(8)]
idxM = C.hnf_det(prods, 8)
log("(d) [L : Z-span(L.L)] =", idxM)

# (e) Lsbar and Ls ------------------------------------------------------------
gen_sb = [C.to_L(C.omul(C.L_BASIS[k], sb)) for k in range(8)]
gen_s = [C.to_L(C.omul(C.L_BASIS[k], s)) for k in range(8)]
log("(e) L_k.sbar in L for all k:", all(C.is_integral(v) for v in gen_sb),
    "| L_k.s in L for all k:", all(C.is_integral(v) for v in gen_s))
H_sb = C.hnf([[int(t) for t in v] for v in gen_sb])
H_s = C.hnf([[int(t) for t in v] for v in gen_s])
idx_sb = C.hnf_det(H_sb, 8)
idx_s = C.hnf_det(H_s, 8)
log("    [L:Lsbar] =", idx_sb, " [L:Ls] =", idx_s)


def in_lattice(c, H):
    sol = C.vecmat([Fr(t) for t in c], C.mat_inv(H))
    return C.is_integral(sol)


twoL = [[2 * int(i == j) for j in range(8)] for i in range(8)]
log("    2L in Lsbar:", all(in_lattice(v, H_sb) for v in twoL), "| 2L in Ls:", all(in_lattice(v, H_s) for v in twoL))
log("    [L : Lsbar + Ls] =", C.hnf_det(H_sb + H_s, 8),
    " => [L : Lsbar cap Ls] = 16*16/1 = 256 = [L:2L], so Lsbar cap Ls = 2L")

# (f) V0, W0 ------------------------------------------------------------------
V0 = C.span_bits([C.to_bits(v) for v in H_sb])
W0 = C.span_bits([C.to_bits(v) for v in H_s])
log("(f) dim V0 =", C.dim_bits(V0), " dim W0 =", C.dim_bits(W0), " V0 cap W0 = {0}:", V0 & W0 == {0})
log("    q = 0 on V0:", all(C.q_form(v) == 0 for v in V0), " q = 0 on W0:", all(C.q_form(v) == 0 for v in W0))
log("    V0 basis (bitmask, bit k <-> L_{k+1}):", C.basis_bits(V0), " W0 basis:", C.basis_bits(W0))

# (g) Wilson's Lambda ---------------------------------------------------------
U = C.leech_basis_from_pair(V0, W0)
detU = C.hnf_det(U, 24)
log("(g) Lambda basis: 24 x 24 HNF in L^3-coords, [L^3 : Lambda] =", detU)
# Wilson's conditions by exact lattice membership (not via the mod-2 shortcut)
okW = True
for u in U:
    x, y, z = u[0:8], u[8:16], u[16:24]
    add = lambda a, b: [p + q for p, q in zip(a, b)]
    okW &= in_lattice(add(x, y), H_sb) and in_lattice(add(x, z), H_sb) and in_lattice(add(y, z), H_sb)
    okW &= in_lattice(add(add(x, y), z), H_s)
log("    every basis vector satisfies Wilson's conditions (exact membership):", okW)
# the conditions are 12 independent F_2-linear conditions on (L/2L)^3 ------------
# annihilators of V0 and W0 (functionals f in F_2^8 with f.v = 0 on the subspace)
def annihilator(S):
    return [f for f in range(256) if all(bin(f & v).count('1') % 2 == 0 for v in S)]
annV, annW = annihilator(V0), annihilator(W0)
fb = []
for f in C.basis_bits(frozenset(annV)):
    fb.append(f | (f << 8))            # f(x + y)
    fb.append(f | (f << 16))           # f(x + z)
for f in C.basis_bits(frozenset(annW)):
    fb.append(f | (f << 8) | (f << 16))  # f(x + y + z)
rank = 0
rows = list(fb)
for bit in range(24):
    piv = next((r for r in rows if (r >> bit) & 1), None)
    if piv is None:
        continue
    rows.remove(piv)
    rows = [r ^ piv if (r >> bit) & 1 else r for r in rows]
    rank += 1
log("    rank of the 12 mod-2 conditions on (L/2L)^3 =", rank,
    "=> [L^3 : Lambda] = 2^%d = %d, equal to the index of the HNF basis: so the basis spans Lambda" % (rank, 2 ** rank))
assert detU == 2 ** rank

GL3 = np.kron(np.eye(3, dtype=np.int64), G)
U_np = np.array(U, dtype=np.int64)
GramLam = U_np @ GL3 @ U_np.T
detLam = C.mat_det(GramLam.tolist())
log("    Gram(Lambda) det =", detLam, "= 2^%d" % (int(detLam).bit_length() - 1))
half = GramLam // 2
log("    Gram/2 integral:", bool(np.all(GramLam % 2 == 0)), " diag(Gram/2) even:", bool(np.all(np.diag(half) % 2 == 0)),
    " det(Gram/2) =", C.mat_det(half.tolist()))
log("    => Lambda/sqrt2 is even unimodular of rank 24 (min norm checked with PARI in s1_leech_min.gp)")

# (h) mod-2 algebra -----------------------------------------------------------
tab = C.mul_table_mod2(C.T_STD)
comm = all(tab[u, v] == tab[v, u] for u in range(256) for v in range(256))
units = [e for e in range(256) if all(tab[e, v] == v and tab[v, e] == v for v in range(256))]
log("(h) A = L/2L: commutative:", comm, "| two-sided identity elements:", units)
AA = C.span_bits([int(tab[u, v]) for u in range(256) for v in range(256)])
log("    dim span(A.A) =", C.dim_bits(AA))

# outputs ---------------------------------------------------------------------
with open('lambda_basis_Lcoords.txt', 'w') as f:
    f.write("# Z-basis of Wilson's Lambda, rows in L^3-coordinates (L-basis of leech_core.py)\n")
    for u in U:
        f.write(" ".join(map(str, u)) + "\n")
with open('lambda_basis_ecoords_doubled.txt', 'w') as f:
    f.write("# same basis, rows = 2 * (e-coordinates of x | y | z), 24 integers\n")
    for u in U:
        e = []
        for blk in range(3):
            e += [int(2 * t) for t in C.from_L(u[8 * blk:8 * blk + 8])]
        f.write(" ".join(map(str, e)) + "\n")
with open('lambda_gram.gp', 'w') as f:
    f.write("G = [" + ";".join(",".join(str(int(t)) for t in row) for row in GramLam) + "];\n")

with open('s1_basics.out', 'w') as f:
    f.write("\n".join(out_lines) + "\n")
