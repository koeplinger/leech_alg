"""Baez-Egan projected Jordan product phi on Lambda: exact span index
[Lambda : S_phi] of the sublattice S_phi := Z-span{ phi(u, v) : u, v
in Lambda }, companion to verify_product_span_index.py (which found
[Lambda : S] = 65536 = 2^16 for the sigma-twisted product star).

Paper Section 6, "Anatomy of the Baez--Egan closure": representing
u = (x,y,z), v = (x',y',z') in Lambda subset O^3 as off-diagonal
Hermitian matrices X, Y in h3(O), direct computation gives
    2(XY + YX) = D + R,
with D in Z^3 diagonal,
    D = 4(<x,x'> + <y,y'>, <x,x'> + <z,z'>, <y,y'> + <z,z'>)
(so tr D = 8 <u,v>), and R the off-diagonal residue, lying in Lambda.
The projected Jordan product is phi(u, v) := R.

Parametrization (Hermitian, zero diagonal; chosen so the diagonal of
2(XY + YX) matches the stated D verbatim -- the span index is invariant
under block-permutation conventions):
    X_12 = x, X_21 = xbar, X_13 = ybar, X_31 = y, X_23 = z, X_32 = zbar.
Expanding 2(XY + YX) entry-wise and reading the off-diagonal blocks
back through the same parametrization gives the closed form used here:
    phi(u, v) = ( 2(ybar zpbar + ypbar zbar),
                  2(zbar xpbar + zpbar xbar),
                  2(xbar ypbar + xpbar ybar) ),
which is manifestly symmetric in u <-> v (Jordan product).

Validations (all must pass before the index is computed):
  (1) D-formula: on 20 random pairs of Z-combinations of the Lambda
      basis, all nine entries of 2(XY + YX) computed from the octonion
      matrix product agree with D + R: diagonal entries are scalars in
      Z matching the stated D, off-diagonal entries are Hermitian and
      match the closed-form phi blocks.
  (2) Lambda-closure: phi(B_i, B_j) in Lambda for all 24 x 24 = 576
      basis pairs (the paper's "an order in the sense of Theorem 1.1").
  (3) Commutativity: phi(u, v) = phi(v, u) on 20 random pairs.

Then, exactly as in verify_product_span_index.py: the Hermite normal
form of the 576 basis-pair products in L^3 coordinates (integers)
gives |det| = [L^3 : S_phi]; dividing by [L^3 : Lambda] = 4096 gives
[Lambda : S_phi].  Finally, as in verify_product_span_structure.py,
checks whether 2*Lambda is contained in S_phi (membership of the 24
doubled basis vectors via the exact HNF inverse).  All arithmetic exact.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
from verify_consecutive_twists_exact import B, omul, l3, detB, inLambda, add

def conj(a):                              # octonion conjugation
    return [a[0]] + [-t for t in a[1:]]

def scl(a, c):                            # scalar multiple of an octonion
    return [c * t for t in a]

def ip(a, b):                             # Euclidean R^8 inner product
    return sum(p * q for p, q in zip(a, b))

def phi(u, v):
    """Off-diagonal residue R of 2(XY + YX), read back into O^3 blocks."""
    x, y, z = u
    xp, yp, zp = v
    rx = add(omul(conj(y), conj(zp)), omul(conj(yp), conj(z)))
    ry = add(omul(conj(z), conj(xp)), omul(conj(zp), conj(x)))
    rz = add(omul(conj(x), conj(yp)), omul(conj(xp), conj(y)))
    return (scl(rx, 2), scl(ry, 2), scl(rz, 2))

def scale(u, c):
    return tuple([c * t for t in blk] for blk in u)

def lincomb(coeffs):
    acc = None
    for c, b in zip(coeffs, B):
        t = scale(b, c)
        acc = t if acc is None else tuple(add(p, q) for p, q in zip(acc, t))
    return acc

def jordan_matrix_entries(u, v):
    """All nine entries of M = 2(XY + YX) from the octonion matrix
    product, with X, Y the parametrized Hermitian matrices of u, v."""
    x, y, z = u
    xp, yp, zp = v
    Xm = [[None, x, conj(y)], [conj(x), None, z], [y, conj(z), None]]
    Ym = [[None, xp, conj(yp)], [conj(xp), None, zp], [yp, conj(zp), None]]
    zero = [0] * 8
    M = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            acc = zero
            for k in range(3):
                if k != i and k != j:     # zero diagonals kill k=i, k=j terms
                    acc = add(acc, omul(Xm[i][k], Ym[k][j]),
                                   omul(Ym[i][k], Xm[k][j]))
            M[i][j] = scl(acc, 2)
    return M

if __name__ == "__main__":
    print("=" * 70)
    print("Baez-Egan phi: index of S_phi = Z-span(phi(Lambda, Lambda)) in Lambda")
    print("=" * 70)
    print(f"[L^3 : Lambda] = {detB}  (expect 4096 = 2^12)")

    # --- validation (1): D-formula and R identification, 20 random pairs
    random.seed(20260710)
    N_RAND = 20
    pairs = []
    for _ in range(N_RAND):
        u = lincomb([random.randint(-2, 2) for _ in range(24)])
        v = lincomb([random.randint(-2, 2) for _ in range(24)])
        pairs.append((u, v))
    ok_D = 0
    for u, v in pairs:
        x, y, z = u
        xp, yp, zp = v
        M = jordan_matrix_entries(u, v)
        D = [4 * (ip(x, xp) + ip(y, yp)),
             4 * (ip(x, xp) + ip(z, zp)),
             4 * (ip(y, yp) + ip(z, zp))]
        diag_ok = all(M[i][i] == [D[i]] + [0] * 7 and D[i].denominator == 1
                      for i in range(3))
        herm_ok = all(M[j][i] == conj(M[i][j])
                      for i in range(3) for j in range(3) if i != j)
        px, py, pz = phi(u, v)
        resid_ok = (M[0][1] == px and M[2][0] == py and M[1][2] == pz)
        if diag_ok and herm_ok and resid_ok:
            ok_D += 1
    print(f"(1) D-formula match (diag = stated D in Z, off-diag Hermitian = phi):"
          f" {ok_D}/{N_RAND}")

    # --- validation (2): Lambda-closure on all 576 basis pairs
    prods = []
    n_in = 0
    for bi in B:
        for bj in B:
            p = phi(bi, bj)
            if inLambda(p):
                n_in += 1
            prods.append(l3(p))
    print(f"(2) Lambda-closure: {n_in}/576 basis-pair products in Lambda")

    # --- validation (3): commutativity on 20 random pairs
    ok_comm = sum(1 for u, v in pairs if phi(u, v) == phi(v, u))
    print(f"(3) commutativity phi(u,v) = phi(v,u): {ok_comm}/{N_RAND}")

    if not (ok_D == N_RAND and n_in == 576 and ok_comm == N_RAND):
        print("VALIDATION FAILED -- not computing the index.")
        sys.exit(1)

    # --- span index via HNF, as in verify_product_span_index.py
    G = Matrix(prods).T                     # 24 x 576, columns = products
    H = hermite_normal_form(G)              # 24 x 24
    det_S = abs(H.det())                    # [L^3 : S_phi]
    print(f"[L^3 : S_phi]    = {det_S}")
    idx = det_S // detB
    print(f"[Lambda : S_phi] = {det_S} / {detB} = {idx}")
    if idx == 1:
        print("=> S_phi = Lambda: the phi products Z-SPAN the full Leech lattice.")
    else:
        print(f"=> S_phi is a proper sublattice of index {idx}.")

    # --- is 2*Lambda contained in S_phi?  (verify_product_span_structure.py)
    Hinv = H.inv()                           # exact rational inverse

    def in_S(vec24):
        c = Hinv * Matrix(vec24)
        return all(t.q == 1 for t in c)

    two_in = sum(1 for b in B if in_S(l3(scale(b, 2))))
    print(f"2*basis vectors of Lambda in S_phi: {two_in}/24")
    if two_in == 24:
        print("=> 2*Lambda is contained in S_phi; S_phi is the full preimage of")
        print("   an F_2-subspace of Lambda/2Lambda.")
