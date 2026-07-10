"""Validation + structure of S = Z-span(Lambda star Lambda), following
verify_product_span_index.py (which found [Lambda : S] = 65536 = 2^16).

Checks:
  (a) Consistency: random products u star v (u, v random Z-combinations
      of the Lambda basis) lie in S.  (Must hold by bilinearity; failure
      would indicate a bug in the span computation.)
  (b) Is 2*Lambda contained in S?  If yes, S/2Lambda is an F_2-subspace
      of Lambda/2Lambda of dimension 24 - 16 = 8, i.e., S is the full
      preimage of an 8-dimensional subspace of Lambda/2Lambda -- a
      2-primary structural statement in the same vein as Section A.1.
  (c) The HNF diagonal of S (in Lambda coordinates), showing where the
      index 2^16 sits.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
from verify_consecutive_twists_exact import B, twisted, l3, detB, inLambda, add

SIGMA = (0, 2, 1, 3, 4, 5, 6, 7)
star = twisted(SIGMA)

def scale(u, c):
    return tuple([c * x for x in blk] for blk in u)

def lincomb(coeffs):
    acc = None
    for c, b in zip(coeffs, B):
        t = scale(b, c)
        acc = t if acc is None else tuple(add(x, y) for x, y in zip(acc, t))
    return acc

if __name__ == "__main__":
    print("=" * 70)
    print("Validation + structure of S = Z-span(Lambda star Lambda)")
    print("=" * 70)

    prods = [l3(star(bi, bj)) for bi in B for bj in B]
    G = Matrix(prods).T
    H = hermite_normal_form(G)              # 24 x 24 HNF basis of S (L^3 coords)
    det_S = abs(H.det())
    print(f"[L^3 : S] = {det_S}   [Lambda : S] = {det_S // detB}")

    Hinv = H.inv()                           # exact rational inverse

    def in_S(vec24):
        c = Hinv * Matrix(vec24)
        return all(x.q == 1 for x in c)

    # (a) random products lie in S
    random.seed(20260607)
    ok = 0
    N_RAND = 20
    for _ in range(N_RAND):
        u = lincomb([random.randint(-2, 2) for _ in range(24)])
        v = lincomb([random.randint(-2, 2) for _ in range(24)])
        p = star(u, v)
        assert inLambda(p), "product left Lambda?!"
        if in_S(l3(p)):
            ok += 1
    print(f"(a) random products in S: {ok}/{N_RAND}  (expect {N_RAND}/{N_RAND})")

    # (b) 2*Lambda subset of S?
    two_in = sum(1 for b in B if in_S(l3(scale(b, 2))))
    print(f"(b) 2*basis vectors of Lambda in S: {two_in}/24")
    if two_in == 24:
        print("    => 2*Lambda is contained in S; S is the full preimage of an")
        print("       F_2-subspace of Lambda/2Lambda of dimension 24 - 16 = 8.")

    # (c) HNF diagonal in Lambda coordinates: express S-basis in B-coords
    Bmat = Matrix([l3(b) for b in B]).T      # 24 x 24, columns = Lambda basis
    C = Bmat.inv() * H                       # S-basis in Lambda coordinates
    assert all(x.q == 1 for x in C), "S-basis not integral in Lambda coords?!"
    HC = hermite_normal_form(C.applyfunc(lambda x: int(x)))
    diag = [HC[i, i] for i in range(24)]
    print(f"(c) HNF diagonal of S in Lambda coordinates:")
    print(f"    {diag}")
    from math import prod
    print(f"    product = {prod(diag)}  (expect 65536 = 2^16)")
