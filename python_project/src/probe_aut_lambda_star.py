"""First-step probe of Aut(Lambda, +, star): test specific candidate
automorphisms in the Conway group Co_0 = Aut(Lambda, +) for whether they
preserve the Z_3-symmetric triple-octonion product star.

Each candidate g is tested on all 24 * 24 = 576 basis-pair products: g is
in Aut(Lambda, +, star) iff g(Lambda) = Lambda and g(B_i star B_j) =
g(B_i) star g(B_j) for every pair of basis vectors B_i, B_j of Lambda.

Uses the exact Lambda Z-basis and star machinery from
verify_consecutive_twists_exact.py (loaded as a module).
"""
from fractions import Fraction as F
from itertools import combinations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_consecutive_twists_exact import (
    B, inLambda, twisted, omul, add, vec, h,
)

# The paper's sigma = (1 2)
SIGMA = (0, 2, 1, 3, 4, 5, 6, 7)
star = twisted(SIGMA)

def neg(u):
    return tuple([-c for c in blk] for blk in u)

def block_cycle(u):                     # (x,y,z) -> (y,z,x)
    x, y, z = u
    return (y, z, x)

def block_swap_12(u):                   # (x,y,z) -> (y,x,z)
    x, y, z = u
    return (y, x, z)

def apply_sigma_block(u):               # (sigma(x), sigma(y), sigma(z))
    def s(v):
        r = [F(0)]*8
        for j in range(8):
            r[SIGMA[j]] = v[j]
        return r
    return tuple(s(blk) for blk in u)

def test(name, g):
    """Test whether g is in Aut(Lambda, +, star).
    Returns: (preserves_lambda, preserves_star, witness)."""
    # 1. g preserves Lambda?
    pres_l = all(inLambda(g(b)) for b in B)
    if not pres_l:
        return (False, None, "g moves Lambda (a basis vector leaves Lambda)")
    # 2. g preserves star?  Check on all 576 basis-pair products.
    for i, bi in enumerate(B):
        for j, bj in enumerate(B):
            lhs = g(star(bi, bj))
            rhs = star(g(bi), g(bj))
            if lhs != rhs:
                return (True, False, f"basis pair (B_{i+1}, B_{j+1}) violates star: lhs != rhs")
    return (True, True, None)

if __name__ == "__main__":
    print("="*70)
    print("Probe of Aut(Lambda, +, star): test specific candidates")
    print("="*70)
    candidates = [
        ("Z_3 cyclic block permutation (x,y,z) -> (y,z,x)", block_cycle),
        ("S_3 transposition of blocks 1<->2", block_swap_12),
        ("componentwise negation -I_24", neg),
        ("sigma applied to each block", apply_sigma_block),
    ]
    in_group = []
    out_group = []
    for name, g in candidates:
        pl, ps, w = test(name, g)
        if pl is False:
            print(f"  {name}: g does NOT preserve Lambda -> NOT in Co_0")
            out_group.append(name)
        elif ps is False:
            print(f"  {name}: preserves Lambda, does NOT preserve star -> in Co_0, NOT in Aut(Lambda,+,star)")
            out_group.append(name)
            print(f"      witness: {w}")
        else:
            print(f"  {name}: preserves Lambda AND star -> IN Aut(Lambda,+,star)")
            in_group.append(name)

    print()
    print("Lower-bound subgroup of Aut(Lambda, +, star) so far:")
    for n in in_group: print("  +", n)
    print("Excluded:")
    for n in out_group: print("  -", n)
