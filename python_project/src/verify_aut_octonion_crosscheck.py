"""Task B, part 3: independent cross-checks of the three octonion-automorphism
stabilisers that drive the Aut(Lambda,+,star) enumeration.

verify_aut_lambda_star.py enumerates, completely,
        Stab(L)      := Aut(O,.s) ^ {A : A(L)     = L}      order  1344
        Stab(Ls)     := Aut(O,.s) ^ {A : A(Ls)    = Ls}     order   168
        Stab(Lsbar)  := Aut(O,.s) ^ {A : A(Lsbar) = Lsbar}  order 12096
by a backtracking search over the images of three minimal vectors that
generate (O,.s) as an algebra.  The final answer rests on those enumerations
being COMPLETE, so they are cross-checked here in two independent ways.

  (i)  SAME SEARCH, DIFFERENT SEED.  The search is re-run from a different
       first generator m1, which yields a different generating triple, different
       Gram filters and a different candidate list.  A complete search must
       return the same SET of matrices.

  (ii) ASSUMPTION-FREE BRUTE FORCE for Stab(L).  e_0 + e_i and e_0 - e_i are
       roots of L for i = 1..7; an automorphism of (O,.s) fixes e_0 and is
       orthogonal, so A(e_0 +- e_i) = e_0 +- A(e_i) must again be a root of L;
       and the only roots of L of the form e_0 + (imaginary unit) are e_0 +- e_j.
       Hence EVERY A in Stab(L) is a signed permutation of e_1..e_7 fixing e_0.
       All 2^7 * 7! = 645120 signed permutations are therefore tested directly.
       (The forcing argument uses roots of norm 2 and does NOT apply to Ls or
       Lsbar, whose minimal vectors have norm 4 -- and indeed those two
       stabilisers do contain non-monomial elements.)

Also exports 8x8 generators of the three groups for GAP identification
(expected: 1344 = AGL(3,2) = 2^3:L_3(2);  168 = L_3(2);  12096 = G_2(2)).
"""
import os, sys, time
from fractions import Fraction as F
from itertools import permutations, product as iproduct

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_consecutive_twists_exact import Lb, Lsb, Lsbarb, inL, inLs, inLsbar
from verify_aut_lambda_star import (
    octonion_autos_stabilising, MIN_L, MIN_LS, MIN_LSBAR, odot, E, _app, _mul8)

def key(A): return tuple(tuple(r) for r in A)
I8 = [[F(1) if r == c else F(0) for c in range(8)] for r in range(8)]

# the .s multiplication table on imaginary units:  e_a .s e_b = SGN[a][b] e_TGT[a][b]
SGN = [[0]*8 for _ in range(8)]; TGT = [[0]*8 for _ in range(8)]
for a in range(1,8):
    for b in range(1,8):
        p = odot(E[a], E[b])
        nz = [k for k in range(8) if p[k] != 0]
        assert len(nz) == 1
        TGT[a][b] = nz[0]; SGN[a][b] = int(p[nz[0]])

def brute_force_monomial_stab_L():
    """all signed permutations of e_1..e_7 fixing e_0 that are .s-automorphisms
    and preserve L: an exhaustive pass over all 645120 of them."""
    monos = []
    for perm in permutations(range(1,8)):
        pi = [0]*8
        for c in range(1,8): pi[c] = perm[c-1]
        # first the multiplication-table condition on the index map (sign-free part)
        if any(0 in (TGT[a][b],) or pi[TGT[a][b]] != TGT[pi[a]][pi[b]]
               for a in range(1,8) for b in range(1,8) if a != b):
            continue
        for signs in iproduct((1,-1), repeat=7):
            eps = [1] + list(signs)
            ok = True
            for a in range(1,8):
                for b in range(1,8):
                    if a == b: continue
                    c = TGT[a][b]
                    if SGN[a][b]*eps[c] != eps[a]*eps[b]*SGN[pi[a]][pi[b]]:
                        ok = False; break
                if not ok: break
            if not ok: continue
            A = [[F(0)]*8 for _ in range(8)]
            A[0][0] = F(1)
            for c in range(1,8): A[pi[c]][c] = F(eps[c])
            monos.append(A)
    return [A for A in monos if all(inL(_app(A, v)) for v in Lb)], len(monos)

def closure_gens(G):
    """a small generating set of the group G (list of matrices), by BFS closure."""
    target = {key(a) for a in G}
    gens = []
    cur = {key(I8): I8}
    for A in G:
        if key(A) in cur: continue
        gens.append(A)
        frontier = list(cur.values())
        while frontier:
            nxt = []
            for X in frontier:
                for g in gens:
                    P = _mul8(X, g)
                    if key(P) not in cur:
                        cur[key(P)] = P; nxt.append(P)
            frontier = nxt
        if len(cur) == len(target): break
    return gens, len(cur)

def main():
    t0 = time.time()
    print("="*72)
    print("cross-checks on the octonion-automorphism lattice stabilisers")
    print("="*72)
    groups = {}
    for name, minv, inG, basis in [("Stab(L)",     MIN_L,     inL,     Lb),
                                   ("Stab(Ls)",    MIN_LS,    inLs,    Lsb),
                                   ("Stab(Lsbar)", MIN_LSBAR, inLsbar, Lsbarb)]:
        G0 = octonion_autos_stabilising(minv, inG, basis, start=0)
        G1 = octonion_autos_stabilising(minv, inG, basis, start=37)
        groups[name] = G0
        print(f"(i)  {name:<12} |G| = {len(G0):>5}   re-run from a different"
              f" generating triple gives the same set: "
              f"{ {key(a) for a in G0} == {key(a) for a in G1} }")
        nonmono = sum(1 for A in G0
                      if any(sum(1 for r in range(8) if A[r][c] != 0) != 1
                             for c in range(8)))
        print(f"     non-monomial elements: {nonmono}")

    print()
    Gm, nmono = brute_force_monomial_stab_L()
    print("(ii) brute force over all 645120 signed permutations of e_1..e_7:")
    print("     monomial .s-automorphisms (before the L-condition):", nmono,
          "(= |2^3 : L_3(2)|)")
    print("     of these, preserving L:", len(Gm))
    print("     equal as a SET to the backtracking Stab(L):",
          {key(a) for a in Gm} == {key(a) for a in groups["Stab(L)"]})

    path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "..", "gap_project", "octonion_stabilisers_gens.g"))
    with open(path, "w") as f:
        f.write("# generators of the three octonion-automorphism lattice stabilisers\n")
        f.write("# (8x8 rational matrices; see verify_aut_octonion_crosscheck.py)\n")
        for name, gname in [("Stab(L)","StabL"), ("Stab(Ls)","StabLs"),
                            ("Stab(Lsbar)","StabLsbar")]:
            gg, sz = closure_gens(groups[name])
            print(f"     GAP export {gname}: {len(gg)} generators, closure {sz}"
                  f" = |G| : {sz == len(groups[name])}")
            f.write(f"{gname} := [\n")
            for A in gg:
                f.write("[" + ",".join("[" + ",".join(
                    (str(A[r][c].numerator) if A[r][c].denominator == 1
                     else f"{A[r][c].numerator}/{A[r][c].denominator}")
                    for c in range(8)) + "]" for r in range(8)) + "],\n")
            f.write("];\n")
    print("     written to gap_project/octonion_stabilisers_gens.g")
    print("\nruntime %.0f s" % (time.time()-t0))

if __name__ == "__main__":
    main()
