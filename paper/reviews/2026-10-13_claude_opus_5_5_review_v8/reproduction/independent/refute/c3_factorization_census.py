#!/usr/bin/env python3
"""
c3_factorization_census.py -- validates, on every permutation twist of the
seven imaginary units, the structural formula

     [Lambda : S_pi] = [B : B.B] * [A : A.A]^2 ,  A = pi(L sbar), B = pi(L s)

(standard product on the right-hand side), which explains the span indices
2^16, 2^10 and 1.  By-products: the closure census per cycle type (paper
Section 5.6), the index per parity (claim C2) and the subalgebra criterion
(claim C5) on all 5,039 permutations.

Run:   python3 c3_factorization_census.py > c3_factorization_census.out
Needs: leech_twist_lib.py.  Exact integer arithmetic.
Convention (paper Section 5.6): pi(e_a) = e_{pi(a)}, e0 fixed, and
x ._pi y := pi^{-1}(pi(x) . pi(y)).  Closure is decided DIRECTLY from the 576
basis-pair products (Wilson's three conditions), not from the criterion.
"""
import itertools
import sys
import time
from collections import Counter, defaultdict

import leech_twist_lib as T

t0 = time.time()


def say(*a):
    print(*a)
    sys.stdout.flush()


def cycle_type(p):
    seen, lens = set(), []
    for a in range(1, 8):
        if a in seen:
            continue
        n, b = 0, a
        while b not in seen:
            seen.add(b)
            b = p[b]
            n += 1
        if n > 1:
            lens.append(n)
    return tuple(sorted(lens, reverse=True))


def parity(p):
    return sum(l - 1 for l in cycle_type(p)) % 2   # 0 even, 1 odd


def square_index(rows):
    X = T.Lattice(rows)
    prods = [T.omul2(a, b) for a in rows for b in rows]
    if not all(q in X for q in prods):
        return None          # not closed under the product
    return T.sub_index(prods, X)


def subalgebra_mod2(rows):
    """X.X inside X (X contains 2L, so this is 'X/2L is a subalgebra of L/2L')."""
    X = T.Lattice(rows)
    return all(T.omul2(a, b) in X for a in rows for b in rows)


def closes_direct(tw):
    for bi in T.LAMBDA_BASIS2:
        for bj in T.LAMBDA_BASIS2:
            if not T.in_Lambda2(T.tprod2(tw, bi, bj)):
                return False
    return True


census = defaultdict(Counter)
by_parity = defaultdict(Counter)
crit = Counter()
fact_ok = fact_bad = 0
n = 0
for img in itertools.permutations(range(1, 8)):
    p = {a: img[a - 1] for a in range(1, 8)}
    if all(p[a] == a for a in p):
        continue
    n += 1
    tw = T.Twist(T.perm_matrix(p))
    ct = cycle_type(p)
    A_rows = [tw.apply(b) for b in T.LSBAR_LAT.rows]
    B_rows = [tw.apply(b) for b in T.LS_LAT.rows]
    critV, critW = subalgebra_mod2(A_rows), subalgebra_mod2(B_rows)
    cl = closes_direct(tw)
    crit[(cl, critV and critW)] += 1
    if not cl:
        census[ct]["no"] += 1
        continue
    o = T.closure_and_index(tw)
    idx = o["index"]
    census[ct]["yes"] += 1
    by_parity["even" if parity(p) == 0 else "odd"][idx] += 1
    iA, iB = square_index(A_rows), square_index(B_rows)
    if iB * iA * iA == idx:
        fact_ok += 1
    else:
        fact_bad += 1
        say("FACTORISATION MISMATCH", p, idx, iA, iB)
    if n % 500 == 0:
        say(f"... {n} permutations done ({time.time() - t0:.0f} s)")

say("\nclosure per cycle type (direct 576-product test):")
tot_yes = 0
for ct in sorted(census, key=lambda c: (sum(c), c)):
    y, no = census[ct]["yes"], census[ct]["no"]
    tot_yes += y
    say(f"   {str(ct):12s} {y:4d} / {y + no:4d}")
say(f"   total        {tot_yes:4d} / {n}")
say("\nspan index of the closing twists, by parity of pi:")
for par in ("even", "odd"):
    say(f"   {par:5s}: " + ", ".join(f"index {k} (=2^{k.bit_length() - 1}): {v}" for k, v in sorted(by_parity[par].items())))
say("\nfactorisation [Lambda:S] = [B:B.B][A:A.A]^2 holds for", fact_ok, "closing twists; fails for", fact_bad)
say("\ncriterion check (closes directly, both V and W subalgebras): counts", dict(crit))
say(f"\nruntime {time.time() - t0:.0f} s")
