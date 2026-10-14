#!/usr/bin/env python3
"""
c3_other_claims.py -- independent checks of claims C4 and C5 (context for
the C3 review).

C4: signed permutations of the imaginary units (e0 fixed,
    e_a -> eps_a e_{pi(a)}) give 480 distinct twisted multiplication tables;
    168 of them close on Lambda, with span indices 2^10 (112) and 2^16 (56).
C5: L/2L carries 270 Lagrangians (maximal totally singular subspaces of
    q(v) = N(v)/2 mod 2); 142 are subalgebras and 8 two-sided ideals of
    (L/2L, mu-bar); 4,704 of the 17,280 ordered complementary pairs consist
    of two subalgebras.

Run:   python3 c3_other_claims.py > c3_other_claims.out
Needs: leech_twist_lib.py.  Exact integer arithmetic.
Conventions: main.tex v8 (see leech_twist_lib.py); twisted product
x ._g y = g^{-1}(g x . g y) with g the signed permutation matrix.
"""
import itertools
import sys
import time
from collections import Counter

import leech_twist_lib as T

t0 = time.time()


def say(*a):
    print(*a)
    sys.stdout.flush()


# ------------------------------------------------------------------ C4
say("=" * 78)
say("C4: signed permutations of e1..e7 as twists")
say("=" * 78)
MT = T.MT
tables = {}
for img in itertools.permutations(range(1, 8)):
    pi = [0] + list(img)                         # pi[a] for a = 0..7, pi[0] = 0
    inv = [0] * 8
    for a in range(8):
        inv[pi[a]] = a
    for signs in range(128):
        eps = [1] + [(-1 if (signs >> (a - 1)) & 1 else 1) for a in range(1, 8)]
        # table entry for (a,b): g^{-1}(g e_a . g e_b) = sign * e_c
        key = []
        for a in range(8):
            for b in range(8):
                sg, c = MT[pi[a]][pi[b]]
                c0 = inv[c]                       # g^{-1} e_c = eps[c0] e_{c0}
                key.append(eps[a] * eps[b] * sg * eps[c0] * (c0 + 1))
        key = tuple(key)
        even = sum(1 for a in range(1, 8) if eps[a] < 0) % 2 == 0     # g(L) = L iff even
        if key not in tables:
            tables[key] = [{a: pi[a] for a in range(1, 8)}, {a: eps[a] for a in range(1, 8)}, even]
        elif even and not tables[key][2]:
            # prefer a representative that preserves L
            tables[key] = [{a: pi[a] for a in range(1, 8)}, {a: eps[a] for a in range(1, 8)}, True]
say("signed permutations:", 5040 * 128, "| distinct twisted tables:", len(tables),
    f"({time.time() - t0:.0f} s)")
# NOTE: a signed permutation with an odd number of sign changes maps L = D8+
# onto the other coset lattice L' (it flips the parity of half-integer
# coordinate sums).  The twisted product still maps L x L into L (L' is also
# closed), so such twists are admissible; they are counted separately.
c4 = Counter()
c4_split = Counter()
for key, (pi, eps, even) in tables.items():
    tw = T.Twist(T.signed_perm_matrix(pi, eps), require_L=False)
    assert tw.maps_L_to_L == even
    o = T.closure_and_index(tw)
    outcome = "no closure" if not o["closes"] else f"index {o['index']}"
    c4[outcome] += 1
    c4_split[("some representative has g(L)=L" if even else "every representative has g(L)=L'", outcome)] += 1
say("closure / index over the distinct tables:", dict(c4), f"({time.time() - t0:.0f} s)")
say("split by whether the table class contains a twist g with g(L) = L:")
for k, v in sorted(c4_split.items()):
    say(f"   {k[0]:34s} {k[1]:14s} {v}")

# ------------------------------------------------------------------ C5
say("\n" + "=" * 78)
say("C5: Lagrangians of L/2L and their algebra properties")
say("=" * 78)
# vectors of L/2L as 8-bit masks of L-basis coordinates mod 2
G = [[int(x) for x in row] for row in T.L_LAT.gram()]      # integral, even diagonal


def q(m):
    c = [(m >> k) & 1 for k in range(8)]
    n = sum(c[i] * c[j] * G[i][j] for i in range(8) for j in range(8))
    return (n // 2) % 2


def b(m1, m2):
    c1 = [(m1 >> k) & 1 for k in range(8)]
    c2 = [(m2 >> k) & 1 for k in range(8)]
    return sum(c1[i] * c2[j] * G[i][j] for i in range(8) for j in range(8)) % 2


# mod-2 product from the products of L-basis vectors
prod_tab = [[0] * 8 for _ in range(8)]
for i in range(8):
    for j in range(8):
        c = T.L_LAT.coords(T.omul2(T.L_BASIS2[i], T.L_BASIS2[j]))
        prod_tab[i][j] = sum(((c[k] % 2) << k) for k in range(8))


def mu(m1, m2):
    out = 0
    for i in range(8):
        if (m1 >> i) & 1:
            for j in range(8):
                if (m2 >> j) & 1:
                    out ^= prod_tab[i][j]
    return out


singular = [m for m in range(1, 256) if q(m) == 0]
say("nonzero singular vectors:", len(singular), "(O+(8,2): 135)")


def span(vs):
    S = {0}
    for v in vs:
        S |= {s ^ v for s in S}
    return frozenset(S)


lag = set()
frontier = {span([v]) for v in singular}
for dim in range(2, 5):
    nxt = set()
    for Sp in frontier:
        for w in singular:
            if w in Sp:
                continue
            if all(b(w, x) == 0 for x in Sp):
                nxt.add(span(list(Sp) + [w]))
    frontier = nxt
lag = frontier
say("Lagrangians (4-dim totally singular subspaces):", len(lag))
full = range(256)
def is_sub(V):
    return all(mu(x, y) in V for x in V for y in V)
def is_ideal(V):
    return all(mu(l, x) in V and mu(x, l) in V for l in range(256) for x in V)
subs = [V for V in lag if is_sub(V)]
ideals = [V for V in lag if is_ideal(V)]
say("subalgebras:", len(subs), "| two-sided ideals:", len(ideals))
pairs = [(V, W) for V in lag for W in lag if V & W == frozenset({0})]
both = sum(1 for V, W in pairs if V in subs and W in subs) if False else \
    sum(1 for V, W in pairs if is_sub(V) and is_sub(W))
say("ordered complementary pairs:", len(pairs), "| with both members subalgebras:", both)
# the pairs realised by the paper's sigma and by the C3 witness
def lag_of(rows):
    return span([sum(((c % 2) << k) for k, c in enumerate(T.L_LAT.coords(r))) for r in rows])
for name, tw in (("identity", T.Twist([[int(i == j) for j in range(8)] for i in range(8)])),
                 ("sigma", T.Twist(T.perm_matrix({1: 2, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}))),
                 ("s_r witness", T.Twist(T.reflection([1, 1, 1, 1, 1, 1, 1, -1])))):
    V = lag_of([tw.apply(x) for x in T.LSBAR_LAT.rows])
    W = lag_of([tw.apply(x) for x in T.LS_LAT.rows])
    say(f"{name:12s} V_g Lagrangian: {V in lag}, W_g Lagrangian: {W in lag}, complementary: {V & W == frozenset({0})}, "
        f"V sub: {is_sub(V)}, V ideal: {is_ideal(V)}, W sub: {is_sub(W)}, W ideal: {is_ideal(W)}")
say(f"\nruntime {time.time() - t0:.0f} s")
