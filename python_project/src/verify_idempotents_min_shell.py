"""Exhaustive idempotent / square-nilpotent search over the full minimal
shell Min(Lambda) of the Leech lattice under the paper's product star.

For every one of the 196,560 minimal vectors u of Lambda (norm 8) in
Wilson's octonionic representation, compute u star u under the paper's
sigma-twisted triple-octonion product and check
  (a) u star u = u   (idempotent),
  (b) u star u = 0   (square-nilpotent),
and tally the histogram of norms N(u star u) (actual norm convention:
minimal vectors have N = 8, so norm multiplicativity gives N(u*u) = 64).

All arithmetic is exact, in doubled-integer coordinates (Nsq = 4*N;
minimal vectors have Nsq = 32).  The three Wilson families are
enumerated exhaustively, not sampled:
  Type 1    (720):  (2 lambda, 0, 0)                        + cyclic shifts,
  Type 2 (11,520):  (lambda sbar, (lambda sbar) j, 0)       + cyclic shifts,
  Type 3 (184,320): ((lambda s) j, lambda k, (lambda j) k)  + cyclic shifts,
with lambda over the 240 E8 roots and j, k over the 16 unit octonions
(same construction formulas as rand_min in verify_section5_properties.py).
Before the search runs, the enumeration is deduplicated and asserted to
give exactly 196,560 distinct vectors, each of norm Nsq = 32.

Reuses the exact star machinery from verify_section5_properties.py.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_section5_properties import (star, eq, Nsq, omul2,
                                        droots, dunits, ds, dsbar, Z, cyc)

ZERO = (Z, Z, Z)


def enumerate_min_shell():
    """All 196,560 minimal vectors of Lambda, doubled-int coords."""
    seen = set()
    vecs = []
    counts = [0, 0, 0]                       # per-type tallies
    for lam in droots:
        # Type 1: (2 lambda, 0, 0) and cyclic shifts
        blk = [2 * c for c in lam]
        for sh in range(3):
            u = cyc([blk, Z, Z], sh)
            key = tuple(u[0]) + tuple(u[1]) + tuple(u[2])
            if key not in seen:
                seen.add(key); vecs.append(u); counts[0] += 1
        # Type 2: (lambda sbar, (lambda sbar) j, 0) and cyclic shifts
        lsb = omul2(lam, dsbar)
        for j in dunits:
            lsbj = omul2(lsb, j)
            for sh in range(3):
                u = cyc([lsb, lsbj, Z], sh)
                key = tuple(u[0]) + tuple(u[1]) + tuple(u[2])
                if key not in seen:
                    seen.add(key); vecs.append(u); counts[1] += 1
        # Type 3: ((lambda s) j, lambda k, (lambda j) k) and cyclic shifts
        ls = omul2(lam, ds)
        for j in dunits:
            lsj = omul2(ls, j)
            lj = omul2(lam, j)
            for k in dunits:
                lk = omul2(lam, k)
                ljk = omul2(lj, k)
                for sh in range(3):
                    u = cyc([lsj, lk, ljk], sh)
                    key = tuple(u[0]) + tuple(u[1]) + tuple(u[2])
                    if key not in seen:
                        seen.add(key); vecs.append(u); counts[2] += 1
    return vecs, counts


if __name__ == "__main__":
    t0 = time.time()
    vecs, counts = enumerate_min_shell()
    assert counts == [720, 11520, 184320], f"per-type counts off: {counts}"
    assert len(vecs) == 196560, f"expected 196560 distinct, got {len(vecs)}"
    assert all(Nsq(u) == 32 for u in vecs), "a listed vector is not minimal"
    print(f"enumerated {len(vecs)} distinct minimal vectors "
          f"(types: {counts[0]} + {counts[1]} + {counts[2]}), "
          f"all of norm 8   [{time.time()-t0:.0f} s]")

    idem = []                                # u with u*u = u
    nilp = []                                # u with u*u = 0
    hist = {}                                # N(u*u) -> count (actual norm)
    for u in vecs:
        s2 = star(u, u)
        if eq(s2, u):
            idem.append(u)
        if eq(s2, ZERO):
            nilp.append(u)
        n = Nsq(s2) // 4
        hist[n] = hist.get(n, 0) + 1
    el = time.time() - t0

    print(f"idempotents (u*u = u):        {len(idem)}")
    print(f"square-nilpotents (u*u = 0):  {len(nilp)}")
    print(f"min N(u*u): {min(hist)}    max N(u*u): {max(hist)}")
    print("histogram of N(u*u):")
    for n in sorted(hist):
        print(f"  {n:>4}: {hist[n]}")
    assert sum(hist.values()) == 196560
    print(f"runtime: {el:.0f} s")
