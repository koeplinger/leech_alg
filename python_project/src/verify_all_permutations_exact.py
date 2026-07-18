"""Closure census over ALL of S_7: for every permutation pi of the seven
imaginary octonion basis units, does the pi-twisted triple product close on
the Leech lattice?

Extends verify_all_cycles_exact.py (single cycles) and
verify_consecutive_twists_exact.py (transpositions, 3-cycles, (2,2)-doubles)
to the full symmetric group: all 5,040 permutations, grouped by cycle type
(the 15 partitions of 7, with fixed points suppressed in the labels).

For each pi the test is exact and complete: closure of the pi-twisted
Z_3-symmetric triple product is checked on all 24 x 24 = 576 basis-pair
products of a Z-basis of Lambda, in exact rational arithmetic (by
bilinearity this settles closure on all of Lambda).  The per-permutation
test is closes() from verify_consecutive_twists_exact.py, reused unchanged.

Runtime: roughly an hour single-threaded; parallelized over CPU cores with
multiprocessing (chunked), typically 10-20 minutes.
"""
import os
import sys
from collections import defaultdict
from itertools import permutations
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from verify_consecutive_twists_exact import closes


def cycle_type(perm):
    """Cycle type of a permutation of {1..7} given as a tuple of images
    (perm[i-1] = image of i), as a sorted tuple of cycle lengths >= 2."""
    seen = [False] * 8
    lengths = []
    for start in range(1, 8):
        if seen[start]:
            continue
        n = 0
        j = start
        while not seen[j]:
            seen[j] = True
            j = perm[j - 1]
            n += 1
        if n >= 2:
            lengths.append(n)
    return tuple(sorted(lengths, reverse=True))


def label(ct):
    if not ct:
        return "identity"
    if len(ct) == 1:
        return f"{ct[0]}-cycle"
    return "(" + ",".join(str(n) for n in ct) + ")"


def _test(perm):
    # closes() expects a length-8 map over coordinates 0..7 fixing 0
    return perm, closes((0,) + perm)


def main():
    print("=" * 70)
    print("Closure census over ALL of S_7 (5,040 permutations, exact)")
    print("=" * 70)

    perms = list(permutations(range(1, 8)))
    assert len(perms) == 5040

    tested = defaultdict(int)
    closed = defaultdict(int)

    with Pool() as pool:
        for i, (perm, ok) in enumerate(pool.imap_unordered(_test, perms, chunksize=32)):
            ct = cycle_type(perm)
            tested[ct] += 1
            if ok:
                closed[ct] += 1
            if (i + 1) % 500 == 0:
                print(f"  ... {i + 1}/5040 tested", flush=True)

    print()
    print(f"{'cycle type':<12} {'tested':>7} {'close':>7} {'rate':>8}")
    print("-" * 38)
    # sort: identity, then by (total cycle weight, type)
    order = sorted(tested, key=lambda ct: (len(ct) == 0, sum(ct), ct))
    total_t = total_c = 0
    for ct in order:
        t, c = tested[ct], closed[ct]
        if ct:
            total_t += t
            total_c += c
        print(f"{label(ct):<12} {t:>7} {c:>7} {100*c/t:>7.2f}%")
    print("-" * 38)
    print(f"{'non-identity':<12} {total_t:>7} {total_c:>7} {100*total_c/total_t:>7.2f}%")

    # sanity: totals per type match the class sizes of S_7
    expected = {(): 1, (2,): 21, (3,): 70, (2, 2): 105, (4,): 210,
                (3, 2): 420, (2, 2, 2): 105, (5,): 504, (4, 2): 630,
                (3, 3): 280, (6,): 840, (5, 2): 504, (4, 3): 420,
                (3, 2, 2): 210, (7,): 720}
    ok_sizes = all(tested[ct] == n for ct, n in expected.items())
    print(f"\nclass sizes match S_7: {ok_sizes}")
    # sanity: the previously published single-cycle and (2,2) rows
    prev = {(2,): 21, (3,): 35, (4,): 0, (5,): 252, (6,): 210, (7,): 336,
            (2, 2): 42, (): 0}
    ok_prev = all(closed[ct] == n for ct, n in prev.items())
    print(f"reproduces the published single-cycle and (2,2) counts: {ok_prev}")
    print()
    print("ALL PASS" if ok_sizes and ok_prev else "MISMATCH -- INVESTIGATE")
    return 0 if ok_sizes and ok_prev else 1


if __name__ == "__main__":
    raise SystemExit(main())
