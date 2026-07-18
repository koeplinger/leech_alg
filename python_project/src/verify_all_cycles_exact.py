"""Exact-arithmetic survey: closure of pi-twisted product on Lambda for
all pure n-cycles pi in S_7, n in {3,4,5,6,7}.

Re-uses the 24-vector Z-basis of Lambda and the closes() test from
verify_consecutive_twists_exact.py (each call: 24*24 basis-pair products
all required to lie in Lambda, exact arithmetic, no sampling).

Counts of n-cycles in S_7:
  3-cycles:  C(7,3) * 2!  =  35 *   2 =  70
  4-cycles:  C(7,4) * 3!  =  35 *   6 = 210
  5-cycles:  C(7,5) * 4!  =  21 *  24 = 504
  6-cycles:  C(7,6) * 5!  =   7 * 120 = 840
  7-cycles:  C(7,7) * 6!  =   1 * 720 = 720


RUNTIME NOTE (2026-07-17): about 15-20 minutes single-threaded on the
project's reference machine (2,344 permutations x 576 exact basis-pair
products).  For the full-S_7 census, including all mixed cycle types and
multiprocessing, see verify_all_permutations_exact.py.
"""
from itertools import combinations, permutations
from verify_consecutive_twists_exact import closes


def n_cycles(n):
    """All n-cycles in S_7 (acting on {1,...,7}, fixing 0)."""
    out = []
    for support in combinations(range(1, 8), n):
        # All cyclic arrangements: fix support[0], permute the rest
        for tail in permutations(support[1:]):
            arrangement = (support[0],) + tail
            p = list(range(8))
            for i in range(n):
                p[arrangement[i]] = arrangement[(i + 1) % n]
            out.append(tuple(p))
    return out


if __name__ == "__main__":
    print("Exact n-cycle closure survey (each permutation tested on the")
    print("24x24 = 576 basis-pair products of a Z-basis of Lambda).")
    print()
    print(f"{'cycle type':<12}{'tested':>10}{'close':>10}{'pct':>10}")
    print("-" * 42)
    grand_tested = 0
    grand_closed = 0
    for n in (3, 4, 5, 6, 7):
        cycles = n_cycles(n)
        n_close = sum(closes(p) for p in cycles)
        grand_tested += len(cycles)
        grand_closed += n_close
        pct = 100.0 * n_close / len(cycles) if cycles else 0.0
        print(f"{n}-cycles    {len(cycles):>10}{n_close:>10}{pct:>9.2f}%")
    print("-" * 42)
    pct = 100.0 * grand_closed / grand_tested
    print(f"{'total':<12}{grand_tested:>10}{grand_closed:>10}{pct:>9.2f}%")
