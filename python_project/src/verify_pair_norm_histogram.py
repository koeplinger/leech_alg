"""Sampled histogram of pair-product norms N(u star v) on Min(Lambda) x
Min(Lambda), cited in the last paragraph of paper Appendix D
(`app:samples`); Remark 5.1 (`rem:product-norms`) refers to that paragraph.

Created: 2026-07-17.  Updated: 2026-10-24.

Reuses the doubled-integer star and the family-weighted minimal-vector
sampler of verify_section5_properties.py, unchanged; exact integer arithmetic
throughout (doubled coordinates, so norms appear x4).

Reported: the set of observed values of N(u star v), the mode and its share,
and whether the values 0 and 8 occur.  The paper's claim: on 10^6 sampled
pairs the norms took only the values {16, 32, 48, 64, 80, 96, 112, 128}, all
multiples of 16, with mode 64.

Invocation for the paper's figure (the two seeds of the identity table of
Appendix D, verify_section5_properties.py):
    python3 verify_pair_norm_histogram.py 20260601 500000
    python3 verify_pair_norm_histogram.py 20260602 500000
and the counts are summed.  Defaults: seed 20260601, 500,000 samples.
"""
import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from verify_section5_properties import Nsq, rand_min, star


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20260601
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 500_000
    random.seed(seed)

    hist = Counter()
    for _ in range(n):
        u, v = rand_min(), rand_min()
        hist[Nsq(star(u, v)) // 4] += 1     # Nsq is 4x the actual norm

    print(f"seed {seed}, {n} sampled pairs of Min(Lambda) x Min(Lambda)")
    print(f"{'N(u*v)':>7} {'count':>9} {'share':>8}")
    for val in sorted(hist):
        print(f"{val:>7} {hist[val]:>9} {100*hist[val]/n:>7.2f}%")
    values = sorted(hist)
    mode = max(hist, key=hist.get)
    ok = (all(v % 16 == 0 for v in values)
          and 0 not in hist and 8 not in hist and mode == 64)
    print(f"\nobserved values: {values}")
    print(f"all multiples of 16: {all(v % 16 == 0 for v in values)}")
    print(f"0 or 8 observed: {0 in hist or 8 in hist}")
    print(f"mode: {mode} ({100*hist[mode]/n:.2f}%)")
    print("ALL PASS" if ok else "UNEXPECTED -- INVESTIGATE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
