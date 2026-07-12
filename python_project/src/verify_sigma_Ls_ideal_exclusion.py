"""Exclude the "left ideal" explanation of Lemma 4.4.

Section 7 of the paper asks for a structural reason why sigma(Ls) is closed
under the standard octonion product (Lemma 4.4), where Ls is not.  One
candidate shape for such a reason is that sigma(Ls) might be an ideal of the
octonion order L.  This script refutes that candidate at the integer level.

Counted here, for each pair of Z-bases (8 x 8 = 64 basis-pair products), is
how many products land back in the target lattice:

  L . sigma(Ls)          <= sigma(Ls)       -- left ideal?
  sigma(Ls) . L          <= sigma(Ls)       -- right ideal?
  sigma(Ls) . sigma(Ls)  <= sigma(Ls)       -- Lemma 4.4 (subalgebra)
  L . sigma(Ls_bar)      <= sigma(Ls_bar)   -- Lemma 4.3 (left ideal)
  sigma(Ls_bar) . L      <= sigma(Ls_bar)   -- two-sided?
  L . Ls                 <= Ls              -- control (Ls is not closed)

Membership is decided in exact rational arithmetic: a product lies in the
target lattice iff its coordinates in the target Z-basis are integers.

Expected (and obtained): 32/64, 21/64, 64/64, 64/64, 64/64, 24/64.  So
sigma(Ls) is neither a left nor a right ideal of L; its companion
sigma(Ls_bar) is a two-sided ideal.  That asymmetry -- ideal on one
sublattice, mere subalgebra on the other -- is exactly what the mod-2
reduction records (V an ideal, W only a subalgebra; see
verify_mod2_quotient.py), and it is what makes Lemma 4.4 the hard one.

The mod-2 picture already implies the exclusion (an integer-level ideal
would reduce to a mod-2 ideal), but the integer-level counts are recorded
here directly.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from verify_mod2_quotient import L, Ls, sLs, sLsbar, coeffs, omul


def count_products_inside(left, right, target):
    """Count how many of the |left| x |right| products lie in Z-span(target)."""
    inside = 0
    total = 0
    for a in left:
        for b in right:
            total += 1
            c = coeffs(omul(a, b), target)
            if c is not None and all(x.denominator == 1 for x in c):
                inside += 1
    return inside, total


TESTS = [
    ("L . sigma(Ls)         <= sigma(Ls)      [left ideal?]", L, sLs, sLs, 32),
    ("sigma(Ls) . L         <= sigma(Ls)      [right ideal?]", sLs, L, sLs, 21),
    ("sigma(Ls) . sigma(Ls) <= sigma(Ls)      [Lemma 4.4]", sLs, sLs, sLs, 64),
    ("L . sigma(Lsbar)      <= sigma(Lsbar)   [Lemma 4.3]", L, sLsbar, sLsbar, 64),
    ("sigma(Lsbar) . L      <= sigma(Lsbar)   [two-sided?]", sLsbar, L, sLsbar, 64),
    ("L . Ls                <= Ls             [control]", L, Ls, Ls, 24),
]


if __name__ == "__main__":
    print("=" * 72)
    print("IDEAL EXCLUSION FOR sigma(Ls)  (Section 7, open question 1)")
    print("=" * 72)
    ok = True
    for label, left, right, target, expected in TESTS:
        inside, total = count_products_inside(left, right, target)
        match = inside == expected
        ok &= match
        print(f"  {label}: {inside:2d}/{total}"
              f"   (expected {expected:2d}) {'OK' if match else 'MISMATCH'}")
    print()
    print("  sigma(Ls) is NOT an ideal of L (neither side); sigma(Lsbar) is a")
    print("  two-sided ideal.  The 'left ideal' candidate for a structural")
    print("  reason for Lemma 4.4 is therefore excluded.")
    print()
    print("ALL COUNTS AS EXPECTED" if ok else "UNEXPECTED COUNTS -- INVESTIGATE")
    sys.exit(0 if ok else 1)
