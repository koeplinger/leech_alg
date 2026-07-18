"""Reproduce Egan's lower bound 244,035,421 on the number of copies of the
Leech lattice in O^3 = (E_8)^3.

Source: Greg Egan, in J. C. Baez, "Integral octonions (Part 10)",
The n-Category Cafe, 1 December 2014,
https://golem.ph.utexas.edu/category/2014/12/integral_octonions_part_10.html
(Correction, 18 July 2026: the part-10 date fixed from 12 December to
1 December 2014; 12 December is part 11's date.)

This is Egan's DERIVATION, an orbit-stabiliser lower bound; it is NOT a
re-run of his enumeration (which counts the embeddings one by one, and which
this project does not reproduce). It reproduces the arithmetic that turns
two standard group orders into the figure quoted in Section 6 of the paper.

The argument.
  The isometry group of the lattice O^3 = (E_8)^3 acts on the set of Leech
  sublattices of the type Egan considers. By orbit-stabiliser, the orbit of
  one such lattice Lambda has size
        |Isom(O^3)| / |Stab(Lambda)|.
  The stabiliser of Lambda inside Isom(O^3) is a subgroup of the full
  automorphism group of the Leech lattice, Aut(Lambda) = Co_0, so
        |Stab(Lambda)| <= |Co_0|.
  Hence the orbit -- and therefore the total number of such Leech lattices --
  is at least
        |Isom(O^3)| / |Co_0|,
  and, being an integer, at least the ceiling of that ratio.

The inputs (both standard).
  |Aut(E_8)| = |W(E_8)| = 696,729,600, the order of the E_8 Weyl group. It is
  reproduced here as the product of the degrees of the basic invariants of
  W(E_8): 2, 8, 12, 14, 18, 20, 24, 30 (a standard fact: for a finite Coxeter
  group, |W| equals the product of its invariant degrees).
  |Isom(O^3)| = 3! * |Aut(E_8)|^3, the automorphisms of (E_8)^3 as a lattice:
  an automorphism of each of the three E_8 factors, composed with a permutation
  of the three factors (the wreath product Aut(E_8) wr S_3).
  |Co_0| = 2^22 * 3^9 * 5^4 * 7^2 * 11 * 13 * 23 = 8,315,553,613,086,720,000,
  the order of Conway's group Co_0 = Aut(Leech) (Conway 1969; ATLAS).
"""
from math import prod


# |W(E_8)| as the product of the degrees of the basic invariants.
E8_INVARIANT_DEGREES = (2, 8, 12, 14, 18, 20, 24, 30)
W_E8 = prod(E8_INVARIANT_DEGREES)

# |Co_0| from its prime factorisation.
CO0_FACTORISATION = {2: 22, 3: 9, 5: 4, 7: 2, 11: 1, 13: 1, 23: 1}
CO0 = prod(p**e for p, e in CO0_FACTORISATION.items())

EGAN_LOWER_BOUND = 244_035_421


def ceil_div(a, b):
    """Exact integer ceiling of a / b for positive integers."""
    return -(-a // b)


def main():
    print("=" * 68)
    print("Egan's lower bound on Leech lattices in O^3 = (E_8)^3")
    print("=" * 68)

    print(f"\n|W(E_8)| = product of invariant degrees "
          f"{E8_INVARIANT_DEGREES}")
    print(f"        = {W_E8:,}")
    assert W_E8 == 696_729_600, "W(E_8) order mismatch"

    print(f"\n|Co_0|   = 2^22 * 3^9 * 5^4 * 7^2 * 11 * 13 * 23")
    print(f"        = {CO0:,}")
    assert CO0 == 8_315_553_613_086_720_000, "Co_0 order mismatch"

    isom_O3 = 6 * W_E8**3          # 3! * |Aut(E_8)|^3 = Aut(E_8) wr S_3
    print(f"\n|Isom(O^3)| = 3! * |W(E_8)|^3")
    print(f"           = {isom_O3:,}")

    ratio_num, ratio_den = isom_O3, CO0
    floor = ratio_num // ratio_den
    bound = ceil_div(ratio_num, ratio_den)
    print(f"\n|Isom(O^3)| / |Co_0| = {ratio_num / ratio_den:.6f}")
    print(f"   floor = {floor:,}")
    print(f"   ceil  = {bound:,}   <-- lower bound (orbit is an integer)")

    ok = bound == EGAN_LOWER_BOUND
    print(f"\nreproduces Egan's 244,035,421 : {ok}")
    print("=" * 68)
    print("ALL PASS" if ok else "MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
