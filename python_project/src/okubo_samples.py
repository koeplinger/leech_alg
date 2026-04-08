"""
okubo_samples.py — Worked examples of Petersson isotope constructions.

This file constructs several Okubo algebras from the standard octonion
algebra, each using a different order-3 automorphism τ.  It is meant to
be read and run as an instructive companion to okubo.py.

Run:
    cd python_project && PYTHONPATH=src python3 src/okubo_samples.py

How an order-3 automorphism τ is built
========================================

Every order-3 automorphism of the octonions has the following anatomy:

1. Pick a FANO LINE — three imaginary basis elements {e_a, e_b, e_c}
   satisfying the quaternion relations  e_a * e_b = e_c  (and cyclic
   rotations).  Together with e_0, they span a quaternion subalgebra H
   that τ will fix pointwise.

2. The remaining four imaginary elements must be split into two ROTATION
   PAIRS.  Each pair is rotated by 2π/3 in its plane.  The pairing is
   not arbitrary: each pair must come from a Fano triple that also
   contains a common fixed element (the "mediator").

3. Not every mediator works — only certain ones yield a valid algebra
   automorphism.  Each valid mediator gives rise to a valid τ, and
   hence to a valid Okubo algebra.

This file demonstrates the construction for three different Fano lines,
and shows one invalid pairing to illustrate the constraint.

References
----------
[MarraniCorradettiZucconi2025]
    A. Marrani, D. Corradetti, F. Zucconi, "Physics with non-unital
    algebras? An invitation to the Okubo algebra", J. Phys. A 58 (2025)
    075202.  eq. (1.5)–(1.6).
"""

import sys
import numpy as np
from octonions import STANDARD_ALGEBRA, STANDARD_FANO_TRIPLES
from okubo import OkuboAlgebra, from_octonions, _validate_order3, _validate_automorphism


# ---------------------------------------------------------------------------
# Verification helper
# ---------------------------------------------------------------------------

def verify_okubo(alg: OkuboAlgebra) -> bool:
    """Check the essential properties that uniquely characterise the Okubo
    algebra and print a summary.  Returns True if all checks pass."""
    C = alg._struct
    basis = alg.basis
    ok = True

    # 1. Composition norm:  n(x*y) = n(x) n(y)  on all 64 basis pairs
    comp_ok = True
    for ei in basis:
        for ej in basis:
            if abs((ei * ej).norm_sq() - ei.norm_sq() * ej.norm_sq()) > 1e-10:
                comp_ok = False
    status = "PASS" if comp_ok else "FAIL"
    print(f"    Composition norm  n(x*y) = n(x)n(y) ............ {status}")
    ok = ok and comp_ok

    # 2. Symmetric composition:  (x*y)*x = n(x)·y  on all 64 basis pairs
    sym_ok = True
    for ei in basis:
        nx = ei.norm_sq()
        for ej in basis:
            if not ((ei * ej) * ei).is_close(nx * ej, atol=1e-10):
                sym_ok = False
            if not (ei * (ej * ei)).is_close(nx * ej, atol=1e-10):
                sym_ok = False
    status = "PASS" if sym_ok else "FAIL"
    print(f"    Symmetric composition  (x*y)*x = n(x)y ........ {status}")
    ok = ok and sym_ok

    # 3. Non-unital:  no left identity exists
    A = np.zeros((64, 8))
    b = np.zeros(64)
    for j in range(8):
        for k in range(8):
            A[j * 8 + k, :] = C[:, j, k]
            b[j * 8 + k] = 1.0 if j == k else 0.0
    u, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    nonunital = np.linalg.norm(A @ u - b) > 0.1
    status = "PASS" if nonunital else "FAIL"
    print(f"    Non-unital (no identity element) ............... {status}")
    ok = ok and nonunital

    # 4. Division algebra:  L_{e_i} invertible for all basis elements
    div_ok = True
    for i in range(8):
        a = np.zeros(8); a[i] = 1.0
        L = np.einsum('i,ijk->kj', a, C)
        if abs(np.linalg.det(L)) < 1e-10:
            div_ok = False
    status = "PASS" if div_ok else "FAIL"
    print(f"    Division algebra (L_a invertible) .............. {status}")
    ok = ok and div_ok

    # 5. Derivation algebra dimension = 8  (Aut = SU(3))
    n_eq = 512
    Ad = np.zeros((n_eq, 64))
    for i in range(8):
        for j in range(8):
            for b_ in range(8):
                row = i * 64 + j * 8 + b_
                for k in range(8):
                    Ad[row, k * 8 + b_] += C[i, j, k]
                for a_ in range(8):
                    Ad[row, i * 8 + a_] -= C[a_, j, b_]
                    Ad[row, j * 8 + a_] -= C[i, a_, b_]
    dim_der = 64 - np.linalg.matrix_rank(Ad, tol=1e-8)
    der_ok = (dim_der == 8)
    status = "PASS" if der_ok else "FAIL"
    print(f"    dim Der = {dim_der}  (expect 8 for SU(3)) .............. {status}")
    ok = ok and der_ok

    return ok


def build_tau(pair1, pair2):
    """Build an 8×8 τ matrix that fixes all indices not in pair1 or pair2,
    and rotates each pair by 2π/3."""
    T = np.eye(8)
    c, s = -0.5, np.sqrt(3.0) / 2.0   # cos 2π/3, sin 2π/3
    for (a, b) in [pair1, pair2]:
        T[a, a] = c;  T[a, b] = -s
        T[b, a] = s;  T[b, b] = c
    return T


# ===================================================================
#  Example 1:  Fix the Fano line {e_1, e_3, e_7}
#              (this is the "standard" τ from the literature)
# ===================================================================

print("""
========================================================================
EXAMPLE 1:  Fix the quaternion subalgebra {e_0, e_1, e_3, e_7}
========================================================================

Fano line L7:  e_7 * e_1 = e_3,  e_1 * e_3 = e_7,  e_3 * e_7 = e_1.

The complementary elements are {e_2, e_4, e_5, e_6}.  The mediating
fixed element is e_3, which appears in two cross-triples:

    (2, 3, 5):  e_2 and e_5 are connected through fixed e_3
    (3, 4, 6):  e_4 and e_6 are connected through fixed e_3

So the rotation pairs are (e_2, e_5) and (e_4, e_6), both at angle 2π/3.
This is the automorphism τ from [MarraniCorradettiZucconi2025] eq. (1.5).
""")

tau1 = build_tau((2, 5), (4, 6))
alg1 = from_octonions(STANDARD_ALGEBRA, tau1)
print(f"  Algebra: {alg1}")
print(f"  Idempotent e_0:  e_0 * e_0 = {alg1.basis_element(0) * alg1.basis_element(0)}")
print()
all_pass = verify_okubo(alg1)
print(f"\n  Conclusion: {'This is an Okubo algebra.' if all_pass else 'NOT an Okubo algebra!'}")


# ===================================================================
#  Example 2:  Fix the Fano line {e_4, e_5, e_7}
# ===================================================================

print("""
========================================================================
EXAMPLE 2:  Fix the quaternion subalgebra {e_0, e_4, e_5, e_7}
========================================================================

Fano line L4:  e_4 * e_5 = e_7,  e_5 * e_7 = e_4,  e_7 * e_4 = e_5.

The complementary elements are {e_1, e_2, e_3, e_6}.  The mediating
fixed element is e_4, which appears in two cross-triples:

    (1, 2, 4):  e_1 and e_2 are connected through fixed e_4
    (3, 4, 6):  e_3 and e_6 are connected through fixed e_4

So the rotation pairs are (e_1, e_2) and (e_3, e_6), both at angle 2π/3.
Note: this line also admits a second valid mediator e_5 with pairs
(e_2, e_3) and (e_1, e_6), giving a different but equally valid τ.
""")

tau2 = build_tau((1, 2), (3, 6))
alg2 = from_octonions(STANDARD_ALGEBRA, tau2)
print(f"  Algebra: {alg2}")
print(f"  Idempotent e_0:  e_0 * e_0 = {alg2.basis_element(0) * alg2.basis_element(0)}")
print()
all_pass = verify_okubo(alg2)
print(f"\n  Conclusion: {'This is an Okubo algebra.' if all_pass else 'NOT an Okubo algebra!'}")


# ===================================================================
#  Example 3:  Fix the Fano line {e_1, e_2, e_4}
# ===================================================================

print("""
========================================================================
EXAMPLE 3:  Fix the quaternion subalgebra {e_0, e_1, e_2, e_4}
========================================================================

Fano line L1:  e_1 * e_2 = e_4,  e_2 * e_4 = e_1,  e_4 * e_1 = e_2.

The complementary elements are {e_3, e_5, e_6, e_7}.  The mediating
fixed element is e_4, which appears in two cross-triples:

    (3, 4, 6):  e_3 and e_6 are connected through fixed e_4
    (4, 5, 7):  e_5 and e_7 are connected through fixed e_4

So the rotation pairs are (e_3, e_6) and (e_5, e_7), both at angle 2π/3.
""")

tau3 = build_tau((3, 6), (5, 7))
alg3 = from_octonions(STANDARD_ALGEBRA, tau3)
print(f"  Algebra: {alg3}")
print(f"  Idempotent e_0:  e_0 * e_0 = {alg3.basis_element(0) * alg3.basis_element(0)}")
print()
all_pass = verify_okubo(alg3)
print(f"\n  Conclusion: {'This is an Okubo algebra.' if all_pass else 'NOT an Okubo algebra!'}")


# ===================================================================
#  Example 4:  INVALID pairing — same line, wrong mediator
# ===================================================================

print("""
========================================================================
EXAMPLE 4 (NEGATIVE):  Fix {e_0, e_1, e_3, e_7} with wrong pairing
========================================================================

Same Fano line L7 as Example 1, but now we try to mediate through e_1
instead of e_3.  This gives pairs (e_2, e_4) and (e_5, e_6):

    (1, 2, 4):  e_2 and e_4 are connected through fixed e_1
    (5, 6, 1):  e_5 and e_6 are connected through fixed e_1

Although each pair individually sits in a valid Fano triple with e_1,
the resulting map is NOT an algebra automorphism.  The cross products
between elements of different pairs are incompatible with this pairing.
""")

tau_bad = build_tau((2, 4), (5, 6))
_validate_order3(tau_bad)
print("  Order-3 check:  PASS  (τ³ = I)")
try:
    _validate_automorphism(STANDARD_ALGEBRA, tau_bad)
    print("  Automorphism check:  PASS")
    print("\n  (Unexpected — this should have failed!)")
except ValueError as e:
    print(f"  Automorphism check:  FAIL")
    print(f"    Reason: {e}")
    print("\n  Conclusion: Not every Fano-compatible pairing gives a valid τ.")
    print("  The mediating element must be chosen correctly.")
