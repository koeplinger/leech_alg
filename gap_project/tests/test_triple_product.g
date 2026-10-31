#############################################################################
##
##  test_triple_product.g — The Z_3-routed triple product on R^24 (paper
##  Definition 3.3, `def:triple`) and the main theorem (Theorem 1.1,
##  `thm:main`: closure of Lambda under *_sigma).
##
##  Created: 2026-04-28.  Updated: 2026-10-28.
##
##  This is the paper's main result, verified here on a representative set of
##  Leech minimal vectors:
##    - all 720 type-1 vectors  (one block of 2*lambda, two zero blocks);
##    - a sample of type-2 vectors  (two nonzero blocks);
##    - a sample of type-3 vectors  (three nonzero blocks).
##
##  Theorem 1.1 says  Lambda *_sigma Lambda  subset Lambda.  By bilinearity,
##  closure on a generating set of Lambda over Z is sufficient.  The 196,560
##  minimal vectors generate Lambda; this file spot-checks the products of
##  all 720 type-1 vectors with six right operands, and type-2 / type-3
##  sample grids.  (The 12M+ random-pair tests are in
##  python_project/src/trial_007_{scaled_test,fast,exhaust}.py; here each
##  product is checked exactly against Wilson's three conditions.)
##
##  Also checks two instances of paper Remark 5.2 (Section 5.1,
##  `rem:no-identity`), that (R^24, *_sigma) has no multiplicative identity:
##    - (e_0, e_0, e_0) is not one;
##    - no (a, b, c) with components in J is a left identity.
##
#############################################################################

Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
LoadAllSrc(GAP_PROJECT_ROOT);
ResetTestCounters();

Print("\n=== test_triple_product.g  (paper Theorem 1.1) ===\n");

J := J_UNITS();
ROOTS := WilsonE8Roots();
ZEROVEC  := [0,0,0,0,0,0,0,0];

InLeechTriple := function(uv)
    return IsInLeech(uv[1], uv[2], uv[3]);
end;

# ---------------------------------------------------------------------------
# Theorem 1.1 spot-check: u *_sigma v  in Lambda  for sample u, v in Lambda.
# We use the 720 type-1 vectors against six right operands, then sample
# grids of type-2 x {type-2, type-3} and type-3 x type-3 to exercise all
# three Wilson conditions nontrivially.
# ---------------------------------------------------------------------------

# Type-1 closure: full check over the 720 type-1 minimal vectors as left
# operand, paired with a representative right operand from each family.
T1  := LeechType1Vectors();

# A few sample type-2 / type-3 left and right operands, built from generators.
type2_sample := [];
for ridx in [1, 60, 120, 180, 240] do
    for jidx in [1, 5, 9, 13] do
        lam  := ROOTS[ridx];
        j    := J[jidx];
        a    := OctMult(lam, WILSON_S_BAR);
        b    := OctMult(a, j);
        Add(type2_sample, [a, b, ZEROVEC]);
        Add(type2_sample, [ZEROVEC, a, b]);
        Add(type2_sample, [b, ZEROVEC, a]);
    od;
od;

type3_sample := [];
for ridx in [1, 60, 120, 180, 240] do
    for jidx in [1, 5] do
        for kidx in [1, 9] do
            lam := ROOTS[ridx];
            j   := J[jidx];
            k   := J[kidx];
            a   := OctMult(OctMult(lam, WILSON_S), j);
            b   := OctMult(lam, k);
            c   := OctMult(OctMult(lam, j), k);
            Add(type3_sample, [a, b, c]);
        od;
    od;
od;

# Theorem 1.1:  Lambda *_sigma Lambda  subset Lambda  (sample-checked).
nfail := 0;
ntot  := 0;

# Type-1 x {type-1, type-2, type-3}: exhaustive over the 720 left operands times
# 6 right operands (3 type-1, 2 type-2, 1 type-3, the type-3 vector to exercise
# condition 3 hardest).
right_ops := [
    T1[1], T1[241], T1[481],                       # three type-1 vectors, each (2r, 0, 0)
    type2_sample[1], type2_sample[10],
    type3_sample[1]
];
for u in T1 do
    for v in right_ops do
        if not InLeechTriple(TripleProductSigma(u, v)) then nfail := nfail + 1; fi;
        ntot := ntot + 1;
    od;
od;
CHECK(Concatenation("Theorem 1.1: type-1 x { 3 type-1, 2 type-2, 1 type-3 } closure (",
                    String(ntot), " products)"),
      nfail = 0);

# Type-2 x type-2 and type-2 x type-3 sample.
nfail := 0;
ntot  := 0;
for u in type2_sample do
    for v in type2_sample do
        if not InLeechTriple(TripleProductSigma(u, v)) then nfail := nfail + 1; fi;
        ntot := ntot + 1;
    od;
    for v in type3_sample do
        if not InLeechTriple(TripleProductSigma(u, v)) then nfail := nfail + 1; fi;
        ntot := ntot + 1;
    od;
od;
CHECK(Concatenation("Theorem 1.1: type-2 x { type-2 u type-3 } closure (",
                    String(ntot), " products)"),
      nfail = 0);

# Type-3 x type-3 sample.
nfail := 0;
ntot  := 0;
for u in type3_sample do
    for v in type3_sample do
        if not InLeechTriple(TripleProductSigma(u, v)) then nfail := nfail + 1; fi;
        ntot := ntot + 1;
    od;
od;
CHECK(Concatenation("Theorem 1.1: type-3 x type-3 closure (", String(ntot), " products)"),
      nfail = 0);

# ---------------------------------------------------------------------------
# The standard (untwisted) triple product fails on Lambda (companion
# Example 8.1).
#
# With the standard product in all nine blocks, P_1+P_2+P_3 =
# (x+y+z)*(x'+y'+z'), so Wilson's condition (3) for the product asks that
# the product of two block sums, both in Ls, lie in Ls; but Ls*Ls is not
# contained in Ls (companion Example 4.2; paper Remark 4.5).
# ---------------------------------------------------------------------------

# The check searches the type-3 sample for a pair (u, v) whose untwisted
# product is not in Lambda.
found_untwisted_fail := false;
for u in type3_sample do
    for v in type3_sample do
        if not InLeechTriple(TripleProduct(u, v)) then
            found_untwisted_fail := true;
            break;
        fi;
    od;
    if found_untwisted_fail then break; fi;
od;
CHECK("Companion Example 8.1: standard (untwisted) triple product fails on at least one (u, v) in Lambda x Lambda",
      found_untwisted_fail);

# ---------------------------------------------------------------------------
# Paper Section 5.1, Remark 5.2: no multiplicative identity in (R^24, *_sigma).
#
# Concretely, for u = (e_0, e_0, e_0) and any v = (x', y', z'):
#   block-1 of  u *_sigma v  =  e_0*x' + e_0*y' + e_0*z'  =  x' + y' + z',
# which is NOT v[1] in general.  So (e_0, e_0, e_0) is NOT a left identity.
# (And no element is, by the matrix argument in the paper's Remark.)
# ---------------------------------------------------------------------------

e0 := OctBasis(0);
u_putative := [e0, e0, e0];
v_test     := [ OctBasis(1) + 2 * OctBasis(2),
                OctBasis(3),
                OctBasis(4) + OctBasis(5) ];
prod := TripleProductSigma(u_putative, v_test);
expected_first_block := v_test[1] + v_test[2] + v_test[3];
CHECK("Paper Remark 5.2: (e_0,e_0,e_0) *_sigma v has block-1 = x'+y'+z', NOT x'  (so it is not an identity)",
      prod[1] = expected_first_block and prod[1] <> v_test[1]);

# Demonstrate that NO triple is a left identity, on the same v_test, by trying
# each candidate (a,b,c) with all components in J: a left identity would need
# to satisfy block-1 = x' for every v, including v_test, which constrains
# (a,b,c) inconsistently.  We take a representative search:
no_identity := true;
for a in J do
    for b in J do
        for c in J do
            uu := [a, b, c];
            if TripleProductSigma(uu, v_test) = v_test then
                # Check this candidate against another v also.
                v2 := [ OctBasis(2), OctBasis(0) - OctBasis(3), OctBasis(7) ];
                if TripleProductSigma(uu, v2) = v2 then
                    no_identity := false;
                fi;
            fi;
        od;
    od;
od;
CHECK("No (a,b,c) with components in J = +-{e_0,...,e_7} is a multiplicative identity",
      no_identity);

PrintTestSummary();
