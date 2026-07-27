#############################################################################
##
##  CORRECTION (2026-07-19): the paper's numbering changed after this file
##  was written.  Every "Theorem 1.2" below is Theorem 1.1 of the published
##  paper (doi:10.13140/RG.2.2.22093.19686, frozen source
##  paper/main_2026-07-19.tex), and "Definition 3.6" is Definition 3.3 there.
##  The mathematics and the tests are unaffected.
##
##  test_triple_product.g — The Z_3-routed triple product on R^24 (paper
##  Definition 3.6) and the main theorem (Theorem 1.2: closure of Lambda
##  under *_sigma).
##
##  This is the paper's main result, verified here on a representative set of
##  Leech minimal vectors:
##    - all 720 type-1 vectors  (one block of 2*lambda, two zero blocks);
##    - a sample of type-2 vectors  (two non-zero blocks);
##    - a sample of type-3 vectors  (three non-zero blocks).
##
##  Theorem 1.2 says  Lambda *_sigma Lambda  subset Lambda.  By bilinearity,
##  closure on a generating set of Lambda over Z is sufficient.  The 196,560
##  minimal vectors generate Lambda; we verify on type-1 (exhaustive) and on
##  type-2 / type-3 samples.  (The 12M+ random-pair test is in
##  python_project/src/trial_007_kirmse_twist.py; we replicate the symbolic
##  verification per Wilson's three conditions here.)
##
##  Also verifies (paper Section 4 "Algebraic properties"):
##    - The ambient algebra (R^24, *) has NO multiplicative identity
##      (paper Section 4 fix; (e_0, e_0, e_0) is NOT one).
##
#############################################################################

Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
LoadAllSrc(GAP_PROJECT_ROOT);
ResetTestCounters();

Print("\n=== test_triple_product.g  (paper Theorem 1.2) ===\n");

J := J_UNITS();
ROOTS := WilsonE8Roots();
ZEROVEC  := [0,0,0,0,0,0,0,0];

InLeechTriple := function(uv)
    return IsInLeech(uv[1], uv[2], uv[3]);
end;

# ---------------------------------------------------------------------------
# Theorem 1.2 spot-check: u *_sigma v  in Lambda  for sample u, v in Lambda.
# We use type-1 x type-1, then a small grid of type-2 x type-2 and
# type-3 x type-3 to exercise all three Wilson conditions non-trivially.
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

# Theorem 1.2:  Lambda *_sigma Lambda  subset Lambda  (sample-checked).
nfail := 0;
ntot  := 0;

# Type-1 x type-1: exhaustive over the 720 left operands times 6 right operands
# (2 from each type, including a type-3 vector to exercise condition 3 hardest).
right_ops := [
    T1[1], T1[241], T1[481],                       # one from each cyclic perm family
    type2_sample[1], type2_sample[10],
    type3_sample[1]
];
for u in T1 do
    for v in right_ops do
        if not InLeechTriple(TripleProductSigma(u, v)) then nfail := nfail + 1; fi;
        ntot := ntot + 1;
    od;
od;
CHECK(Concatenation("Theorem 1.2: type-1 x { 3 type-1, 2 type-2, 1 type-3 } closure (",
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
CHECK(Concatenation("Theorem 1.2: type-2 x { type-2 u type-3 } closure (",
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
CHECK(Concatenation("Theorem 1.2: type-3 x type-3 closure (", String(ntot), " products)"),
      nfail = 0);

# ---------------------------------------------------------------------------
# The standard (untwisted) triple product fails on Lambda — paper's pre-twist
# obstruction, and the companion's Example 4.5.
#
# By Definition 3.6, P_1+P_2+P_3 = (x+y+z)*(x'+y'+z'), and this is in Ls
# only if Ls*Ls subset Ls, which is false.  We exhibit a concrete failing
# (u, v) pair: take u and v to be type-3 Leech vectors whose triple sums
# are the witness (a, b) of Companion Example 4.5.
# ---------------------------------------------------------------------------

# Standard triple product fails on a witness pair drawn from Lambda.
# Construction: pick u, v in Lambda whose triple sums fall outside Ls *
# Ls's image in Ls.  We just iterate over the type-3 sample and report any
# untwisted-failure.
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
CHECK("Companion Example 4.7: standard (untwisted) triple product fails on at least one (u, v) in Lambda x Lambda",
      found_untwisted_fail);

# ---------------------------------------------------------------------------
# Section 4 of paper: no multiplicative identity in (R^24, *_sigma).
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
CHECK("Paper Section 4 fix: (e_0,e_0,e_0) *_sigma v has block-1 = x'+y'+z', NOT x'  (so it is not an identity)",
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
