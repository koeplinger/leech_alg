#############################################################################
##
##  test_companion_examples.g — The explicit hand-checkable numerical
##  examples in the companion note (paper/companion.tex), recomputed in GAP
##  to give the reader a fully independent re-derivation.
##
##  Verifies (companion Sections 2-8):
##    1.  e_1 * e_2  =  e_4   (companion Section 3, Example after the
##        Fano-triple list).
##    2.  e_1 * e_3  =  e_7   (same example, second multiplication).
##    3.  e_2 * e_7  = -e_6   (same example, third multiplication).
##    4.  s * s  =  -3/2 e_0  -  (1/2)(e_1 + e_2 + ... + e_7)   in L.
##        (Companion Section 3, "A sample product".)
##    5.  Companion Example 4.5:
##           a = (-1, 0, 0, 0, 1, 0, 1, 1)   and
##           b = (-1, 1, 0, 0, 0, 1, 0, 1)
##        both lie in Ls;  their octonion product is
##           a * b = (0, -2, 2, 1, -2, -1, -1, -1) in L,
##        and a * b is NOT in Ls.
##    6.  Companion Example 4.6 / Section 5: under  sigma = (1 2),
##           e_1 *_sigma e_2  =  -e_4    (vs.  e_1 * e_2 = +e_4).
##    7.  Companion Section 8: sigma fixes L  (sigma(L) = L) but moves Ls
##        (sigma(Ls) ne Ls), with explicit witness  (-1, 0, 1, 0, 0, 1, 0, 1).
##    8.  Companion final box ("Summary of concrete data"):
##         (i)   sigma(Ls) ne Ls
##         (ii)  Ls * Ls  not subset Ls
##         (iii) sigma(Ls) * sigma(Ls)  subset sigma(Ls).
##
##  All arithmetic is exact (rational); no floating-point.
##
#############################################################################

Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
LoadAllSrc(GAP_PROJECT_ROOT);
ResetTestCounters();

Print("\n=== test_companion_examples.g ===\n");

# 1-3. Specific multiplication examples (companion Section 3, after the
#      Fano-triple list).
e1 := OctBasis(1);
e2 := OctBasis(2);
e3 := OctBasis(3);
e4 := OctBasis(4);
e6 := OctBasis(6);
e7 := OctBasis(7);
CHECK("Companion Example (Section 3):  e_1 * e_2  =  e_4",
      OctMult(e1, e2) =  e4);
CHECK("Companion Example (Section 3):  e_1 * e_3  =  e_7",
      OctMult(e1, e3) =  e7);
CHECK("Companion Example (Section 3):  e_2 * e_7  = -e_6",
      OctMult(e2, e7) = -e6);

# 4. s*s = -3/2 e_0 - (1/2) sum_{k=1..7} e_k, in L.
expected_ss := [-3/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2];
ss := OctMult(WILSON_S, WILSON_S);
CHECK("Companion Section 3: s * s  =  -3/2 e_0  -  (1/2)(e_1 + ... + e_7)",
      ss = expected_ss);
CHECK("Companion Section 3: s * s  is in L",
      IsInL(ss));

# 5. Companion Example 4.5: explicit a, b in Ls with a*b not in Ls.
a := [-1, 0, 0, 0, 1, 0, 1, 1];
b := [-1, 1, 0, 0, 0, 1, 0, 1];
ab_expected := [ 0, -2, 2, 1, -2, -1, -1, -1 ];

CHECK("Companion Example 4.5: a in Ls",
      IsInLs(a));
CHECK("Companion Example 4.5: b in Ls",
      IsInLs(b));
CHECK("Companion Example 4.5: a * b = (0, -2, 2, 1, -2, -1, -1, -1)",
      OctMult(a, b) = ab_expected);
CHECK("Companion Example 4.5: a * b is in L",
      IsInL(OctMult(a, b)));
CHECK("Companion Example 4.5: a * b is NOT in Ls   (so Ls * Ls not subset Ls)",
      not IsInLs(OctMult(a, b)));

# 6. Companion Example (Section 5): under sigma = (1 2), e_1 *_sigma e_2 = -e_4.
CHECK("Companion Section 5: e_1 *_sigma e_2  =  -e_4   (vs  e_1 * e_2 = +e_4)",
      OctMultSigma(e1, e2) = -e4);

# 7. Companion Section 8: sigma fixes L (Lemma A) but moves Ls (sigma(Ls) ne Ls).
CHECK("Companion Section 8: sigma(L) = L   (sigma(b) is a Z-comb of L_BASIS for all 8 basis b)",
      ForAll(L_BASIS, b -> IsIntegerCombination(ApplySigma(b), L_BASIS)));

witness := [-1, 0, 1, 0, 0, 1, 0, 1];
CHECK("Companion Section 8 witness  (-1, 0, 1, 0, 0, 1, 0, 1)  in sigma(Ls)",
      IsInSigmaLs(witness));
CHECK("Companion Section 8 witness  (-1, 0, 1, 0, 0, 1, 0, 1)  NOT in Ls",
      not IsInLs(witness));

# 8. Companion final box: three concrete data points.
sigma_Ls_neq_Ls := false;
for v in SIGMA_LS_BASIS do
    if not IsInLs(v) then sigma_Ls_neq_Ls := true; fi;
od;
CHECK("Companion final summary (i):   sigma(Ls) ne Ls",
      sigma_Ls_neq_Ls);

ls_not_closed := false;
for v in LS_BASIS do for w in LS_BASIS do
    if not IsInLs(OctMult(v, w)) then ls_not_closed := true; fi;
od; od;
CHECK("Companion final summary (ii):  Ls * Ls  not subset Ls",
      ls_not_closed);

sigma_ls_closed := true;
for v in SIGMA_LS_BASIS do for w in SIGMA_LS_BASIS do
    if not IsInSigmaLs(OctMult(v, w)) then sigma_ls_closed := false; fi;
od; od;
CHECK("Companion final summary (iii): sigma(Ls) * sigma(Ls)  subset sigma(Ls)",
      sigma_ls_closed);

PrintTestSummary();
