#############################################################################
##
##  test_e8_wilson.g — Wilson's E_8 lattice L = D_8^+ (paper Section 2.2,
##  Wilson [Wilson2009] Section 2).  Mirrors python_project/tests/test_e8_wilson.py.
##
##  Verifies:
##    1.  Type-1 root count = 112.
##    2.  Type-2 root count = 128.
##    3.  Total root count  = 240, all distinct.
##    4.  Type-1 roots have exactly two non-zero coords, each +-1.
##    5.  Type-2 roots have all coords +-1/2.
##    6.  All 240 roots have squared norm 2.
##    7.  All 240 roots pass IsInL.
##    8.  Single basis elements e_i are NOT in L (integer coords, sum 1, odd).
##    9.  Doubled basis elements 2*e_i ARE in L.
##    10. (1/2)(e_0 + e_1 + ... + e_7) (zero minus signs) is NOT in L.
##    11. Full 240x240 Gram matrix has integer entries (integer lattice).
##    12. The 240 roots span Q^8 (rank 8).
##    13. A Z-basis built from WILSON_S + 7 type-1 roots has |det| = 1
##        (L is unimodular / self-dual).
##    14. Wilson's distinguished elements:
##        s = (-1/2, 1/2, ..., 1/2) is in L (norm 2, type-2 root).
##        sbar = (-1/2, -1/2, ..., -1/2) has norm 2 but is NOT in L
##        (paper Section 2.3: even sum -4 falls in the excluded branch).
##    15. Even lattice: pairwise sums and differences of roots have even norm.
##    16. L * L  subset L  (all 57600 root products lie in L).
##        This is the Coxeter / paper Lemma 4.2 closure.
##    17. r * sbar in L for all 240 roots  (Lsbar subset L; Wilson Section 2).
##    18. r * s    in L for all 240 roots  (Ls    subset L; Wilson Section 2).
##
#############################################################################

Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
LoadAllSrc(GAP_PROJECT_ROOT);
ResetTestCounters();

Print("\n=== test_e8_wilson.g ===\n");

T1 := WilsonType1Roots();
T2 := WilsonType2Roots();
ALL_ROOTS := WilsonE8Roots();

# 1-3. Counts.
CHECK("Type-1 root count = 112",  Length(T1) = 112);
CHECK("Type-2 root count = 128",  Length(T2) = 128);
CHECK("Total root count  = 240",  Length(ALL_ROOTS) = 240);
CHECK("All 240 roots are distinct",
      Length(Set(ALL_ROOTS)) = 240);

# 4-5. Coordinate structure.
CHECK("Type-1 roots have exactly two non-zero coords, each +-1",
      ForAll(T1, r -> Number(r, x -> x <> 0) = 2 and
                       ForAll(r, x -> x in [-1, 0, 1])));
CHECK("Type-2 roots have all coords +-1/2",
      ForAll(T2, r -> ForAll(r, x -> x in [-1/2, 1/2])));

# 6. Norm.
CHECK("All 240 roots have squared norm 2",
      ForAll(ALL_ROOTS, r -> OctNormSq(r) = 2));

# 7-10. Membership predicate.
CHECK("All 240 roots pass IsInL",
      ForAll(ALL_ROOTS, IsInL));
CHECK("Single basis elements e_i are NOT in L",
      ForAll([0..7], i -> not IsInL(OctBasis(i))));
CHECK("Doubled basis elements 2*e_i ARE in L",
      ForAll([0..7], i -> IsInL(2 * OctBasis(i))));
CHECK("(1/2)(e_0+...+e_7) is NOT in L  (zero minus signs, even sum)",
      not IsInL([1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2]));

# 11. Integer Gram matrix.
gram := List(ALL_ROOTS, r -> List(ALL_ROOTS, s -> r * s));   # 240x240
CHECK("All 240*240 root inner products are integers",
      ForAll(gram, row -> ForAll(row, x -> IsInt(x))));

# 12. Rank.
CHECK("The 240 roots span Q^8 (rank 8)",
      RankMat(ALL_ROOTS) = 8);

# 13. Unimodular Z-basis.
CHECK("Z-basis WILSON_S + 7 type-1 roots is unimodular  (|det| = 1)",
      AbsoluteValue(DeterminantMat(L_BASIS)) = 1);

# 14. Wilson's s and sbar.
CHECK("WILSON_S has squared norm 2",
      OctNormSq(WILSON_S) = 2);
CHECK("WILSON_S is in L  (type-2 root, sum = 3 = odd)",
      IsInL(WILSON_S));
CHECK("WILSON_S_BAR has squared norm 2",
      OctNormSq(WILSON_S_BAR) = 2);
CHECK("WILSON_S_BAR is NOT in L  (paper Section 2.3, sum = -4 = even, excluded branch)",
      not IsInL(WILSON_S_BAR));

# 15. Even lattice: sums and differences of roots have even squared norm.
ok := true;
for r1 in ALL_ROOTS{[1..20]} do
    for r2 in ALL_ROOTS{[1..20]} do
        if not (OctNormSq(r1 + r2) mod 2 = 0) then ok := false; fi;
        if not (OctNormSq(r1 - r2) mod 2 = 0) then ok := false; fi;
    od;
od;
CHECK("Sample sums and differences of roots have even squared norm  (even lattice)", ok);

# 16. Coxeter closure  L * L  subset L  on all 240^2 root pairs (paper Lemma 4.2).
nfail := 0;
for r1 in ALL_ROOTS do
    for r2 in ALL_ROOTS do
        if not IsInL(OctMult(r1, r2)) then nfail := nfail + 1; fi;
    od;
od;
CHECK("L * L  subset L  on all 240*240 = 57,600 root pairs  (Coxeter / Lemma 4.2)",
      nfail = 0);

# 17-18. Lsbar and Ls subset L  (Wilson Section 2; main paper, fix to error 2.1).
nfail := 0;
for r in ALL_ROOTS do
    if not IsInL(OctMult(r, WILSON_S_BAR)) then nfail := nfail + 1; fi;
od;
CHECK("r * sbar  in L  for all 240 roots  (so Lsbar subset L despite sbar notin L)",
      nfail = 0);

nfail := 0;
for r in ALL_ROOTS do
    if not IsInL(OctMult(r, WILSON_S)) then nfail := nfail + 1; fi;
od;
CHECK("r * s     in L  for all 240 roots  (Ls subset L)",
      nfail = 0);

PrintTestSummary();
