#############################################################################
##
##  test_sublattices.g — The sublattices Ls, Lsbar, sigma(Ls), sigma(Lsbar).
##
##  Verifies (paper Section 2.3, Wilson [Wilson2009] Section 2):
##    1.  LS_BASIS and LS_BAR_BASIS each have rank 8.
##    2.  Their Gram determinants equal 256 = 16^2 = [L : Ls]^2 = [L : Lsbar]^2.
##    3.  IsInLs and IsInLsBar accept all 240 root products r*s, r*sbar.
##    4.  IsInLs and IsInLsBar reject every E8 root r itself (min norm in Ls
##        is 2*2 = 4, while N(r) = 2).
##    5.  2L subset Ls and 2L subset Lsbar  (Wilson:  2L = Ls cap Lsbar).
##    6.  Ls + Lsbar generates L (joint Z-span has rank 8).
##    7.  sigma(L) = L  (paper Lemma 4.1).
##    8.  sigma(Ls) ne Ls  (paper Remark 4.5; companion Section 7).
##        Explicit witness: a sigma-image of an Ls basis vector that is not in Ls.
##    9.  Ls * Ls  not subset Ls  (companion Example 4.5; main paper Remark 4.5).
##        Explicit witness: an Ls-basis pair whose product is not in Ls.
##    10. The companion's specific witness vector
##           (-1, 0, 1, 0, 0, 1, 0, 1)  in sigma(Ls) but not in Ls.
##
#############################################################################

Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
LoadAllSrc(GAP_PROJECT_ROOT);
ResetTestCounters();

Print("\n=== test_sublattices.g ===\n");

ALL_ROOTS := WilsonE8Roots();

# 1. Rank 8.
CHECK("LS_BASIS has rank 8",
      RankMat(LS_BASIS) = 8);
CHECK("LS_BAR_BASIS has rank 8",
      RankMat(LS_BAR_BASIS) = 8);

# 2. Gram determinant = 256 (index 16 in unimodular L).
gram_Ls    := List(LS_BASIS,     r -> List(LS_BASIS,     s -> r * s));
gram_Lsbar := List(LS_BAR_BASIS, r -> List(LS_BAR_BASIS, s -> r * s));
CHECK("det Gram(Ls)     = 256  (so [L : Ls]    = 16)",
      DeterminantMat(gram_Ls) = 256);
CHECK("det Gram(Lsbar)  = 256  (so [L : Lsbar] = 16)",
      DeterminantMat(gram_Lsbar) = 256);

# 3. IsInLs / IsInLsBar accept all 240 root products.
CHECK("IsInLs accepts r * s for all 240 roots r",
      ForAll(ALL_ROOTS, r -> IsInLs(OctMult(r, WILSON_S))));
CHECK("IsInLsBar accepts r * sbar for all 240 roots r",
      ForAll(ALL_ROOTS, r -> IsInLsBar(OctMult(r, WILSON_S_BAR))));

# 4. IsInLs / IsInLsBar reject every E8 root  (min norm in Ls / Lsbar is 4).
CHECK("IsInLs rejects every E8 root  (norm 2 < min norm 4 in Ls)",
      ForAll(ALL_ROOTS, r -> not IsInLs(r)));
CHECK("IsInLsBar rejects every E8 root  (norm 2 < min norm 4 in Lsbar)",
      ForAll(ALL_ROOTS, r -> not IsInLsBar(r)));

# 5. 2L subset Ls and 2L subset Lsbar.
CHECK("2L subset Ls    (2*r in Ls    for all 240 roots r)",
      ForAll(ALL_ROOTS, r -> IsInLs(2 * r)));
CHECK("2L subset Lsbar (2*r in Lsbar for all 240 roots r)",
      ForAll(ALL_ROOTS, r -> IsInLsBar(2 * r)));

# 6. Ls + Lsbar generates L (joint span has rank 8).
joint := Concatenation(LS_BASIS, LS_BAR_BASIS);
CHECK("Z-span of LS_BASIS u LS_BAR_BASIS has rank 8  (Ls + Lsbar = L)",
      RankMat(joint) = 8);

# 7. sigma(L) = L  (paper Lemma 4.1).
CHECK("Lemma 4.1: sigma(b) is a Z-combination of L_BASIS for all 8 basis b",
      ForAll(L_BASIS, b -> IsIntegerCombination(ApplySigma(b), L_BASIS)));

# 8. sigma(Ls) ne Ls  (Remark 4.5).
witness_index := 0;
for i in [1..8] do
    if not IsIntegerCombination(SIGMA_LS_BASIS[i], LS_BASIS) then
        witness_index := i;
        break;
    fi;
od;
CHECK("Remark 4.5:  sigma(Ls) ne Ls    (some sigma(Ls)-basis vector is not in Ls)",
      witness_index > 0);

# 9. Ls * Ls  not subset Ls  (companion Example 4.5).
counter_pair := fail;
for i in [1..8] do
    for j in [1..8] do
        if not IsInLs(OctMult(LS_BASIS[i], LS_BASIS[j])) then
            counter_pair := [i, j];
            break;
        fi;
    od;
    if counter_pair <> fail then break; fi;
od;
CHECK("Companion Example 4.5:  Ls * Ls  NOT subset Ls   (explicit witness pair found)",
      counter_pair <> fail);

# 10. Companion's specific witness vector.
witness := [-1, 0, 1, 0, 0, 1, 0, 1];
CHECK("Companion witness (-1,0,1,0,0,1,0,1)  is in sigma(Ls)",
      IsInSigmaLs(witness));
CHECK("Companion witness (-1,0,1,0,0,1,0,1)  is NOT in Ls",
      not IsInLs(witness));

PrintTestSummary();
