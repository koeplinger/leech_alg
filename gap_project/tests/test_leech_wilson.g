#############################################################################
##
##  test_leech_wilson.g — Wilson's Leech lattice Lambda inside L^3
##  (paper Definition 2.4, Wilson [Wilson2009] Section 3).
##
##  Mirrors python_project/tests/test_leech_wilson.py (the parts that don't
##  depend on enumerating the 184,320 type-3 vectors, which would be slow in
##  GAP).  We verify all type-1 vectors directly and use small samples of
##  type-2 / type-3 generators.
##
##  Verifies:
##    1.  J = { +-e_t : t = 0..7 } has 16 elements, each of squared norm 1.
##    2.  All 240 root products r * sbar lie in L  (Lsbar subset L).
##    3.  All 240 root products r * s    lie in L  (Ls    subset L).
##    4.  IsInLsBar / IsInLs reject every E8 root  (norm 2 < min 4 in Ls(*)).
##    5.  Type-1 minimal vector count = 720 = 3 * 240.
##    6.  Every type-1 minimal vector has ambient squared norm 8.
##    7.  Every type-1 minimal vector satisfies Wilson's three conditions.
##    8.  All 720 type-1 vectors are distinct.
##    9.  A small sample of type-2 vectors (lambda*sbar, (lambda*sbar)*j, 0)
##        satisfies all three conditions.
##    10. A small sample of type-3 vectors ((lambda*s)*j, lambda*k, (lambda*j)*k)
##        satisfies all three conditions.
##    11. (root, 0, 0) is NOT in Lambda  (root has norm 2; Lsbar min norm = 4
##        so condition 2 fails).
##    12. (root, root, 0) is NOT in Lambda  (x+z = root not in Lsbar).
##
#############################################################################

Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
LoadAllSrc(GAP_PROJECT_ROOT);
ResetTestCounters();

Print("\n=== test_leech_wilson.g ===\n");

J := J_UNITS();
ROOTS := WilsonE8Roots();

# 1. J basics.
CHECK("|J| = 16",
      Length(J) = 16);
CHECK("Every j in J has squared norm 1",
      ForAll(J, j -> OctNormSq(j) = 1));

# 2. Lsbar subset L on roots.
CHECK("Lsbar subset L:  r * sbar in L for all 240 roots",
      ForAll(ROOTS, r -> IsInL(OctMult(r, WILSON_S_BAR))));

# 3. Ls subset L on roots.
CHECK("Ls subset L:  r * s in L for all 240 roots",
      ForAll(ROOTS, r -> IsInL(OctMult(r, WILSON_S))));

# 4. Sublattices reject roots.
CHECK("IsInLsBar rejects every E8 root  (min norm in Lsbar is 4)",
      ForAll(ROOTS, r -> not IsInLsBar(r)));
CHECK("IsInLs    rejects every E8 root  (min norm in Ls    is 4)",
      ForAll(ROOTS, r -> not IsInLs(r)));

# 5-8. Type-1 minimal vectors.
T1 := LeechType1Vectors();
CHECK("Type-1 minimal vector count = 720  (3 * 240)",
      Length(T1) = 720);
CHECK("Every type-1 minimal vector has ambient squared norm 8",
      ForAll(T1, t -> Sum(t, c -> OctNormSq(c)) = 8));
CHECK("Every type-1 minimal vector satisfies Wilson's three conditions",
      ForAll(T1, t -> IsInLeech(t[1], t[2], t[3])));
CHECK("All 720 type-1 vectors are distinct",
      Length(Set(T1)) = 720);

# 9. Sample type-2 vectors:  ( lambda*sbar,  (lambda*sbar)*j,  0 )  and cyclic perms.
type2_sample_ok := true;
sampled := 0;
for ridx in [1, 60, 120, 180, 240] do          # 5 lambdas spread across 240 roots
    for jidx in [1, 5, 9, 13] do               # 4 j's spread across 16 unit octonions
        lam   := ROOTS[ridx];
        j     := J[jidx];
        a     := OctMult(lam, WILSON_S_BAR);
        b     := OctMult(a, j);
        zero  := [0,0,0,0,0,0,0,0];
        if not IsInLeech(a, b, zero) then type2_sample_ok := false; fi;
        if not IsInLeech(zero, a, b) then type2_sample_ok := false; fi;
        if not IsInLeech(b, zero, a) then type2_sample_ok := false; fi;
        sampled := sampled + 3;
    od;
od;
CHECK(Concatenation("Type-2 sample (", String(sampled), " vectors) satisfies Wilson's three conditions"),
      type2_sample_ok);

# 10. Sample type-3 vectors:  ( (lambda*s)*j,  lambda*k,  (lambda*j)*k )  and perms.
type3_sample_ok := true;
sampled := 0;
for ridx in [1, 60, 120, 180, 240] do
    for jidx in [1, 5] do
        for kidx in [1, 9] do
            lam := ROOTS[ridx];
            j   := J[jidx];
            k   := J[kidx];
            a   := OctMult(OctMult(lam, WILSON_S), j);
            b   := OctMult(lam, k);
            c   := OctMult(OctMult(lam, j), k);
            if not IsInLeech(a, b, c) then type3_sample_ok := false; fi;
            sampled := sampled + 1;
        od;
    od;
od;
CHECK(Concatenation("Type-3 sample (", String(sampled), " vectors) satisfies Wilson's three conditions"),
      type3_sample_ok);

# 11. Non-members.
zero := [0,0,0,0,0,0,0,0];
CHECK("(root, 0, 0) is NOT in Lambda for any of the 240 roots",
      ForAll(ROOTS, r -> not IsInLeech(r, zero, zero)));
CHECK("(root, root, 0) is NOT in Lambda for any of the first 20 roots",
      ForAll(ROOTS{[1..20]}, r -> not IsInLeech(r, r, zero)));

PrintTestSummary();
