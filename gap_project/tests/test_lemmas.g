#############################################################################
##
##  test_lemmas.g — The four lemmas of the symbolic proof in
##  paper Section 4 (and the accompanying non-triviality remark).
##
##  Mirrors python_project/src/symbolic_proof_checks.py.
##
##  All checks are exact (rational arithmetic in GAP); no floating point.
##
##  Lemma 4.1 (paper):  sigma(L) = L.
##  Lemma 4.2 (paper):  L * L  subset L          (Coxeter / maximal order).
##  Lemma 4.3 (paper):  L * sigma(Lsbar)  subset sigma(Lsbar).
##  Lemma 4.4 (paper):  sigma(Ls) * sigma(Ls)  subset sigma(Ls).
##
##  Remark 4.5:  sigma(Ls) ne Ls   AND   Ls * Ls  not subset Ls.
##  Bonus (untwisted condition 2):   L * Lsbar  subset Lsbar.
##
##  Each lemma reduces, by Z-bilinearity, to a finite check on basis pairs.
##
#############################################################################

Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
LoadAllSrc(GAP_PROJECT_ROOT);
ResetTestCounters();

Print("\n=== test_lemmas.g  (paper Section 4) ===\n");

# ---------------------------------------------------------------------------
# Lemma 4.1: sigma(L) = L.
# ---------------------------------------------------------------------------

CHECK("Lemma 4.1:  sigma(b) in Z-span(L_BASIS)  for all 8 basis vectors of L",
      ForAll(L_BASIS, b -> IsIntegerCombination(ApplySigma(b), L_BASIS)));

# ---------------------------------------------------------------------------
# Lemma 4.2: L * L subset L  (already checked in test_e8_wilson over all 240^2
# root pairs; here we re-check on the 8x8 = 64 Z-basis products as the
# minimal-Z-bilinear certification).
# ---------------------------------------------------------------------------

nfail := 0;
for i in [1..8] do
    for j in [1..8] do
        if not IsIntegerCombination(OctMult(L_BASIS[i], L_BASIS[j]), L_BASIS) then
            nfail := nfail + 1;
        fi;
    od;
od;
CHECK("Lemma 4.2:  L * L  subset L   (all 64 basis-basis products in Z-span(L_BASIS))",
      nfail = 0);

# ---------------------------------------------------------------------------
# Lemma 4.3:  L * sigma(Lsbar) subset sigma(Lsbar).
# ---------------------------------------------------------------------------

nfail := 0;
for i in [1..8] do
    for j in [1..8] do
        if not IsIntegerCombination(OctMult(L_BASIS[i], SIGMA_LS_BAR_BASIS[j]),
                                    SIGMA_LS_BAR_BASIS) then
            nfail := nfail + 1;
        fi;
    od;
od;
CHECK("Lemma 4.3:  L * sigma(Lsbar)  subset sigma(Lsbar)   (all 64 basis-basis products)",
      nfail = 0);

# ---------------------------------------------------------------------------
# Lemma 4.4:  sigma(Ls) * sigma(Ls) subset sigma(Ls).
# ---------------------------------------------------------------------------

nfail := 0;
for i in [1..8] do
    for j in [1..8] do
        if not IsIntegerCombination(OctMult(SIGMA_LS_BASIS[i], SIGMA_LS_BASIS[j]),
                                    SIGMA_LS_BASIS) then
            nfail := nfail + 1;
        fi;
    od;
od;
CHECK("Lemma 4.4:  sigma(Ls) * sigma(Ls)  subset sigma(Ls)   (all 64 basis-basis products)",
      nfail = 0);

# ---------------------------------------------------------------------------
# Remark 4.5:  sigma(Ls) ne Ls   AND   Ls * Ls  not subset Ls.
# ---------------------------------------------------------------------------

# sigma(Ls) ne Ls.
witness := false;
for i in [1..8] do
    if not IsIntegerCombination(SIGMA_LS_BASIS[i], LS_BASIS) then
        witness := true; break;
    fi;
od;
CHECK("Remark 4.5(a):  sigma(Ls) ne Ls   (a sigma(Ls)-basis vector is not in Ls)",
      witness);

# Ls * Ls  not subset Ls.
counter := false;
for i in [1..8] do
    for j in [1..8] do
        if not IsInLs(OctMult(LS_BASIS[i], LS_BASIS[j])) then
            counter := true; break;
        fi;
    od;
    if counter then break; fi;
od;
CHECK("Remark 4.5(b):  Ls * Ls  NOT subset Ls   (basis pair with product outside Ls)",
      counter);

# ---------------------------------------------------------------------------
# Bonus:  L * Lsbar  subset Lsbar.   This is the untwisted condition 2.
# ---------------------------------------------------------------------------

nfail := 0;
for i in [1..8] do
    for j in [1..8] do
        if not IsIntegerCombination(OctMult(L_BASIS[i], LS_BAR_BASIS[j]), LS_BAR_BASIS) then
            nfail := nfail + 1;
        fi;
    od;
od;
CHECK("Bonus:  L * Lsbar  subset Lsbar  (untwisted condition 2 also closes)",
      nfail = 0);

PrintTestSummary();
