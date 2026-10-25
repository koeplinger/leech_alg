#############################################################################
##
##  test_twist.g — The sigma-twisted octonion product (paper Definitions
##  3.1 and 3.2, `def:twist` and `def:twisted-product`), for
##  sigma = (1 2).  All transposition twists are equivalent up to
##  relabeling (paper Section 4, `sec:proof`, footnote to sigma = (1 2)).
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
##
##  Verifies:
##    1.  ApplySigma is an involution: sigma o sigma = id  on R^8.
##    2.  ApplySigma is an isometry: it preserves the Euclidean norm.
##    3.  sigma( x * y )  =  sigma(x) *_sigma sigma(y): sigma is an algebra
##        isomorphism from (O, *) to (O, *_sigma), as the paper states after
##        Definition 3.2.
##    4.  *_sigma satisfies the composition law N(x *_sigma y) = N(x) N(y) on
##        all 64 basis pairs.
##    5.  *_sigma differs from * on at least one basis pair: it is a
##        different bilinear product on R^8, isomorphic to * via sigma
##        (item 3).
##    6.  Remark 4.5 obstruction restated through the twist:
##         Ls * Ls  not subset Ls            (standard product fails on Ls)
##         Ls *_sigma Ls   subset Ls         (twisted product succeeds on Ls)
##    7.  Companion Section 10: the twisted product preserves L:
##           L *_sigma L  subset L
##        (explicit on the 64 L-basis products).  Hence Coxeter's E_8 closure
##        is not undone by the sigma-twist.
##
#############################################################################

Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
LoadAllSrc(GAP_PROJECT_ROOT);
ResetTestCounters();

Print("\n=== test_twist.g  (paper Section 3) ===\n");

# 1. Involution.
CHECK("ApplySigma is an involution  (sigma o sigma = id  on the 8 basis vectors)",
      ForAll([0..7], i -> ApplySigma(ApplySigma(OctBasis(i))) = OctBasis(i)));

# 2. Isometry.
CHECK("ApplySigma is an isometry  (preserves squared norm on basis vectors)",
      ForAll([0..7], i -> OctNormSq(ApplySigma(OctBasis(i))) = OctNormSq(OctBasis(i))));

# 3. sigma(x*y) = sigma(x) *_sigma sigma(y)  (paper, after Definition 3.2).
ok := true;
for i in [0..7] do for j in [0..7] do
    x := OctBasis(i); y := OctBasis(j);
    if ApplySigma(OctMult(x, y)) <> OctMultSigma(ApplySigma(x), ApplySigma(y)) then
        ok := false; break;
    fi;
od; if not ok then break; fi; od;
CHECK("sigma is an algebra isomorphism:  sigma(x*y) = sigma(x) *_sigma sigma(y)  on all 64 basis pairs",
      ok);

# 4. Composition law for *_sigma.
ok := true;
for i in [0..7] do for j in [0..7] do
    x := OctBasis(i); y := OctBasis(j);
    if OctNormSq(OctMultSigma(x, y)) <> OctNormSq(x) * OctNormSq(y) then
        ok := false; break;
    fi;
od; if not ok then break; fi; od;
CHECK("Composition law for *_sigma:  N(x *_sigma y) = N(x) N(y)  on all 64 basis pairs",
      ok);

# 5. *_sigma differs from * on at least one basis pair.
diff := false;
for i in [0..7] do for j in [0..7] do
    if OctMultSigma(OctBasis(i), OctBasis(j)) <> OctMult(OctBasis(i), OctBasis(j)) then
        diff := true;
    fi;
od; od;
CHECK("*_sigma is a genuinely different bilinear product on R^8  (some basis pair disagrees)",
      diff);

# 6. Twist resolves the Ls obstruction.
ls_fails := false;
for i in [1..8] do for j in [1..8] do
    if not IsInLs(OctMult(LS_BASIS[i], LS_BASIS[j])) then
        ls_fails := true;
    fi;
od; od;
CHECK("Standard *  fails on Ls    (Ls * Ls    not subset Ls)",
      ls_fails);

twisted_ok := true;
for i in [1..8] do for j in [1..8] do
    if not IsInLs(OctMultSigma(LS_BASIS[i], LS_BASIS[j])) then
        twisted_ok := false;
    fi;
od; od;
CHECK("Twisted *_sigma succeeds on Ls   (Ls *_sigma Ls   subset Ls)",
      twisted_ok);

# 7. L *_sigma L  subset L  (companion Section 10: the twist keeps L closed).
nfail := 0;
for i in [1..8] do for j in [1..8] do
    if not IsIntegerCombination(OctMultSigma(L_BASIS[i], L_BASIS[j]), L_BASIS) then
        nfail := nfail + 1;
    fi;
od; od;
CHECK("Companion Section 10:  L *_sigma L  subset L   (twist preserves L)",
      nfail = 0);

PrintTestSummary();
