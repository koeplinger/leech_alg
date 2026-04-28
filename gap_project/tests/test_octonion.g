#############################################################################
##
##  test_octonion.g — Octonion algebra basics (paper Section 2.1).
##
##  Verifies (mirrors python_project/tests/test_octonions.py):
##    1.  e_0 is left and right identity.
##    2.  e_i^2 = -e_0 for i = 1..7.
##    3.  Imaginary basis pairs anticommute.
##    4.  Algebra is non-commutative.
##    5.  Algebra is non-associative.
##    6.  Composition law N(xy) = N(x) N(y) on all 64 basis pairs.
##    7.  x * conj(x) = N(x) e_0 (and conj(x) * x).
##    8.  Conjugate is an involution.
##    9.  Left alternativity on all 64 basis pairs.
##    10. Right alternativity on all 64 basis pairs.
##    11. Flexibility on all 64 basis pairs.
##    12. Total antisymmetry of the associator on all 512 basis triples
##        (paper-strength result; certifies alternativity for all of R^8).
##    13. (Left) Moufang identity z*(x*(z*y)) = ((z*x)*z)*y on basis triples.
##    14. Division algebra: each basis element has a two-sided inverse.
##    15. Power-associativity on basis elements (cubes).
##    16. The seven Fano triples encode Dixon's rule e_a * e_{a+1} = e_{a+3}.
##
##  Plus, via the LOOPS package (independent verification):
##    L1. The 16 signed basis units form a loop of order 16.
##    L2. That loop is a Moufang loop.
##    L3. It is non-associative.
##    L4. It is isomorphic to MoufangLoop(16, 3) (the standard octonion loop).
##
#############################################################################

LoadPackage("loops");
Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
LoadAllSrc(GAP_PROJECT_ROOT);
ResetTestCounters();

Print("\n=== test_octonion.g ===\n");

# 1. Identity.
e0 := OctBasis(0);
CHECK("Left  identity:  e_0 * e_i = e_i for all i",
     ForAll([0..7], i -> OctMult(e0, OctBasis(i)) = OctBasis(i)));
CHECK("Right identity:  e_i * e_0 = e_i for all i",
     ForAll([0..7], i -> OctMult(OctBasis(i), e0) = OctBasis(i)));

# 2. Imaginary squares.
CHECK("e_i^2 = -e_0  for i = 1..7",
     ForAll([1..7], i -> OctMult(OctBasis(i), OctBasis(i)) = -e0));

# 3. Anticommutativity of distinct imaginary pairs.
ok := true;
for i in [1..7] do for j in [1..7] do
    if i <> j and OctMult(OctBasis(i), OctBasis(j)) <> -OctMult(OctBasis(j), OctBasis(i)) then
        ok := false; break;
    fi;
od; if not ok then break; fi; od;
CHECK("Distinct imaginary basis units anticommute  (e_i e_j = - e_j e_i)", ok);

# 4. Non-commutativity (sanity).
CHECK("Algebra is NOT commutative  (e_1 * e_2 <> e_2 * e_1)",
     OctMult(OctBasis(1), OctBasis(2)) <> OctMult(OctBasis(2), OctBasis(1)));

# 5. Non-associativity (sanity).
CHECK("Algebra is NOT associative  ((e_1*e_2)*e_3 <> e_1*(e_2*e_3))",
     OctMult(OctMult(OctBasis(1), OctBasis(2)), OctBasis(3))
       <> OctMult(OctBasis(1), OctMult(OctBasis(2), OctBasis(3))));

# 6. Composition law on basis pairs.
ok := true;
for i in [0..7] do for j in [0..7] do
    if OctNormSq(OctMult(OctBasis(i), OctBasis(j))) <> OctNormSq(OctBasis(i)) * OctNormSq(OctBasis(j)) then
        ok := false; break;
    fi;
od; if not ok then break; fi; od;
CHECK("Composition law  N(e_i e_j) = N(e_i) N(e_j)  on all 64 basis pairs", ok);

# 7. Conjugate formula and (8) involution.
ok := true;
for i in [0..7] do
    e := OctBasis(i);
    if OctMult(e, OctConjugate(e)) <> OctNormSq(e) * e0 then ok := false; break; fi;
    if OctMult(OctConjugate(e), e) <> OctNormSq(e) * e0 then ok := false; break; fi;
    if OctConjugate(OctConjugate(e)) <> e then ok := false; break; fi;
od;
CHECK("Conjugate: x*conj(x) = conj(x)*x = N(x) e_0  AND  conj is an involution", ok);

# 9. Left alternativity (x*x)*y = x*(x*y) on all 64 basis pairs.
ok := true;
for i in [0..7] do for j in [0..7] do
    x := OctBasis(i); y := OctBasis(j);
    if OctMult(OctMult(x, x), y) <> OctMult(x, OctMult(x, y)) then
        ok := false; break;
    fi;
od; if not ok then break; fi; od;
CHECK("Left  alternativity  (x*x)*y = x*(x*y)  on all 64 basis pairs", ok);

# 10. Right alternativity x*(y*y) = (x*y)*y on all 64 basis pairs.
ok := true;
for i in [0..7] do for j in [0..7] do
    x := OctBasis(i); y := OctBasis(j);
    if OctMult(x, OctMult(y, y)) <> OctMult(OctMult(x, y), y) then
        ok := false; break;
    fi;
od; if not ok then break; fi; od;
CHECK("Right alternativity  x*(y*y) = (x*y)*y  on all 64 basis pairs", ok);

# 11. Flexibility (x*y)*x = x*(y*x).
ok := true;
for i in [0..7] do for j in [0..7] do
    x := OctBasis(i); y := OctBasis(j);
    if OctMult(OctMult(x, y), x) <> OctMult(x, OctMult(y, x)) then
        ok := false; break;
    fi;
od; if not ok then break; fi; od;
CHECK("Flexibility  (x*y)*x = x*(y*x)  on all 64 basis pairs", ok);

# 12. Antisymmetry of the associator on all 8^3 = 512 basis triples.
ok := true;
for i in [0..7] do for j in [0..7] do for k in [0..7] do
    x := OctBasis(i); y := OctBasis(j); z := OctBasis(k);
    a_xyz := OctMult(OctMult(x, y), z) - OctMult(x, OctMult(y, z));
    a_yxz := OctMult(OctMult(y, x), z) - OctMult(y, OctMult(x, z));
    a_xzy := OctMult(OctMult(x, z), y) - OctMult(x, OctMult(z, y));
    a_yzx := OctMult(OctMult(y, z), x) - OctMult(y, OctMult(z, x));
    if a_xyz <> -a_yxz or a_xyz <> -a_xzy or a_xyz <> a_yzx then
        ok := false;
    fi;
od; od; od;
CHECK("Associator alternates on all 512 basis triples  (certifies alternativity on all of R^8)", ok);

# 13. (Left) Moufang identity z*(x*(z*y)) = ((z*x)*z)*y on basis triples.
ok := true;
for i in [0..7] do for j in [0..7] do for k in [0..7] do
    x := OctBasis(i); y := OctBasis(j); z := OctBasis(k);
    lhs := OctMult(z, OctMult(x, OctMult(z, y)));
    rhs := OctMult(OctMult(OctMult(z, x), z), y);
    if lhs <> rhs then ok := false; fi;
od; od; od;
CHECK("Moufang identity  z*(x*(z*y)) = ((z*x)*z)*y  on all 512 basis triples", ok);

# 14. Division-algebra: every nonzero basis element has two-sided inverse.
ok := true;
for i in [0..7] do
    e := OctBasis(i);
    inv_e := (1 / OctNormSq(e)) * OctConjugate(e);
    if OctMult(e, inv_e) <> e0 or OctMult(inv_e, e) <> e0 then
        ok := false; break;
    fi;
od;
CHECK("Division algebra: every basis e_i has two-sided inverse e_i^{-1} = conj(e_i)/N(e_i)", ok);

# 15. Power-associativity on basis: (e_i^2)*e_i = e_i*(e_i^2).
ok := true;
for i in [0..7] do
    x := OctBasis(i);
    x2 := OctMult(x, x);
    if OctMult(x2, x) <> OctMult(x, x2) then ok := false; break; fi;
od;
CHECK("Power-associativity  (x^2)*x = x*(x^2)  on basis", ok);

# 16. Standard Fano-plane multiplication rule e_a * e_{a+1} = e_{a+3}.
ok := true;
for a in [1..7] do
    b := ((a - 1 + 1) mod 7) + 1;             # a + 1 in {1..7}
    c := ((a - 1 + 3) mod 7) + 1;             # a + 3 in {1..7}
    if OctMult(OctBasis(a), OctBasis(b)) <> OctBasis(c) then
        ok := false; break;
    fi;
od;
CHECK("Dixon's standard rule  e_a * e_{a+1} = e_{a+3}  for a = 1..7  (paper eq. (1))", ok);

# ---------------------------------------------------------------------------
# LOOPS independent verification.
# ---------------------------------------------------------------------------

cayley := BuildSignedOctonionCayleyTable();
M := LoopByCayleyTable(cayley);

CHECK("LOOPS: the 16 signed octonion units form a loop of order 16",
     Size(M) = 16);
CHECK("LOOPS: the loop is a Moufang loop",
     IsMoufangLoop(M));
CHECK("LOOPS: the loop is NOT associative",
     not IsAssociative(M));
# Compare to the LOOPS-library octonion loop  MoufangLoop(16, 3):
M_ref := MoufangLoop(16, 3);
CHECK("LOOPS: our signed-units loop is isomorphic to MoufangLoop(16, 3)",
     IsomorphismLoops(M, M_ref) <> fail);

PrintTestSummary();
