#############################################################################
##
##  run_01_setup_checks.g                                          (GAP 4.12)
##
##  Purpose: validate every ingredient before any twist is tested.
##    A. octonion product: composition, alternativity, the paper's own
##       worked example (Remark 4.5) and all 192 rows of the paper's
##       Appendix-A coefficient tables (LL, LN, MM), which pins down the
##       Fano-sign convention, the left-multiple convention for L s and
##       L sbar, and the action of sigma = (1 2).
##    B. L = D8+ : 240 roots, a Z-basis, L.L in L.
##    C. L sbar, L s : indices 16, sum = L, 2L inside both (=> intersection 2L).
##    D. Lambda : a Z-basis from the mod-2 description, index 2^12 in L^3,
##       Gram determinant 2^24, even (norms in 4Z, inner products in 2Z).
##    E. minimum 8 and 196,560 minimal vectors (GAP ShortestVectors on
##       Gram/2, i.e. the unimodular scaling: no vectors of norm <= 3 and
##       98,280 +- pairs of norm 4), which identifies Lambda as sqrt(2) Leech.
##    F. the paper's three families of minimal vectors (Sec. 2.3) are exactly
##       the minimal shell found in E.
##
##  Conventions: see the header of leech_lib.g (row vectors, index k <->
##  e_{k-1}, g(x) = x*G, pi(e_a) = e_{a^p}).
##
##  Run:  gap -q -o 4g run_01_setup_checks.g < /dev/null > run_01_setup_checks.out
##
#############################################################################

Read("leech_lib.g");
Read("paper_tables.g");     # produced by extract_paper_tables.sh from main.tex

Reset(GlobalMersenneTwister, 20261013);   # fixed seed: reproducible random tests

ok := function(name, cond)
  if cond = true then Print("  [ok]   ", name, "\n");
  else Print("  [FAIL] ", name, "\n"); fi;
end;

Print("=== run_01_setup_checks.g ===\n");
Print("GAP ", GAPInfo.Version, "\n\n");

#############################################################################
Print("A. Octonion product\n");
#############################################################################
N8 := x -> x * x;           # norm form N(x) = sum x_i^2
randv := function() return List([1..8], i -> Random([-6..6]) / 2); end;
pairs := List([1..500], i -> [randv(), randv()]);;
ok("composition N(xy) = N(x)N(y) on 500 random rational pairs",
   ForAll(pairs, p -> N8(OctMul(p[1], p[2])) = N8(p[1]) * N8(p[2])));
ok("left alternativity (xx)y = x(xy) on 500 random pairs",
   ForAll(pairs, p -> OctMul(OctMul(p[1], p[1]), p[2]) = OctMul(p[1], OctMul(p[1], p[2]))));
ok("right alternativity (yx)x = y(xx) on 500 random pairs",
   ForAll(pairs, p -> OctMul(OctMul(p[2], p[1]), p[1]) = OctMul(p[2], OctMul(p[1], p[1]))));
ok("e_i e_j = -e_j e_i for distinct imaginary units",
   ForAll(Combinations([2..8], 2), c -> OctMul(UnitVec8(c[1]), UnitVec8(c[2])) = -OctMul(UnitVec8(c[2]), UnitVec8(c[1]))));
ok("e1 e2 = e4 (first Fano triple)", OctMul(UnitVec8(2), UnitVec8(3)) = UnitVec8(5));

s    := [-1, 1, 1, 1, 1, 1, 1, 1] / 2;
sbar := [-1,-1,-1,-1,-1,-1,-1,-1] / 2;
SIG  := PermTwistMat((1,2));           # sigma = (1 2): swaps e1 and e2
sig  := x -> x * SIG;

##  Paper basis {L_k} (Appendix A) and M_k = sigma(L_k . s), N_k = sigma(L_k . sbar)
PL := [ s, [1,1,0,0,0,0,0,0], [1,-1,0,0,0,0,0,0], [1,0,1,0,0,0,0,0],
        [1,0,0,1,0,0,0,0], [1,0,0,0,1,0,0,0], [1,0,0,0,0,1,0,0], [1,0,0,0,0,0,1,0] ];
PM := List(PL, l -> sig(OctMul(l, s)));
PN := List(PL, l -> sig(OctMul(l, sbar)));
##  M_k and N_k as printed in the paper (transcribed from Appendix A)
PM_printed := [ [-3,-1,-1,-1,-1,-1,-1,-1]/2, [-1,0,0,0,1,0,1,1], [0,1,1,1,0,1,0,0],
                [-1,0,1,0,0,1,0,1], [-1,1,1,0,0,0,1,0], [-1,1,0,1,0,0,0,1],
                [-1,0,1,1,1,0,0,0], [-1,1,0,0,1,1,0,0] ];
PN_printed := [ [2,0,0,0,0,0,0,0], [0,0,-1,0,-1,0,-1,-1], [-1,-1,0,-1,0,-1,0,0],
                [0,-1,-1,0,0,-1,0,-1], [0,-1,-1,-1,0,0,-1,0], [0,-1,0,-1,-1,0,0,-1],
                [0,0,-1,-1,-1,-1,0,0], [0,-1,0,0,-1,-1,-1,0] ];
ok("M_k = sigma(L_k . s) agree with the paper's printed M_1..M_8", PM = PM_printed);
ok("N_k = sigma(L_k . sbar) agree with the paper's printed N_1..N_8", PN = PN_printed);

checkTable := function(tbl, A, Bm, C)     # rows [i,j,c1..c8]: A_i . B_j = sum c_k C_k
  return ForAll(tbl, row -> OctMul(A[row[1]], Bm[row[2]]) = row{[3..10]} * C);
end;
ok("paper Table LL (64 rows): L_i . L_j = sum c_k L_k", checkTable(PAPER_TABLE_LL, PL, PL, PL));
ok("paper Table LN (64 rows): L_i . N_j = sum c_k N_k", checkTable(PAPER_TABLE_LN, PL, PN, PN));
ok("paper Table MM (64 rows): M_i . M_j = sum c_k M_k", checkTable(PAPER_TABLE_MM, PM, PM, PM));

#############################################################################
Print("\nB. L = D8+\n");
#############################################################################
LAT := BuildLattices();;
BL := LAT.BL;;
Print("  Z-basis of L used below (rows; Hermite basis of the root span):\n");
for b in BL do Print("    ", b, "\n"); od;
ok("240 roots, all in L, all of norm 2",
   Length(LAT.roots) = 240 and ForAll(LAT.roots, InL) and ForAll(LAT.roots, r -> r * r = 2));
ok("|det(basis of L)| = 1 (covolume of D8+)", AbsInt(DeterminantMat(BL)) = 1);
ok("Gram matrix of L is even unimodular",
   DeterminantMat(BL * TransposedMat(BL)) = 1 and ForAll([1..8], i -> IsEvenInt(BL[i]*BL[i])));
ok("paper basis {L_k} is also a Z-basis of L",
   ForAll(PL, InL) and AbsInt(DeterminantMat(PL)) = 1);
ok("L . L in L (all 64 basis products in L)",
   ForAll(BL, a -> ForAll(BL, b -> InL(OctMul(a, b)))));
ok("e0 not in L, s in L, sbar not in L", not InL(UnitVec8(1)) and InL(s) and not InL(sbar));

#############################################################################
Print("\nC. Wilson's sublattices L sbar and L s (left multiples)\n");
#############################################################################
ok("L sbar inside L and L s inside L",
   ForAll(LAT.BLsbar, InL) and ForAll(LAT.BLs, InL));
Print("  |det| of bases: L sbar ", AbsInt(DeterminantMat(LAT.BLsbar)),
      ", L s ", AbsInt(DeterminantMat(LAT.BLs)), "  (= index in L)\n");
ok("[L : L sbar] = [L : L s] = 16",
   AbsInt(DeterminantMat(LAT.BLsbar)) = 16 and AbsInt(DeterminantMat(LAT.BLs)) = 16);
ok("2L inside L sbar and inside L s",
   ForAll(2 * BL, x -> InLatticeInv(x, LAT.BLsbarInv) and InLatticeInv(x, LAT.BLsInv)));
SumBasis := BaseIntMat(2 * Concatenation(LAT.BLsbar, LAT.BLs)) / 2;;
ok("L sbar + L s = L (so L sbar cap L s has index 16*16/1 = 256 = [L:2L], i.e. equals 2L)",
   AbsInt(DeterminantMat(SumBasis)) = 1);
u1 := [-1,0,0,0,1,0,1,1];; u2 := [-1,1,0,0,0,1,0,1];;
v  := [-1,0,1,0,0,1,0,1];; w  := [0,0,-1,0,-1,0,-1,-1];;
ok("paper Remark 4.5: u1, u2 in L s", InLatticeInv(u1, LAT.BLsInv) and InLatticeInv(u2, LAT.BLsInv));
ok("paper Remark 4.5: u1.u2 = (0,-2,2,1,-2,-1,-1,-1), in L sbar but not in L s",
   OctMul(u1, u2) = [0,-2,2,1,-2,-1,-1,-1] and InLatticeInv(OctMul(u1,u2), LAT.BLsbarInv)
   and not InLatticeInv(OctMul(u1,u2), LAT.BLsInv));
ok("paper Remark 4.5: v in sigma(L s) but not in L s",
   InLatticeInv(v, Inverse(PM)) and not InLatticeInv(v, LAT.BLsInv));
ok("paper Remark 4.5: w in sigma(L sbar) but not in L sbar",
   InLatticeInv(w, Inverse(PN)) and not InLatticeInv(w, LAT.BLsbarInv));

#############################################################################
Print("\nD. Lambda = {(x,y,z) in L^3 : x+y,x+z,y+z in L sbar; x+y+z in L s}\n");
#############################################################################
B := LAT.B;;
Print("  dim over GF(2) of Lambda/2L^3 inside L^3/2L^3: ", LAT.dimK, "\n");
ok("[L^3 : Lambda] = |det(coordinate matrix)| = 2^12", AbsInt(DeterminantMat(LAT.Clam)) = 2^12);
ok("each of the 24 basis vectors satisfies Wilson's definition", ForAll(B, u -> InLambdaDef(LAT, u)));
ok("2L^3 inside Lambda (24 generators 2*(basis of L in each block))",
   ForAll(2 * DirectSumMat(BL, BL, BL), u -> InLambdaDef(LAT, u)));
Gram := B * TransposedMat(B);;
Print("  Gram determinant = ", DeterminantMat(Gram), " = 2^", Log2Exact(DeterminantMat(Gram)), "\n");
ok("Gram determinant = 2^24", DeterminantMat(Gram) = 2^24);
ok("all inner products even, all norms divisible by 4 (sqrt(2) x even unimodular)",
   ForAll(Flat(Gram), IsEvenInt) and ForAll([1..24], i -> Gram[i][i] mod 4 = 0));
##  definition-based membership agrees with basis-based membership on random tests
BL3 := DirectSumMat(BL, BL, BL);;
rndA := List([1..200], i -> List([1..24], k -> Random([-3..3])) * B);;          # in Lambda
rndB := List([1..200], i -> List([1..24], k -> Random([-3..3])) * BL3);;        # in L^3
rndC := List([1..200], i -> List([1..24], k -> Random([-3..3])) * B + BL3[Random([1..24])]);;  # near misses
rnd := Concatenation(rndA, rndB, rndC);;
ok("on 600 test vectors: in-Lambda by Wilson's definition <=> integral coordinates in the basis",
   ForAll(rnd, u -> InLambdaDef(LAT, u) = ForAll(u * LAT.Binv, IsInt)));
Print("  of these, in Lambda: random Lambda-combinations ", Number(rndA, u -> InLambdaDef(LAT, u)),
      "/200, random L^3 vectors ", Number(rndB, u -> InLambdaDef(LAT, u)),
      "/200, Lambda + one basis vector of L^3 ", Number(rndC, u -> InLambdaDef(LAT, u)), "/200\n");
Print("  Z-basis of Lambda (rows, ambient coordinates, blocks of 8):\n");
for b in B do Print("    ", b, "\n"); od;

#############################################################################
Print("\nE. Minimum and minimal vectors (ShortestVectors on Gram/2)\n");
#############################################################################
U := Gram / 2;;
t0 := Runtime();;
sv3 := ShortestVectors(U, 3);;
Print("  vectors with Gram/2-norm <= 3: ", Length(sv3.vectors), "  (", Runtime() - t0, " ms)\n");
t0 := Runtime();;
sv4 := ShortestVectors(U, 4);;
Print("  vectors with Gram/2-norm <= 4 (one per +- pair): ", Length(sv4.vectors),
      ", norms ", Collected(sv4.norms), "  (", Runtime() - t0, " ms)\n");
ok("minimum norm of Lambda is 8 (no nonzero vector of norm <= 6)", Length(sv3.vectors) = 0);
ok("exactly 2 * 98,280 = 196,560 vectors of norm 8",
   Length(sv4.vectors) = 98280 and ForAll(sv4.norms, n -> n = 4));
Print("  => Gram/2 is even unimodular of rank 24 without roots: the Leech lattice.\n");

#############################################################################
Print("\nF. The paper's three families of minimal vectors (Sec. 2.3)\n");
#############################################################################
J := Concatenation(List([1..8], t -> UnitVec8(t)), List([1..8], t -> -UnitVec8(t)));;
cyc := function(x, y, z) return [Concatenation(x,y,z), Concatenation(z,x,y), Concatenation(y,z,x)]; end;
zero8 := ListWithIdenticalEntries(8, 0);;
fam := [];;
for lam in LAT.roots do Append(fam, cyc(2*lam, zero8, zero8)); od;                       # type 1
n1 := Length(fam);;
for lam in LAT.roots do
  ls := OctMul(lam, sbar);
  for j in J do Append(fam, cyc(ls, OctMul(ls, j), zero8)); od;                          # type 2
od;
n2 := Length(fam) - n1;;
for lam in LAT.roots do
  lss := OctMul(lam, s);
  for j in J do
    lsj := OctMul(lss, j);  lj := OctMul(lam, j);
    for k in J do Append(fam, cyc(lsj, OctMul(lam, k), OctMul(lj, k))); od;              # type 3
  od;
od;
n3 := Length(fam) - n1 - n2;;
Print("  generated: type1 ", n1, ", type2 ", n2, ", type3 ", n3, ", total ", Length(fam), "\n");
ok("all generated vectors have norm 8", ForAll(fam, u -> u * u = 8));
ok("all generated vectors satisfy Wilson's definition of Lambda", ForAll(fam, u -> InLambdaDef(LAT, u)));
famset := Set(fam);;
Print("  distinct vectors among them: ", Length(famset), "\n");
minset := Set(Concatenation(List(sv4.vectors, c -> c * B), List(sv4.vectors, c -> -c * B)));;
ok("the 196,560 family vectors are exactly the minimal shell found in E", famset = minset);

Print("\n=== end run_01 ===\n");
QUIT;
