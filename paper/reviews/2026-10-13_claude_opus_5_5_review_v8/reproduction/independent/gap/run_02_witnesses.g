#############################################################################
##
##  run_02_witnesses.g                                             (GAP 4.12)
##
##  Detailed, cross-checked analysis of five individual twists:
##    id          the untwisted product (paper: does not close)
##    sigma       the transposition (1,2), i.e. e1 <-> e2       (claim C1)
##    (5,6,7)     pi: e5 -> e6 -> e7 -> e5                       (claim C2 witness)
##    (5,7,6)     its inverse pi^{-1}: e5 -> e7 -> e6 -> e5      (convention check)
##    s_r         reflection in r = (e0+e1+e2+e3+e4+e5+e6-e7)/2  (claim C3 witness)
##
##  For each twist g (matrix G, row convention g(x) = x*G) the script reports:
##   1. G preserves L, G is orthogonal; identity element of ._g.
##   2. Closure from the 576 basis products (integrality of Lambda-coordinates).
##      If it fails: a basis pair whose product leaves Lambda and the violated
##      Wilson condition(s).
##   3. Independent code path: the 576 products and 300 random pairs of
##      Lambda-vectors recomputed with the naive formula P1,P2,P3 of
##      Def. 3.3 (TripleProdDirect) and tested with Wilson's definition
##      (InLambdaDef), with no use of the Lambda basis or of TripleLeftMat.
##   4. Span S: elementary divisors of the 576x24 coordinate matrix, index,
##      whether 2*Lambda <= S, and the basis-free check
##      det Gram(Z-basis of S) = 2^24 * [Lambda:S]^2.
##   5. Second construction: Lambda_g := {(x,y,z) in L^3 : x+y,.. in g(L sbar),
##      x+y+z in g(L s)} built from scratch (different Z-basis), with the
##      UNTWISTED product.  (Lambda, ._g) ~ (Lambda_g, .) via g^3, so the
##      closure verdict and the elementary divisors must agree.
##   6. Mod-2 reading: V_g = g(L sbar)/2L, W_g = g(L s)/2L in L/2L with the
##      product inherited from L.L in L; subalgebra / ideal properties and
##      dim span(V.V), dim span(W.W).  Prediction (proved in the report):
##        dim (S + 2L^3)/2L^3 = dim span(W.W) + 2 dim span(V.V),
##      hence [Lambda:S] = 2^(12 - that) * [2L^3 : S cap 2L^3].
##   7. For index 1: two sets A, B of 24 basis products with |det| = 2^a and
##      |det| = m odd (coordinate matrices in the Lambda basis); then
##      2^a*Lambda <= S and m*Lambda <= S, so S = Lambda.  Written in exact
##      integer form (doubled coordinates) to cert_index1_reflection.txt.
##
##  Conventions: see leech_lib.g.  Permutation convention: GAP cycle (5,6,7)
##  means pi(5)=6, pi(6)=7, pi(7)=5, i.e. pi(e_a) = e_{pi(a)}, and the twisted
##  product is x ._pi y = pi^{-1}(pi(x).pi(y)) as in the paper (Sec. 5.6).
##
##  Run:  gap -q -o 4g run_02_witnesses.g < /dev/null > run_02_witnesses.out
##
#############################################################################

Read("leech_lib.g");
Reset(GlobalMersenneTwister, 20261013);

LAT := BuildLattices();;
M2  := Mod2Setup(LAT);;
B := LAT.B;;  Binv := LAT.Binv;;
GramDetLambda := DeterminantMat(B * TransposedMat(B));;   # = 2^24

ok := function(name, cond)
  if cond = true then Print("    [ok]   ", name, "\n");
  else Print("    [FAIL] ", name, "\n"); fi;
end;

##  Which of Wilson's conditions does a vector violate?
WilsonViolations := function(u)
  local x, y, z, out;
  x := u{[1..8]}; y := u{[9..16]}; z := u{[17..24]};
  out := [];
  if not (InL(x) and InL(y) and InL(z)) then Add(out, "(1) some block not in L"); fi;
  if not InLatticeInv(x + y, LAT.BLsbarInv) then Add(out, "(2) x+y not in L sbar"); fi;
  if not InLatticeInv(x + z, LAT.BLsbarInv) then Add(out, "(2) x+z not in L sbar"); fi;
  if not InLatticeInv(y + z, LAT.BLsbarInv) then Add(out, "(2) y+z not in L sbar"); fi;
  if not InLatticeInv(x + y + z, LAT.BLsInv) then Add(out, "(3) x+y+z not in L s"); fi;
  return out;
end;

PairOfRow := idx -> [QuoInt(idx - 1, 24) + 1, RemInt(idx - 1, 24) + 1];

RESULTS := rec();;

AnalyzeTwist := function(label, G)
  local Gi, Tg, sd, idx, pr, a, c, prod, viol, direct576, rndpairs, rndok,
        Sb, gdet, k, BV, BW, W2, sd2, U, tl, pv, pw, dimS2, pred, e, x, y, ident;
  Print("\n==================================================================\n");
  Print("Twist ", label, "\n");
  Print("==================================================================\n");
  Print("  G (row convention, g(x) = x*G):\n");
  for x in G do Print("    ", x, "\n"); od;
  Gi := Inverse(G);
  ## 1. basic properties
  ok("g(L) = L  (G maps a basis of L into L and |det G| = 1)",
     ForAll(LAT.BL, b -> InL(b * G)) and AbsInt(DeterminantMat(G)) = 1);
  ok("G is orthogonal (G G^T = I), so ._g is norm-multiplicative", G * TransposedMat(G) = IdentityMat(8));
  ident := ShallowCopy(IdentityMat(8)[1]) * Gi;          # g^{-1}(e0)
  Print("    identity element of ._g is g^{-1}(e0) = ", ident, "\n");
  e := List([1..20], i -> List([1..8], k -> Random([-4..4]) / 2));
  ok("g^{-1}(e0) is a two-sided identity of ._g (20 random tests)",
     ForAll(e, x -> TwMulDirect(G, Gi, ident, x) = x and TwMulDirect(G, Gi, x, ident) = x));
  ## 2. closure from the basis products
  Tg := TwistedTensor(G);
  sd := SpanData(Tg, B, Binv);
  Print("  CLOSURE (576 basis products have integral Lambda-coordinates): ", sd.closes, "\n");
  if not sd.closes then
    idx := First([1..576], i -> not ForAll(sd.rows[i], IsInt));
    pr := PairOfRow(idx);
    prod := TripleProdDirect(G, B[pr[1]], B[pr[2]]);
    viol := WilsonViolations(prod);
    Print("    witness: b_", pr[1], " * b_", pr[2], " = ", prod, "\n");
    Print("    violated: ", viol, "\n");
    Print("    number of the 576 basis products outside Lambda: ",
          Number(sd.rows, r -> not ForAll(r, IsInt)), "\n");
  fi;
  ## 3. independent code path (naive formula + Wilson's definition)
  direct576 := List([1..576], i -> TripleProdDirect(G, B[PairOfRow(i)[1]], B[PairOfRow(i)[2]]));
  ok("naive P1,P2,P3 formula reproduces all 576 products (vs. coordinates * basis)",
     ForAll([1..576], i -> direct576[i] = sd.rows[i] * B));
  ok("naive products in Lambda by Wilson's definition <=> closure verdict",
     ForAll(direct576, p -> InLambdaDef(LAT, p)) = sd.closes);
  rndpairs := List([1..300], i -> [List([1..24], k -> Random([-2..2])) * B,
                                   List([1..24], k -> Random([-2..2])) * B]);
  rndok := Number(rndpairs, p -> InLambdaDef(LAT, TripleProdDirect(G, p[1], p[2])));
  Print("    300 random pairs of Lambda-vectors: product in Lambda for ", rndok, " of 300\n");
  ## 5. second construction (untwisted product on Lambda_g)
  BV := List(LAT.BLsbar, v -> v * G);
  BW := List(LAT.BLs, v -> v * G);
  W2 := BuildWilsonTypeBasis(LAT.BL, BV, BW);
  ok("Lambda_g = g^3(Lambda) as lattices (bases differ by a unimodular matrix)",
     AbsInt(DeterminantMat((B * DirectSumMat(G, G, G)) * W2.Binv)) = 1
     and ForAll(Flat((B * DirectSumMat(G, G, G)) * W2.Binv), IsInt));
  sd2 := SpanData(TwistedTensor(IdentityMat(8)), W2.B, W2.Binv);
  ok("untwisted product on Lambda_g: same closure verdict", sd2.closes = sd.closes);
  if sd.closes then
    ok("untwisted product on Lambda_g: same elementary divisors", sd2.ed = sd.ed);
  fi;
  ## 6. mod-2 reading
  tl := TwistLagrangians(LAT, G);
  pv := SubspaceProps(LAT, M2, tl.V);
  pw := SubspaceProps(LAT, M2, tl.W);
  Print("  MOD-2 READING (L/2L with the product inherited from L.L in L):\n");
  Print("    V_g = g(L sbar)/2L: dim ", pv.dim, ", totally singular ", pv.totallySingular,
        ", subalgebra ", pv.subalgebra, ", left ideal ", pv.leftIdeal,
        ", right ideal ", pv.rightIdeal, ", dim span(V.V) = ", pv.dimSquare, "\n");
  Print("    W_g = g(L s)/2L:    dim ", pw.dim, ", totally singular ", pw.totallySingular,
        ", subalgebra ", pw.subalgebra, ", left ideal ", pw.leftIdeal,
        ", right ideal ", pw.rightIdeal, ", dim span(W.W) = ", pw.dimSquare, "\n");
  ok("criterion: closes <=> V_g and W_g both subalgebras", sd.closes = (pv.subalgebra and pw.subalgebra));
  if sd.closes then
    ## 4. span
    Print("  SPAN S = Z-span of all products:\n");
    Print("    elementary divisors of the 576x24 coordinate matrix: ", EDString(sd.ed), "\n");
    Print("    Lambda/S = ", JoinStringsWithSeparator(List(Filtered(sd.ed, d -> d > 1), d -> Concatenation("Z/", String(d))), " + "), "\n");
    Print("    [Lambda : S] = ", sd.index, " = 2^", Log2Exact(sd.index), "\n");
    ok("2*Lambda <= S (every elementary divisor divides 2)", ForAll(sd.ed, d -> 2 mod d = 0));
    Print("    dim_F2 S/2Lambda = ", Number(sd.ed, d -> d = 1), "\n");
    gdet := SpanGramDet(LAT, sd.rows);
    Print("    det Gram(Z-basis of S) = 2^", Log2Exact(gdet), "  (= 2^24 * index^2 = 2^", 24 + 2*Log2Exact(sd.index), ")\n");
    ok("basis-free index check det Gram(S) = det Gram(Lambda) * index^2", gdet = GramDetLambda * sd.index^2);
    dimS2 := DimSmod2L3(LAT, sd.rows);
    pred := pw.dimSquare + 2 * pv.dimSquare;
    Print("    dim (S + 2L^3)/2L^3 = ", dimS2, ";  predicted dim span(W.W) + 2 dim span(V.V) = ", pred, "\n");
    ok("prediction holds", dimS2 = pred);
    Print("    index factorization: [Lambda : S + 2L^3] = 2^", 12 - dimS2,
          ",  [2L^3 : S cap 2L^3] = 2^", Log2Exact(sd.index) - (12 - dimS2), "\n");
    if pv.dimSquare = 0 then
      ##  V.V = 0 forces a = b = 0 in the mod-2 coordinates: x = y = z mod 2L
      ok("V.V = 0, so every product has its three blocks congruent mod 2L (all 576 basis products)",
         ForAll(direct576, p -> InLatticeInv(p{[1..8]} - p{[9..16]}, Inverse(2 * LAT.BL))
                            and InLatticeInv(p{[1..8]} - p{[17..24]}, Inverse(2 * LAT.BL))));
    fi;
  fi;
  return rec(G := G, sd := sd, pv := pv, pw := pw);
end;

RESULTS.id    := AnalyzeTwist("id (untwisted)", IdentityMat(8));;
RESULTS.sigma := AnalyzeTwist("sigma = (1,2): e1 <-> e2", PermTwistMat((1,2)));;
RESULTS.c567  := AnalyzeTwist("pi = GAP (5,6,7): e5 -> e6 -> e7 -> e5", PermTwistMat((5,6,7)));;
RESULTS.c576  := AnalyzeTwist("pi = GAP (5,7,6): e5 -> e7 -> e6 -> e5 (inverse of the above)", PermTwistMat((5,7,6)));;
RSTAR := [1,1,1,1,1,1,1,-1] / 2;;
RESULTS.refl  := AnalyzeTwist("s_r, r = (e0+e1+e2+e3+e4+e5+e6-e7)/2", RootReflectionMat(RSTAR));;

#############################################################################
##  7. Certificate that S = Lambda for the reflection witness
#############################################################################
Print("\n==================================================================\n");
Print("Certificate for [Lambda : S] = 1 (reflection witness)\n");
Print("==================================================================\n");
sd := RESULTS.refl.sd;;
G := RESULTS.refl.G;;
cert := CoprimeDetCertificate(sd.rows, 200);;
if cert = fail then
  Print("  no certificate found\n");
else
  pairsA := List(cert.selA, PairOfRow);;
  pairsB := List(cert.selB, PairOfRow);;
  prodsA := List(pairsA, p -> TripleProdDirect(G, B[p[1]], B[p[2]]));;
  prodsB := List(pairsB, p -> TripleProdDirect(G, B[p[1]], B[p[2]]));;
  Print("  subset A: basis pairs (a,c) of products b_a * b_c:\n    ", pairsA, "\n");
  Print("    det of the Lambda-coordinate matrix of A: ", cert.detA, " = 2^", Log2Exact(AbsInt(cert.detA)), "\n");
  Print("  subset B: basis pairs (a,c):\n    ", pairsB, "\n");
  Print("    det of the Lambda-coordinate matrix of B: ", cert.detB, " = ", Factors(AbsInt(cert.detB)), " (odd)\n");
  ok("all 48 products lie in Lambda (Wilson's definition, naive product formula)",
     ForAll(Concatenation(prodsA, prodsB), p -> InLambdaDef(LAT, p)));
  gA := DeterminantMat(prodsA * TransposedMat(prodsA));;
  gB := DeterminantMat(prodsB * TransposedMat(prodsB));;
  ok("basis-free: det Gram(A) = 2^24 * detA^2 and det Gram(B) = 2^24 * detB^2",
     gA = 2^24 * cert.detA^2 and gB = 2^24 * cert.detB^2);
  ok("gcd(detA, detB) = 1, hence 2^a*Lambda + m*Lambda = Lambda <= S, i.e. S = Lambda",
     GcdInt(cert.detA, cert.detB) = 1);
  ## write the exact certificate (doubled coordinates: all integers)
  f := "cert_index1_reflection.txt";;
  PrintTo(f, "# Certificate: for the twist g = s_r, r = (e0+e1+e2+e3+e4+e5+e6-e7)/2, the Z-span S of all\n",
    "# products u *_g v (u, v in Lambda) is all of Lambda.  Generated by run_02_witnesses.g (GAP ", GAPInfo.Version, ").\n#\n",
    "# Conventions: coordinates w.r.t. e0..e7 per block; R^24 = blocks (x,y,z).\n",
    "# ALL VECTORS ARE GIVEN IN DOUBLED COORDINATES (2*v), so every entry is an integer.\n",
    "# Octonion product: Fano triples (1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3);\n",
    "#   e_a e_b = e_c cyclically, anticyclic negative, e_i^2 = -e0, e0 = identity.\n",
    "# Twist: g(x) = x - <x,r> r (a reflection, g = g^{-1}); x ._g y = g^{-1}(g(x) . g(y)).\n",
    "# Triple product: P1 = x.x'+z.y'+y.z', P2 = y.y'+x.z'+z.x', P3 = z.z'+y.x'+x.y' (all ._g).\n",
    "# Lambda: (x,y,z) in L^3 with x+y, x+z, y+z in L sbar and x+y+z in L s; L = D8+ (integer\n",
    "#   coordinates with even sum, or half-integer coordinates with odd sum); L s = {l.s : l in L},\n",
    "#   L sbar = {l.sbar : l in L}, s = (-e0+e1+...+e7)/2, sbar = (-e0-e1-...-e7)/2.\n#\n",
    "# How to check:\n",
    "#  (i)   each BASIS row lies in Lambda and det Gram(BASIS) = 2^24 = covolume(Lambda)^2, so BASIS\n",
    "#        is a Z-basis of Lambda;\n",
    "#  (ii)  closure: all 576 products BASIS[a] *_g BASIS[c] lie in Lambda (so S <= Lambda);\n",
    "#  (iii) the 24 products listed under A lie in Lambda and det Gram(A) = 2^24 * DET_A^2, so\n",
    "#        span(A) has index |DET_A| = 2^", Log2Exact(AbsInt(cert.detA)), " in Lambda and 2^",
    Log2Exact(AbsInt(cert.detA)), " * Lambda <= S;\n",
    "#  (iv)  likewise for B with |DET_B| odd, so DET_B * Lambda <= S;\n",
    "#  (v)   gcd(DET_A, DET_B) = 1, so Lambda <= S.  With (ii): S = Lambda.\n#\n");
  AppendTo(f, "r2 = ", 2 * RSTAR, "\n\n");
  AppendTo(f, "BASIS2 = [  # 24 rows, doubled coordinates\n");
  for k in [1..24] do AppendTo(f, "  ", 2 * B[k], ",\n"); od;
  AppendTo(f, "]\n\n");
  AppendTo(f, "PAIRS_A = ", pairsA, "  # 1-based (a,c): product BASIS[a] *_g BASIS[c]\n");
  AppendTo(f, "DET_A = ", cert.detA, "  # det of the Lambda-coordinate matrix of the A products (= +-2^", Log2Exact(AbsInt(cert.detA)), ")\n\n");
  AppendTo(f, "PRODUCTS_A2 = [  # doubled coordinates of the A products\n");
  for k in [1..24] do AppendTo(f, "  ", 2 * prodsA[k], ",\n"); od;
  AppendTo(f, "]\n\n");
  AppendTo(f, "PAIRS_B = ", pairsB, "\n");
  AppendTo(f, "DET_B = ", cert.detB, "  # odd: ", Factors(AbsInt(cert.detB)), "\n\n");
  AppendTo(f, "PRODUCTS_B2 = [  # doubled coordinates of the B products\n");
  for k in [1..24] do AppendTo(f, "  ", 2 * prodsB[k], ",\n"); od;
  AppendTo(f, "]\n");
  Print("  certificate written to ", f, "\n");
  ## the same data in PARI/GP syntax, for verify_cert_index1.gp
  fg := "cert_index1_reflection.gp";;
  gpvec := v -> Concatenation("[", JoinStringsWithSeparator(List(v, String), ","), "]");;
  gpmat := m -> Concatenation("[", JoinStringsWithSeparator(List(m, gpvec), ", \\\n  "), "]");;
  PrintTo(fg, "\\\\ certificate data written by run_02_witnesses.g; doubled coordinates (2*v)\n");
  AppendTo(fg, "r2 = ", gpvec(2 * RSTAR), ";\n");
  AppendTo(fg, "BASIS2 = ", gpmat(2 * B), ";\n");
  AppendTo(fg, "PAIRS_A = ", gpmat(pairsA), ";\n");
  AppendTo(fg, "DET_A = ", cert.detA, ";\n");
  AppendTo(fg, "PRODUCTS_A2 = ", gpmat(2 * prodsA), ";\n");
  AppendTo(fg, "PAIRS_B = ", gpmat(pairsB), ";\n");
  AppendTo(fg, "DET_B = ", cert.detB, ";\n");
  AppendTo(fg, "PRODUCTS_B2 = ", gpmat(2 * prodsB), ";\n");
  Print("  PARI/GP copy written to ", fg, "\n");
fi;

Print("\n=== end run_02 ===\n");
QUIT;
