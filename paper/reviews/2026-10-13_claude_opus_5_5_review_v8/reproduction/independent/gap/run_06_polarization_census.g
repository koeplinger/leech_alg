#############################################################################
##
##  run_06_polarization_census.g                                   (GAP 4.12)
##
##  Claim C5 and a complete census over all polarizations of L/2L.
##
##  Lbar = L/2L = GF(2)^8 (coordinates in the basis of L, mod 2) carries
##    * the quadratic form q(c) = N(lift)/2 mod 2 (plus type, Witt index 4);
##    * the product mu-bar inherited from L.L in L (untwisted octonions).
##  A Lagrangian is a 4-dim totally singular subspace.  For an ORDERED pair
##  (V, W) of complementary Lagrangians let Vhat, What be their preimages in L
##  and define the Wilson-type lattice
##      Lambda(V,W) = {(x,y,z) in L^3 : x+y, x+z, y+z in Vhat ; x+y+z in What}.
##  (Wilson's Lambda is Lambda(V0,W0) with V0 = L sbar/2L, W0 = L s/2L; for a
##  twist g, (Lambda, ._g) is isomorphic via g^3 to (Lambda(gV0, gW0), .).)
##
##  The script
##   1. enumerates all Lagrangians (claim: 270), and counts subalgebras
##      (claim: 142) and two-sided ideals (claim: 8);
##   2. enumerates the ordered complementary pairs (claim: 17,280) and those
##      with both members subalgebras (claim: 4,704);
##   3. for EVERY ordered complementary pair builds a Z-basis of Lambda(V,W)
##      from scratch, checks Gram determinant 2^24 and evenness, tests
##      closure under the UNTWISTED triple product (576 basis products) and,
##      if closed, computes the elementary divisors of the span S;
##   4. checks "closed <=> V and W both subalgebras" on all 17,280 pairs;
##   5. spot-checks minimum 8 (ShortestVectors) on a few lattices of each class;
##   6. locates the pairs realized by sigma, (5,6,7) and the reflection witness;
##   7. tests the mod-2 index formula on every closed pair:
##        log2 [Lambda : S] = 24 - ( dim W.W + 2 dim V.V + 2 delta(V) + delta(W) ),
##      delta(U) = dim (U + U.Lbar + Lbar.U)/U  (0 iff U is a two-sided ideal).
##      Proved parts (see the report): dim (S+2L^3)/2L^3 = dim W.W + 2 dim V.V
##      exactly, and if 2*Lambda <= S then dim S/2Lambda >= that + 2 delta(V)
##      + delta(W); equality is the observation tested here;
##   8. checks that the Weyl group W(E8) = <root reflections> acts transitively
##      on the 17,280 ordered complementary pairs (orbit of (V0, W0) under the
##      mod-2 images of the 120 reflections), so every pair is (g V0, g W0) for
##      some isometry g of L, and (Lambda, *_g) = (Lambda(gV0, gW0), *) via g^3;
##   9. splits the census by the two families of Lagrangians (U and U' are in
##      the same family iff dim(U cap U') is even) and checks, on all 5040
##      permutations and 120 reflections, that g(L sbar)/2L lies in the family
##      of L sbar/2L exactly when det g = +1.
##
##  Run:  gap -q -o 6g run_06_polarization_census.g < /dev/null > run_06_polarization_census.out
##
#############################################################################

Read("leech_lib.g");
Reset(GlobalMersenneTwister, 20261013);

LAT := BuildLattices();;
M2  := Mod2Setup(LAT);;
BL := LAT.BL;;
T0 := TwistedTensor(IdentityMat(8));;        # untwisted structure tensor

Print("=== run_06_polarization_census.g ===\n");

ALL := List(Tuples([0, 1], 8), t -> t * Z(2));;
NONZERO := Filtered(ALL, v -> not IsZero(v));;
q := v -> QForm(LAT, v);;
bil := function(u, v) return (q(u + v) - q(u) - q(v)) mod 2; end;;
SING := Filtered(NONZERO, v -> q(v) = 0);;
Print("nonzero singular vectors of (L/2L, q): ", Length(SING), " (plus type: 135)\n");

## ------------------------------------------------------------ 1. Lagrangians
level := List(SING, v -> SpanMod2([v]));;
for k in [2..4] do
  keys := [];  new := [];
  for U in level do
    for v in SING do
      if not v in U and ForAll(BasisVectors(Basis(U)), b -> bil(b, v) = 0) then
        Wsp := SpanMod2(Concatenation(BasisVectors(Basis(U)), [v]));
        key := AsSSortedList(Wsp);
        pos := PositionSorted(keys, key);
        if pos > Length(keys) or keys[pos] <> key then
          Add(keys, key, pos);  Add(new, Wsp, pos);
        fi;
      fi;
    od;
  od;
  level := new;
  Print("totally singular subspaces of dimension ", k, ": ", Length(level), "\n");
od;
LAGR := level;;
LKEYS := List(LAGR, AsSSortedList);;      # sorted, aligned with LAGR
PROPS := List(LAGR, U -> SubspaceProps(LAT, M2, U));;
Print("\nLagrangians: ", Length(LAGR), "\n");
Print("  subalgebras: ", Number(PROPS, p -> p.subalgebra), "\n");
Print("  two-sided ideals: ", Number(PROPS, p -> p.leftIdeal and p.rightIdeal),
      ", left ideals: ", Number(PROPS, p -> p.leftIdeal), ", right ideals: ", Number(PROPS, p -> p.rightIdeal), "\n");
Print("  subalgebras by (left ideal, right ideal, dim U.U): ",
      Collected(List(Filtered(PROPS, p -> p.subalgebra), p -> [p.leftIdeal, p.rightIdeal, p.dimSquare])), "\n");

## ------------------------------------------------------------ 2. pairs
PAIRS := [];;
for i in [1..Length(LAGR)] do
  for j in [1..Length(LAGR)] do
    if RankMat(Concatenation(BasisVectors(Basis(LAGR[i])), BasisVectors(Basis(LAGR[j])))) = 8 then
      Add(PAIRS, [i, j]);
    fi;
  od;
od;
Print("\nordered complementary pairs: ", Length(PAIRS),
      "  (per Lagrangian: ", Collected(List([1..Length(LAGR)], i -> Number(PAIRS, p -> p[1] = i))), ")\n");
Print("  with both members subalgebras: ", Number(PAIRS, p -> PROPS[p[1]].subalgebra and PROPS[p[2]].subalgebra), "\n");

## ------------------------------------------------------------ 3. census
Preimage := function(U)        # Z-basis (ambient coordinates) of the preimage of U in L
  local lifts;
  lifts := List(BasisVectors(Basis(U)), v -> List(v, IntFFE));
  return BaseIntMat(Concatenation(lifts, 2 * IdentityMat(8))) * BL;
end;;
PRE := List(LAGR, Preimage);;

DEFECT := List(LAGR, U -> IdealDefect(M2, U));;
t0 := Runtime();;
CENSUS := [];;
Gm := fail;; W2 := fail;; sd := fail;; edstr := fail;;   # (declared before the loop: avoids parser warnings)
for pr in PAIRS do
  W2 := BuildWilsonTypeBasis(BL, PRE[pr[1]], PRE[pr[2]]);
  Gm := W2.B * TransposedMat(W2.B);
  if W2.dimK <> 12 or DeterminantMat(Gm) <> 2^24 or not ForAll(Flat(Gm), IsEvenInt)
     or not ForAll([1..24], i -> Gm[i][i] mod 4 = 0) then
    Error("pair ", pr, ": lattice is not sqrt(2) * even unimodular");
  fi;
  sd := SpanData(T0, W2.B, W2.Binv);
  if sd.closes then edstr := EDString(sd.ed); else edstr := "-"; fi;
  Add(CENSUS, rec(pair := pr, closes := sd.closes, ed := edstr,
                  crit := PROPS[pr[1]].subalgebra and PROPS[pr[2]].subalgebra));
  if sd.closes then
    CENSUS[Length(CENSUS)].edlist := sd.ed;
    CENSUS[Length(CENSUS)].dimS2L3 := DimSmod2L3(W2, sd.rows);
  fi;
od;
Print("\ncensus of all ", Length(PAIRS), " ordered pairs: ", Runtime() - t0, " ms\n");
Print("  every Lambda(V,W) has a Z-basis with Gram determinant 2^24, even, norms in 4Z: true\n");
Print("  closed pairs: ", Number(CENSUS, c -> c.closes), "\n");
Print("  criterion closed <=> (V and W subalgebras) on all pairs: ", ForAll(CENSUS, c -> c.closes = c.crit), "\n");
Print("  closed pairs by elementary divisors of S: ", Collected(List(Filtered(CENSUS, c -> c.closes), c -> c.ed)), "\n");
Print("  closed pairs by (V: left ideal, right ideal, dim V.V | W: left, right, dim W.W | elementary divisors):\n");
for c in Collected(List(Filtered(CENSUS, c -> c.closes), c ->
      [PROPS[c.pair[1]].leftIdeal, PROPS[c.pair[1]].rightIdeal, PROPS[c.pair[1]].dimSquare,
       PROPS[c.pair[2]].leftIdeal, PROPS[c.pair[2]].rightIdeal, PROPS[c.pair[2]].dimSquare, c.ed])) do
  Print("    ", c[1], " : ", c[2], "\n");
od;
Print("  non-closed pairs by (V subalgebra, W subalgebra): ",
      Collected(List(Filtered(CENSUS, c -> not c.closes), c -> [PROPS[c.pair[1]].subalgebra, PROPS[c.pair[2]].subalgebra])), "\n");

## ------------------------------------------------------------ 5. minimum spot checks
Print("\nspot check of the minimum (ShortestVectors on Gram/2, bound 3; 0 vectors = minimum 8):\n");
for e in Set(List(CENSUS, c -> c.ed)) do
  sel := Filtered(CENSUS, c -> c.ed = e);
  for c in sel{[1, Random([1..Length(sel)])]} do
    W2 := BuildWilsonTypeBasis(BL, PRE[c.pair[1]], PRE[c.pair[2]]);
    Print("  class ", String(e, -10), " pair ", c.pair, ": vectors of norm <= 6: ",
          Length(ShortestVectors(W2.B * TransposedMat(W2.B) / 2, 3).vectors), "\n");
  od;
od;

## ------------------------------------------------------------ 6. known twists
PairOfTwist := function(G)
  local tl;
  tl := TwistLagrangians(LAT, G);
  return [Position(LKEYS, AsSSortedList(tl.V)), Position(LKEYS, AsSSortedList(tl.W))];
end;;
Print("\nknown twists located in the census (pair = [index of V, index of W]):\n");
for tw in [ ["untwisted", IdentityMat(8)], ["sigma (1,2)", PermTwistMat((1,2))],
            ["(5,6,7)", PermTwistMat((5,6,7))], ["(5,7,6)", PermTwistMat((5,7,6))],
            ["reflection r*", RootReflectionMat([1,1,1,1,1,1,1,-1]/2)] ] do
  pr := PairOfTwist(tw[2]);
  c := First(CENSUS, c -> c.pair = pr);
  Print("  ", String(tw[1], -16), " pair ", pr, "  closes ", c.closes, "  ", c.ed, "\n");
od;
## pairs realized by the 5040 permutations and the 120 reflections
permPairs := Set(List(AsList(SymmetricGroup(7)), p -> PairOfTwist(PermTwistMat(p))));;
reflPairs := Set(List(Filtered(LAT.roots, r -> First(r, x -> x <> 0) > 0), r -> PairOfTwist(RootReflectionMat(r))));;
Print("  distinct pairs realized by the 5040 permutations: ", Length(permPairs),
      ", of which closed: ", Number(permPairs, p -> First(CENSUS, c -> c.pair = p).closes), "\n");
Print("  distinct pairs realized by the 120 reflections: ", Length(reflPairs),
      ", of which closed: ", Number(reflPairs, p -> First(CENSUS, c -> c.pair = p).closes), "\n");

## ------------------------------------------------------------ 7. index formula
Print("\nmod-2 index formula on the ", Number(CENSUS, c -> c.closes), " closed pairs:\n");
cl := Filtered(CENSUS, c -> c.closes);;
Print("  2*Lambda <= S for every closed pair: ", ForAll(cl, c -> ForAll(c.edlist, e -> 2 mod e = 0)), "\n");
Print("  dim (S+2L^3)/2L^3 = dim W.W + 2 dim V.V for every closed pair: ",
      ForAll(cl, c -> c.dimS2L3 = PROPS[c.pair[2]].dimSquare + 2 * PROPS[c.pair[1]].dimSquare), "\n");
pred := c -> 24 - (PROPS[c.pair[2]].dimSquare + 2 * PROPS[c.pair[1]].dimSquare
                   + 2 * DEFECT[c.pair[1]] + DEFECT[c.pair[2]]);;
Print("  log2 [Lambda:S] = 24 - (dim W.W + 2 dim V.V + 2 delta(V) + delta(W)) for every closed pair: ",
      ForAll(cl, c -> Log2Exact(Product(c.edlist)) = pred(c)), "\n");
Print("  closed pairs by (dim V.V, delta(V), dim W.W, delta(W), elementary divisors):\n");
for x in Collected(List(cl, c -> [PROPS[c.pair[1]].dimSquare, DEFECT[c.pair[1]],
                                  PROPS[c.pair[2]].dimSquare, DEFECT[c.pair[2]], c.ed])) do
  Print("    ", x[1], " : ", x[2], "\n");
od;
## ------------------------------------------------------------ 8. transitivity
## induced action of a linear map g with g(L) = L on L/2L (row convention, coordinates in BL)
Mod2Mat := G -> (BL * G * LAT.BLinv) * Z(2);;
gens := List(Filtered(LAT.roots, r -> First(r, x -> x <> 0) > 0), r -> Mod2Mat(RootReflectionMat(r)));;
canon := U -> OnSubspacesByCanonicalBasis(BasisVectors(Basis(U)), One(gens[1]));;
actPair := function(pr, g) return [OnSubspacesByCanonicalBasis(pr[1], g), OnSubspacesByCanonicalBasis(pr[2], g)]; end;;
V0 := SpanMod2(List(LAT.BLsbar, x -> ClassMod2(LAT, x)));;
W0 := SpanMod2(List(LAT.BLs, x -> ClassMod2(LAT, x)));;
orb := Orbit(Group(gens), [canon(V0), canon(W0)], actPair);;
allPairs := Set(List(PAIRS, p -> [canon(LAGR[p[1]]), canon(LAGR[p[2]])]));;
Print("\norbit of (L sbar/2L, L s/2L) under W(E8) (generated by the 120 reflections, acting on L/2L): ",
      Length(orb), " ordered pairs; equals the set of all ordered complementary pairs: ", Set(orb) = allPairs, "\n");
Print("order of the image of W(E8) in GL(8,2): ", Size(Group(gens)), "  (|O+(8,2)| = 348364800)\n");
## ------------------------------------------------------------ 9. families
FAMILY := List(LAGR, U -> Dimension(Intersection(U, V0)) mod 2);;   # 0: family of V0 = L sbar/2L
Print("\nLagrangians per family (0 = family of L sbar/2L): ", Collected(FAMILY), "\n");
Print("complementary pairs lie in one family: ", ForAll(PAIRS, p -> FAMILY[p[1]] = FAMILY[p[2]]), "\n");
Print("closed pairs by (family, elementary divisors): ",
      Collected(List(cl, c -> [FAMILY[c.pair[1]], c.ed])), "\n");
famOfTwist := G -> FAMILY[PairOfTwist(G)[1]];;
Print("family of g(L sbar)/2L is that of L sbar/2L iff det g = +1, on all 5040 permutations: ",
      ForAll(AsList(SymmetricGroup(7)), p -> (famOfTwist(PermTwistMat(p)) = 0) = (SignPerm(p) = 1)), "\n");
Print("  ... and on all 120 reflections (det = -1): ",
      ForAll(Filtered(LAT.roots, r -> First(r, x -> x <> 0) > 0), r -> famOfTwist(RootReflectionMat(r)) = 1), "\n");
Print("=== end run_06 ===\n");
QUIT;
