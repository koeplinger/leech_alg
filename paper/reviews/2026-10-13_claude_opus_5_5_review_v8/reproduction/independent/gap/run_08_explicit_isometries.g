#############################################################################
##
##  run_08_explicit_isometries.g                                   (GAP 4.12)
##
##  Explicit twists for every span-index class of the polarization census
##  (run_06): for each class, an isometry g of L, written as a short product
##  of root reflections, such that the twisted product *_g on WILSON'S Lambda
##  closes with that span index.  The verdict is computed directly (576 basis
##  products, elementary divisors), independently of the census machinery.
##
##  Method: breadth-first search over the ordered complementary pairs
##  (V, W) = (V0 g, W0 g) reachable from (V0, W0) = (L sbar/2L, L s/2L) by the
##  mod-2 images of the 120 root reflections (right action on canonical
##  bases).  Each pair is classified by the mod-2 invariants
##  (dim V.V, delta(V), dim W.W, delta(W)) of run_06; the first pair found in
##  each class gives a word k_1 ... k_m, and g has matrix
##  G = G_{k_1} G_{k_2} ... G_{k_m} (row convention: s_{k_1} is applied first).
##
##  Run:  gap -q -o 4g run_08_explicit_isometries.g < /dev/null > run_08_explicit_isometries.out
##
#############################################################################

Read("leech_lib.g");

LAT := BuildLattices();;
M2  := Mod2Setup(LAT);;
B := LAT.B;;  Binv := LAT.Binv;;
BL := LAT.BL;;

Print("=== run_08_explicit_isometries.g ===\n");

reps := Filtered(LAT.roots, r -> First(r, c -> c <> 0) > 0);;
GR := List(reps, RootReflectionMat);;
Mod2Mat := G -> (BL * G * LAT.BLinv) * Z(2);;
GR2 := List(GR, Mod2Mat);;
canonM := M -> OnSubspacesByCanonicalBasis(M, One(GR2[1]));;
V0 := canonM(List(LAT.BLsbar, x -> ClassMod2(LAT, x)));;
W0 := canonM(List(LAT.BLs, x -> ClassMod2(LAT, x)));;
act := function(pr, g) return [OnSubspacesByCanonicalBasis(pr[1], g), OnSubspacesByCanonicalBasis(pr[2], g)]; end;;

ClassOf := function(pr)       # (dim V.V, delta(V), dim W.W, delta(W)) and the closure criterion
  local V, W, pv, pw;
  V := SpanMod2(pr[1]);  W := SpanMod2(pr[2]);
  pv := SubspaceProps(LAT, M2, V);  pw := SubspaceProps(LAT, M2, W);
  return [pv.subalgebra and pw.subalgebra, pv.dimSquare, IdealDefect(M2, V), pw.dimSquare, IdealDefect(M2, W)];
end;

TARGETS := [ [true, 0, 0, 4, 4], [true, 1, 2, 4, 4], [true, 4, 4, 0, 0], [true, 4, 4, 1, 2], [true, 4, 4, 4, 4] ];;
FOUND := [];;   # entries [class, word]

## breadth-first search with words
start := [V0, W0];;
KEYS := [start];;  WORDS := [[]];;  frontier := [1];;  SORTED := [start];;
c0 := ClassOf(start);;
if c0 in TARGETS then Add(FOUND, [c0, []]); fi;
depth := 0;;  G := fail;;  sd := fail;;  pr := fail;;  pos := fail;;  cl := fail;;
while Length(FOUND) < Length(TARGETS) and frontier <> [] do
  depth := depth + 1;
  newfrontier := [];
  for idx in frontier do
    for k in [1..Length(GR2)] do
      pr := act(KEYS[idx], GR2[k]);
      pos := PositionSorted(SORTED, pr);
      if pos > Length(SORTED) or SORTED[pos] <> pr then
        Add(SORTED, pr, pos);
        Add(KEYS, pr);  Add(WORDS, Concatenation(WORDS[idx], [k]));
        Add(newfrontier, Length(KEYS));
        cl := ClassOf(pr);
        if cl in TARGETS and not cl in List(FOUND, f -> f[1]) then
          Add(FOUND, [cl, WORDS[Length(WORDS)]]);
        fi;
      fi;
    od;
    if Length(FOUND) = Length(TARGETS) then break; fi;
  od;
  frontier := newfrontier;
  Print("  depth ", depth, ": ", Length(KEYS), " pairs reached, classes found ", Length(FOUND), "\n");
od;

Print("\nexplicit twists (g = product of reflections s_r, applied left to right), verified directly on Wilson's Lambda:\n");
for f in FOUND do
  G := IdentityMat(8);
  for k in f[2] do G := G * GR[k]; od;
  sd := SpanData(TwistedTensor(G), B, Binv);
  Print("  class (dim V.V, delta V, dim W.W, delta W) = ", f[1]{[2..5]}, ": roots ", List(f[2], k -> reps[k]), "\n");
  Print("    det g = ", DeterminantMat(G), ", g(L) = L: ", ForAll(BL, b -> InL(b * G)),
        ", closes: ", sd.closes);
  if sd.closes then
    Print(", elementary divisors ", EDString(sd.ed), ", predicted index 2^",
          24 - (f[1][4] + 2 * f[1][2] + 2 * f[1][3] + f[1][5]), "\n");
  else
    Print("\n");
  fi;
od;
Print("=== end run_08 ===\n");
QUIT;
