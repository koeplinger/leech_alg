#############################################################################
##
##  run_04_reflections.g                                           (GAP 4.12)
##
##  All 120 root reflections of L = D8+ as twists (claim C3):
##    s_r(x) = x - <x,r> r, r a root (r.r = 2), one r per pair +-r
##    (representative: first nonzero coordinate positive).
##  s_r is an isometry with s_r(L) = L and s_r = s_r^{-1}, so the twisted
##  product x ._g y = g^{-1}(g(x).g(y)) does not depend on the g / g^{-1}
##  convention.  It is an octonion product (isomorphic to . via g) whose
##  identity element is s_r(e0); s_r fixes e0 iff <r,e0> = 0.
##
##  For each reflection: closure on Lambda (576 basis products), elementary
##  divisors of the span S, the mod-2 data of V = s_r(L sbar)/2L and
##  W = s_r(L s)/2L, and the decomposition of the class of r in
##  L/2L = V0 + W0 (V0 = L sbar/2L, W0 = L s/2L).
##
##  Output: run_04_reflections.out (summary) and refl_census_table.txt.
##
##  Run:  gap -q -o 4g run_04_reflections.g < /dev/null > run_04_reflections.out
##
#############################################################################

Read("leech_lib.g");

LAT := BuildLattices();;
M2  := Mod2Setup(LAT);;
B := LAT.B;;  Binv := LAT.Binv;;

Print("=== run_04_reflections.g ===\n");

reps := Filtered(LAT.roots, r -> First(r, c -> c <> 0) > 0);;
Print("roots: ", Length(LAT.roots), ", reflections (one root per +- pair): ", Length(reps), "\n");

RootKind := function(r)
  if not ForAll(r, IsInt) then return "half-integer"; fi;
  if r[1] = 0 then return "integer, orthogonal to e0"; fi;
  return "integer, involves e0";
end;

##  bases of V0 = L sbar/2L and W0 = L s/2L, to split classes as v + w
V0 := SpanMod2(List(LAT.BLsbar, x -> ClassMod2(LAT, x)));;
W0 := SpanMod2(List(LAT.BLs, x -> ClassMod2(LAT, x)));;
BV0 := BasisVectors(Basis(V0));;  BW0 := BasisVectors(Basis(W0));;
SplitVW := function(c)          # c = v + w, v in V0, w in W0
  local sol;
  sol := SolutionMat(Concatenation(BV0, BW0), c);
  return [sol{[1..4]} * BV0, sol{[5..8]} * BW0];
end;

t0 := Runtime();;
DATA := [];;
G := fail;;          # (declared before the loop: avoids a parser warning)
for r in reps do
  G := RootReflectionMat(r);
  if not (ForAll(LAT.BL, b -> InL(b * G)) and G * G = IdentityMat(8)) then Error("bad reflection"); fi;
  sd := SpanData(TwistedTensor(G), B, Binv);
  tl := TwistLagrangians(LAT, G);
  pv := SubspaceProps(LAT, M2, tl.V);
  pw := SubspaceProps(LAT, M2, tl.W);
  vw := SplitVW(ClassMod2(LAT, r));
  d := rec( r := r, kind := RootKind(r), closes := sd.closes, fixesE0 := (r[1] = 0),
            unit := ShallowCopy(IdentityMat(8)[1]) * G,
            Vsub := pv.subalgebra, Wsub := pw.subalgebra, dimVV := pv.dimSquare, dimWW := pw.dimSquare,
            Videal := [pv.leftIdeal, pv.rightIdeal], Wideal := [pw.leftIdeal, pw.rightIdeal],
            v := vw[1], w := vw[2],
            rInLsbar := InLatticeInv(r, LAT.BLsbarInv), rInLs := InLatticeInv(r, LAT.BLsInv) );
  if sd.closes then
    d.ed := sd.ed; d.index := sd.index; d.dimS2L3 := DimSmod2L3(LAT, sd.rows);
    d.twoLambdaInS := ForAll(sd.ed, e -> 2 mod e = 0);
  fi;
  Add(DATA, d);
od;
Print("computation time: ", Runtime() - t0, " ms\n\n");

VerdictStr := function(d)        # elementary divisors, or "-" if the twist does not close
  if d.closes then return EDString(d.ed); fi;
  return "-";
end;

Print("Claim C3, counts by verdict (closes, elementary divisors of S):\n");
for c in Collected(List(DATA, d -> [d.closes, VerdictStr(d)])) do
  Print("  ", c[1], " : ", c[2], "\n");
od;
Print("\nby root kind:\n");
for k in ["integer, orthogonal to e0", "integer, involves e0", "half-integer"] do
  sel := Filtered(DATA, d -> d.kind = k);
  Print("  ", k, " (", Length(sel), " reflections): ",
        Collected(List(sel, d -> [d.closes, VerdictStr(d)])), "\n");
od;

## sigma = (1 2) is the reflection in e1 - e2
rs := [0, 1, -1, 0, 0, 0, 0, 0];;
Print("\nsigma = (1,2) equals the reflection in e1 - e2: ", RootReflectionMat(rs) = PermTwistMat((1,2)), "\n");
dsig := First(DATA, d -> d.r = rs);;
Print("  its verdict: closes ", dsig.closes, ", ", EDString(dsig.ed), "\n");
RSTAR := [1,1,1,1,1,1,1,-1] / 2;;
dstar := First(DATA, d -> d.r = RSTAR);;
Print("witness r = (e0+e1+...+e6-e7)/2: closes ", dstar.closes, ", ", EDString(dstar.ed), "\n\n");

Print("Index-1 reflections (S = Lambda), roots r:\n");
for d in Filtered(DATA, d -> d.closes and d.index = 1) do Print("  ", d.r, "   (", d.kind, ")\n"); od;
Print("\nIndex-2^16 reflections, roots r:\n");
for d in Filtered(DATA, d -> d.closes and d.index = 2^16) do Print("  ", d.r, "   (", d.kind, ")\n"); od;
others := Filtered(DATA, d -> d.closes and not d.index in [1, 2^16]);;
Print("\nclosing reflections with other indices: ", Length(others), "\n");

Print("\nAll closing reflections have 2*Lambda <= S: ", ForAll(Filtered(DATA, d -> d.closes), d -> d.twoLambdaInS), "\n");
Print("Mod-2 criterion (C5) closes <=> V and W subalgebras: ", ForAll(DATA, d -> d.closes = (d.Vsub and d.Wsub)), "\n");
Print("Mod-2 prediction dim (S+2L^3)/2L^3 = dim W.W + 2 dim V.V: ",
      ForAll(Filtered(DATA, d -> d.closes), d -> d.dimS2L3 = d.dimWW + 2 * d.dimVV), "\n");
Print("Closing reflections by (dim V.V, dim W.W, V left ideal, V right ideal, W left, W right, elementary divisors):\n");
for c in Collected(List(Filtered(DATA, d -> d.closes),
          d -> [d.dimVV, d.dimWW, d.Videal[1], d.Videal[2], d.Wideal[1], d.Wideal[2], EDString(d.ed)])) do
  Print("  ", c[1], " : ", c[2], "\n");
od;
Print("Non-closing reflections by (V subalgebra, W subalgebra): ",
      Collected(List(Filtered(DATA, d -> not d.closes), d -> [d.Vsub, d.Wsub])), "\n");
Print("\nclass of r split as v + w (v in L sbar/2L, w in L s/2L): every r has v <> 0 and w <> 0: ",
      ForAll(DATA, d -> not IsZero(d.v) and not IsZero(d.w)), "\n");
Print("  number of distinct v: ", Length(Set(List(DATA, d -> d.v))), ", distinct (v,w): ", Length(Set(List(DATA, d -> [d.v, d.w]))), "\n");
Print("  verdict by v (15 nonzero v, 8 reflections each):\n");
for v in Set(List(DATA, d -> d.v)) do
  Print("    v = ", List(v, IntFFE), " : ",
        Collected(List(Filtered(DATA, d -> d.v = v), d -> VerdictStr(d))), "\n");
od;
Print("  verdict by w:\n");
for w in Set(List(DATA, d -> d.w)) do
  Print("    w = ", List(w, IntFFE), " : ",
        Collected(List(Filtered(DATA, d -> d.w = w), d -> VerdictStr(d))), "\n");
od;

f := "refl_census_table.txt";;
PrintTo(f, "# All 120 root reflections s_r(x) = x - <x,r> r of L = D8+ as twists x ._g y = g^{-1}(g(x).g(y)).\n",
  "# columns: r | kind | closes | elementary divisors of S | identity element s_r(e0) of ._g |\n",
  "#          dim span(V.V), dim span(W.W) | V left/right ideal | class of r = v + w (v in Lsbar/2L, w in Ls/2L; coordinates in the basis of L)\n");
for d in DATA do
  AppendTo(f, d.r, " | ", d.kind, " | ", d.closes, " | ");
  if d.closes then AppendTo(f, EDString(d.ed)); else AppendTo(f, "-"); fi;
  AppendTo(f, " | ", d.unit, " | ", d.dimVV, ",", d.dimWW, " | ", d.Videal, " | ", List(d.v, IntFFE), " + ", List(d.w, IntFFE), "\n");
od;
Print("\nper-reflection table written to ", f, "\n=== end run_04 ===\n");
QUIT;
