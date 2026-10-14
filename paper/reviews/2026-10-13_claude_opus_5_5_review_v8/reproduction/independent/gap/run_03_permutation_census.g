#############################################################################
##
##  run_03_permutation_census.g                                    (GAP 4.12)
##
##  All 5,040 permutations pi of the seven imaginary units (identity
##  included), twist x ._pi y = pi^{-1}(pi(x).pi(y)), triple product of
##  Def. 3.3.  For each: closure on Lambda (576 basis products) and, if it
##  closes, the elementary divisors of the span S.  Also the mod-2 data of
##  V = pi(L sbar)/2L and W = pi(L s)/2L (subalgebra flags, dim of squares).
##
##  Checks:
##   * per cycle type, the closing counts against the table of paper Sec. 5.6;
##   * claim C2: every closing EVEN permutation has elementary divisors
##     1^14 2^10 and every closing ODD one has 1^8 2^16;
##   * 3-cycles: for each 3-subset exactly one orientation closes (paper);
##   * the closure criterion "V and W both subalgebras of L/2L" (claim C5)
##     agrees with the computed closure on every permutation;
##   * the mod-2 prediction dim (S+2L^3)/2L^3 = dim W.W + 2 dim V.V.
##
##  Conventions (see leech_lib.g): a GAP permutation p stands for pi with
##  pi(e_a) = e_{a^p}; the per-permutation table lists pi as the image list
##  [pi(1),...,pi(7)].  Counts per cycle type do not depend on the pi versus
##  pi^{-1} convention (inversion preserves cycle type and parity), but the
##  verdict for an individual permutation does: the table also reports the
##  verdict for pi^{-1}.
##
##  Output: run_03_permutation_census.out (summary) and
##          perm_census_table.txt (one line per permutation).
##
##  Run:  gap -q -o 4g run_03_permutation_census.g < /dev/null > run_03_permutation_census.out
##
#############################################################################

Read("leech_lib.g");

LAT := BuildLattices();;
M2  := Mod2Setup(LAT);;
B := LAT.B;;  Binv := LAT.Binv;;

CycleTypeOf := function(p)    # partition of the non-trivial cycles, descending
  local cs, out, i;
  cs := CycleStructurePerm(p);
  out := [];
  for i in [1..Length(cs)] do
    if IsBound(cs[i]) then Append(out, ListWithIdenticalEntries(cs[i], i + 1)); fi;
  od;
  return Reversed(SortedList(out));
end;

TypeName := function(ct)
  if ct = [] then return "id"; fi;
  if Length(ct) = 1 then return Concatenation(String(ct[1]), "-cycle"); fi;
  return Concatenation("(", JoinStringsWithSeparator(List(ct, String), ","), ")");
end;

##  paper Sec. 5.6 table: [cycle type, number of permutations, number closing]
PAPER := [ [[2],21,21], [[3],70,35], [[2,2],105,42], [[4],210,0], [[3,2],420,126],
           [[5],504,252], [[2,2,2],105,21], [[3,3],280,112], [[4,2],630,294],
           [[6],840,210], [[3,2,2],210,105], [[4,3],420,168], [[5,2],504,42],
           [[7],720,336] ];

perms := AsSSortedList(SymmetricGroup(7));;
Print("=== run_03_permutation_census.g ===\n");
Print("permutations to test (identity included): ", Length(perms), "\n");

t0 := Runtime();;
DATA := [];;
for p in perms do
  G := PermTwistMat(p);
  sd := SpanData(TwistedTensor(G), B, Binv);
  tl := TwistLagrangians(LAT, G);
  pv := SubspaceProps(LAT, M2, tl.V);
  pw := SubspaceProps(LAT, M2, tl.W);
  rec1 := rec( p := p, ct := CycleTypeOf(p), sign := SignPerm(p), closes := sd.closes,
               Vsub := pv.subalgebra, Wsub := pw.subalgebra, dimVV := pv.dimSquare,
               dimWW := pw.dimSquare, Videal := [pv.leftIdeal, pv.rightIdeal],
               Wideal := [pw.leftIdeal, pw.rightIdeal] );
  if sd.closes then
    rec1.ed := sd.ed;
    rec1.index := sd.index;
    rec1.dimS2L3 := DimSmod2L3(LAT, sd.rows);
  fi;
  Add(DATA, rec1);
od;
Print("computation time: ", Runtime() - t0, " ms\n\n");

PosOf := p -> Position(perms, p);
ParityName := function(sgn) if sgn = 1 then return "even"; fi; return "odd"; end;

## ---------------------------------------------------------------- table
Print("Closing counts per cycle type (this computation vs. paper Sec. 5.6):\n");
Print("  type       perms  close  paper  match  parity  elementary divisors of closing twists\n");
allmatch := true;;
for row in PAPER do
  sel := Filtered(DATA, d -> d.ct = row[1]);
  cl := Filtered(sel, d -> d.closes);
  m := (Length(sel) = row[2] and Length(cl) = row[3]);
  allmatch := allmatch and m;
  Print("  ", String(TypeName(row[1]), -9), String(Length(sel), 6), String(Length(cl), 7),
        String(row[3], 7), "  ", String(m, -6), String(ParityName(sel[1].sign), -6), "  ",
        Collected(List(cl, d -> EDString(d.ed))), "\n");
od;
idd := DATA[PosOf(())];;
Print("  identity: closes = ", idd.closes, "\n");
closing := Filtered(DATA, d -> d.closes and d.p <> ());;
Print("  total non-identity closing: ", Length(closing), " (paper: 1764); all rows match the paper: ", allmatch, "\n\n");

## ---------------------------------------------------------------- C2
ev := Filtered(closing, d -> d.sign = 1);;
od_ := Filtered(closing, d -> d.sign = -1);;
Print("Claim C2 (span index by parity):\n");
Print("  closing even permutations: ", Length(ev), ", elementary divisors ", Collected(List(ev, d -> EDString(d.ed))), "\n");
Print("  closing odd permutations:  ", Length(od_), ", elementary divisors ", Collected(List(od_, d -> EDString(d.ed))), "\n");
Print("  C2 holds: ", Length(ev) = 1176 and Length(od_) = 588
      and ForAll(ev, d -> d.ed = Concatenation(ListWithIdenticalEntries(14, 1), ListWithIdenticalEntries(10, 2)))
      and ForAll(od_, d -> d.ed = Concatenation(ListWithIdenticalEntries(8, 1), ListWithIdenticalEntries(16, 2))), "\n\n");

## ---------------------------------------------------------------- 3-cycles
three := Filtered(DATA, d -> d.ct = [3]);;
subsets := Set(List(three, d -> Set(MovedPoints(d.p))));;
Print("3-cycles: ", Length(subsets), " three-element subsets; exactly one orientation closes for each: ",
      ForAll(subsets, S -> Number(three, d -> Set(MovedPoints(d.p)) = S and d.closes) = 1), "\n");
Print("  closing 3-cycles (GAP notation, pi(a) = a^p): ", List(Filtered(three, d -> d.closes), d -> d.p), "\n\n");

## ---------------------------------------------------------------- convention
nInvAgree := Number(DATA, d -> d.closes = DATA[PosOf(d.p^-1)].closes);;
Print("Convention dependence: pi and pi^{-1} have the same closure verdict for ", nInvAgree,
      " of 5040 permutations\n");
Print("  (involutions always agree; e.g. (5,6,7) closes but (5,7,6) does not)\n\n");

## ---------------------------------------------------------------- C5 / mod-2
Print("Mod-2 criterion (claim C5): closes <=> V and W both subalgebras of L/2L: ",
      ForAll(DATA, d -> d.closes = (d.Vsub and d.Wsub)), "\n");
Print("Mod-2 prediction dim (S+2L^3)/2L^3 = dim W.W + 2 dim V.V on all closing twists: ",
      ForAll(Filtered(DATA, d -> d.closes), d -> d.dimS2L3 = d.dimWW + 2 * d.dimVV), "\n");
Print("Closing twists by (dim V.V, dim W.W, V left ideal, V right ideal, elementary divisors):\n");
for c in Collected(List(Filtered(DATA, d -> d.closes),
                        d -> [d.dimVV, d.dimWW, d.Videal[1], d.Videal[2], EDString(d.ed)])) do
  Print("  ", c[1], " : ", c[2], "\n");
od;
Print("Non-closing twists by (V subalgebra, W subalgebra): ",
      Collected(List(Filtered(DATA, d -> not d.closes), d -> [d.Vsub, d.Wsub])), "\n");

## ---------------------------------------------------------------- table file
f := "perm_census_table.txt";;
PrintTo(f, "# All 5040 permutations pi of e1..e7 (pi(e0) = e0); twist x ._pi y = pi^{-1}(pi(x).pi(y)).\n",
  "# columns: GAP cycle notation (pi(a) = a^p) | images [pi(1),...,pi(7)] | cycle type | sign |\n",
  "#          closes | elementary divisors of S (if closes) | closes for pi^{-1} |\n",
  "#          dim span(V.V), dim span(W.W) for V = pi(L sbar)/2L, W = pi(L s)/2L | V sub, W sub\n");
for d in DATA do
  AppendTo(f, String(d.p), " | ", List([1..7], a -> a^d.p), " | ", TypeName(d.ct), " | ", d.sign, " | ",
           d.closes, " | ");
  if d.closes then AppendTo(f, EDString(d.ed)); else AppendTo(f, "-"); fi;
  AppendTo(f, " | ", DATA[PosOf(d.p^-1)].closes, " | ", d.dimVV, ",", d.dimWW, " | ", d.Vsub, ",", d.Wsub, "\n");
od;
Print("\nper-permutation table written to ", f, "\n");
Print("=== end run_03 ===\n");
QUIT;
