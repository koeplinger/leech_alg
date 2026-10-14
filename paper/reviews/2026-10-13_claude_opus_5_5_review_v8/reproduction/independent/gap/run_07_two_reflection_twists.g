#############################################################################
##
##  run_07_two_reflection_twists.g                                 (GAP 4.12)
##
##  Twists by products of two root reflections, g = s_{r2} o s_{r1}
##  (row convention: G = G_{r1} * G_{r2}, g(x) = x*G), r1, r2 roots of L
##  with r2 <> +-r1 (one root per +- pair, ordered pairs: 120*119 = 14,280).
##  Purpose: explicit witnesses for every span-index class found in the
##  polarization census (run_06), in particular 2^8 and 2^5, which no
##  coordinate permutation and no single reflection realizes.
##
##  For every g: the mod-2 data of V = g(L sbar)/2L, W = g(L s)/2L; for every
##  g with V and W subalgebras (the closure criterion), the full closure test
##  on the 576 basis products and the elementary divisors of S.  A sample of
##  the remaining g is also run through the full test to confirm non-closure.
##
##  Run:  gap -q -o 4g run_07_two_reflection_twists.g < /dev/null > run_07_two_reflection_twists.out
##
#############################################################################

Read("leech_lib.g");
Reset(GlobalMersenneTwister, 20261013);

LAT := BuildLattices();;
M2  := Mod2Setup(LAT);;
B := LAT.B;;  Binv := LAT.Binv;;

Print("=== run_07_two_reflection_twists.g ===\n");
reps := Filtered(LAT.roots, r -> First(r, c -> c <> 0) > 0);;
GR := List(reps, RootReflectionMat);;

t0 := Runtime();;
DATA := [];;
G := fail;; tl := fail;; pv := fail;; pw := fail;; sd := fail;;
for i in [1..120] do
  for j in [1..120] do
    if i <> j then
      G := GR[i] * GR[j];
      tl := TwistLagrangians(LAT, G);
      pv := SubspaceProps(LAT, M2, tl.V);
      pw := SubspaceProps(LAT, M2, tl.W);
      Add(DATA, rec(i := i, j := j, crit := pv.subalgebra and pw.subalgebra,
                    dimVV := pv.dimSquare, dimWW := pw.dimSquare));
    fi;
  od;
od;
Print("mod-2 data for ", Length(DATA), " twists: ", Runtime() - t0, " ms\n");

t0 := Runtime();;
for d in Filtered(DATA, d -> d.crit) do
  sd := SpanData(TwistedTensor(GR[d.i] * GR[d.j]), B, Binv);
  d.closes := sd.closes;
  if sd.closes then d.ed := EDString(sd.ed); fi;
od;
sample := Filtered(DATA, d -> not d.crit);;
sample := sample{List([1..300], k -> Random([1..Length(sample)]))};;
for d in sample do
  d.closes := SpanData(TwistedTensor(GR[d.i] * GR[d.j]), B, Binv).closes;
od;
Print("full closure/span tests: ", Runtime() - t0, " ms\n\n");

Print("twists passing the mod-2 criterion: ", Number(DATA, d -> d.crit), " of ", Length(DATA), "\n");
Print("  all of them close (576 basis products): ", ForAll(Filtered(DATA, d -> d.crit), d -> d.closes), "\n");
Print("  300 random twists failing the criterion: none closes: ", ForAll(sample, d -> not d.closes), "\n");
Print("closing twists by (dim V.V, dim W.W, elementary divisors): ",
      Collected(List(Filtered(DATA, d -> d.crit), d -> [d.dimVV, d.dimWW, d.ed])), "\n\n");

Print("first witness (in the order of the loops) for each span-index class:\n");
for e in Set(List(Filtered(DATA, d -> d.crit), d -> d.ed)) do
  d := First(DATA, d -> d.crit and d.ed = e);
  Print("  ", String(e, -10), "  g = s_r2 o s_r1 with r1 = ", reps[d.i], ", r2 = ", reps[d.j], "\n");
od;
Print("=== end run_07 ===\n");
QUIT;
