#############################################################################
##
##  run_05_signed_permutations.g                                   (GAP 4.12)
##
##  Claim C4: signed permutations of the imaginary units (e0 fixed),
##     g(e_a) = eps_a e_{a^p}   (a = 1..7, eps_a = +-1, p in S_7),
##  2^7 * 7! = 645,120 maps.  Two maps give the same twisted product
##  x ._g y = g^{-1}(g(x).g(y)) iff their structure constants agree; we
##  compute the structure constants of every one of the 645,120 maps
##  (no group theory assumed) and keep one representative per distinct
##  table.  For each distinct table: closure on Lambda and the elementary
##  divisors of the span S.
##
##  Structure constants of a signed permutation, i <> j imaginary:
##    g(e_i) g(e_j) = eps_i eps_j e_{p(i)} e_{p(j)} = eps_i eps_j s e_k
##    and g^{-1}(e_k) = eps_m e_m with m = p^{-1}(k), so
##    e_i ._g e_j = eps_i eps_j s eps_m e_m.
##  (e0 rows/columns and the squares e_i ._g e_i = -e0 never change.)
##  The 42 signed indices (eps_i eps_j s eps_m) * m form the table signature.
##
##  NOTE: g(L) = L only when the number of minus signs is even (negating an
##  odd number of coordinates flips the parity of the coordinate sum of the
##  half-integer vectors of D8+); this is checked below for all 128 sign
##  vectors.  The paper's twists require g(L) = L, so tables are reported
##  separately according to whether they come from L-preserving maps.  The
##  closure test itself (576 basis products of Lambda) is applied to all.
##
##  Run:  gap -q -o 4g run_05_signed_permutations.g < /dev/null > run_05_signed_permutations.out
##
#############################################################################

Read("leech_lib.g");

LAT := BuildLattices();;
M2  := Mod2Setup(LAT);;
B := LAT.B;;  Binv := LAT.Binv;;

Print("=== run_05_signed_permutations.g ===\n");

OFFDIAG := Filtered(Cartesian([1..7], [1..7]), ij -> ij[1] <> ij[2]);;   # 42 ordered pairs

Signature := function(p, eps)
  local pinv, sig, ij, i, j, t, k, m;
  pinv := p^-1;
  sig := [];
  for ij in OFFDIAG do
    i := ij[1]; j := ij[2];
    t := OCT_TABLE[(i^p) + 1][(j^p) + 1];     # e_{p(i)} e_{p(j)} = t[1] e_{t[2]-1}
    k := t[2] - 1;
    m := k^pinv;
    Add(sig, eps[i] * eps[j] * t[1] * eps[m] * m);
  od;
  return sig;
end;

##  sanity: the signature agrees with the full twisted tensor on random maps
Reset(GlobalMersenneTwister, 20261013);
checkSig := function(p, eps)
  local Tg, sig, n, ij, row;
  Tg := TwistedTensor(SignedPermTwistMat(p, eps));
  sig := Signature(p, eps);
  for n in [1..42] do
    ij := OFFDIAG[n];
    row := Tg[ij[1] + 1][ij[2] + 1];
    if row <> SignInt(sig[n]) * ShallowCopy(IdentityMat(8)[AbsInt(sig[n]) + 1]) then return false; fi;
  od;
  return true;
end;
Print("signature formula agrees with TwistedTensor on 200 random signed permutations: ",
      ForAll([1..200], i -> checkSig(Random(SymmetricGroup(7)), List([1..7], k -> Random([1,-1])))), "\n");

perms := AsSSortedList(SymmetricGroup(7));;
signs := Tuples([1, -1], 7);;
t0 := Runtime();;
SIGS := [];;   REPS := [];;   MULT := [];;   HASPERM := [];;   HASEVEN := [];;   HASODD := [];;
for p in perms do
  for eps in signs do
    sig := Signature(p, eps);
    ev := IsEvenInt(Number(eps, e -> e = -1));
    pos := PositionSorted(SIGS, sig);
    if pos > Length(SIGS) or SIGS[pos] <> sig then
      Add(SIGS, sig, pos);  Add(REPS, [p, eps], pos);  Add(MULT, 1, pos);
      Add(HASPERM, ForAll(eps, e -> e = 1), pos);
      Add(HASEVEN, ev, pos);  Add(HASODD, not ev, pos);
    else
      MULT[pos] := MULT[pos] + 1;
      if ForAll(eps, e -> e = 1) then HASPERM[pos] := true; fi;
      if ev then HASEVEN[pos] := true; else HASODD[pos] := true; fi;
    fi;
  od;
od;
Print("enumerated ", Length(perms) * Length(signs), " signed permutations in ", Runtime() - t0, " ms\n");
Print("distinct twisted tables: ", Length(SIGS), "\n");
Print("maps per table: ", Collected(MULT), "  (= order of the group of signed-permutation automorphisms of O)\n");
Print("tables realized by an unsigned permutation: ", Number(HASPERM, x -> x), "\n");
Print("g(L) = L  <=>  even number of minus signs (checked on all 128 sign vectors): ",
      ForAll(signs, eps -> ForAll(LAT.BL, b -> InL(b * SignedPermTwistMat((), eps)))
                           = IsEvenInt(Number(eps, e -> e = -1))), "\n");
Print("tables with members of both sign parities: ", Number([1..Length(SIGS)], n -> HASEVEN[n] and HASODD[n]), "\n");
Print("tables from L-preserving maps (even sign count): ", Number(HASEVEN, x -> x),
      "; of these realized by an unsigned permutation: ", Number([1..Length(SIGS)], n -> HASEVEN[n] and HASPERM[n]), "\n");
Print("tables from maps with g(L) <> L (odd sign count): ", Number(HASODD, x -> x), "\n\n");

## prefer an unsigned permutation as representative where one exists
for n in [1..Length(SIGS)] do
  if HASPERM[n] then
    REPS[n] := [First(perms, p -> Signature(p, ListWithIdenticalEntries(7, 1)) = SIGS[n]), ListWithIdenticalEntries(7, 1)];
  fi;
od;
Print("representatives preserve L exactly for the even-sign tables: ",
      ForAll([1..Length(SIGS)], n -> HASEVEN[n] = IsEvenInt(Number(REPS[n][2], e -> e = -1))), "\n");

t0 := Runtime();;
RES := [];;
for n in [1..Length(SIGS)] do
  G := SignedPermTwistMat(REPS[n][1], REPS[n][2]);
  sd := SpanData(TwistedTensor(G), B, Binv);
  if sd.closes then edstr := EDString(sd.ed); else edstr := "-"; fi;
  x := rec(closes := sd.closes, ed := edstr, hasperm := HASPERM[n], preservesL := HASEVEN[n],
           det := DeterminantMat(G));
  if HASEVEN[n] then
    tl := TwistLagrangians(LAT, G);
    pv := SubspaceProps(LAT, M2, tl.V);  pw := SubspaceProps(LAT, M2, tl.W);
    x.crit := (pv.subalgebra and pw.subalgebra);  x.dimVV := pv.dimSquare;  x.dimWW := pw.dimSquare;
  fi;
  Add(RES, x);
od;
Print("closure/span computation for ", Length(SIGS), " tables: ", Runtime() - t0, " ms\n\n");

Print("Claim C4 (one entry per distinct twisted table; closure tested on Wilson's Lambda for all):\n");
Print("  closing tables: ", Number(RES, x -> x.closes), " of ", Length(RES), "\n");
Print("  by elementary divisors: ", Collected(List(Filtered(RES, x -> x.closes), x -> x.ed)), "\n");
Print("  tables from L-preserving maps (= the 240 unsigned-permutation tables): ",
      Collected(List(Filtered(RES, x -> x.preservesL), x -> x.ed)), "\n");
Print("  tables from maps with g(L) <> L: ",
      Collected(List(Filtered(RES, x -> not x.preservesL), x -> x.ed)), "\n");
Print("  split by det(g) = +-1 (closing tables):\n");
for dt in [1, -1] do
  Print("    det ", dt, ": ", Collected(List(Filtered(RES, x -> x.closes and x.det = dt), x -> x.ed)), "\n");
od;
Print("  mod-2 criterion (C5) closes <=> V,W subalgebras on all L-preserving tables: ",
      ForAll(Filtered(RES, x -> x.preservesL), x -> x.closes = x.crit), "\n");
Print("  closing L-preserving tables by (dim V.V, dim W.W, elementary divisors): ",
      Collected(List(Filtered(RES, x -> x.closes and x.preservesL), x -> [x.dimVV, x.dimWW, x.ed])), "\n");
## Odd sign counts are conjugation composed with even ones: for c(x) = xbar
## (all seven imaginary signs flipped), c(a).c(b) = c(b.a), hence
## x ._{c g} y = y ._g x, and since the Z3 routing is symmetric in the two
## block indices, u *_{c g} v = v *_g u: the opposite algebra, with the same
## closure verdict and the same span.
OPP := List(OFFDIAG, ij -> Position(OFFDIAG, [ij[2], ij[1]]));;
Print("\nflipping all seven signs gives the opposite table (T'(i,j) = T(j,i)) for every table: ",
      ForAll(REPS, rp -> Signature(rp[1], -rp[2]) = Signature(rp[1], rp[2]){OPP}), "\n");
Print("and the same verdict/elementary divisors: ",
      ForAll([1..Length(SIGS)], n -> RES[n].ed = RES[PositionSorted(SIGS, Signature(REPS[n][1], -REPS[n][2]))].ed), "\n");
Print("=== end run_05 ===\n");
QUIT;
