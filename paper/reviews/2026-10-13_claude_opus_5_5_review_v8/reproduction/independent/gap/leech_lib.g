#############################################################################
##
##  leech_lib.g                                                    (GAP 4.12)
##
##  Independent GAP implementation for checking the span-index claims about
##  the Z3-symmetric triple-octonion product on Wilson's Leech lattice
##  (paper v8: "An order on the Leech lattice from a Z3-symmetric
##  triple-octonion product").  Written from the paper's definitions only.
##
##  All arithmetic is exact (GAP integers and rationals; GF(2) only for the
##  mod-2 linear algebra that cuts Lambda out of L^3).
##
##  CONVENTIONS
##  -----------
##  (1) Vectors are ROW vectors (GAP lists).  In R^8 the list index k = 1..8
##      holds the coefficient of e_{k-1}: index 1 <-> e0 (the identity),
##      index 2 <-> e1, ..., index 8 <-> e7.  A vector of R^24 = O+O+O is the
##      concatenation of three 8-blocks: block 1 = indices 1..8, block 2 =
##      9..16, block 3 = 17..24.
##  (2) Octonion product: Fano triples (1,2,4),(2,3,5),(3,4,6),(4,5,7),
##      (5,6,1),(6,7,2),(7,1,3).  For a triple (a,b,c): e_a e_b = e_c,
##      e_b e_c = e_a, e_c e_a = e_b, and the three reversed products are the
##      negatives; e_i e_i = -e0 for i = 1..7; e0 is the two-sided identity.
##  (3) L = D8+ : all coordinates integers with even sum, or all coordinates
##      in Z+1/2 with ODD sum.  s = (-e0+e1+...+e7)/2, sbar = (-e0-...-e7)/2.
##      L sbar = { l . sbar : l in L },  L s = { l . s : l in L }  (LEFT
##      multiples, untwisted product).
##  (4) Lambda = { (x,y,z) in L^3 : x+y, x+z, y+z in L sbar ; x+y+z in L s }.
##  (5) A linear map g of R^8 acts on row vectors from the right:
##      g(x) = x * G for an 8x8 rational matrix G ("row convention").
##      Twisted product:  x ._g y := g^{-1}( g(x) . g(y) ) = ((x*G).(y*G))*G^{-1}.
##  (6) Coordinate permutations of the imaginary units: a GAP permutation p
##      of [1..7] stands for the map pi with pi(e0) = e0 and
##      pi(e_a) = e_{a^p}  (a^p = GAP image of the point a under p).
##      Its matrix in row convention has G[a+1][a^p+1] = 1, G[1][1] = 1.
##      Example: the GAP cycle (5,6,7) is pi: e5 -> e6 -> e7 -> e5, i.e.
##      pi(5)=6, pi(6)=7, pi(7)=5.  The resulting twisted product is the
##      paper's  x ._pi y = pi^{-1}( pi(x) . pi(y) )  (paper Sec. 5.6).
##      For sigma = (1,2) (an involution) the convention does not matter.
##  (7) Root reflections: for a root r of L (r.r = 2),
##      s_r(x) = x - <x,r> r, matrix G = I - r^T r (symmetric, G^2 = I),
##      so g = g^{-1} and the g / g^{-1} convention does not matter.
##  (8) Triple product (paper Def. 3.3), every block product twisted:
##        P1 = x.x' + z.y' + y.z',  P2 = y.y' + x.z' + z.x',
##        P3 = z.z' + y.x' + x.y'.
##      With blocks indexed 0,1,2: u_i . v_j lands in block k = -(i+j) mod 3.
##  (9) Closure: u * v in Lambda for all u,v in Lambda; by bilinearity it is
##      enough that the 576 products b_a * b_c of a Z-basis (b_1..b_24) of
##      Lambda have INTEGER coordinates in that basis.
## (10) Span: S = Z-span of all products = Z-span of the 576 basis products.
##      With C the 576 x 24 integer matrix of their coordinates, Lambda/S is
##      isomorphic to (+)_i Z/d_i where d_1 | d_2 | ... | d_24 are the
##      elementary divisors of C (GAP: ElementaryDivisorsMat), and the index
##      [Lambda : S] = d_1 * ... * d_24 (when all d_i are nonzero).
##
#############################################################################

SizeScreen([200, 50]);    # avoid GAP's 80-column line wrapping in outputs

#############################################################################
##  Octonion multiplication table
#############################################################################

FANO_TRIPLES := [ [1,2,4], [2,3,5], [3,4,6], [4,5,7], [5,6,1], [6,7,2], [7,1,3] ];

##  OCT_TABLE[i+1][j+1] = [ sign, k+1 ]  means  e_i * e_j = sign * e_k.
BuildOctTable := function()
  local T, i, j, t, a, b, c, set;
  T := List([1..8], i -> List([1..8], j -> fail));
  set := function(i, j, sgn, k)      # i,j,k are octonion indices 0..7
    if T[i+1][j+1] <> fail then Error("table entry assigned twice"); fi;
    T[i+1][j+1] := [sgn, k+1];
  end;
  for i in [0..7] do set(0, i, 1, i); od;            # e0 e_i = e_i
  for i in [1..7] do set(i, 0, 1, i); od;            # e_i e0 = e_i
  for i in [1..7] do set(i, i, -1, 0); od;           # e_i e_i = -e0
  for t in FANO_TRIPLES do
    a := t[1]; b := t[2]; c := t[3];
    set(a, b,  1, c);  set(b, c,  1, a);  set(c, a,  1, b);   # cyclic
    set(b, a, -1, c);  set(c, b, -1, a);  set(a, c, -1, b);   # anticyclic
  od;
  for i in [1..8] do for j in [1..8] do
    if T[i][j] = fail then Error("table incomplete"); fi;
  od; od;
  return T;
end;
OCT_TABLE := BuildOctTable();

##  Untwisted octonion product of two coefficient vectors of length 8.
OctMul := function(x, y)
  local z, i, j, t;
  z := ListWithIdenticalEntries(8, 0);
  for i in [1..8] do
    if x[i] <> 0 then
      for j in [1..8] do
        if y[j] <> 0 then
          t := OCT_TABLE[i][j];
          z[t[2]] := z[t[2]] + t[1] * x[i] * y[j];
        fi;
      od;
    fi;
  od;
  return z;
end;

UnitVec8 := i -> ShallowCopy(IdentityMat(8)[i]);   # e_{i-1} as a list (i = 1..8)
IdentityMat8 := IdentityMat(8);

#############################################################################
##  Twisted products
#############################################################################

##  Structure tensor of x ._g y = ((x*G).(y*G)) * G^{-1}:
##  Tg[i] is the 8x8 matrix whose row j is  e_{i-1} ._g e_{j-1}.
##  Hence the left multiplication matrix of x is  Sum_i x[i]*Tg[i]  and
##  x ._g y = y * (Sum_i x[i]*Tg[i]).
TwistedTensor := function(G)
  local Gi;
  Gi := Inverse(G);
  return List([1..8], i -> List([1..8], j -> OctMul(G[i], G[j]) * Gi));
end;

##  Direct (naive) twisted product, used only for cross-checks.
TwMulDirect := function(G, Gi, x, y)
  return OctMul(x * G, y * G) * Gi;
end;

##  Left multiplication matrix (8x8) of x under the twisted product Tg.
TwLeftMat := function(Tg, x)
  local M, i;
  M := NullMat(8, 8);
  for i in [1..8] do
    if x[i] <> 0 then M := M + x[i] * Tg[i]; fi;
  od;
  return M;
end;

##  24x24 left-multiplication matrix of u for the triple product:
##  v * TripleLeftMat(Tg,u) = u * v.  Block (jb,kb) (v-block jb -> output
##  block kb, blocks 0,1,2) is the left multiplication by u-block
##  ib = -(jb+kb) mod 3.
TripleLeftMat := function(Tg, u)
  local M, Lb, ib, jb, kb, r, c;
  Lb := List([0..2], b -> TwLeftMat(Tg, u{[8*b+1 .. 8*b+8]}));
  M := NullMat(24, 24);
  for jb in [0..2] do
    for kb in [0..2] do
      ib := (-(jb + kb)) mod 3;
      for r in [1..8] do
        for c in [1..8] do
          M[8*jb + r][8*kb + c] := Lb[ib+1][r][c];
        od;
      od;
    od;
  od;
  return M;
end;

##  Direct (naive) triple product straight from the paper's formulas
##  P1 = x.x' + z.y' + y.z', P2 = y.y' + x.z' + z.x', P3 = z.z' + y.x' + x.y'.
TripleProdDirect := function(G, u, v)
  local Gi, m, x, y, z, xp, yp, zp;
  Gi := Inverse(G);
  m := function(a, b) return TwMulDirect(G, Gi, a, b); end;
  x := u{[1..8]};  y := u{[9..16]};  z := u{[17..24]};
  xp := v{[1..8]}; yp := v{[9..16]}; zp := v{[17..24]};
  return Concatenation( m(x,xp) + m(z,yp) + m(y,zp),
                        m(y,yp) + m(x,zp) + m(z,xp),
                        m(z,zp) + m(y,xp) + m(x,yp) );
end;

#############################################################################
##  Twist matrices
#############################################################################

##  Row-convention matrix of the coordinate permutation pi(e_a) = e_{a^p}.
PermTwistMat := function(p)
  local G, a;
  G := NullMat(8, 8);
  G[1][1] := 1;
  for a in [1..7] do G[a+1][(a^p)+1] := 1; od;
  return G;
end;

##  Row-convention matrix of the signed permutation e_a -> sgn[a] e_{a^p}.
SignedPermTwistMat := function(p, sgn)
  local G, a;
  G := NullMat(8, 8);
  G[1][1] := 1;
  for a in [1..7] do G[a+1][(a^p)+1] := sgn[a]; od;
  return G;
end;

##  Row-convention matrix of the reflection s_r(x) = x - <x,r> r (r.r = 2).
RootReflectionMat := function(r)
  if r * r <> 2 then Error("r must have norm 2"); fi;
  return IdentityMat8 - TransposedMat([r]) * [r];
end;

#############################################################################
##  Lattices
#############################################################################

##  Membership in L = D8+.
InL := function(x)
  if ForAll(x, IsInt) then return IsEvenInt(Sum(x)); fi;
  if ForAll(x, c -> IsInt(c - 1/2)) then return IsOddInt(Sum(x)); fi;
  return false;
end;

##  Membership of x in the full-rank lattice with basis matrix A, given A^{-1}.
InLatticeInv := function(x, Ainv)
  return ForAll(x * Ainv, IsInt);
end;

##  The 240 roots of L = D8+.
RootsOfL := function()
  local R, i, j, si, sj, v, signs, m;
  R := [];
  for i in [1..8] do for j in [i+1..8] do
    for si in [1,-1] do for sj in [1,-1] do
      v := ListWithIdenticalEntries(8, 0); v[i] := si; v[j] := sj;
      Add(R, v);
    od; od;
  od; od;
  for signs in Tuples([1,-1], 8) do
    m := Number(signs, t -> t = -1);
    if IsOddInt(m) then Add(R, signs / 2); fi;   # odd number of minus signs <=> odd sum
  od;
  return R;
end;

##  Z-basis of the Wilson-type lattice
##     { (x,y,z) in L^3 : x+y, x+z, y+z in Vhat ; x+y+z in What }
##  where BL is a Z-basis of L and BV, BW are Z-bases of sublattices Vhat,
##  What of L that contain 2L (for Lambda itself: Vhat = L sbar, What = L s).
##  Method: L^3 = Z^24 in the coordinates of BL (blockwise); since 2L^3 lies in
##  the lattice, it is the full preimage of the GF(2)-solution space K of the
##  16 linear conditions "class of x+y in Vhat/2L" etc.; a Z-basis is the
##  Hermite basis of lift(K) together with 2*I_24.
BuildWilsonTypeBasis := function(BL, BV, BW)
  local BLinv, CV2, CW2, QV, QW, A, put, K, Kint, Clam, Bmat;
  BLinv := Inverse(BL);
  ##  coordinates of the sublattice bases in the basis of L, reduced mod 2:
  ##  their row spaces are Vhat/2L and What/2L inside L/2L = GF(2)^8
  CV2 := (BV * BLinv) * Z(2);
  CW2 := (BW * BLinv) * Z(2);
  ##  Q-matrices: a class c in L/2L lies in the subspace iff c * Q = 0
  QV := TransposedMat(NullspaceMat(TransposedMat(CV2)));
  QW := TransposedMat(NullspaceMat(TransposedMat(CW2)));
  ##  condition matrix A (24 x 16 over GF(2)) acting on c = (cx, cy, cz):
  ##  columns 1-4: x+y in Vhat; 5-8: x+z in Vhat; 9-12: y+z in Vhat;
  ##  13-16: x+y+z in What
  A := NullMat(24, 16, GF(2));
  put := function(rowb, colb, Q)
    local r, c;
    for r in [1..8] do for c in [1..4] do
      A[8*(rowb-1) + r][4*(colb-1) + c] := A[8*(rowb-1) + r][4*(colb-1) + c] + Q[r][c];
    od; od;
  end;
  put(1, 1, QV); put(2, 1, QV);                 # x + y
  put(1, 2, QV); put(3, 2, QV);                 # x + z
  put(2, 3, QV); put(3, 3, QV);                 # y + z
  put(1, 4, QW); put(2, 4, QW); put(3, 4, QW);  # x + y + z
  K := NullspaceMat(A);
  Kint := List(K, v -> List(v, IntFFE));
  Clam := BaseIntMat(Concatenation(Kint, 2 * IdentityMat(24)));
  Bmat := Clam * DirectSumMat(BL, BL, BL);
  return rec(dimK := Length(K), Clam := Clam, B := Bmat, Binv := Inverse(Bmat),
             rankV := RankMat(CV2), rankW := RankMat(CW2));
end;

##  Build L, L sbar, L s and a Z-basis of Lambda.  Returns a record.
BuildLattices := function()
  local R, BL, s, sbar, BLsbar, BLs, W;
  R := RootsOfL();
  if Length(R) <> 240 or not ForAll(R, InL) or not ForAll(R, r -> r*r = 2) then
    Error("roots of L wrong");
  fi;
  ##  Z-basis of L: Hermite basis of the Z-span of the 240 roots (E8 is
  ##  generated by its roots), computed on 2*roots to stay integral.
  BL := BaseIntMat(2 * R) / 2;
  if Length(BL) <> 8 or AbsInt(DeterminantMat(BL)) <> 1 or not ForAll(BL, InL) then
    Error("basis of L wrong");
  fi;
  s    := [-1, 1, 1, 1, 1, 1, 1, 1] / 2;
  sbar := [-1,-1,-1,-1,-1,-1,-1,-1] / 2;
  ##  Z-bases of L sbar and L s: images of the basis of L under the injective
  ##  linear maps l -> l . sbar and l -> l . s (left multiples).
  BLsbar := List(BL, l -> OctMul(l, sbar));
  BLs    := List(BL, l -> OctMul(l, s));
  W := BuildWilsonTypeBasis(BL, BLsbar, BLs);
  return rec( roots := R, BL := BL, BLinv := Inverse(BL), s := s, sbar := sbar,
              BLsbar := BLsbar, BLs := BLs, BLsbarInv := Inverse(BLsbar),
              BLsInv := Inverse(BLs), dimK := W.dimK, Clam := W.Clam,
              B := W.B, Binv := W.Binv );
end;

##  Membership in Lambda straight from Wilson's definition (no Lambda basis used).
InLambdaDef := function(lat, u)
  local x, y, z;
  x := u{[1..8]}; y := u{[9..16]}; z := u{[17..24]};
  return InL(x) and InL(y) and InL(z)
     and InLatticeInv(x + y, lat.BLsbarInv)
     and InLatticeInv(x + z, lat.BLsbarInv)
     and InLatticeInv(y + z, lat.BLsbarInv)
     and InLatticeInv(x + y + z, lat.BLsInv);
end;

#############################################################################
##  Closure test and span index
#############################################################################

##  For the twist with structure tensor Tg, form the 576 x 24 matrix of
##  coordinates (in the basis B of Lambda) of all products b_a * b_c,
##  row index 24*(a-1)+c.  Returns a record with
##    closes : all 576*24 coordinates are integers
##    ed     : elementary divisors (only if closes)
##    index  : product of the elementary divisors (0 if rank < 24)
SpanData := function(Tg, B, Binv)
  local rows, a, Ca, closes, ed;
  rows := [];
  for a in [1..24] do
    Ca := B * TripleLeftMat(Tg, B[a]) * Binv;   # row c = coords of b_a * b_c
    Append(rows, Ca);
  od;
  closes := ForAll(rows, r -> ForAll(r, IsInt));
  if not closes then
    return rec(closes := false, rows := rows);
  fi;
  ed := ElementaryDivisorsMat(Integers, rows);
  return rec(closes := true, rows := rows, ed := ed, index := Product(ed));
end;

##  Compact text form of an elementary-divisor list, e.g. "1^8 2^16".
EDString := function(ed)
  return JoinStringsWithSeparator(
           List(Collected(ed), p -> Concatenation(String(p[1]), "^", String(p[2]))), " ");
end;

##  Exponent of 2 in a positive integer (for printing indices as 2^k).
Log2Exact := function(n)
  local k;
  k := 0;
  while IsEvenInt(n) and n > 0 do n := n / 2; k := k + 1; od;
  if n <> 1 then return fail; fi;
  return k;
end;

##  Certificate that the rows (coordinates of products in a basis of Lambda)
##  generate Z^24, i.e. S = Lambda, in a form that is easy to re-check:
##  two 24-subsets A and B of the rows with |det A| = 2^a and |det B| = m odd.
##  Since span(A) has index 2^a in Lambda, 2^a*Lambda <= span(A) <= S; likewise
##  m*Lambda <= S; gcd(2^a, m) = 1 then gives Lambda <= S.
##  (A single 24-subset with det +-1 was searched for and not found: local
##  minima of |det| are either highly 2-divisible or large and odd.)
##  B: greedy selection of rows that are independent over GF(2) (odd det).
##  A: local search from a greedy Q-independent selection.  Replacing selected
##     row k by another row r multiplies the determinant d by c_k, the k-th
##     coordinate of r in the current selection, so the new determinant is the
##     integer d*c_k; the swap is made whenever the odd part of |d*c_k| is
##     smaller than the odd part of |d|.  When no such swap exists the search
##     restarts from a shuffled row order; it stops when the odd part is 1.
OddPart := function(n)
  while n <> 0 and IsEvenInt(n) do n := n / 2; od;
  return n;
end;

CoprimeDetCertificate := function(rows, maxRestarts)
  local n, selB, i, order, attempt, sel, M, Minv, d, improved, r, coeffs, cand, k;
  n := Length(rows[1]);
  ##  B: odd determinant
  selB := [];
  for i in [1..Length(rows)] do
    if RankMat(rows{Concatenation(selB, [i])} * Z(2)) = Length(selB) + 1 then Add(selB, i); fi;
    if Length(selB) = n then break; fi;
  od;
  if Length(selB) < n then return fail; fi;
  ##  A: determinant a power of 2
  for attempt in [1..maxRestarts] do
    if attempt = 1 then order := List([1..Length(rows)]);
    else order := Shuffle(List([1..Length(rows)])); fi;
    sel := [];
    for i in order do
      if RankMat(rows{Concatenation(sel, [i])}) = Length(sel) + 1 then Add(sel, i); fi;
      if Length(sel) = n then break; fi;
    od;
    repeat
      M := rows{sel};
      d := AbsInt(DeterminantMat(M));
      if OddPart(d) = 1 then
        return rec(selA := sel, detA := DeterminantMat(M),
                   selB := selB, detB := DeterminantMat(rows{selB}));
      fi;
      Minv := Inverse(M);
      improved := false;
      for r in order do
        if not r in sel then
          coeffs := rows[r] * Minv;     # rows[r] = sum_j coeffs[j] * M[j]
          cand := Filtered([1..n], j -> coeffs[j] <> 0 and
                     OddPart(AbsInt(NumeratorRat(d * coeffs[j]))) < OddPart(d));
          if cand <> [] then
            k := cand[PositionMinimum(List(cand, j -> OddPart(AbsInt(NumeratorRat(d * coeffs[j])))))];
            sel[k] := r;
            improved := true;
            break;
          fi;
        fi;
      od;
    until not improved;
  od;
  return fail;
end;

#############################################################################
##  The mod-2 quotient Lbar = L/2L
##  (classes are GF(2)-vectors of coordinates in the basis BL of L)
#############################################################################

ZERO8_GF2 := ListWithIdenticalEntries(8, 0 * Z(2));

##  Class in L/2L of a vector x of L.
ClassMod2 := function(lat, x)
  return (x * lat.BLinv) * Z(2);
end;

##  Product that L/2L inherits from L.L in L: structure constants
##  MU[i][j] = class of BL[i].BL[j]  (untwisted octonion product).
Mod2Setup := function(lat)
  return rec( MU := List([1..8], i -> List([1..8], j ->
                       ClassMod2(lat, OctMul(lat.BL[i], lat.BL[j])))),
              full := VectorSpace(GF(2), IdentityMat(8, GF(2))) );
end;

##  mu-bar(a, b) for classes a, b.
MuBar := function(M2, a, b)
  local z, i, j;
  z := ZERO8_GF2;
  for i in [1..8] do
    if not IsZero(a[i]) then
      for j in [1..8] do
        if not IsZero(b[j]) then z := z + M2.MU[i][j]; fi;
      od;
    fi;
  od;
  return z;
end;

SpanMod2 := vecs -> VectorSpace(GF(2), vecs, ZERO8_GF2);

##  Span of mu-bar(U1, U2) for subspaces U1, U2 of L/2L.
ProductSpanMod2 := function(M2, U1, U2)
  return SpanMod2(Concatenation(List(Basis(U1), a -> List(Basis(U2), b -> MuBar(M2, a, b)))));
end;

##  Quadratic form q(c) = N(lift)/2 mod 2 on L/2L.
QForm := function(lat, c)
  local x;
  x := List(c, IntFFE) * lat.BL;
  return (x * x / 2) mod 2;
end;

##  Algebraic properties of a subspace U of L/2L.
SubspaceProps := function(lat, M2, U)
  local UU, LU, UL;
  UU := ProductSpanMod2(M2, U, U);
  LU := ProductSpanMod2(M2, M2.full, U);
  UL := ProductSpanMod2(M2, U, M2.full);
  return rec( dim := Dimension(U),
              totallySingular := ForAll(AsList(U), c -> QForm(lat, c) = 0),
              subalgebra := IsSubset(U, UU),
              leftIdeal := IsSubset(U, LU),       # Lbar . U in U
              rightIdeal := IsSubset(U, UL),      # U . Lbar in U
              dimSquare := Dimension(UU) );       # dim span(U.U)
end;

##  The two Lagrangians moved by a twist g (matrix G):
##  V_g = g(L sbar)/2L and W_g = g(L s)/2L.
TwistLagrangians := function(lat, G)
  local gV, gW;
  gV := List(lat.BLsbar, x -> x * G);
  gW := List(lat.BLs, x -> x * G);
  if not ForAll(Concatenation(gV, gW), InL) then Error("g does not preserve L"); fi;
  return rec( V := SpanMod2(List(gV, x -> ClassMod2(lat, x))),
              W := SpanMod2(List(gW, x -> ClassMod2(lat, x))) );
end;

##  dim over GF(2) of (S + 2L^3)/2L^3 for the span S given by its 576 rows of
##  Lambda-coordinates (L^3-coordinates of the rows are rows * Clam).
DimSmod2L3 := function(lat, rows)
  return RankMat((rows * lat.Clam) * Z(2));
end;

##  Gram determinant of a Z-basis of S (basis-free index check:
##  det Gram(S) = det Gram(Lambda) * [Lambda:S]^2).
SpanGramDet := function(lat, rows)
  local Sb, A;
  Sb := BaseIntMat(rows);
  if Length(Sb) < 24 then return 0; fi;
  A := Sb * lat.B;
  return DeterminantMat(A * TransposedMat(A));
end;

##  dim (U + span(U.Lbar) + span(Lbar.U)) / U : how far U is from a
##  two-sided ideal of Lbar (0 iff U is a two-sided ideal).
IdealDefect := function(M2, U)
  local gens;
  gens := Concatenation(BasisVectors(Basis(U)),
                        BasisVectors(Basis(ProductSpanMod2(M2, U, M2.full))),
                        BasisVectors(Basis(ProductSpanMod2(M2, M2.full, U))));
  return Dimension(SpanMod2(gens)) - Dimension(U);
end;
