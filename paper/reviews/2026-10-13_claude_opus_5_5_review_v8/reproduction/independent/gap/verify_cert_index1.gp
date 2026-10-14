\\ ===========================================================================
\\  verify_cert_index1.gp                                         (PARI/GP)
\\
\\  Second-system check, written from the paper's definitions and sharing no
\\  code with the GAP library:
\\   (a) re-implements the octonion product, L = D8+, L s, L sbar, Wilson's
\\       membership test for Lambda, twisted products and the Z3-routed
\\       triple product (Def. 3.3);
\\   (b) verifies the index-1 certificate cert_index1_reflection.gp written by
\\       run_02_witnesses.g (all vectors there are doubled, i.e. integral);
\\   (c) recomputes, with PARI's own mathnf/matsnf, the elementary divisors of
\\       the span S for sigma = (1 2), for pi: e5 -> e6 -> e7 -> e5, for its
\\       inverse, and for the reflection s_r, r = (e0+...+e6-e7)/2, using the
\\       certificate's Z-basis of Lambda (re-validated here first).
\\
\\  Conventions: vectors are row vectors, x[k] = coefficient of e_{k-1};
\\  a permutation pi is given by its image list P = [pi(1),...,pi(7)] and acts
\\  by pi(e_a) = e_{pi(a)}; twisted product x ._g y = g^{-1}(g(x) . g(y)).
\\
\\  Run:  gp -q verify_cert_index1.gp < /dev/null > verify_cert_index1.out
\\ ===========================================================================

read("cert_index1_reflection.gp");

\\ ---- octonion product: e_i e_j = SG[i+1,j+1] * e_{IX[i+1,j+1]}
FANO = [[1,2,4],[2,3,5],[3,4,6],[4,5,7],[5,6,1],[6,7,2],[7,1,3]];
SG = matrix(8,8); IX = matrix(8,8);
{
  for(i = 0, 7, SG[1,i+1] = 1; IX[1,i+1] = i; SG[i+1,1] = 1; IX[i+1,1] = i);
  for(i = 1, 7, SG[i+1,i+1] = -1; IX[i+1,i+1] = 0);
  for(t = 1, 7,
    my(a = FANO[t][1], b = FANO[t][2], c = FANO[t][3]);
    SG[a+1,b+1] =  1; IX[a+1,b+1] = c;  SG[b+1,a+1] = -1; IX[b+1,a+1] = c;
    SG[b+1,c+1] =  1; IX[b+1,c+1] = a;  SG[c+1,b+1] = -1; IX[c+1,b+1] = a;
    SG[c+1,a+1] =  1; IX[c+1,a+1] = b;  SG[a+1,c+1] = -1; IX[a+1,c+1] = b);
}
omul(x, y) =
{
  my(z = vector(8));
  for(i = 1, 8, if(x[i], for(j = 1, 8, if(y[j], z[IX[i,j]+1] += SG[i,j] * x[i] * y[j]))));
  z;
}

\\ ---- L = D8+ (odd-sum half-integer coset), its sublattices, Lambda
isint(t) = type(t) == "t_INT";
inL(x) =
{
  if(#select(t -> !isint(t), x) == 0, return(vecsum(x) % 2 == 0));
  if(#select(t -> denominator(t) != 2, x) == 0, return(vecsum(x) % 2 == 1));
  0;
}
s  = [-1,1,1,1,1,1,1,1] / 2;
sb = [-1,-1,-1,-1,-1,-1,-1,-1] / 2;
\\ basis of L: L1 = s, L2 = e0+e1, L3 = e0-e1, L_{k+2} = e0+e_k (k = 2..6)
{
  PL = [s, [1,1,0,0,0,0,0,0], [1,-1,0,0,0,0,0,0], [1,0,1,0,0,0,0,0], [1,0,0,1,0,0,0,0],
        [1,0,0,0,1,0,0,0], [1,0,0,0,0,1,0,0], [1,0,0,0,0,0,1,0]];
}
rowsmat(v) = matrix(#v, #v[1], i, j, v[i][j]);
BL = rowsmat(PL);
BLs  = rowsmat(vector(8, k, omul(PL[k], s)));     \\ L s    = { l . s }
BLsb = rowsmat(vector(8, k, omul(PL[k], sb)));    \\ L sbar = { l . sbar }
BLsI = BLs^(-1); BLsbI = BLsb^(-1);
inLat(x, Minv) = #select(t -> !isint(t), x * Minv) == 0;
inLam(u) =
{
  my(x = u[1..8], y = u[9..16], z = u[17..24]);
  inL(x) && inL(y) && inL(z) && inLat(x+y, BLsbI) && inLat(x+z, BLsbI)
         && inLat(y+z, BLsbI) && inLat(x+y+z, BLsI);
}
print("basis of L valid (all in L, |det| = 1): ", #select(v -> !inL(v), PL) == 0 && abs(matdet(BL)) == 1);
print("[L : L s] = ", abs(matdet(BLs)), ",  [L : L sbar] = ", abs(matdet(BLsb)));

\\ ---- twists and the triple product
R = r2 / 2;
refl(x) = x - (x * R~) * R;                       \\ s_r, an involution
pfwd(P, x) = {my(y = vector(8)); y[1] = x[1]; for(a = 1, 7, y[P[a]+1] = x[a+1]); y}   \\ pi(x)
pinv(P, y) = {my(x = vector(8)); x[1] = y[1]; for(a = 1, 7, x[a+1] = y[P[a]+1]); x}   \\ pi^{-1}(y)
mkperm(P) = [x -> pfwd(P, x), y -> pinv(P, y)];
triple(u, v, g, gi) =
{
  my(m = (a, b) -> gi(omul(g(a), g(b))),
     x = u[1..8], y = u[9..16], z = u[17..24], xp = v[1..8], yp = v[9..16], zp = v[17..24]);
  concat([m(x,xp) + m(z,yp) + m(y,zp), m(y,yp) + m(x,zp) + m(z,xp), m(z,zp) + m(y,xp) + m(x,yp)]);
}

\\ ---- [L^3 : Lambda] = 2^12, computed here independently:
\\ Lambda is the kernel of phi: L^3 -> L/Lsbar x L/Lsbar x L/Ls,
\\ (x,y,z) -> (x+y, x+z, x+y+z) (the condition y+z in Lsbar follows from the
\\ other two because 2L <= Lsbar).  Since 2L lies in both sublattices, the
\\ class of v in L/M is encoded injectively by 2*frac(v*M^-1) mod 2, and
\\ |image(phi)| = 2^(F2-rank of the images of the 24 basis vectors of L^3).
cls(v, Minv) = apply(t -> (2 * frac(t)) % 2, v * Minv);
{
  my(Z8 = vector(8), imgs = matrix(24, 24), k = 0, xs, rk);
  for(blk = 1, 3, for(i = 1, 8,
    my(x = Z8, y = Z8, z = Z8, row);
    if(blk == 1, x = PL[i], if(blk == 2, y = PL[i], z = PL[i]));
    row = concat([cls(x+y, BLsbI), cls(x+z, BLsbI), cls(x+y+z, BLsI)]);
    k++; for(j = 1, 24, imgs[k,j] = row[j])));
  rk = matrank(Mod(imgs, 2));
  print("[L^3 : Lambda] = 2^", rk, "  (each of L/Lsbar, L/Ls has order 16 = 2^4; image rank ", rk, " of 12)");
}

\\ ---- the Z-basis of Lambda from the certificate, re-validated
\\ (24 vectors of Lambda with Gram determinant 2^24 = ([L^3:Lambda] * covol(L^3))^2
\\  span a sublattice of index 1, i.e. form a Z-basis of Lambda)
B = BASIS2 / 2;
BM = rowsmat(B); BI = BM^(-1);
{
  print("certificate basis: all 24 rows in Lambda: ", #select(u -> !inLam(u), B) == 0,
        ",  det Gram = 2^", valuation(matdet(BM * BM~), 2), " (exact: ", matdet(BM * BM~) == 2^24, ")");
}

\\ ---- span data for a twist: closure and elementary divisors
spandata(g, gi, label) =
{
  my(C = matrix(576, 24), k = 0, closes = 1, defin = 1, H, d);
  for(a = 1, 24, for(c = 1, 24,
    my(p = triple(B[a], B[c], g, gi), co = p * BI);
    k++; for(j = 1, 24, C[k,j] = co[j]);
    if(#select(t -> !isint(t), co), closes = 0);
    if(!inLam(p), defin = 0)));
  if(!closes,
    print(label, ": does not close (Wilson-definition test agrees: ", defin == closes, ")"); return(0));
  H = mathnf(C~);                \\ columns = the 576 coordinate vectors
  d = matsnf(H);
  print(label, ": closes (Wilson-definition test agrees: ", defin == closes, "); #HNF columns = ", #H,
        "; elementary divisors: ", d, "; index = 2^", valuation(vecprod(d), 2));
  d;
}

idt = x -> x;
spandata(idt, idt, "untwisted");
SIG = mkperm([2,1,3,4,5,6,7]);   spandata(SIG[1], SIG[2], "sigma = (1 2)");
C567 = mkperm([1,2,3,4,6,7,5]);  spandata(C567[1], C567[2], "pi: e5->e6->e7->e5");
C576 = mkperm([1,2,3,4,7,5,6]);  spandata(C576[1], C576[2], "pi: e5->e7->e6->e5");
spandata(refl, refl, "reflection s_r, r = (e0+..+e6-e7)/2");

\\ ---- the certificate itself
checkset(pairs, prods2, detc, label) =
{
  my(P = vector(24, k, triple(B[pairs[k][1]], B[pairs[k][2]], refl, refl)), CM, ok1, ok2);
  ok1 = (P == prods2 / 2);
  ok2 = (#select(p -> !inLam(p), P) == 0);
  CM = rowsmat(vector(24, k, P[k] * BI));
  print(label, ": products match the listed ones: ", ok1, ", all in Lambda: ", ok2,
        ", det(coordinates) = ", matdet(CM), " (listed ", detc, "), det Gram = 2^24 * det^2: ",
        matdet(rowsmat(P) * rowsmat(P)~) == 2^24 * detc^2);
  matdet(CM) == detc && ok1 && ok2;
}
okA = checkset(PAIRS_A, PRODUCTS_A2, DET_A, "set A");
okB = checkset(PAIRS_B, PRODUCTS_B2, DET_B, "set B");
{
  print("DET_A = 2^", valuation(DET_A, 2), " * ", DET_A / 2^valuation(DET_A, 2), ";  DET_B odd: ", DET_B % 2 == 1,
        ";  gcd(DET_A, DET_B) = ", gcd(DET_A, DET_B));
}
print("CERTIFICATE VALID (S = Lambda for the reflection twist): ", okA && okB && gcd(DET_A, DET_B) == 1);
quit;
