\\ f6_pari_independent.gp -- an independent second implementation in PARI/GP
\\ (no Python data used except where stated in part 5 and 6).
\\
\\ Conventions (as in leech_final.py, restated):
\\  * octonions: basis e0..e7, e0 = identity; Fano triples (1,2,4),(2,3,5),(3,4,6),
\\    (4,5,7),(5,6,1),(6,7,2),(7,1,3); e_a e_b = e_c cyclically, anticyclic negative,
\\    e_i^2 = -e0.  Vectors are COLUMN vectors of ACTUAL (rational) coordinates.
\\  * L = D8+ (integer coordinates with even sum, or half-integers with ODD sum);
\\    its basis here is the HNF of the 240 roots (different from the Python basis).
\\  * Ls = {l.s}, Lsbar = {l.sbar} (left multiples); Lambda by Wilson's 3 conditions,
\\    basis = HNF of the generators (w,w,w), (v,v,0), (v,0,v), (2L)^3.
\\  * twist g: x ._g y = g^{-1}(g(x).g(y)); permutation pi = [pi(1),...,pi(7)] acts by
\\    g(e_a) = e_{pi(a)} (convention A, the paper's); reflection s_r(x) = x - <x,r> r.
\\  * triple product: P1 = xx'+zy'+yz', P2 = yy'+xz'+zx', P3 = zz'+yx'+xy' (u = left).
\\
\\ Parts: 1. build L, Ls, Lsbar, Lambda; Gram det and qfminim (minimum 8, 196560).
\\        2. closure + index (mathnf, matsnf) for identity, sigma, the 3-cycle and
\\           its inverse, and the C3 witness reflection.
\\        3. all 120 root reflections.
\\        4. index formula factors [Ls:<Ls._g Ls>], [Lsbar:<Lsbar._g Lsbar>].
\\        5. check of the Python certificate cert_index1_final.gp from scratch.
\\        6. mathnf/matsnf of the Python coordinate matrices f2_matrices.gp.

default(parisizemax, 4000000000);
FANO = [[1,2,4],[2,3,5],[3,4,6],[4,5,7],[5,6,1],[6,7,2],[7,1,3]];
TS = matrix(8,8); TC = matrix(8,8);
for(a=0,7, TS[1,a+1]=1; TC[1,a+1]=a; TS[a+1,1]=1; TC[a+1,1]=a);
for(a=1,7, TS[a+1,a+1]=-1; TC[a+1,a+1]=0);
{
for(t=1,7, my(a=FANO[t][1], b=FANO[t][2], c=FANO[t][3]);
  foreach([[a,b,c],[b,c,a],[c,a,b]], u,
    TS[u[1]+1,u[2]+1] = 1;  TC[u[1]+1,u[2]+1] = u[3];
    TS[u[2]+1,u[1]+1] = -1; TC[u[2]+1,u[1]+1] = u[3]));
}
omul(x,y) = my(z = vectorv(8)); for(a=1,8, if(x[a], for(b=1,8, if(y[b], z[TC[a,b]+1] += TS[a,b]*x[a]*y[b])))); z;

inL(x) = {
  my(d = 2*x, ev);
  if (denominator(d) != 1, return(0));
  ev = sum(i=1,8, d[i]%2==0);
  if (ev == 8, return(denominator(x) == 1 && vecsum(Vec(x)) % 2 == 0));
  if (ev == 0, return(vecsum(Vec(x)) % 2 == 1));
  0
};
inlat(B, x) = denominator(matsolve(B, x)) == 1;

\\ ---------------------------------------------------------------- 1. lattices
ROOTS = List();
for(i=1,8, for(j=i+1,8, foreach([[1,1],[1,-1],[-1,1],[-1,-1]], sg, my(v=vectorv(8)); v[i]=sg[1]; v[j]=sg[2]; listput(ROOTS, v))));
forvec(sg = vector(8,i,[0,1]), my(v = vectorv(8,i,(-1)^sg[i]/2)); if(inL(v), listput(ROOTS, v)));
ROOTS = Vec(ROOTS);
print("1. roots of L: ", #ROOTS, ";  all in L: ", prod(i=1,#ROOTS, inL(ROOTS[i])));
LB = mathnf(matconcat(vector(#ROOTS, i, 2*ROOTS[i]))) / 2;     \\ columns: basis of L
print("   det(L basis) = ", matdet(LB), "  (covolume of D8+ is 1)");
s = vectorv(8, i, if(i==1, -1/2, 1/2));  sb = vectorv(8, i, -1/2);
LSB = matconcat(vector(8, k, omul(LB[,k], sb)));   \\ basis of Lsbar
LS  = matconcat(vector(8, k, omul(LB[,k], s)));    \\ basis of Ls
print("   [L:Lsbar] = ", abs(matdet(LSB)/matdet(LB)), "   [L:Ls] = ", abs(matdet(LS)/matdet(LB)));
blk(a,b,c) = concat(concat(a,b),c);
Z8 = vectorv(8);
GEN = List();
{
for(k=1,8, listput(GEN, blk(LS[,k],LS[,k],LS[,k])); listput(GEN, blk(LSB[,k],LSB[,k],Z8)); listput(GEN, blk(LSB[,k],Z8,LSB[,k]));
           listput(GEN, blk(2*LB[,k],Z8,Z8)); listput(GEN, blk(Z8,2*LB[,k],Z8)); listput(GEN, blk(Z8,Z8,2*LB[,k])));
}
LAM = mathnf(2*matconcat(Vec(GEN))) / 2;           \\ columns: basis of Lambda
inLam(u) = {
  my(x = u[1..8], y = u[9..16], z = u[17..24]);
  inL(x) && inL(y) && inL(z) && inlat(LSB, x+y) && inlat(LSB, x+z) && inlat(LSB, y+z) && inlat(LS, x+y+z)
};
print("   Lambda basis: 24 columns, all satisfy Wilson's definition: ", prod(k=1,24, inLam(LAM[,k])));
print("   [L^3 : Lambda] = ", abs(matdet(LAM) / matdet(LB)^3));
G = LAM~ * LAM;
print("   Gram det = ", matdet(G), " = 2^", valuation(matdet(G), 2));
Gr = qflllgram(G); G2 = Gr~ * G * Gr;
q7 = qfminim(G2, 7, 0)[1]; q8 = qfminim(G2, 8, 0)[1];
cal = qfminim(LB~*LB, 2, 0)[1];
print("   qfminim calibration: E8 = L at norm <= 2 gives count ", cal, " (so the count includes both v and -v)");
print("   qfminim on Lambda: vectors of norm <= 7: ", q7, ";  of norm <= 8: ", q8, "  -> minimum 8 with ", q8, " minimal vectors");

\\ ---------------------------------------------------------------- 2. twists
permg(pi) = (x -> my(y = vectorv(8)); y[1] = x[1]; for(a=1,7, y[pi[a]+1] = x[a+1]); y);
permginv(pi) = (x -> my(y = vectorv(8)); y[1] = x[1]; for(a=1,7, y[a+1] = x[pi[a]+1]); y);
reflg(r) = (x -> x - (x~*r) * r);
\\ twisted constants K[a,b] = g^{-1}(g e_a . g e_b), so x ._g y = sum x_a y_b K[a,b]
twK(g, gi) = { my(E = matid(8), K = vector(8, a, vector(8)));
  for(a=1,8, for(b=1,8, K[a][b] = gi(omul(g(E[,a]), g(E[,b]))))); K };
tw(K, x, y) = my(z = vectorv(8)); for(a=1,8, if(x[a], for(b=1,8, if(y[b], z += x[a]*y[b]*K[a][b])))); z;
triple(K, u, v) = {
  my(x=u[1..8], y=u[9..16], z=u[17..24], X=v[1..8], Y=v[9..16], W=v[17..24]);
  blk(tw(K,x,X) + tw(K,z,Y) + tw(K,y,W), tw(K,y,Y) + tw(K,x,W) + tw(K,z,X), tw(K,z,W) + tw(K,y,X) + tw(K,x,Y))
};
\\ returns [number of failing basis products, index, elementary divisors]
analyse(K) = {
  my(C = List(), nf = 0, c, H);
  for(p=1,24, for(q=1,24,
     my(P = triple(K, LAM[,p], LAM[,q]));
     if(!inLam(P), nf++, c = matsolve(LAM, P); if(denominator(c) != 1, error("inLam/matsolve mismatch")); listput(C, c))));
  if(nf, return([nf, 0, []]));
  H = mathnf(matconcat(Vec(C)));
  [0, abs(matdet(H)), select(d -> d > 1, matsnf(H))]
};
{
TW = [["identity", [1,2,3,4,5,6,7], 0], ["sigma = (1 2)", [2,1,3,4,5,6,7], 0],
      ["3-cycle e5->e6->e7 (C2 witness)", [1,2,3,4,6,7,5], 0], ["inverse 3-cycle e5->e7->e6", [1,2,3,4,7,5,6], 0],
      ["reflection r=(e0+..+e6-e7)/2 (C3 witness)", [1,1,1,1,1,1,1,-1]~/2, 1]];
}
KW = vector(#TW);
print("2. witnesses (failing products, [Lambda:S], elementary divisors > 1):");
{
for(i=1,#TW, my(t = TW[i], K);
  K = if(t[3], twK(reflg(t[2]), reflg(t[2])), twK(permg(t[2]), permginv(t[2])));
  KW[i] = K;
  my(res = analyse(K));
  print("   ", t[1], ": failing = ", res[1], if(res[1], "  (does not close)",
        Str("  index = ", res[2], " = 2^", valuation(res[2],2), "  divisors: ", res[3]))));
}

\\ ---------------------------------------------------------------- 3. reflections
print("3. all 120 root reflections (one per +-pair, first nonzero coordinate > 0):");
cnt = Map();
{
for(i=1,#ROOTS, my(r = ROOTS[i], f = 0, res, key);
  for(j=1,8, if(r[j], f = r[j]; break));
  if(f > 0,
    res = analyse(twK(reflg(r), reflg(r)));
    key = if(res[1], "no closure", Str("index ", res[2]));
    mapput(cnt, key, if(mapisdefined(cnt, key), mapget(cnt, key), 0) + 1)));
print("   ", Mat(cnt));
}

\\ ---------------------------------------------------------------- 4. 8-dim factors
\\ h(X) = [X : <X ._g X>] for X = Ls, Lsbar
hfac(K, X) = {
  my(P = List(), H);
  for(i=1,8, for(j=1,8, listput(P, tw(K, X[,i], X[,j]))));
  H = mathnf(2*matconcat(Vec(P))) / 2;
  if(matrank(H) < 8, return(0));
  for(k=1,#P, if(!inlat(X, P[k]), return(-1)));
  abs(matdet(H) / matdet(X))
};
print("4. index formula factors [Ls:<Ls._g Ls>], [Lsbar:<Lsbar._g Lsbar>] and predicted index:");
{
for(i=2,#TW, if(i != 4, my(a = hfac(KW[i], LS), b = hfac(KW[i], LSB)); print("   ", TW[i][1], ": ", a, ", ", b, "  -> ", a*b^2)));
}

\\ ---------------------------------------------------------------- 5. Python certificate
print("5. checking cert_index1_final.gp (written by f3_certificates.py) from scratch:");
read("cert_index1_final.gp");
BAS = BASIS2~ / 2;                                       \\ columns, actual coordinates
{
print("   its 24 basis vectors lie in Lambda: ", prod(k=1,24, inLam(BAS[,k])),
      ";  det of their coordinates in this script's Lambda basis: ", matdet(matsolve(LAM, BAS)));
}
K5 = KW[5];
{
my(g = 0);
for(t=1,NSETS, my(PA = eval(Str("PAIRS", t)), DT = eval(Str("DET", t)), M);
  M = matconcat(vector(24, k, my(P = triple(K5, BAS[,PA[k,1]], BAS[,PA[k,2]])); if(!inLam(P), error("product not in Lambda")); matsolve(BAS, P)));
  print("   set ", t, ": 24 products recomputed, all in Lambda; det = ", matdet(M), "  (listed ", DT, ")");
  g = gcd(g, matdet(M)));
print("   gcd of the determinants = ", g, "  =>  [Lambda : S] = 1");
}

\\ ---------------------------------------------------------------- 6. Python matrices
print("6. PARI mathnf/matsnf of the Python coordinate matrices (f2_matrices.gp):");
read("f2_matrices.gp");
{
foreach([["sigma", t1], ["3-cycle", t3], ["reflection", t5]], m,
  my(H = mathnf(m[2]~)); print("   ", m[1], ": index ", abs(matdet(H)), ", divisors > 1: ", select(d -> d > 1, matsnf(H))));
}
