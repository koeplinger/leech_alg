\\ s6_independent.gp -- independent PARI/GP reimplementation (no Python data used).
\\
\\ Conventions (as in the paper): e0 identity, Fano triples
\\ (1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3), cyclic +, anticyclic -,
\\ e_i^2 = -e0.  Vectors are COLUMN vectors of e-coordinates; GP index k+1 <-> e_k.
\\ L = D8+ (odd half-integer coset) = Z-span of {e_i +- e_j} and s.
\\ Ls = {l.s}, Lsbar = {l.sbar} (left multiples).  Lambda = Wilson's Leech lattice.
\\ Twist g (matrix on column vectors): x ._g y = g^{-1}(g x . g y);
\\ triple product P1 = x.x'+z.y'+y.z', P2 = y.y'+x.z'+z.x', P3 = z.z'+y.x'+x.y'.
\\ Index [Lambda:S] = det of the HNF of the 24 x 576 matrix of product coordinates
\\ (PARI mathnf), elementary divisors by matsnf.
\\ Output: s6_independent.out (via tee).

F = [[1,2,4],[2,3,5],[3,4,6],[4,5,7],[5,6,1],[6,7,2],[7,1,3]];
SG = matrix(8,8); IX = matrix(8,8);
for(a=0,7, SG[1,a+1]=1; IX[1,a+1]=a; SG[a+1,1]=1; IX[a+1,1]=a);
for(a=1,7, SG[a+1,a+1]=-1; IX[a+1,a+1]=0);
for(t=1,7, my(a=F[t][1], b=F[t][2], c=F[t][3]); \
  foreach([[a,b,c],[b,c,a],[c,a,b]], q, \
    SG[q[1]+1,q[2]+1]=1; IX[q[1]+1,q[2]+1]=q[3]; \
    SG[q[2]+1,q[1]+1]=-1; IX[q[2]+1,q[1]+1]=q[3]));
om(x,y) = {my(z=vectorv(8)); for(a=1,8, if(x[a], for(b=1,8, if(y[b], z[IX[a,b]+1] += SG[a,b]*x[a]*y[b])))); z};

\\ ---------------------------------------------------------------- L, Ls, Lsbar
s = vectorv(8, k, if(k==1, -1/2, 1/2));
sbar = vectorv(8, k, -1/2);
gens = List();
for(i=1,8, for(j=i+1,8, foreach([[1,1],[1,-1],[-1,1],[-1,-1]], e, \
  my(v=vectorv(8)); v[i]=e[1]; v[j]=e[2]; listput(gens, v))));
listput(gens, s);
M = matrix(8, #gens, i, j, gens[j][i]);
BL = mathnf(2*M)/2;
inlat(B, v) = denominator(matsolve(B, v)) == 1;
print("L: det(basis) = ", matdet(BL), ", Gram det = ", matdet(BL~*BL), ", e0 in L: ", inlat(BL, vectorv(8,k,k==1)));
\\ L.L in L on a basis
okLL = 1; for(i=1,8, for(j=1,8, okLL = okLL && inlat(BL, om(BL[,i], BL[,j]))));
print("L.L in L (64 basis products): ", okLL);
BLsb = mathnf(2*matrix(8,8,i,j, om(BL[,j], sbar)[i]))/2;
BLs  = mathnf(2*matrix(8,8,i,j, om(BL[,j], s)[i]))/2;
print("Lsbar in L: ", vecprod(vector(8,k, inlat(BL, BLsb[,k]))), "  Ls in L: ", vecprod(vector(8,k, inlat(BL, BLs[,k]))), \
      "  [L:Lsbar] = ", abs(matdet(BLsb)/matdet(BL)), "  [L:Ls] = ", abs(matdet(BLs)/matdet(BL)));

\\ ---------------------------------------------------------------- Lambda
Z8 = vectorv(8);
cat3(a,b,c) = concat(concat(a,b),c);
LG = List();
for(k=1,8, listput(LG, cat3(2*BL[,k], Z8, Z8)); listput(LG, cat3(Z8, 2*BL[,k], Z8)); listput(LG, cat3(Z8, Z8, 2*BL[,k])));
for(k=1,8, my(a=BLsb[,k]); listput(LG, cat3(a, Z8, a)); listput(LG, cat3(Z8, a, a)));
for(k=1,8, my(w=BLs[,k]); listput(LG, cat3(w, w, w)));
BLam = mathnf(2*matrix(24, #LG, i, j, LG[j][i]))/2;
wil(u) = {my(x=u[1..8], y=u[9..16], z=u[17..24]); \
  inlat(BL,x) && inlat(BL,y) && inlat(BL,z) && inlat(BLsb,x+y) && inlat(BLsb,x+z) && inlat(BLsb,y+z) && inlat(BLs,x+y+z)};
GLam = BLam~*BLam;
print("Lambda: 24 basis vectors satisfy Wilson's conditions: ", vecprod(vector(24,k, wil(BLam[,k]))), \
      "  Gram det = ", matdet(GLam), "  [L^3:Lambda] = ", abs(matdet(BLam)) / abs(matdet(BL))^3);
r = qfminim(GLam, 7, 0); print("vectors of norm <= 7: ", r[1]);
r = qfminim(GLam, 8, 0); print("vectors of norm <= 8: ", r[1]);
\\ (the generating set contains 2L^3, so Lambda/2L^3 has dimension 24 - log2[L^3:Lambda]; together with the
\\  12 independent mod-2 conditions this shows the HNF basis spans all of Lambda; see s1_basics.py)

\\ ---------------------------------------------------------------- twisted products
BLinv = BLam^(-1);
tw(a, b, g, gi) = gi*om(g*a, g*b);
tp(u, v, g, gi) = {my(x=u[1..8], y=u[9..16], z=u[17..24], xp=v[1..8], yp=v[9..16], zp=v[17..24]); \
  cat3(tw(x,xp,g,gi)+tw(z,yp,g,gi)+tw(y,zp,g,gi), tw(y,yp,g,gi)+tw(x,zp,g,gi)+tw(z,xp,g,gi), \
       tw(z,zp,g,gi)+tw(y,xp,g,gi)+tw(x,yp,g,gi))};
analyze(g) = {my(gi = g^(-1), C = matrix(24, 576), n = 0, nfail = 0, H);
  for(i=1,24, for(j=1,24, my(c = BLinv*tp(BLam[,i], BLam[,j], g, gi)); n++; \
     if(denominator(c) != 1, nfail++); C[,n] = c));
  if(nfail, return([0, nfail, 0, []]));
  H = mathnf(C); [1, 0, abs(matdet(H)), matsnf(H)]};
permmat(p) = {my(P = matrix(8,8)); P[1,1] = 1; for(a=1,7, P[p[a]+1, a+1] = 1); P};   \\ g(e_a) = e_{p[a]}
refl(rv) = matid(8) - rv*rv~;
show(name, res) = {if(res[1], my(ed = select(t->t>1, res[4])); \
   print(name, ": closes, [Lambda:S] = ", res[3], " = 2^", valuation(res[3],2), \
         ", Lambda/S has ", #ed, " cyclic factors ", Set(ed)), \
   print(name, ": does NOT close (", res[2], " of 576 products outside Lambda)"))};

show("identity                      ", analyze(matid(8)));
show("sigma = (1 2)                 ", analyze(permmat([2,1,3,4,5,6,7])));
show("3-cycle conv A: g(e5)=e6,g(e6)=e7,g(e7)=e5", analyze(permmat([1,2,3,4,6,7,5])));
show("3-cycle conv B: inverse       ", analyze(permmat([1,2,3,4,7,5,6])));
rw = vectorv(8, k, if(k==8, -1/2, 1/2));
show("reflection in r=(e0+...+e6-e7)/2", analyze(refl(rw)));
show("reflection in e1-e2 (= sigma) ", analyze(refl(vectorv(8,k, (k==2)-(k==3)))));

\\ all 120 root reflections (one root per +-pair: first nonzero coordinate positive)
roots = List();
for(i=1,8, for(j=i+1,8, foreach([[1,1],[1,-1]], e, my(v=vectorv(8)); v[i]=e[1]; v[j]=e[2]; listput(roots, v))));
forvec(e = vector(7, k, [0,1]), my(sg = concat([1], vector(7,k, 1-2*e[k]))); \
  if(#select(t->t<0, sg) % 2 == 1, listput(roots, vectorv(8, k, sg[k]/2))));
print("number of roots up to sign: ", #roots);
tally = Map();
{for(k=1, #roots, my(rv = roots[k], res = analyze(refl(rv)), key, typ);
   typ = if(denominator(rv) == 2, "half-integer", if(rv[1] != 0, "int +-e0+-ej", "int e0 not involved"));
   key = [typ, if(res[1], valuation(res[3],2), "no close")];
   mapput(tally, key, if(mapisdefined(tally, key), mapget(tally, key), 0) + 1));}
print("120 reflections: [root type, log2 index] -> count:");
M2 = Mat(tally); for(i=1, matsize(M2)[1], print("   ", M2[i,1], " : ", M2[i,2]));
