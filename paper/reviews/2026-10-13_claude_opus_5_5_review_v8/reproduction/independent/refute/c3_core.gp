\\ =========================================================================
\\ c3_core.gp -- shared PARI/GP definitions (no output) for the scripts
\\ c3_reflection_census.gp, c3_mult_dets.gp, c3_certificate.gp.
\\ c3_witness.gp is self-contained and does not use this file.
\\ The caller must set default(parisizemax, ...) BEFORE read("c3_core.gp")
\\ (changing the stack inside read() aborts the read).
\\
\\ CONVENTIONS (main.tex v8), identical to c3_witness.gp:
\\  * e_a e_b = e_c for Fano triples (1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),
\\    (6,7,2),(7,1,3) and cyclic shifts, reversed products negative, e_i^2=-e0;
\\    column vector entry k+1 = coefficient of e_k.
\\  * L = D8+ with the ODD-sum half-integer coset; s=(-e0+e1+..+e7)/2,
\\    sbar=(-e0-..-e7)/2; L sbar = {l.sbar}, L s = {l.s} (LEFT multiples).
\\  * Lambda = {(x,y,z) in L^3: x+y,x+z,y+z in L sbar; x+y+z in L s}.
\\  * x ._g y = g^{-1}(g(x).g(y)); s_r(x) = x - <x,r> r.
\\  * P1 = x.x'+z.y'+y.z', P2 = y.y'+x.z'+z.x', P3 = z.z'+y.x'+x.y'.
\\ =========================================================================
fano = [[1,2,4],[2,3,5],[3,4,6],[4,5,7],[5,6,1],[6,7,2],[7,1,3]];
SG = matrix(8,8); IX = matrix(8,8);
for(a=0,7, SG[1,a+1]=1; IX[1,a+1]=a; SG[a+1,1]=1; IX[a+1,1]=a);
for(a=1,7, SG[a+1,a+1]=-1; IX[a+1,a+1]=0);
{
for(t=1,7, my(a=fano[t][1], b=fano[t][2], c=fano[t][3]);
  SG[a+1,b+1]= 1; IX[a+1,b+1]=c;  SG[b+1,c+1]= 1; IX[b+1,c+1]=a;
  SG[c+1,a+1]= 1; IX[c+1,a+1]=b;
  SG[b+1,a+1]=-1; IX[b+1,a+1]=c;  SG[c+1,b+1]=-1; IX[c+1,b+1]=a;
  SG[a+1,c+1]=-1; IX[a+1,c+1]=b);
}
omul(x,y) = {
  my(z=vectorv(8));
  for(a=1,8, if(x[a], for(b=1,8, if(y[b], z[IX[a,b]+1] += SG[a,b]*x[a]*y[b]))));
  z;
}
E(k) = { my(v=vectorv(8)); v[k+1]=1; v; }
nrm(x) = sum(i=1,#x, x[i]^2);
inL(x) = {
  my(d=denominator(x));
  if(d==1, return(sum(i=1,8,x[i])%2==0));
  if(d==2, for(i=1,8, if(denominator(x[i])!=2, return(0)));
           return( (2*sum(i=1,8,x[i]))%4==2 ));
  0;
}
inlat(B, x) = denominator(matsolve(B,x)) == 1;
s    = vectorv(8,i, if(i==1,-1/2, 1/2));
sbar = vectorv(8,i, -1/2);
LB = matconcat([s, E(0)+E(1), E(0)-E(1), E(0)+E(2), E(0)+E(3), E(0)+E(4), E(0)+E(5), E(0)+E(6)]);
AB0 = matconcat(vector(8,k, omul(LB[,k], sbar)));   \\ L sbar
BB0 = matconcat(vector(8,k, omul(LB[,k], s)));      \\ L s
IL(B) = abs(matdet(matsolve(LB,B)));
PA = matsolve(LB, AB0); PB = matsolve(LB, BB0);
QA = 2*PA^(-1); QB = 2*PB^(-1);
Z8 = matrix(8,8);
K = matconcat([QA, QA, Z8; QA, Z8, QA; QB, QB, QB]);
KER = matkermod(K, 2);
LamC = mathnf(matconcat([lift(KER), 2*matid(24)]));
LB3 = matconcat([LB, Z8, Z8; Z8, LB, Z8; Z8, Z8, LB]);
LamB = LB3 * LamC;
GramLam = LamB~ * LamB;
LamBi = LamB^(-1);
inLambda(P) = {
  my(x=P[1..8], y=P[9..16], z=P[17..24]);
  inL(x) && inL(y) && inL(z) && inlat(AB0,x+y) && inlat(AB0,x+z) && inlat(AB0,y+z) && inlat(BB0,x+y+z);
}
tw(G, Gi, x, y) = Gi * omul(G*x, G*y);
tprod(G, Gi, u, v) = {
  my(x=u[1..8], y=u[9..16], z=u[17..24], xp=v[1..8], yp=v[9..16], zp=v[17..24]);
  concat([ tw(G,Gi,x,xp) + tw(G,Gi,z,yp) + tw(G,Gi,y,zp),
           tw(G,Gi,y,yp) + tw(G,Gi,x,zp) + tw(G,Gi,z,xp),
           tw(G,Gi,z,zp) + tw(G,Gi,y,xp) + tw(G,Gi,x,yp) ]);
}
\\ [closes?, #failing basis pairs, 24x576 coordinate matrix]; column (i-1)*24+j = coords(b_i *_g b_j)
products(G) = {
  my(Gi=G^(-1), C=matrix(24,576), bad=0, P, c);
  for(i=1,24, for(j=1,24,
    P = tprod(G, Gi, LamB[,i], LamB[,j]);
    if(!inLambda(P), bad++; next);
    c = LamBi * P;
    if(denominator(c)!=1, error("Wilson test and Lambda basis disagree"));
    C[,(i-1)*24+j] = c));
  [bad==0, bad, C];
}
spanindex(C) = abs(matdet(mathnf(C)));
reflmat(r) = matid(8) - r*r~;
