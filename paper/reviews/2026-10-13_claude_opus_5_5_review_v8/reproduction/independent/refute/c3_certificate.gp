\\ =========================================================================
\\ c3_certificate.gp -- 24 single products u_k *_g v_k that form a Z-basis
\\ of Lambda, for the twist g = s_r, r = (e0+...+e6-e7)/2 (claim C3's
\\ witness).  Such a set is a compact certificate of S_g = Lambda: it is
\\ checked with 24 products, Wilson's membership test, and one determinant
\\ (or, basis-free, one Gram determinant = 2^24 = det Lambda).
\\
\\ Why a construction is needed: greedy/random selection among the 576
\\ basis-pair products always stops at 23 primitive vectors (1,001 seeded
\\ runs with the HNF basis, 300 with an LLL basis of minimal vectors).  The
\\ 24th vector would have to make the determinant exactly +-1, and products
\\ (norm >= 16, about 64 typically) are far longer than a typical basis vector
\\ of Lambda (det Gram = 2^24), so a random product essentially never does.
\\
\\ Construction:
\\  * For a factor u, P_u = L_u(Lambda) = {u *_g v : v in Lambda} (and
\\    R_u(Lambda) = {v *_g u}) is a full-rank sublattice ALL of whose elements
\\    are single products.  The pool consists of the LLL-reduced bases of P_u
\\    for u in a fixed list of factors (the 24 basis vectors and the first
\\    200 minimal vectors in qfminim's deterministic order), both sides.
\\  * Greedy primitive selection over the pool (fixed seeds): keep a pool
\\    vector if the chosen set stays primitive in Lambda.
\\  * Last vector(s): with the chosen set E primitive of rank 23, let phi be
\\    the primitive functional vanishing on E.  For a factor u with
\\    gcd(phi(P_u)) = 1, Bezout gives v with phi(u * v) = 1, and v is reduced
\\    modulo {w : phi(u*w) = 0} by LLL + Babai rounding (this changes u*v only
\\    by elements of span(E), so the determinant stays +-1).
\\  * Finally the 24 products are RECOMPUTED from the printed factors with
\\    the octonion code, tested against Wilson's conditions, and their
\\    determinant and Gram determinant are evaluated.
\\ Run: gp -q c3_certificate.gp > c3_certificate.gp.out
\\ Conventions: see c3_core.gp (identical to main.tex v8).
\\ =========================================================================
default(parisizemax, 4*10^9);
read("c3_core.gp");

r = vectorv(8,i,1/2); r[8] = -1/2;
G = reflmat(r); Gi = G^(-1);
res = products(G);
if(!res[1], error("witness does not close"));
C = res[3];
printf("witness closes; [Lambda:S] = %s\n", spanindex(C));
mkT(C) = {
  my(TL = matrix(576,24), TR = matrix(576,24));
  for(i=1,24, for(j=1,24, for(k=1,24,
    TL[(j-1)*24+k, i] = C[k, (i-1)*24+j];      \\ L_{b_i} column j = coords(b_i * b_j)
    TR[(j-1)*24+k, i] = C[k, (j-1)*24+i])));   \\ R_{b_i} column j = coords(b_j * b_i)
  [TL, TR];
}
reshape(w) = matconcat(vector(24, j, w[(j-1)*24+1 .. j*24]));
T = mkT(C); TL = T[1]; TR = T[2];
Mop(side, c) = if(side==1, reshape(TL*c), reshape(TR*c));   \\ side 1: v -> u*v ; side 2: v -> v*u

mv = qfminim(GramLam, 8, , 2)[3];
NMIN = 200;
facs = concat(vector(24,i, matid(24)[,i]), vector(NMIN, n, mv[,n]));

\\ pool of short single products: entries [side, factor index, beta (coords of v), product coords]
pool = List();
{
for(n=1, #facs, for(side=1,2,
  my(M = Mop(side, facs[n]), U);
  if(matdet(M)==0, next);
  U = qflllgram(M~*GramLam*M);
  for(j=1,24, listput(pool, [side, n, U[,j], M*U[,j]]))));
}
printf("pool: %d short single products from %d factors x 2 sides\n", #pool, #facs);

primitive(M) = { my(d=matsnf(M)); for(i=1,#d, if(d[i]!=0 && d[i]!=1, return(0))); matrank(M)==#M; }
shuffled(n) = { my(v=vector(n,i,i), k, tmp); forstep(i=n,2,-1, k=random(i)+1; tmp=v[i]; v[i]=v[k]; v[k]=tmp); v; }
\\ Babai rounding of beta modulo the integer lattice spanned by the columns of Kb, in the metric of p = M*beta
babai(beta, Kb, M) = {
  if(#Kb==0, return(beta));
  my(Kr = Kb*qflll(M*Kb), Gr, t);                 \\ LLL in product coordinates, exact integers
  Gr = (M*Kr)~*GramLam*(M*Kr);
  t = round(Gr^(-1) * ((M*Kr)~*GramLam*(M*beta)));
  beta - Kr*t;
}

{
my(found = 0, seed = 0, sel, Ecols);
while(!found && seed < 50,
  seed++; setrand(seed);
  sel = List(); Ecols = matrix(24,0);
  my(ord = shuffled(#pool));
  for(t=1, #ord, if(#sel==23, break);
    my(cand = pool[ord[t]][4], Mn = matconcat([Ecols, cand]));
    if(matrank(Mn)==#sel+1 && primitive(Mn), Ecols = Mn; listput(sel, pool[ord[t]])));
  if(#sel < 23, next);
  \\ last vector by Bezout on phi(P_u)
  my(phi = matker(Ecols~)[,1]); phi = phi/content(phi);
  for(n=1, #facs, for(side=1,2,
    my(M = Mop(side, facs[n]), w = phi~*M, H, beta, Kb);
    if(content(w)!=1, next);
    H = mathnf(w, 1);                       \\ w * H[2] = [0,...,0,1]
    beta = H[2][,24]; Kb = H[2][,1..23];     \\ phi(u * Kb) = 0
    if(w*beta != 1, error("bezout"));
    beta = babai(beta, Kb, M);
    listput(sel, [side, n, beta, M*beta]);
    found = 1; break(2)));
  printf("seed %d: greedy reached 23, completion %s\n", seed, if(found, "found", "failed"));
  if(found, CERT = sel));
if(!found, error("no certificate found"));
}

\\ ---------- independent re-verification from the printed factors -----------
{
my(P = matrix(24,24), ok = 1, maxv = 0, side, u, v, p);
printf("\nCERTIFICATE: 24 products p_k = (left factor) *_g (right factor); vectors in R^24 as (x|y|z):\n");
for(k=1,24,
  side = CERT[k][1]; u = LamB*facs[CERT[k][2]]; v = LamB*CERT[k][3];
  if(!inLambda(u) || !inLambda(v), ok = 0);
  p = if(side==1, tprod(G, Gi, u, v), tprod(G, Gi, v, u));
  if(!inLambda(p), ok = 0);
  P[,k] = LamBi*p;
  maxv = max(maxv, nrm(v));
  printf("k=%2d  left  = %s\n      right = %s\n      norms: left %s, right %s, product %s\n",
     k, if(side==1,u,v)~, if(side==1,v,u)~, nrm(if(side==1,u,v)), nrm(if(side==1,v,u)), nrm(p)));
printf("\nall factors and all 24 products pass Wilson's test: %d\n", ok);
printf("det of the 24 products in Lambda-coordinates = %s\n", matdet(P));
printf("Gram determinant of the 24 product vectors = %s ; det Gram(Lambda) = %s\n", matdet(P~*GramLam*P), matdet(GramLam));
printf("largest norm of a factor: %s\n", maxv);
}
