\\ =========================================================================
\\ c3_mult_dets.gp -- is S = Lambda "trivially" explained by a unit?
\\ For the twist g = s_r, r = (e0+...+e6-e7)/2 (claim C3's witness) and,
\\ for comparison, sigma = (1 2):
\\   (a) exact det of left multiplication L_u: v -> u *_g v and right
\\       multiplication R_u: v -> v *_g u on Lambda, for ALL 196,560 minimal
\\       vectors u (98,280 up to sign; det(L_{-u}) = det(L_u) in rank 24);
\\   (b) test of the predicted divisibility N(x+y+z)^4 | det(L_u), det(R_u)
\\       for u = (x,y,z) (proof in the report: L_u is block diagonal on
\\       R^24 = D + T and acts on D = {(a,a,a)} with determinant N(x+y+z)^4);
\\   (c) one-sided identities: is L_e = I or R_e = I solvable for e in
\\       Lambda (x) Q?  (exact linear algebra over Q);
\\   (d) the multiples of the idempotents (eps,0,0) built from the identity
\\       element eps = g^{-1}(e0) of ._g.
\\ Run: gp -q c3_mult_dets.gp > c3_mult_dets.gp.out
\\ Conventions: see c3_core.gp (identical to main.tex v8).
\\ =========================================================================
default(parisizemax, 4*10^9);
read("c3_core.gp");

mv = qfminim(GramLam, 8, , 2);
printf("minimal vectors: %d (stored %d, one per +-pair)\n", mv[1], #mv[3]);
MV = mv[3];

\\ vec(L_u) = TL*c, vec(R_u) = TR*c for u = LamB*c; vec stacks columns
mkT(C) = {
  my(TL = matrix(576,24), TR = matrix(576,24));
  for(i=1,24, for(j=1,24, for(k=1,24,
    TL[(j-1)*24+k, i] = C[k, (i-1)*24+j];      \\ L_{b_i} column j = coords(b_i * b_j)
    TR[(j-1)*24+k, i] = C[k, (j-1)*24+i])));   \\ R_{b_i} column j = coords(b_j * b_i)
  [TL, TR];
}
reshape(w) = matconcat(vector(24, j, w[(j-1)*24+1 .. j*24]));
blocksum(u) = u[1..8] + u[9..16] + u[17..24];

analyse(name, G) = {
  my(res = products(G), T, TL, TR, hl = Map(), hr = Map(), unimod = 0, divfail = 0, cnt = 0,
     dl, dr, u, n4, minL = -1, minR = -1, zeroL = 0, zeroR = 0, oddL = 0, oddR = 0, v2l = Map(), v2r = Map());
  printf("\n== %s ==  closes: %d, [Lambda:S] = %s\n", name, res[1], if(res[1], spanindex(res[3]), "n/a"));
  if(!res[1], return);
  T = mkT(res[3]); TL = T[1]; TR = T[2];
  for(n=1, #MV,
    my(c = MV[,n]);
    dl = matdet(reshape(TL*c)); dr = matdet(reshape(TR*c));
    u = LamB*c; n4 = nrm(blocksum(u))^4;
    if(abs(dl)==1 || abs(dr)==1, unimod++);
    if(n4==0, if(dl!=0 || dr!=0, divfail++), if(dl % n4 || dr % n4, divfail++));
    if(dl==0, zeroL++, if(minL<0 || abs(dl)<minL, minL=abs(dl)));
    if(dr==0, zeroR++, if(minR<0 || abs(dr)<minR, minR=abs(dr)));
    if(dl, oddL = max(oddL, abs(dl) >> valuation(dl,2)));
    if(dr, oddR = max(oddR, abs(dr) >> valuation(dr,2)));
    my(a = if(dl, valuation(dl,2), -1), b = if(dr, valuation(dr,2), -1), t=0);
    if(mapisdefined(v2l, a, &t), mapput(v2l, a, t+1), mapput(v2l, a, 1));
    if(mapisdefined(v2r, b, &t), mapput(v2r, b, t+1), mapput(v2r, b, 1)));
  printf("   minimal vectors with |det L_u| = 1 or |det R_u| = 1: %d\n", unimod);
  printf("   violations of N(x+y+z)^4 | det: %d\n", divfail);
  printf("   singular L_u: %d ; singular R_u: %d (these are the u with x+y+z = 0)\n", zeroL, zeroR);
  printf("   smallest nonzero |det L_u| = %s = 2^%d ; smallest nonzero |det R_u| = %s = 2^%d\n",
     minL, valuation(minL,2), minR, valuation(minR,2));
  printf("   largest odd part of a nonzero det L_u: %d ; of a nonzero det R_u: %d (1 = every nonzero det is a power of 2)\n", oddL, oddR);
  printf("   distribution of v_2(det L_u) (key -1 = singular): %s\n", Mat(v2l)~);
  printf("   distribution of v_2(det R_u) (key -1 = singular): %s\n", Mat(v2r)~);
  \\ (c) one-sided identities in Lambda (x) Q
  my(Iv = concat(Vec(matid(24))), sl = matinverseimage(TL, Iv), sr = matinverseimage(TR, Iv));
  printf("   left identity in Lambda (x) Q exists: %d ; right identity exists: %d\n", #sl>0, #sr>0);
  \\ (d) multiples of the idempotents built from eps = g^{-1}(e0)
  my(eps = G^(-1)*E(0), Z=vectorv(8));
  printf("   eps = g^{-1}(e0) = %s ; 2eps in L: %d ; 4eps in 2L: %d\n", eps~, inL(2*eps), inL(2*eps));
  my(cands = [concat([4*eps,Z,Z]), concat([Z,4*eps,Z]), concat([Z,Z,4*eps]), concat([4*eps,4*eps,4*eps]),
              concat([2*eps,2*eps,Z]), concat([2*eps,2*eps,2*eps]), concat([2*eps,-2*eps,Z])],
     names = ["(4eps,0,0)","(0,4eps,0)","(0,0,4eps)","(4eps,4eps,4eps)","(2eps,2eps,0)","(2eps,2eps,2eps)","(2eps,-2eps,0)"]);
  for(t=1,#cands, my(w=cands[t], cw);
    if(!inLambda(w), printf("   %-18s not in Lambda\n", names[t]); next);
    cw = LamBi*w;
    printf("   %-18s in Lambda: det L = %s, det R = %s\n", names[t], factor(matdet(reshape(TL*cw))), factor(matdet(reshape(TR*cw)))));
}

r = vectorv(8,i,1/2); r[8] = -1/2;
SIG = matid(8); SIG[2,2]=0; SIG[3,3]=0; SIG[2,3]=1; SIG[3,2]=1;
analyse("s_r, r = (e0+...+e6-e7)/2", reflmat(r));
analyse("sigma = (1 2)", SIG);
