\\ =========================================================================
\\ c3_reflection_census.gp -- PARI/GP census of the 120 root reflections of
\\ L = D8+ used as twists (claim C3: 64 do not close, 28 close with index
\\ 2^16, 28 close with index 1).  Independent of the Python implementation.
\\ Run: gp -q c3_reflection_census.gp > c3_reflection_census.gp.out
\\ Conventions: see c3_core.gp (identical to main.tex v8).
\\ =========================================================================
default(parisizemax, 2*10^9);
read("c3_core.gp");

\\ the 240 roots of L
roots = List();
{
for(i=1,8, for(j=i+1,8, for(a=0,1, for(b=0,1,
  my(v=vectorv(8)); v[i]=(-1)^a; v[j]=(-1)^b; listput(roots, v)))));
forvec(sg=vector(8,k,[0,1]), my(v=vectorv(8,k, (-1)^sg[k]/2)); if(inL(v), listput(roots, v)));
}
if(#roots!=240, error("root count"));
\\ one representative per +-r pair: keep r whose first nonzero coordinate is positive
reps = select(r -> my(k=1); while(r[k]==0, k++); r[k]>0, Vec(roots));
printf("roots: %d, representatives: %d\n", #roots, #reps);
rtype(r) = if(denominator(r)==2, "half-integer", if(r[1]!=0, "+-e0+-e_j", "+-e_i+-e_j (i,j>=1)"));

tab = Map(); idx1 = List(); idx16 = List();
{
for(n=1,#reps, my(r=reps[n], G=reflmat(r), res, key, c=0);
  res = products(G);
  key = if(res[1], Str(rtype(r), " | closes, index ", spanindex(res[3])), Str(rtype(r), " | no closure"));
  if(mapisdefined(tab, key, &c), mapput(tab, key, c+1), mapput(tab, key, 1));
  if(res[1] && spanindex(res[3])==1, listput(idx1, r~));
  if(res[1] && spanindex(res[3])==2^16, listput(idx16, r~)));
}
{
my(M=Mat(tab));
for(i=1,matsize(M)[1], printf("  %-58s %d\n", M[i,1], M[i,2]));
printf("total index 1: %d ; total index 2^16: %d ; total no closure: %d\n", #idx1, #idx16, #reps-#idx1-#idx16);
printf("sigma's root (0,1,-1,0,...,0) has index 2^16: %d\n", setsearch(Set(Vec(idx16)), [0,1,-1,0,0,0,0,0])>0);
printf("witness root (1,1,1,1,1,1,1,-1)/2 has index 1: %d\n", setsearch(Set(Vec(idx1)), [1,1,1,1,1,1,1,-1]/2)>0);
printf("index-1 roots (one per +-pair):\n"); for(i=1,#idx1, print("   ", idx1[i]));
printf("index-2^16 roots (one per +-pair):\n"); for(i=1,#idx16, print("   ", idx16[i]));
}
