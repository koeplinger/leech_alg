\\ =========================================================================
\\ c3_variant_min.gp -- the mirror convention: Wilson's construction with
\\ RIGHT multiples sbar.L = {sbar.l}, s.L = {s.l} instead of the paper's
\\ left multiples.  (Equivalently: the paper's construction with the
\\ opposite octonion table, since u *_op v = v * u for the triple product.)
\\ Checks that the resulting lattice Lambda_R is again a Leech lattice at
\\ this scale (det 2^24, minimum 8, 196560 minimal vectors).  The closure
\\ and span-index census of the 120 reflections for Lambda_R is in
\\ c3_refute_analysis.py, section 7 (all 56 closing ones have index 2^10;
\\ sigma has index 2^10 there; the C3 witness does not close).
\\ Run: gp -q c3_variant_min.gp > c3_variant_min.gp.out
\\ =========================================================================
default(parisizemax, 2*10^9);
read("c3_core.gp");
ABR = matconcat(vector(8,k, omul(sbar, LB[,k])));    \\ sbar . L
BBR = matconcat(vector(8,k, omul(s, LB[,k])));       \\ s . L
{
printf("[L : sbar.L] = %s, [L : s.L] = %s ; sbar.L = L.sbar: %d ; s.L = L.s: %d\n",
  IL(ABR), IL(BBR), mathnf(matsolve(LB,ABR))==mathnf(matsolve(LB,AB0)), mathnf(matsolve(LB,BBR))==mathnf(matsolve(LB,BB0)));
}
QAR = 2*matsolve(LB,ABR)^(-1); QBR = 2*matsolve(LB,BBR)^(-1);
KR = matconcat([QAR, QAR, Z8; QAR, Z8, QAR; QBR, QBR, QBR]);
LamR = LB3*mathnf(matconcat([lift(matkermod(KR,2)), 2*matid(24)]));
GR = LamR~*LamR;
{
printf("Lambda_R: Gram det %s, vectors of norm <= 7: %d, vectors of norm 8: %d\n",
  matdet(GR), qfminim(GR,7,0,2)[1], qfminim(GR,8,0,2)[1]);
my(Cc = matsolve(LamB, LamR));
printf("Lambda_R = Lambda (the paper's lattice): %d\n", denominator(Cc)==1 && abs(matdet(Cc))==1);
}
