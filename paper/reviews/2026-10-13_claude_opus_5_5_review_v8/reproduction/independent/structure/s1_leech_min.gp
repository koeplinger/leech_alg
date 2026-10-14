\\ s1_leech_min.gp -- PARI/GP check of the Gram matrix written by s1_basics.py
\\ (Wilson's Lambda at the paper's scale: E8 roots of norm 2).
\\ Expected: det 2^24, minimum 8, exactly 196560 vectors of norm 8.
\\ qfminim with flag 0 (exact integral Fincke-Pohst); the count includes both v and -v.
read("lambda_gram.gp");
print("det(Gram) = ", matdet(G), " = 2^", valuation(matdet(G), 2));
H = G/2;
print("Gram/2 integral: ", denominator(H) == 1, "   det(Gram/2) = ", matdet(H));
print("diagonal of Gram/2 all even: ", vecsum(vector(24, i, H[i,i] % 2)) == 0);
r = qfminim(G, 7, 0);
print("vectors of norm <= 7 (excluding 0): ", r[1]);
r = qfminim(G, 8, 0);
print("vectors of norm <= 8 (excluding 0): ", r[1], "   max norm found: ", r[2]);
