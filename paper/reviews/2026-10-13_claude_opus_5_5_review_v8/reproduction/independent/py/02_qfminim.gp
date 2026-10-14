\\ 02_qfminim.gp -- PARI/GP check of the Gram matrix of Lambda written by
\\ 01_build_lambda.py (lambda_gram.gp; e-coordinates, sqrt(2)-scaled Leech
\\ normalization, norm = squared Euclidean length).
\\ Expected: det 2^24, all norms even, minimum 8, 196560 vectors of norm 8,
\\ no nonzero vector of norm < 8.
\\ PARI convention: qfminim(G, B, 0) returns [N, maxnorm, V], N = number of
\\ vectors of norm <= B counting v and -v separately (checked on E8 first).
\\ Run:  gp -q 02_qfminim.gp > 02_qfminim.out.txt
\\ (default stack size suffices)
\\ calibration on E8 (Cartan matrix): 240 roots expected
E8 = [2,-1,0,0,0,0,0,0;-1,2,-1,0,0,0,0,0;0,-1,2,-1,0,0,0,-1;0,0,-1,2,-1,0,0,0;0,0,0,-1,2,-1,0,0;0,0,0,0,-1,2,-1,0;0,0,0,0,0,-1,2,0;0,0,-1,0,0,0,0,2];
print("calibration: det(E8 Cartan) = ", matdet(E8), ", qfminim(E8,2,0)[1] = ", qfminim(E8,2,0)[1], " (expect 240)");
read("lambda_gram.gp");
print("Gram matrix G of Lambda: size ", matsize(G), ", symmetric: ", G == G~);
print("det G = ", matdet(G), " = 2^", valuation(matdet(G), 2));
print("all entries even: ", content(G) % 2 == 0, "   (content = ", content(G), ")");
H = G/2;
print("G/2 integral with even diagonal: ", type(H) == "t_MAT" && content(H) == 1 && vecsum(vector(24, i, H[i,i] % 2)) == 0);
U = qflllgram(G); G2 = U~ * G * U;
print("LLL change of basis unimodular: ", abs(matdet(U)) == 1);
t0 = getabstime();
r6 = qfminim(G2, 6, 0);
print("number of nonzero vectors of norm <= 6: ", r6[1]);
r8 = qfminim(G2, 8, 0);
print("number of nonzero vectors of norm <= 8: ", r8[1], "   max norm among them: ", r8[2]);
rmin = qfminim(G2, , 0);
print("minimum (qfminim without bound): ", rmin[2], "  with ", rmin[1], " minimal vectors");
print("time (ms): ", getabstime() - t0);
quit;
