"""Item (1c), characterisation: WHICH consecutive sigma-twists close on Lambda.

The exact test (verify_consecutive_twists_exact.py) found that of the
permutations reachable as a product of two transpositions, 35 of the 70
three-cycles and 42 of the 105 (2,2)-double-transpositions close on Lambda.
This script characterises which ones, against the Fano-line structure.
"""
from itertools import combinations
from verify_consecutive_twists_exact import closes

FANO = [(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)]
lines = {frozenset(t) for t in FANO}

# --- 3-cycles: per 3-subset, the two orientations ---
print("THREE-CYCLES  (a 3-cycle = two transpositions sharing one index)")
tally = {0:0, 1:0, 2:0}
line_tally = {0:0, 1:0, 2:0}
nonline_tally = {0:0, 1:0, 2:0}
for a,b,c in combinations(range(1,8),3):
    p1 = list(range(8)); p1[a],p1[b],p1[c] = b,c,a      # (a b c)
    p2 = list(range(8)); p2[a],p2[b],p2[c] = c,a,b      # (a c b)
    n = closes(tuple(p1)) + closes(tuple(p2))
    tally[n] += 1
    if frozenset((a,b,c)) in lines: line_tally[n] += 1
    else:                           nonline_tally[n] += 1
print(f"  3-subsets by #closing orientations (of 2):  {dict(tally)}")
print(f"    among the 7 Fano-line subsets:            {dict(line_tally)}")
print(f"    among the 28 non-line subsets:            {dict(nonline_tally)}")

# --- (2,2) double transpositions ---
print()
print("(2,2) DOUBLE TRANSPOSITIONS  (two disjoint transpositions)")
cl_line = cl_nonline = 0
tot_line = tot_nonline = 0
for a,b in combinations(range(1,8),2):
    rest = [x for x in range(1,8) if x not in (a,b)]
    for c,d in combinations(rest,2):
        if min(c,d) < min(a,b):     # count each (a b)(c d) once
            continue
        p = list(range(8)); p[a],p[b] = b,a; p[c],p[d] = d,c
        fixed = frozenset(set(range(1,8)) - {a,b,c,d})   # the 3 fixed indices
        cl = closes(tuple(p))
        if fixed in lines:
            tot_line += 1;    cl_line += cl
        else:
            tot_nonline += 1; cl_nonline += cl
print(f"  fixed triple IS a Fano line:     {cl_line} of {tot_line} close")
print(f"  fixed triple is NOT a Fano line: {cl_nonline} of {tot_nonline} close")
print(f"  total: {cl_line+cl_nonline} of {tot_line+tot_nonline} close")
