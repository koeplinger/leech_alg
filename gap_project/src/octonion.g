#############################################################################
##
##  octonion.g — Octonion algebra over Q^8 with the standard Fano-plane
##  multiplication used in the paper.
##
##  Conventions (matching the main paper, Section 2.1, and the Python module
##  python_project/src/octonions.py):
##
##    Basis     : {e_0, e_1, ..., e_7}, e_0 the identity.
##    Internal  : 1-indexed; position k holds the e_{k-1} coefficient.
##                So position 1 = e_0, position 2 = e_1, ..., position 8 = e_7.
##
##  Standard Fano triples (paper eq. (1)):
##    (1,2,4), (2,3,5), (3,4,6), (4,5,7), (5,6,1), (6,7,2), (7,1,3).
##  Encoding:  for each triple (a,b,c),  e_a*e_b = +e_c  (and its two cyclic
##  rotations);  e_b*e_a = -e_c  (and its two anti-cyclic rotations).
##
##  All arithmetic is exact (Rationals).
##
##  References
##    [Wilson2009]  R.A. Wilson, "Octonions and the Leech lattice",
##                  J. Algebra 322 (2009) 2186-2190.
##    [Dixon2010]   G.M. Dixon, "Integral Octonions, Octonion XY-Product, and
##                  the Leech Lattice", arXiv:1011.2541.
##    [Baez2002]    J.C. Baez, "The octonions",
##                  Bull. AMS 39 (2002) 145-205.
##
#############################################################################

# ---------------------------------------------------------------------------
# Fano triples (1-indexed in {1..7}).
# ---------------------------------------------------------------------------

STANDARD_FANO_TRIPLES := [ [1,2,4], [2,3,5], [3,4,6], [4,5,7],
                           [5,6,1], [6,7,2], [7,1,3] ];

# ---------------------------------------------------------------------------
# Build the 8x8 multiplication table.
#   table[i][j] = [sign, k]  meaning  e_{i-1} * e_{j-1} = sign * e_{k-1}.
# Rows/columns are 1..8 (so position 1 is e_0).
# ---------------------------------------------------------------------------

BuildOctMultTable := function(triples)
    local table, i, t, a, b, c;
    table := List([1..8], i -> List([1..8], j -> [0, 0]));

    # e_0 is the two-sided identity.  Position 1 = e_0.
    for i in [1..8] do
        table[1][i] := [1, i];
        table[i][1] := [1, i];
    od;

    # e_i * e_i = -e_0  for i = 1..7  (positions 2..8).
    for i in [2..8] do
        table[i][i] := [-1, 1];
    od;

    # Process each Fano triple (a, b, c) in {1..7}.
    # In our 1-indexed coords this is positions (a+1, b+1, c+1).
    for t in triples do
        a := t[1] + 1;
        b := t[2] + 1;
        c := t[3] + 1;
        # Cyclic (positive)
        table[a][b] := [ 1, c ];
        table[b][c] := [ 1, a ];
        table[c][a] := [ 1, b ];
        # Anti-cyclic (negative)
        table[b][a] := [-1, c ];
        table[c][b] := [-1, a ];
        table[a][c] := [-1, b ];
    od;
    return table;
end;

OCT_TABLE := BuildOctMultTable(STANDARD_FANO_TRIPLES);

# ---------------------------------------------------------------------------
# Octonion arithmetic on Q^8.
# ---------------------------------------------------------------------------

# OctMult(x, y) — bilinear product using OCT_TABLE.  x, y are 8-element lists.
OctMult := function(x, y)
    local r, i, j, st;
    r := [0, 0, 0, 0, 0, 0, 0, 0];
    for i in [1..8] do
        if x[i] = 0 then continue; fi;
        for j in [1..8] do
            if y[j] = 0 then continue; fi;
            st := OCT_TABLE[i][j];
            r[st[2]] := r[st[2]] + st[1] * x[i] * y[j];
        od;
    od;
    return r;
end;

# OctConjugate(x) — flip imaginary parts.
OctConjugate := function(x)
    return Concatenation([x[1]], -x{[2..8]});
end;

# OctNormSq(x) — sum of squares; equals x*conj(x) as an e_0-multiple.
OctNormSq := function(x)
    return Sum([1..8], i -> x[i]^2);
end;

# Basis vector e_k for k in 0..7 (returns a length-8 list).
OctBasis := function(k)
    local v;
    v := [0, 0, 0, 0, 0, 0, 0, 0];
    v[k + 1] := 1;
    return v;
end;

# Convenience: zero vector and identity vector.
OCT_ZERO := [0, 0, 0, 0, 0, 0, 0, 0];
OCT_E0   := [1, 0, 0, 0, 0, 0, 0, 0];

# ---------------------------------------------------------------------------
# Bridge to the LOOPS package: build the 16-element signed-units loop from
# OCT_TABLE and let LOOPS verify the Moufang property.
#
# Element labelling 1..16:
#    1..8   ->  +e_0, +e_1, ..., +e_7
#    9..16  ->  -e_0, -e_1, ..., -e_7
# ---------------------------------------------------------------------------

BuildSignedOctonionCayleyTable := function()
    local T, i, j, ai, asign, bj, bsign, st, sgn, idx;
    T := List([1..16], i -> List([1..16], j -> 0));
    for i in [1..16] do
        if i <= 8 then asign := 1; ai := i; else asign := -1; ai := i - 8; fi;
        for j in [1..16] do
            if j <= 8 then bsign := 1; bj := j; else bsign := -1; bj := j - 8; fi;
            st  := OCT_TABLE[ai][bj];
            sgn := asign * bsign * st[1];
            idx := st[2];
            if sgn = 1 then
                T[i][j] := idx;
            else
                T[i][j] := idx + 8;
            fi;
        od;
    od;
    return T;
end;
