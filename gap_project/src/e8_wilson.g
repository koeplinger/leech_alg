#############################################################################
##
##  e8_wilson.g — Wilson's E_8 lattice L = D_8^+ in the octonionic basis,
##  matching python_project/src/e8_wilson.py exactly.
##
##  Wilson [Wilson2009, Section 2] uses an integer span of 240 octonions:
##    Type-1 (112): ±e_a ± e_b for distinct a, b in {0,...,7}.
##    Type-2 (128): (1/2)(±e_0 ± ... ± e_7) with an ODD number of minus signs.
##
##  Equivalently L = D_8^+ under the membership criterion
##    (all integer coords, even sum) OR (all half-integer coords, odd sum).
##
##  Wilson's distinguished element:
##    s     = (1/2)(-e_0 + e_1 + e_2 + ... + e_7)         in L  (norm 2)
##    sbar  = (1/2)(-e_0 - e_1 - e_2 - ... - e_7)         norm 2, NOT in L
##  (sbar fails the half-integer/odd-sum branch: its sum is -4.  See main
##  paper Section 2.3 and the companion note Section 7.)
##
##  Requires:  src/octonion.g
##
##  References
##    [Wilson2009] R.A. Wilson, J. Algebra 322 (2009) 2186-2190.
##
#############################################################################

# ---------------------------------------------------------------------------
# Lattice membership: L = D_8^+ in Wilson's octonionic basis.
#
#   v in L  iff  ((all v[i] in Z and sum even) or
#                 (all v[i] in Z+1/2 and sum odd integer)).
#
# Implementation uses doubled coordinates 2*v[i], which are always integers
# for elements of L; their parities determine the branch and the sum
# condition becomes (sum mod 4) = 0 (integer/even) or 2 (half-integer/odd).
# ---------------------------------------------------------------------------

IsInL := function(v)
    local doubled, parities;
    doubled := 2 * v;
    if not ForAll(doubled, x -> IsInt(x)) then
        return false;                              # not a multiple of 1/2
    fi;
    parities := List(doubled, x -> x mod 2);
    if ForAll(parities, p -> p = 0) then
        return Sum(doubled) mod 4 = 0;             # integer, even sum
    elif ForAll(parities, p -> p = 1) then
        return Sum(doubled) mod 4 = 2;             # half-integer, odd sum
    fi;
    return false;                                  # mixed parities
end;

# ---------------------------------------------------------------------------
# Root enumerators.
# ---------------------------------------------------------------------------

WilsonType1Roots := function()
    local roots, a, b, sa, sb, v;
    roots := [];
    for a in [1..8] do
        for b in [a+1..8] do
            for sa in [-1, 1] do
                for sb in [-1, 1] do
                    v := [0,0,0,0,0,0,0,0];
                    v[a] := sa;
                    v[b] := sb;
                    Add(roots, v);
                od;
            od;
        od;
    od;
    return roots;
end;

WilsonType2Roots := function()
    local roots, mask, k, num_minus, v;
    roots := [];
    for mask in [0..255] do
        num_minus := 0;
        for k in [0..7] do
            if QuoInt(mask, 2^k) mod 2 = 1 then
                num_minus := num_minus + 1;
            fi;
        od;
        if num_minus mod 2 = 1 then
            v := [];
            for k in [0..7] do
                if QuoInt(mask, 2^k) mod 2 = 1 then
                    Add(v, -1/2);
                else
                    Add(v,  1/2);
                fi;
            od;
            Add(roots, v);
        fi;
    od;
    return roots;
end;

WilsonE8Roots := function()
    return Concatenation(WilsonType1Roots(), WilsonType2Roots());
end;

# ---------------------------------------------------------------------------
# Wilson's distinguished elements.
# ---------------------------------------------------------------------------

WILSON_S     := [-1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2];
WILSON_S_BAR := OctConjugate(WILSON_S);            # = (-1/2,...,-1/2), sum -4

# ---------------------------------------------------------------------------
# Z-basis for L.
#
# Construction: take WILSON_S (a type-2 root, providing the half-integer
# coset that D_8 alone misses), then 7 type-1 roots that together with
# WILSON_S form an 8x8 matrix of full rank.  By the determinant argument in
# python_project/src/e8_wilson.py (see test_e8_wilson.py / TestRankAndSelfDuality),
# this basis has |det| = 1, so it is a unimodular Z-basis of L.
# ---------------------------------------------------------------------------

BuildLBasis := function()
    local basis, roots, r;
    basis := [WILSON_S];
    roots := WilsonE8Roots();
    for r in roots do
        if Length(basis) = 8 then break; fi;
        if RankMat(Concatenation(basis, [r])) = Length(basis) + 1 then
            Add(basis, r);
        fi;
    od;
    if Length(basis) <> 8 then
        Error("Could not build a rank-8 basis of L");
    fi;
    return basis;
end;

L_BASIS := BuildLBasis();
