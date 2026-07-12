#############################################################################
##
##  leech_wilson.g — Wilson's Leech lattice Lambda inside L^3, and the
##  three-condition membership test (paper Definition 2.4).
##  CORRECTION (2026-07-12): "Definition 2.4" is the v3 numbering. In v4 and
##  later, Wilson's three conditions are Definition 2.1 (label def:wilson-leech).
##
##  Lambda = { (x, y, z) in L^3  :
##              (1) x, y, z in L;
##              (2) x+y, x+z, y+z in Lsbar;
##              (3) x + y + z in Ls. }
##
##  Generators of the 196,560 minimal vectors (paper Table 1, after Wilson's
##  Section 3).  J = { +-e_t : t = 0..7 } is the set of 16 unit octonions.
##
##  Requires:  src/octonion.g, src/e8_wilson.g, src/sublattices.g
##
#############################################################################

# ---------------------------------------------------------------------------
# Wilson's three conditions.
# ---------------------------------------------------------------------------

IsInLeech := function(x, y, z)
    if not (IsInL(x) and IsInL(y) and IsInL(z)) then return false; fi;
    if not (IsInLsBar(x + y) and IsInLsBar(x + z) and IsInLsBar(y + z)) then
        return false;
    fi;
    if not IsInLs(x + y + z) then return false; fi;
    return true;
end;

# ---------------------------------------------------------------------------
# J = { +- e_t : t = 0..7 } — 16 unit octonions.
# ---------------------------------------------------------------------------

J_UNITS := function()
    local result, t, e;
    result := [];
    for t in [0..7] do
        e := OctBasis(t);
        Add(result,  e);
        Add(result, -e);
    od;
    return result;
end;

# ---------------------------------------------------------------------------
# Type-1 minimal vectors of Lambda:
#   (2*lambda, 0, 0) and the two cyclic permutations,  lambda in 240 roots.
#   Count = 3 * 240 = 720,  each of squared norm 8.
# ---------------------------------------------------------------------------

LeechType1Vectors := function()
    local result, roots, r, two_r, zero;
    result := [];
    zero  := [0,0,0,0,0,0,0,0];
    roots := WilsonE8Roots();
    for r in roots do
        two_r := 2 * r;
        Add(result, [ two_r, zero,  zero ]);
        Add(result, [ zero,  two_r, zero ]);
        Add(result, [ zero,  zero,  two_r]);
    od;
    return result;
end;
