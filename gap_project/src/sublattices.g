#############################################################################
##
##  sublattices.g — The sublattices Ls, Lsbar of L, and integer-combination
##  membership tests.  Mirrors python_project/src/leech_wilson.py (the
##  Z-basis construction and `is_in_Ls`, `is_in_Ls_bar` predicates).
##
##  Definitions (paper Section 2.3):
##    Ls    = { ell * s     : ell in L },     index 16 in L.
##    Lsbar = { ell * sbar  : ell in L },     index 16 in L.
##  Both are sublattices of L (despite sbar not lying in L).
##
##  Z-basis construction.  If {b_1,...,b_8} is a Z-basis of L and c is any
##  octonion, then {b_1*c, ..., b_8*c} is a Z-basis of L*c, because
##  multiplication on the left by a fixed c is Z-linear:  for n_i in Z,
##    (sum n_i b_i) * c = sum n_i (b_i * c).
##
##  Requires:  src/octonion.g, src/e8_wilson.g
##
##  References
##    [Wilson2009]  R.A. Wilson, J. Algebra 322 (2009), Section 2.
##
#############################################################################

# ---------------------------------------------------------------------------
# Right-multiply each row of a basis by a fixed octonion c.
# ---------------------------------------------------------------------------

RightMultiplyBasis := function(basis, c)
    return List(basis, b -> OctMult(b, c));
end;

# ---------------------------------------------------------------------------
# Bases of Ls and Lsbar (and their sigma-images, defined in twist.g).
# ---------------------------------------------------------------------------

LS_BASIS     := RightMultiplyBasis(L_BASIS, WILSON_S);
LS_BAR_BASIS := RightMultiplyBasis(L_BASIS, WILSON_S_BAR);

# ---------------------------------------------------------------------------
# IsIntegerCombination(v, basis) — true iff v lies in the Z-span of basis.
#
# Uses GAP's exact-rational SolutionMat to find c with c * basis = v, then
# checks that every component of c is an integer.
# ---------------------------------------------------------------------------

IsIntegerCombination := function(v, basis)
    local c;
    c := SolutionMat(basis, v);
    if c = fail then return false; fi;
    return ForAll(c, x -> IsInt(x));
end;

# Membership predicates.
IsInLs    := function(v) return IsIntegerCombination(v, LS_BASIS);     end;
IsInLsBar := function(v) return IsIntegerCombination(v, LS_BAR_BASIS); end;
