#############################################################################
##
##  twist.g — The transposition sigma = (1 2) and the resulting twisted
##  octonion product: paper Definition 3.1 (`def:twist`) with the canonical
##  choice i = 1, j = 2, and Definition 3.2 (`def:twisted-product`).  All
##  transposition twists are equivalent up to relabeling (paper Section 4,
##  `sec:proof`, footnote to sigma = (1 2)).
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
##
##  In our 1-indexed coordinates: e_0 is at position 1, e_1 at position 2,
##  e_2 at position 3, etc.  The transposition sigma = (e_1 e_2) swaps
##  positions 2 and 3.
##
##  Twisted product (paper eq. (3.1), `eq:twist-formula`):
##    x  *_sigma  y  :=  sigma(  sigma(x)  *  sigma(y)  ).
##
##  Sigma-images of the sublattices:
##    SIGMA_LS_BASIS     = sigma( LS_BASIS )
##    SIGMA_LS_BAR_BASIS = sigma( LS_BAR_BASIS )
##
##  Requires:  src/octonion.g, src/e8_wilson.g, src/sublattices.g
##
#############################################################################

# Apply sigma = (e_1 e_2) (a coordinate transposition).
ApplySigma := function(v)
    local w;
    w := ShallowCopy(v);
    w[2] := v[3];
    w[3] := v[2];
    return w;
end;

# Twisted product:  x *_sigma y  =  sigma( sigma(x) * sigma(y) ).
OctMultSigma := function(x, y)
    return ApplySigma( OctMult( ApplySigma(x), ApplySigma(y) ) );
end;

# Sigma-images of the sublattice bases.
SIGMA_LS_BASIS     := List(LS_BASIS,     ApplySigma);
SIGMA_LS_BAR_BASIS := List(LS_BAR_BASIS, ApplySigma);

IsInSigmaLs    := function(v) return IsIntegerCombination(v, SIGMA_LS_BASIS);     end;
IsInSigmaLsBar := function(v) return IsIntegerCombination(v, SIGMA_LS_BAR_BASIS); end;
