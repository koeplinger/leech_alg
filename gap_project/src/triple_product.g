#############################################################################
##
##  triple_product.g — The Z_3-symmetric triple product on R^24 of the paper
##  (Definition 3.6), in both untwisted (star_0) and sigma-twisted (star)
##  forms.
##
##  Decomposition:  R^24 = O_1 (+) O_2 (+) O_3.
##  Routing:        block alpha * block beta -> block gamma where
##                  {alpha, beta, gamma} = {1, 2, 3}, same-block products
##                  staying in the block.
##
##  Triple product u * v with u = (x, y, z), v = (x', y', z'):
##      P_1 = x * x' + z * y' + y * z',
##      P_2 = y * y' + x * z' + z * x',
##      P_3 = z * z' + y * x' + x * y'.
##
##  We provide:
##    TripleProduct(u, v)       -- using the standard octonion product.
##    TripleProductSigma(u, v)  -- using the sigma-twisted product.
##
##  Requires:  src/octonion.g, src/twist.g
##
#############################################################################

# Triple product (untwisted), using a user-supplied octonion product `mul`.
# u and v are length-3 lists of length-8 vectors.
TripleProductWith := function(u, v, mul)
    local x, y, z, xp, yp, zp, P1, P2, P3;
    x  := u[1]; y  := u[2]; z  := u[3];
    xp := v[1]; yp := v[2]; zp := v[3];
    P1 := mul(x, xp) + mul(z, yp) + mul(y, zp);
    P2 := mul(y, yp) + mul(x, zp) + mul(z, xp);
    P3 := mul(z, zp) + mul(y, xp) + mul(x, yp);
    return [P1, P2, P3];
end;

TripleProduct      := function(u, v) return TripleProductWith(u, v, OctMult);      end;
TripleProductSigma := function(u, v) return TripleProductWith(u, v, OctMultSigma); end;
