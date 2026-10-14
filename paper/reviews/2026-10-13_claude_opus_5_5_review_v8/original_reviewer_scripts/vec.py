"""Vectorized exact integer helpers (numpy int64 on scaled coordinates)."""
import sys, pickle
sys.path.insert(0, '.')
import numpy as np
from core import *

# LB-coordinate map for scaled inputs: a = x * LBinv; LBinv has denominators dividing 2
den = 1
for row in LBinv:
    for v in row:
        den = den * v.denominator // np.gcd(den, v.denominator)
LBinv_int = np.array([[int(v * den) for v in row] for row in LBinv], dtype=np.int64)
HV = np.array(H_V, dtype=np.int64)
HW = np.array(H_W, dtype=np.int64)

def omul_arr(X, Y, TT=T):
    """X, Y: (...,8) integer arrays; returns (...,8) product (no rescaling)."""
    return np.einsum('...a,...b,abc->...c', X, Y, TT)

def Lcoords_scaled(X, scale):
    """X: (N,8) int array representing X/scale. Returns (ok, A) with A integer LB-coords if X/scale in L."""
    num = X @ LBinv_int            # = scale*den * coords
    q = scale * den
    ok = np.all(num % q == 0, axis=1)
    A = num // q
    return ok, A

def in_L_scaled(X, scale):
    """exact D8+ (odd half-integer coset) membership of X/scale."""
    # all coords integers with even sum, or all half-integers with odd sum
    Y2 = 2 * X
    ok_div = np.all(Y2 % scale == 0, axis=1)          # 2x integer
    Z = Y2 // scale                                   # = 2x (integer where ok_div)
    allint = np.all(Z % 2 == 0, axis=1)
    allhalf = np.all(Z % 2 != 0, axis=1)
    ssum = Z.sum(axis=1)  # = 2*sum
    ev = (ssum // 2) % 2 == 0
    return ok_div & ((allint & (ssum % 2 == 0) & ev) | (allhalf & (ssum % 2 == 0) & ~ev))

def in_Lambda_scaled(Vv, scale):
    X, Y, Z = Vv[:, 0:8], Vv[:, 8:16], Vv[:, 16:24]
    ok = in_L_scaled(X, scale) & in_L_scaled(Y, scale) & in_L_scaled(Z, scale)
    _, Ax = Lcoords_scaled(X, scale)
    _, Ay = Lcoords_scaled(Y, scale)
    _, Az = Lcoords_scaled(Z, scale)
    Ax %= 2; Ay %= 2; Az %= 2
    c1 = np.all(((Ax + Ay) @ HV.T) % 2 == 0, axis=1)
    c2 = np.all(((Ax + Az) @ HV.T) % 2 == 0, axis=1)
    c3 = np.all(((Ay + Az) @ HV.T) % 2 == 0, axis=1)
    c4 = np.all(((Ax + Ay + Az) @ HW.T) % 2 == 0, axis=1)
    return ok & c1 & c2 & c3 & c4

LamB = pickle.load(open('LamB.pkl', 'rb'))
B2 = np.array([[int(2 * v) for v in row] for row in LamB], dtype=np.int64)   # 2*basis, integer
Binv = mat_inv_fr(LamB)
dB = 1
for row in Binv:
    for v in row:
        dB = dB * v.denominator // np.gcd(dB, v.denominator)
Binv_int = np.array([[int(v * dB) for v in row] for row in Binv], dtype=np.int64)

def lam_coords_scaled(Vv, scale):
    """Vv: (N,24) int array representing Vv/scale. Returns (ok, C) integer Lambda-coords when integral."""
    num = Vv @ Binv_int
    q = scale * dB
    ok = np.all(num % q == 0, axis=1)
    return ok, num // q

def star_arr(U, V, TT=TS):
    """Z3-routed triple product, blockwise with structure tensor TT. U,V: (N,24)."""
    x, y, z = U[:, 0:8], U[:, 8:16], U[:, 16:24]
    xp, yp, zp = V[:, 0:8], V[:, 8:16], V[:, 16:24]
    m = lambda a, b: np.einsum('na,nb,abc->nc', a, b, TT)
    P1 = m(x, xp) + m(z, yp) + m(y, zp)
    P2 = m(y, yp) + m(x, zp) + m(z, xp)
    P3 = m(z, zp) + m(y, xp) + m(x, yp)
    return np.concatenate([P1, P2, P3], 1)

def basis_products(BB, TT=TS):
    """all 576 products of rows of BB (24 x 24 int), row-major (i,j) -> i*24+j."""
    n = BB.shape[0]
    I, Jx = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
    return star_arr(BB[I.ravel()], BB[Jx.ravel()], TT)
