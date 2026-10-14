import sys; sys.path.insert(0,'.')
from vec import *
rng = np.random.default_rng(11)
N=200
U = rng.integers(-2,3,size=(N,24)) @ B2; V = rng.integers(-2,3,size=(N,24)) @ B2; W = rng.integers(-2,3,size=(N,24)) @ B2  # scale 2
S = lambda a,b: star_arr(a,b)
def frac_fail(lhs, rhs): return float(np.mean(np.any(lhs!=rhs,axis=1)))
print("commutativity fails on fraction:", frac_fail(S(U,V), S(V,U)))
print("associativity fails on fraction:", frac_fail(S(S(U,V),W), S(U,S(V,W))))
print("left alternativity fails:", frac_fail(S(S(U,U),V), S(U,S(U,V))))
print("right alternativity fails:", frac_fail(S(S(U,V),V), S(U,S(V,V))))
print("flexibility fails:", frac_fail(S(S(U,V),U), S(U,S(V,U))))
print("power-assoc (x^2 x = x x^2) fails:", frac_fail(S(S(U,U),U), S(U,S(U,U))))
