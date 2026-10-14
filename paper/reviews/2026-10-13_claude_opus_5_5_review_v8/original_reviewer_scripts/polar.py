import sys, json, itertools
sys.path.insert(0,'.')
from vec import *
from collections import Counter
res = json.load(open('census_results.json'))
LB2 = np.array([[int(2*v) for v in b] for b in LB],dtype=np.int64)
Lsb2 = np.array([[int(4*v) for v in b] for b in LsbarB],dtype=np.int64)   # scale 4
Ls2 = np.array([[int(4*v) for v in b] for b in LsB],dtype=np.int64)
def f2span_key(rows):
    # rows: integer LB-coords; return frozenset of the F2 span (as tuples)
    R,_ = f2_rref([[int(x)%2 for x in r] for r in rows], 8)
    return tuple(tuple(r) for r in R)
def permute(X, p):   # apply pi: e_a -> e_{p(a)} to rows of X
    Y = np.zeros_like(X)
    for a in range(8): Y[:,p[a]] = X[:,a]
    return Y
pols = Counter(); pol_closing = set(); pol_nonclosing=set()
equiv_ok = 0; tot=0
for p,closed,idx in res:
    p = tuple(p)
    Vp = permute(Lsb2,p); Wp = permute(Ls2,p)
    okv, Av = Lcoords_scaled(Vp,4); okw, Aw = Lcoords_scaled(Wp,4)
    assert okv.all() and okw.all()
    key = (f2span_key(Av.tolist()), f2span_key(Aw.tolist()))
    (pol_closing if closed else pol_nonclosing).add(key)
    # criterion: L.pi(Lsbar) in pi(Lsbar)  and  pi(Ls).pi(Ls) in pi(Ls)   (mod-2 test suffices since both contain 2L)
    Vset = set(itertools.product((0,1),repeat=8))
    def span_set(R):
        S={(0,)*8}
        for r in R: S |= {tuple((a+b)%2 for a,b in zip(s_,r)) for s_ in S}
        return S
    Vs = span_set(key[0]); Ws = span_set(key[1])
    LV = omul_arr(np.repeat(LB2,8,0), np.tile(Vp,(8,1)))        # scale 8
    WW = omul_arr(np.repeat(Wp,8,0), np.tile(Wp,(8,1)))         # scale 16
    ok1, A1 = Lcoords_scaled(LV,8); ok2, A2 = Lcoords_scaled(WW,16)
    crit = ok1.all() and ok2.all() and all(tuple(int(x)%2 for x in r) in Vs for r in A1.tolist()) and all(tuple(int(x)%2 for x in r) in Ws for r in A2.tolist())
    tot += 1; equiv_ok += (crit == closed)
print("permutations where [closure of star_pi] == [L.pi(Lsbar) in pi(Lsbar) and pi(Ls).pi(Ls) in pi(Ls)]:", equiv_ok, "/", tot)
print("distinct polarizations (pi(Lsbar)/2L, pi(Ls)/2L) from closing pi:", len(pol_closing), "; from non-closing pi:", len(pol_nonclosing),
      "; overlap:", len(pol_closing & pol_nonclosing))

print("\n=== Appendix B: 21 transpositions applied to Kirmse's J_1 ===")
# J_1 code: supports {0} u Fano triple, and complements
fano = [(1,2,4),(2,3,5),(3,4,6),(4,5,7),(5,6,1),(6,7,2),(7,1,3)]
w = [tuple(1 if i in ({0}|set(t)) else 0 for i in range(8)) for t in fano]
def span(gs):
    S={(0,)*8}
    for g in gs: S |= {tuple((a+b)%2 for a,b in zip(x,g)) for x in S}
    return frozenset(S)
J1 = span(w + [(1,)*8])
print("J_1 code size:", len(J1), " doubly even:", all(sum(x)%4==0 for x in J1))
from item5_helpers import lattice_closed
imgs = Counter()
for i,j in itertools.combinations(range(1,8),2):
    C = frozenset(tuple(x[{i:j,j:i}.get(k,k)] for k in range(8)) for x in J1)
    imgs[C]+=1
print("distinct images of J_1 under the 21 transpositions:", len(imgs), "; multiplicities:", sorted(imgs.values()),
      "; closed (paper's multiplication):", sum(lattice_closed(C,T) for C in imgs))
