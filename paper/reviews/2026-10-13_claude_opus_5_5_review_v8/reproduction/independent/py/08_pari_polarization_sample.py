r"""
08_pari_polarization_sample.py -- PARI/GP cross-check of the polarization
census (07): for 10 random closing ordered pairs (V, W) from each of the five
index classes found in 07 (seed 2026), rebuild Lambda(V,W), compute the 576
products of the UNTWISTED triple product in Lambda(V,W)-coordinates, and let
PARI/GP compute rank, |det| of mathnf, matsnf and the F2-rank; compare with
polarization_results.json.  The GP code is the same as in 04_pari_crosscheck.gp
and is passed to `gp -q` on stdin; its raw output is saved to
08_pari_polarization_sample.gp.out.txt.
Run:  python3 08_pari_polarization_sample.py > 08_pari_polarization_sample.out.txt
"""
import json
import random
import subprocess
from fractions import Fraction as F

import numpy as np
import leechlib as LL

D = LL.build_all()
C = D['C']
pol = json.load(open('polarization_results.json'))
lag = pol['lagrangians']


def vec(m):
    return [(m >> i) & 1 for i in range(8)]


def products_for_pair(Vs, Ws):
    Vb, _ = LL.f2_rowspace_basis([vec(m) for m in Vs], 8)
    Wb, _ = LL.f2_rowspace_basis([vec(m) for m in Ws], 8)
    Z8 = [0] * 8
    gens = [w + w + w for w in Wb] + [v + Z8 + v for v in Vb] + [v + v + Z8 for v in Vb]
    gens += [[2 * int(i == j) for j in range(24)] for i in range(24)]
    BLam = LL.hnf(gens, 24)
    inv, den = LL.int_inverse_scaled([[F(x) for x in r] for r in BLam])
    N2 = np.array([[2 * x // den for x in row] for row in inv], dtype=np.int64)
    P = LL.all_basis_products(np.array(BLam, dtype=np.int64), C).reshape(576, 24)
    ok, c2 = LL.lambda_coordinates(P, N2)
    assert ok.all()
    return c2 // 2


rng = random.Random(2026)
closing = [p for p in pol['pairs'] if p['closes']]
classes = sorted({p['index'] for p in closing})
sample = []
for idx in classes:
    members = [p for p in closing if p['index'] == idx]
    sample += rng.sample(members, 10)

gp_in = ['NAMES = [' + ','.join(f'"pair {p["a"]},{p["b"]}"' for p in sample) + '];',
         f'MATS = vector({len(sample)});']
for k, p in enumerate(sample, 1):
    M = products_for_pair(lag[p['a']], lag[p['b']])
    gp_in.append(f'MATS[{k}] = [' + ';'.join(','.join(str(int(a)) for a in row) for row in M) + '];')
gp_in.append(r'''nontriv(d) = select(x -> x != 1, d);
for (k = 1, #MATS, M = MATS[k]; H = mathnf(M~); rk = matsize(H)[2]; idx = if (rk == 24, abs(matdet(H)), 0); d1 = nontriv(matsnf(H)); r2 = matrank(M * Mod(1, 2)); print(k, " | ", NAMES[k], " | ", rk, " | ", idx, " | ", d1, " | ", r2));
quit;''')
out = subprocess.run(['gp', '-q'], input='\n'.join(gp_in), capture_output=True, text=True).stdout
open('08_pari_polarization_sample.gp.out.txt', 'w').write(out)
lines = [l for l in out.splitlines() if l.strip() and l[0].isdigit()]
assert len(lines) == len(sample)
agree = 0
for line, p in zip(lines, sample):
    k, nm, rk, idx, d1, r2 = [t.strip() for t in line.split('|')]
    divs = [int(t) for t in d1.strip('[]').split(',') if t.strip()]
    same = int(rk) == 24 and int(idx) == p['index'] and int(r2) == p['f2rank'] and \
        all(d == 2 for d in divs) and len(divs) == 24 - p['f2rank']
    agree += same
    print(f"{nm:14s} python index {p['index']:>6d}  PARI index {int(idx):>6d}  SNF: {len(divs)} x 2  "
          f"F2-rank py/PARI {p['f2rank']}/{r2}  {'agree' if same else 'DISAGREE'}")
print(f"\nPARI/GP agrees with the Python HNF on {agree} of {len(sample)} sampled closing polarizations "
      f"(10 from each index class {[f'2^{c.bit_length() - 1}' for c in classes]}).")
