"""Dixon's 196,560 vectors are the minimal vectors of a Leech lattice: exact checks.

Created: 2026-10-31.  Updated: 2026-10-31.

leech_dixon.py builds the three families of [Dixon2010] Section 6, eqs.
(20)-(22), at the scale of Wilson's Λ (norm 8).  This script checks, exactly:

  1. the 196,560 vectors are distinct and all have norm 8 (integer
     arithmetic on the coordinates times 2);
  2. their componentwise conjugates (x̄, ȳ, z̄) are exactly Wilson's 196,560
     minimal vectors (leech_wilson.py), family by family;
  3. independently of Wilson's construction, with PARI/GP: the Z-span M of
     the vectors has rank 24; rescaled to norm 4 (M/√2), its Gram matrix on
     a Hermite-normal-form basis is integral with an even diagonal and
     determinant 1; it has no vector of norm 2 and exactly 196,560 vectors of
     norm 4.  So M/√2 is even unimodular of rank 24 without roots, the Leech
     lattice (Conway's uniqueness theorem), and the 196,560 vectors are all
     of its minimal vectors.

Check 3 needs PARI/GP (`gp`, version 2.11 or later) on the PATH; without it
the script reports check 3 as skipped, prints PASS for checks 1 and 2 only,
and exits with code 2.  Coordinates are taken as exact half-integers, and the
script stops if one is not.  About 10 seconds; about 1.3 GB of memory at
peak (Python and gp together).  Prints one line per check and ALL PASS or
FAIL; exit code 0 on ALL PASS.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from leech_dixon import dixon_leech_type1, dixon_leech_type2, dixon_leech_type3
from leech_wilson import leech_type1_vectors, leech_type2_vectors, leech_type3_vectors

CONJ = np.array([1] + [-1] * 7)

GP_SCRIPT = r"""
default(nbthreads, 1);
default(parisize, 600000000);
V = readvec("%s");
n = #V;
M = matrix(24, n, r, c, V[c][r]);
H = mathnf(matrix(24, 3000, r, c, M[r, 1 + ((c * 7919) %% n)]));
{
while (1,
  X = H^(-1) * M;
  out = select(c -> denominator(X[, c]) != 1, vector(n, c, c));
  if (#out == 0, break);
  H = mathnf(concat(H, matrix(24, #out, r, c, M[r, out[c]]))));
}
G4 = (H~ * H) / 8;
U = qflllgram(G4); G4r = U~ * G4 * U;
print("rank ", matsize(H)[2]);
print("integral ", denominator(G4) == 1);
print("even ", vector(24, i, G4[i, i] %% 2) == vector(24));
print("det ", matdet(G4));
print("norm<=2 ", qfminim(G4r, 2, 0)[1]);
print("norm<=4 ", qfminim(G4r, 4, 0)[1]);
quit;
"""


def keys(vecs, conjugate=False):
    """Integer keys 2v (or 2v̄), after checking that 2v is integral exactly."""
    sign = np.concatenate([CONJ] * 3) if conjugate else np.ones(24)
    w = 2 * sign * np.array([np.concatenate(t) for t in vecs])
    k = np.rint(w)
    if not np.array_equal(w, k):
        raise ValueError("a coordinate is not an exact half-integer")
    return [tuple(r) for r in k.astype(int).tolist()]


def main() -> int:
    ok = True
    d = [dixon_leech_type1(), dixon_leech_type2(), dixon_leech_type3()]
    w = [leech_type1_vectors(), leech_type2_vectors(), leech_type3_vectors()]

    allkeys = [k for fam in d for k in keys(fam)]
    distinct = len(set(allkeys)) == len(allkeys) == 196_560
    norms = {sum(x * x for x in k) for k in allkeys}          # |2v|^2 = 4 |v|^2
    print(f"1. {len(allkeys)} vectors, distinct: {distinct}; norms {sorted(n // 4 for n in norms)}")
    ok &= distinct and norms == {32}

    same = [set(keys(df, conjugate=True)) == set(keys(wf)) for df, wf in zip(d, w)]
    print(f"2. conjugates equal Wilson's families 1, 2, 3: {same}")
    ok &= all(same)

    gp = shutil.which("gp")
    if gp is None:
        print("3. skipped: PARI/GP (gp) not found on the PATH")
    else:
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "vectors.gp"
            data.write_text("\n".join("[" + ",".join(map(str, k)) + "]" for k in allkeys) + "\n")
            script = Path(tmp) / "test.gp"
            script.write_text(GP_SCRIPT % data.as_posix())
            out = subprocess.run([gp, "-q", script.as_posix()], capture_output=True, text=True).stdout
        res = dict(line.split(" ", 1) for line in out.strip().splitlines() if " " in line)
        expect = {"rank": "24", "integral": "1", "even": "1", "det": "1",
                  "norm<=2": "0", "norm<=4": "196560"}
        good = all(res.get(k) == v for k, v in expect.items())
        print(f"3. PARI/GP: rank {res.get('rank')}, Gram at norm 4 integral {res.get('integral')}, "
              f"even {res.get('even')}, det {res.get('det')}; vectors of norm <= 2: {res.get('norm<=2')}, "
              f"of norm <= 4: {res.get('norm<=4')}")
        ok &= good

    if ok and gp is None:
        print("PASS (checks 1 and 2); check 3 SKIPPED")
        return 2
    print("ALL PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
