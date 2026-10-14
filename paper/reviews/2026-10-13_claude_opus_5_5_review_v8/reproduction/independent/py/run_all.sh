#!/usr/bin/env bash
# run_all.sh -- reproduce every output in this folder, in dependency order.
# Independent Python re-implementation (label "py") for the review of
# "An order on the Leech lattice from a Z3-symmetric triple-octonion product" (v8).
# All conventions are stated in the header of leechlib.py.
# Run:  bash run_all.sh      (about 4 minutes on 14 cores)
#
#   00_selftest.py            primitives vs sympy (octonion table, HNF, SNF, det, F2/Fp ranks, Z/4 spans)
#   01_build_lambda.py        Z-bases of L, L.sbar, L.s, Lambda; Gram det 2^24; exports Gram for gp
#   02_qfminim.gp             PARI/GP: min 8 and 196,560 minimal vectors of Lambda
#   03_twist_census.py        claims C1-C3: sigma, all 5,040 permutations, all 120 root reflections
#   04_pari_crosscheck.gp     PARI/GP mathnf/matsnf on 108 exported product matrices
#   05_compare_pari.py        line-by-line comparison PARI vs Python
#   06_certificates.py        HNF-free certificates (Fractions from the definitions) for the witnesses
#                             and all 28 index-1 reflections
#   07_polarizations.py       claim C5 + all 17,280 complementary Lagrangian pairs (all of Aut(L))
#   08_pari_polarization_sample.py  PARI/GP check of 50 sampled polarizations (10 per index class)
#   09_index_formula.py       F2-level index formula tested on all 4,704 closing polarizations
#   10_signed_perms.py        claim C4: the 645,120 signed permutations / 480 twisted tables
#   11_span_structure.py      explicit description of S for the index-2^16 class (sigma etc.)
set -euo pipefail
cd "$(dirname "$0")"
{
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  python3 --version
  python3 -c "import numpy, sympy; print('numpy', numpy.__version__, ' sympy', sympy.__version__)"
  echo 'print("PARI/GP ", version())' | gp -q
} > versions.txt 2>&1
python3 00_selftest.py                 > 00_selftest.out.txt
python3 01_build_lambda.py             > 01_build_lambda.out.txt
gp -q 02_qfminim.gp                    > 02_qfminim.out.txt
python3 03_twist_census.py             > 03_twist_census.out.txt
gp -q 04_pari_crosscheck.gp            > 04_pari_crosscheck.out.txt
python3 05_compare_pari.py             > 05_compare_pari.out.txt
python3 06_certificates.py             > 06_certificates.out.txt
python3 07_polarizations.py            > 07_polarizations.out.txt
python3 08_pari_polarization_sample.py > 08_pari_polarization_sample.out.txt
python3 09_index_formula.py            > 09_index_formula.out.txt
python3 10_signed_perms.py             > 10_signed_perms.out.txt
python3 11_span_structure.py           > 11_span_structure.out.txt
echo "all done"
