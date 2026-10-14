#!/bin/bash
# run_all.sh -- regenerate every output in this folder from scratch (about 3 minutes
# on 14 cores; Python 3 with numpy and sympy, PARI/GP).  Order matters: f6 reads
# cert_index1_final.gp (written by f3) and f2_matrices.gp (written by f2).
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
{
  echo "python: $(python3 --version 2>&1)"
  python3 -c "import numpy, sympy; print('numpy:', numpy.__version__); print('sympy:', sympy.__version__)"
  echo "PARI/GP: $(echo 'print(version())' | gp -q)"
  echo "run started: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
} > versions.txt
python3 f1_setup.py            > /dev/null   # -> f1_setup.out.txt, lambda_gram_final.gp
python3 f2_witnesses.py        > /dev/null   # -> f2_witnesses.out.txt, f2_matrices.gp
python3 f3_certificates.py     > /dev/null   # -> f3_certificates.out.txt, cert_index1_final.gp
python3 f4_census.py           > /dev/null   # -> f4_census.out.txt
python3 f5_polarizations.py    > /dev/null   # -> f5_polarizations.out.txt
gp -q f6_pari_independent.gp < /dev/null 2>&1 | sed 's/\x1b\[[0-9;]*m//g' \
   | command grep -v "Warning: new maximum stack" > f6_pari_independent.out.txt
python3 f7_structure_checks.py > /dev/null   # -> f7_structure_checks.out.txt
python3 f8_ideal_classes.py    > /dev/null   # -> f8_ideal_classes.out.txt (ideal type of V per closing pair)
echo "run finished: $(date -u '+%Y-%m-%d %H:%M:%S UTC')" >> versions.txt
