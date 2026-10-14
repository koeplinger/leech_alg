#!/bin/sh
# run_all.sh -- regenerate every output of this GAP / PARI-GP review check, in order.
# Requires GAP 4 (tested 4.12.1) and PARI/GP (tested 2.15.4).  Runtime about 20 minutes.
# Reads the paper source only once, read-only, to extract its Appendix-A tables.
set -e
cd "$(dirname "$0")"
./extract_paper_tables.sh ../../../../main.tex paper_tables.g
gap -q -o 4g run_01_setup_checks.g        < /dev/null > run_01_setup_checks.out
gap -q -o 4g run_02_witnesses.g           < /dev/null > run_02_witnesses.out        # also writes cert_index1_reflection.{txt,gp}
gap -q -o 4g run_03_permutation_census.g  < /dev/null > run_03_permutation_census.out  # also writes perm_census_table.txt
gap -q -o 4g run_04_reflections.g         < /dev/null > run_04_reflections.out       # also writes refl_census_table.txt
gap -q -o 4g run_05_signed_permutations.g < /dev/null > run_05_signed_permutations.out
gap -q -o 6g run_06_polarization_census.g < /dev/null > run_06_polarization_census.out
gap -q -o 4g run_07_two_reflection_twists.g < /dev/null > run_07_two_reflection_twists.out
gap -q -o 4g run_08_explicit_isometries.g < /dev/null > run_08_explicit_isometries.out
gp -q verify_cert_index1.gp               < /dev/null > verify_cert_index1.out
# failure markers: any "[FAIL]" or GAP "Error", a missing end marker, or a failed certificate check
grep -l -E '\[FAIL\]|Error' ./*.out && exit 1 || true
for f in run_0*.out; do
  grep -q '=== end run_0' "$f" || { echo "missing end marker in $f"; exit 1; }
done
grep -q 'CERTIFICATE VALID (S = Lambda for the reflection twist): 1' verify_cert_index1.out \
  || { echo "PARI/GP certificate check failed"; exit 1; }
echo "all outputs regenerated"
