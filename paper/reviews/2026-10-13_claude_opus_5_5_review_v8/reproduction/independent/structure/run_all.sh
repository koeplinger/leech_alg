#!/bin/sh
# run_all.sh -- reproduce every output in this folder (about 2 minutes with 14 worker processes).
# Requirements: python3 with numpy and sympy; PARI/GP (gp).  All arithmetic is exact.
set -e
cd "$(dirname "$0")"
python3 s1_basics.py              # basic objects, Lambda basis, Gram matrix          -> s1_basics.out
gp -q -s 200000000 s1_leech_min.gp < /dev/null > s1_leech_min.out   # min 8, 196560 minimal vectors
python3 s2_criterion.py           # proof checks: closure criterion, 2-adic splitting -> s2_criterion.out
python3 s3_lagrangians.py         # 270 Lagrangians, all 17280 pairs, closure + index  -> s3_*.out/json/tsv
python3 s4_twists.py              # permutations (both conventions), 120 reflections,
                                  # 480 signed tables, W(E8) orbit, index realizers   -> s4_*.out/tsv/json
python3 s5_certificates.py        # HNF-free index certificates for 3 witnesses       -> s5_certificates.*
python3 s5_verify_certificates.py # standalone verifier (own octonion code)           -> s5_verify_certificates.out
gp -q -s 400000000 s6_independent.gp < /dev/null > s6_independent.out  # independent PARI/GP reimplementation
python3 s7_eps_rule.py            # the rule deciding the index, lemma for its lower bound -> s7_eps_rule.out
