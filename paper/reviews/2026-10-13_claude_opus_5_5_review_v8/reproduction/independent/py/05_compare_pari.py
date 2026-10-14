"""
05_compare_pari.py -- compare the PARI/GP results (04_pari_crosscheck.out.txt)
with the Python HNF/SNF results stored in twist_results.json, line by line.
For every exported twist: same rank (24), same index, same nontrivial
elementary divisors (PARI's matsnf on the 24x24 HNF and on the full 24x576
matrix), same F2-rank.
Run:  python3 05_compare_pari.py > 05_compare_pari.out.txt
"""
import json

res = json.load(open('twist_results.json'))
names = json.load(open('pari_crosscheck_names.json'))
byperm = {tuple(r['pi']): r for r in res['permutations']}
byroot = {tuple(r['root_doubled_e_coords']): r for r in res['reflections']}

lines = [l for l in open('04_pari_crosscheck.out.txt') if l.strip() and l[0].isdigit()]
assert len(lines) == len(names), (len(lines), len(names))
n_ok = 0
for line, (name, kind, data) in zip(lines, names):
    k, nm, rk, idx, lg, d1, d2, r2 = [t.strip() for t in line.split('|')]
    assert nm == name
    py = byperm[tuple(data)] if kind == 'perm' else byroot[tuple(data)]
    parse = lambda s: [int(t) for t in s.strip('[]').split(',') if t.strip()]
    same = (int(rk) == 24 and int(idx) == py['span_index']
            and sorted(parse(d1)) == sorted(py['elementary_divisors_nontrivial'])
            and sorted(parse(d2)) == sorted(py['elementary_divisors_nontrivial'])
            and int(r2) == py['f2_rank_S_mod_2Lambda'])
    n_ok += same
    print(f"{int(k):3d}  {name:52s} python index {py['span_index']:>6d}  PARI index {int(idx):>6d}  "
          f"F2-rank py/PARI {py['f2_rank_S_mod_2Lambda']}/{r2}  {'agree' if same else 'DISAGREE'}")
print(f"\nPARI/GP and the Python HNF/SNF agree on {n_ok} of {len(names)} twists "
      f"(index, all elementary divisors, F2-rank).")
