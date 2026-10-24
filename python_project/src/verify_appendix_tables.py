"""Recompute the three coefficient tables of Appendix A of the released v7 of
the paper (paper/main_2026-07-19.tex) and compare them with the printed
tables.

Created: 2026-10-14.  Updated: 2026-10-24.

In v7 these tables (Tables 1 to 3) certify Lemmas 4.2 to 4.4, and v7's
Appendix C states that their basis-pair identities were recomputed in Python
and in GAP as a cross-check against transcription or implementation errors.
symbolic_proof_checks.py checks Lemmas 4.1, 4.3 and 4.4 but not Lemma 4.2;
this script covers all three tables and compares every printed entry with an
exact recomputation.

Conventions are those of v7's Appendix A: the multiplication (2.1)
(verify_ideal_decomposition.omul), the basis L_1 = s = (-e0 + e1 + ... + e7)/2,
L_2 = e0 + e1, L_3 = e0 - e1, L_{k+2} = e0 + e_k (k = 2, ..., 6), and
M_k := sigma(L_k . s), N_k := sigma(L_k . sbar) with sigma = (1 2).

  Table tbl:appendix-LL:  L_i . L_j = sum_k c_k L_k   (Lemma 4.2)
  Table tbl:appendix-LN:  L_i . N_j = sum_k c_k N_k   (Lemma 4.3)
  Table tbl:appendix-MM:  M_i . M_j = sum_k c_k M_k   (Lemma 4.4)

For each table: (1) the 64 coefficient vectors are computed exactly
(fractions.Fraction, Gaussian elimination over Q) and checked to be integral,
which proves the lemma's inclusion by Z-bilinearity; (2) the printed table is
read from the LaTeX source MAIN_TEX by its label (not by line numbers), and
its 64 rows (i, j, c_1, ..., c_8) are compared entry by entry with the
recomputation.  Prints one line per table, then ALL PASS or FAILURE.

MAIN_TEX defaults to paper/main.tex.  The paper's current source (v8 draft)
prints no coefficient tables, so on that file the script stops with a
ValueError: the table labels are not found.  verify_appendix_tables_v7.py
sets MAIN_TEX to paper/main_2026-07-19.tex and runs the check.

Usage: python3 verify_appendix_tables_v7.py   (under a second)

Not cited in the paper (v8 draft); ledger claim 022 cites
verify_appendix_tables_v7.py.
"""
import re
from fractions import Fraction as F
from pathlib import Path

from verify_ideal_decomposition import omul

MAIN_TEX = Path(__file__).resolve().parents[2] / "paper" / "main.tex"
h = F(1, 2)


def unit(k):
    return [F(1) if i == k else F(0) for i in range(8)]


def add(a, b, sign=1):
    return [x + sign * y for x, y in zip(a, b)]


def sig(v):                       # sigma = (1 2): swap coordinates 1 and 2
    return [v[0], v[2], v[1]] + list(v[3:])


s = [-h] + [h] * 7
sbar = [-h] * 8
L = [s, add(unit(0), unit(1)), add(unit(0), unit(1), -1)] + [add(unit(0), unit(k)) for k in range(2, 7)]
M = [sig(omul(b, s)) for b in L]
N = [sig(omul(b, sbar)) for b in L]


def coords(v, basis):
    """Coordinates of v in the given basis of R^8 (exact)."""
    A = [[basis[j][i] for j in range(8)] + [v[i]] for i in range(8)]
    for c in range(8):
        p = next(r for r in range(c, 8) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        A[c] = [x / A[c][c] for x in A[c]]
        for r in range(8):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [A[i][8] for i in range(8)]


def computed(left, right, target):
    return {(i + 1, j + 1): coords(omul(left[i], right[j]), target) for i in range(8) for j in range(8)}


def printed(label):
    tex = MAIN_TEX.read_text(encoding="utf-8")
    end = tex.index(r"\label{" + label + "}")
    start = tex.rindex(r"\begin{table}", 0, end)
    rows = {}
    for line in tex[start:end].splitlines():
        m = re.match(r"\s*([1-8])\s*&\s*([1-8])\s*&(.*)\\\\", line)
        if m:
            rows[(int(m.group(1)), int(m.group(2)))] = [F(int(x)) for x in m.group(3).split("&")]
    return rows


def main():
    ok = True
    for label, lemma, left, right, target in [
        ("tbl:appendix-LL", "Lemma 4.2  L . L in L", L, L, L),
        ("tbl:appendix-LN", "Lemma 4.3  L . sigma(L sbar) in sigma(L sbar)", L, N, N),
        ("tbl:appendix-MM", "Lemma 4.4  sigma(Ls) . sigma(Ls) in sigma(Ls)", M, M, M),
    ]:
        comp = computed(left, right, target)
        integral = all(c.denominator == 1 for v in comp.values() for c in v)
        rows = printed(label)
        match = len(rows) == 64 and all(rows.get(key) == comp[key] for key in comp)
        print(f"{lemma}:  64 products integral: {integral} | printed table {label}: {len(rows)} rows, "
              f"identical to the recomputation: {match}")
        ok &= integral and match
    print("ALL PASS" if ok else "FAILURE")


if __name__ == "__main__":
    main()
