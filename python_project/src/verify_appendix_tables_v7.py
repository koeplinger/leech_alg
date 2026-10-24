"""Run verify_appendix_tables.py on the source of the released v7 of the
paper (paper/main_2026-07-19.tex), whose Appendix A prints the three
coefficient tables certifying Lemmas 4.2 to 4.4 in v7 (Tables 1 to 3).

Created: 2026-10-20.  Updated: 2026-10-24.

All 192 printed rows (64 per table, eight coefficients each) are recomputed
exactly and compared.  Prints one line per table (the 64 products integral;
the printed table has 64 rows, identical to the recomputation) and ALL PASS.

Not cited in the paper (v8 draft, paper/main.tex); cited in ledger claim
022 and in evidence_and_reasoning/references/computer_algebra_systems.txt.

Usage: python3 verify_appendix_tables_v7.py   (under a second)
"""
from pathlib import Path

import verify_appendix_tables as vat

vat.MAIN_TEX = Path(__file__).resolve().parents[2] / "paper" / "main_2026-07-19.tex"

if __name__ == "__main__":
    vat.main()
