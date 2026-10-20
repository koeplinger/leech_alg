"""Run verify_appendix_tables.py on the frozen v7 source of the paper
(paper/main_2026-07-19.tex), which carries the three Appendix A tables that
v8 retired at prompt 311 (20 October 2026).  The numeric verification is
kept: all 192 printed entries of v7 are recomputed exactly and compared.

Usage: python3 verify_appendix_tables_v7.py   (under a second)
"""
from pathlib import Path

import verify_appendix_tables as vat

vat.MAIN_TEX = Path(__file__).resolve().parents[2] / "paper" / "main_2026-07-19.tex"

if __name__ == "__main__":
    vat.main()
