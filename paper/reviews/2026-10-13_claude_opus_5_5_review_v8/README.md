# Review of v8 (13 October 2026)

- `2026-10-13_claude_opus_5_5_review_v8.md`: **the review**.  It covers
  the span index of twisted triple products (index $2^{10}$ for the even
  permutation twists, and index 1 for 28 root-reflection twists), with
  two theorems (the closure criterion and the index formula), the full
  spectrum over $W(E_8)$, and the other findings on v8.
- `reproduction/final/`: the consolidated, path-relative pipeline.  Run
  `./run_all.sh`; it needs Python 3 with numpy and sympy, plus PARI/GP,
  and takes about 2 minutes.  Every script writes its result to a
  matching `*.out.txt`.
- `reproduction/independent/`: four independent implementations (Python;
  GAP 4; PARI/GP, with the 24-product certificate for index 1; Python
  with PARI), each written from the paper's definitions alone.  Their
  certificates and bases are re-checked by `final/f3` and `final/f7`.
- `original_reviewer_scripts/`: the working scripts of the review agent
  that first found the $2^{10}$ and index-1 cases.  They are kept for
  provenance and superseded by `reproduction/final/`, and they load
  `LamB.pkl` (a $\mathbb Z$-basis of $\Lambda$) from their own folder.

Not peer-reviewed; machine-obtained.  Section 5 of the review says what
is proved, what is computed exhaustively, and what is open.
