# Tools

Created: 2026-07-12.  Updated: 2026-10-24.

Project scaffolding.  Nothing here is part of the research; these are the
checks that keep the record honest.

| File | Purpose |
|---|---|
| [`lint_docs.py`](lint_docs.py) | Enforces [`DOCUMENT_GENRES.md`](../DOCUMENT_GENRES.md) mechanically. No LLM, no network, no dependencies. |

```bash
python3 tools/lint_docs.py            # check; exit 1 on any error
python3 tools/lint_docs.py --list     # print the genre of every file
```

Run it before handing work back, and before any commit that touches
documentation.

## What it checks

1. **Genre.**  A current-state document (every document except the prompt log,
   the released versions of the paper and material kept verbatim) must not
   narrate its own history: no dated addenda, notes, correction or status
   blocks, no correction notices, no settlement dates, no supersession
   narratives, no change-history footers.  In a review, which compares versions
   of the paper by nature, only the dated-block, deprecation and footer patterns
   apply.  If such a document is wrong, it is corrected in place; git carries
   the history.
2. **Immutable.**  A prompt log that git already tracks must not be modified or
   deleted.  New entries are fine.
3. **Public record.**  A released version of the paper
   (`paper/main_YYYY-MM-DD.tex` and its PDF) must not be modified or deleted.
4. **Dangling paths.**  Every repository path named in an index document must
   exist.
5. **Index coverage.**  Every program in `python_project/src/` and
   `gap_project/` must appear in the README that indexes it.
6. **Dates.**  A changed document carries the line
   `Created: YYYY-MM-DD.  Updated: YYYY-MM-DD.`, and its `Updated` date is the
   day of the change (the file's modification date); for the paper's working
   source, the subtitle date plays this role.
