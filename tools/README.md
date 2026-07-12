# Tools

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

1. **Genre.**  A current-state document must not narrate its own history: no
   dated correction notices, no settlement dates, no supersession narratives,
   no change-history footers.  If such a document is wrong, it is corrected in
   place; git carries the history.
2. **Immutable.**  A prompt log that git already tracks must not be modified or
   deleted.  New entries are fine.
3. **Frozen.**  A frozen artifact (a released paper version, a claim file, a
   review, a referenced program) must not be rewritten in place: a change to
   one must add a dated addendum, or mark it deprecated and superseded.
4. **Dangling paths.**  Every repository path named in a current-state document
   must exist.
5. **Index coverage.**  Every program in `python_project/src/` and
   `gap_project/` must appear in the README that indexes it.
