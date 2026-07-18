# Document Genres

Every artifact in this repository belongs to exactly one of three genres.  The
genre determines how the artifact may be changed.  This is the rule that
prevents documentation from silently turning into a change log of itself.

Enforced mechanically by [`tools/lint_docs.py`](tools/lint_docs.py).

---

## 1. Immutable

**The prompt log.  Never edited, never rewritten, never deleted.**

| Artifact | |
|---|---|
| `prompt_logs/**` | Every prompt, in order, with what was done in response |

A prompt log entry may be *appended to* while the work it records is still in
progress (the closing `Status:` line).  Once closed, it is finished.  If a
statement in it later turns out to be wrong, the correction goes in a *later*
prompt log entry, never in the old one.

---

## 2. Frozen

**A dated record of what was true, known, or decided at a point in time.**

| Artifact | Frozen when |
|---|---|
| `paper/main_YYYY-MM-DD.tex` and the released PDF | the version is declared good to go |
| `paper/companion.tex`, supplemental and standalone notes | on release |
| `evidence_and_reasoning/key_claims/NNN_*.txt` | on writing |
| `paper/reviews/**` | on writing |
| the version changelogs in `paper/` | on release of the later version |
| `evidence_and_reasoning/trial_NNN_results.md` | when the trial is closed |
| research plans (`evidence_and_reasoning/YYYY-MM-DD_plan.md`) | once engaged |
| computer programs (`python_project/src/**`, `gap_project/**`) | once referenced by a paper, note, or claim |

A frozen artifact is **never rewritten in place.**  It changes in exactly two ways:

- **A dated addendum.**  A block headed with its date, placed at the top of the
  body (below the header), stating what is now known.  The original text stays
  as written, beneath it.  Example: the dated `--- FORWARD CORRECTION (date) ---`
  blocks in `key_claims/008_transposition_twist_order.txt`.  For a program, the
  addendum is a dated correction block in the module docstring.
- **Deprecation.**  Marked deprecated and superseded, naming what supersedes it.

The reason is not ceremony.  A frozen artifact is *evidence*: it records what
was believed on the day it was written, and a reader must be able to see both
the belief and its correction.  Rewriting it destroys the thing it exists to
preserve.

---

## 3. Current state

**Describes what IS.  Forward-evolving.  Edited freely, in place.**

| Artifact | |
|---|---|
| `MANIFESTO.md` | the operating rules |
| `DOCUMENT_GENRES.md` | this file |
| `TRIAL_METHODOLOGY.md` | how a trial is written |
| `CURRENT_STATE.md` | the state of the research |
| every indexing `README.md` | what is in this folder, and what each file does |
| `evidence_and_reasoning/terminology.md` | the glossary |
| `evidence_and_reasoning/editorial_standards.md` | the prose standards |
| `evidence_and_reasoning/research_result.md` | the condensed summary of the finding |
| `LICENSE.md`, `LICENSE-CODE` | the licenses |

**Git is the history of these files.**  They carry no history of their own.

### What must never appear in a current-state document

- A dated correction notice: *"(Corrected 12 July 2026; earlier versions of
  this file said …)"*, *"Correction of record"*.
- A settlement date: *"settled 12 July 2026"*, *"was open through v5"*,
  *"now proved"*, *"SETTLED (was OPEN)"*.
- A supersession narrative: *"this supersedes the earlier statement that …"*.
- A change-history footer: *"Last updated: ⟨date⟩ (new entries: … ; previous
  update: …)"*.  A bare date stamp with no narrative is fine.
- A status block or section restating another document's content — for example,
  the state of the manuscript inside a folder index.
- A section that is not about what the document is about.

If a current-state document said something that is now wrong, **just make it
say the right thing.**  Do not announce the change.  The reader wants the truth,
not its provenance; `git log` and the prompt log carry the provenance.

### The research statement is a special case

`evidence_and_reasoning/research_statement.md` records the research question as
posed.  It is updated only on a *significant* deviation of the research goal —
not on every finding.  Its "Corrections applied" section is part of its purpose:
it records how the question itself had to be sharpened.

---

## Denormalization is deliberate

This project keeps indexes, summaries, ledgers and cross-references that repeat
information held elsewhere.  That is not redundancy to be eliminated.  It serves
two readers at once: a human, who needs the material to be navigable, and an AI
assistant, for which a good index is the difference between reading one file and
reading the whole corpus.

So a documentation consistency sweep **may and should add** to these artifacts
when they have fallen behind.  What it must not do is add content of the wrong
*genre* to them.  The failure mode this file exists to prevent is not
"documentation grew"; it is "documentation started narrating its own history."

---

## Rules for AI assistance

1. **Log the prompt first.**  Before any other action.  (Manifesto §8.)
2. **Subagents get no write access to prose.**  A subagent may read, search,
   compute and report.  Edits to any document are made by the main thread,
   which carries the full context of the conversation and the user's standing
   instructions, and are shown to the user.  A parallel fleet of agents with
   write access turns one bad instruction into many bad edits, and the cost of
   that lands on the human.
3. **A sweep reports before it writes.**  Findings first, diffstat second, edits
   third.
4. **Run the linter** (`python3 tools/lint_docs.py`) before handing work back.
