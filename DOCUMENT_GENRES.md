# Document Genres

Created: 2026-07-12.  Updated: 2026-10-24.

Every artifact in this repository belongs to exactly one of three genres.  The
genre determines how the artifact may be changed.  The rule behind them:
documents describe the current state, and git carries the history.

Enforced mechanically by [`tools/lint_docs.py`](tools/lint_docs.py).

---

## 1. Immutable: the prompt log

**Never edited, never rewritten, never deleted.**

| Artifact | |
|---|---|
| `prompt_logs/**` | Every prompt, in order, with what was done in response |

A prompt log entry may be *appended to* while the work it records is still in
progress (the closing `Status:` line).  Once closed, it is finished.  If a
statement in it later turns out to be wrong, the correction goes in a *later*
prompt log entry, never in the old one.

---

## 2. Public record: the released versions of the paper

**A version of the paper that has been made public is a dated record.  It is
never edited.**

| Artifact | |
|---|---|
| `paper/main_YYYY-MM-DD.tex` and its PDF | the source of a released version, as released |

A preprint (v7 is on ResearchGate) is a dated public record, but not a final
one: the latest released version is the record of the work, and when the
finalized v8 is released it becomes the single record.  A formal journal
publication, once there is one, is a public record of the same kind.  A
correction to a public record goes into the next version.

---

## 3. Current state: everything else

**Describes what IS.  Forward-evolving.  Edited in place.**

Every other artifact: the operating rules and the methodology, `CURRENT_STATE.md`,
the README files, the glossary, the research result and statement, the ledger
and the claim files, the trial results, the plans, the changelogs, the reviews
and notes in `paper/reviews/`, the companion and supplemental notes, the
reference registry, and the programs with their docstrings and comments.

Each document carries two bare dates near its top:

```
Created: YYYY-MM-DD.  Updated: YYYY-MM-DD.
```

`Created` is the day the document was first written; `Updated` is the day of
its last change, set in the same edit that makes the change.  The date is
checked (`date`), not assumed.  The paper's working source carries its date in
the subtitle, which is likewise set to the day of the last change.  A document
that does not carry the line yet gets it at its next edit.

### What must never appear in a current-state document

- A dated addendum, note, correction or status block: *"Addendum (date)"*,
  *"Note (date, prompt N)"*, *"CORRECTION (date)"*, *"Updated date (prompt N)"*,
  *"Status (date): …"* in a stack.
- A correction notice: *"(Corrected 12 July 2026; earlier versions of this file
  said …)"*, *"Correction of record"*.
- A settlement date or narrative: *"settled 12 July 2026"*, *"was open through
  v5"*, *"now proved"*, *"SETTLED (was OPEN)"*.
- A supersession narrative: *"this supersedes the earlier statement that …"*.
- History: what the document, a program or the paper used to say, which
  version did what, or at which prompt something changed.
- A change-history footer: *"Last updated: ⟨date⟩ (new entries: … ; previous
  update: …)"*.
- A status block or section restating another document's content, or a
  section that is not about what the document is about.

If a document says something that is now wrong, **make it say the right thing,
and set `Updated`.**  Do not announce the change.  The reader wants the truth,
not its provenance; `git log` and the prompt log carry the provenance.

### Kinds of current-state document

- **A program** is documented as it is: what it computes, what it prints, and
  where the paper cites it, in the paper's current numbering.
- **A plan** states its status (proposed, engaged or executed) and the
  decisions as taken, in one block below its title; its body is the plan.
- **A changelog** states the net change from one released version of the paper
  to the next, part by part, with the reasons.  It does not narrate the steps
  by which a change was reached.
- **A review or note** records an assessment.  Where it is wrong about a fact,
  it is corrected in place.  Material received from others, such as an
  external review incorporated verbatim or a proposal quoted verbatim, is
  reproduced as received; comments on it belong in the project's own text.
- **The research statement** (`evidence_and_reasoning/research_statement.md`)
  records the research question as posed.  It changes only on a significant
  deviation of the research goal, and its "Corrections applied" section records
  how the question itself had to be sharpened.

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
4. **Date every change.**  Check today's date (`date`) and set `Updated`, and
   for the paper its subtitle date, in the same edit as the change.
5. **Run the linter** (`python3 tools/lint_docs.py`) before handing work back.
