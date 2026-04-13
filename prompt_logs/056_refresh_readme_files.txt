Prompt 056 -- Refresh stale README files across the project
==============================================================
Date: 2026-04-12

User prompt (verbatim):
-----------------------
I just noticed an outdated file README.md in python_project. Update
with with the missing info. Check other README files and references
in this project to make sure they're current.

Survey: README files in the repo (9 total)
-------------------------------------------
- README.md                                              up to date
- LICENSE.md                                             current (added in 052)
- prompt_logs/README.md                                  current
- source_documents/README.md                             current (refreshed in 054)
- evidence_and_reasoning/README.md                       STALE
- evidence_and_reasoning/key_claims/README.md            current (claim 008 added in 051)
- evidence_and_reasoning/references/README.md            current
- python_project/README.md                               STUB
- python_project/src/README.md                           STALE
- python_project/tests/README.md                         STUB

Actions taken:
--------------

1. python_project/README.md (was a stub with TBD setup / running
   instructions): rewritten as the real entry point.  Now describes
   the directory structure as it actually exists (no fictional
   data/ folder), names the dependencies (numpy + pytest, all else
   stdlib), gives the pytest command, gives runnable example
   commands for the symbolic-proof checks and the trial 007
   harnesses, and points readers at src/README.md and
   tests/README.md.

2. python_project/src/README.md (had stale test count "157" and was
   missing five trial-007 files plus consistency_checks.py and
   symbolic_proof_checks.py): rewritten.  Trial table extended with
   the 5 trial-007 files, each with a one-line description and the
   verdict column added (PASS for trial 007, FAIL with reason for
   001-006).  New "Verification scripts" section added for
   symbolic_proof_checks.py and consistency_checks.py.  Test count
   updated 157 -> 197.

3. python_project/tests/README.md (was a stub with no test list):
   rewritten.  Added the actual list of six test files
   (test_octonions, test_e8_wilson, test_e8_dixon, test_leech_wilson,
   test_leech_dixon, test_okubo) each with a one-line subject.
   Reproduced the pytest command at the top.  Total count 197 noted.

4. evidence_and_reasoning/README.md (key-claims table missing 008;
   trial-results table missing 007; structure block missing
   research_result.md, research_statement.md, terminology.md):
   - Key-claims table: 007 wording softened ("Untwisted triple-
     octonion algebra comprehensively ruled out"); 008 row added in
     bold ("Transposition-twisted triple octonion product is an
     order on Lambda").
   - Trial-results table: row added for trial 007 with PASS verdict
     and pointer to paper/main.tex for the symbolic proof.
   - Structure block: research_statement.md, research_result.md,
     terminology.md added to the listing so the actual contents of
     the folder are visible at a glance.

5. Spot-checked the remaining README files (root, prompt_logs/,
   source_documents/, key_claims/, references/, LICENSE.md): all
   already up to date from prompts 043, 051, 052, 054, 055.  No
   edits needed.
