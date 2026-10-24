#!/usr/bin/env python3
"""Documentation linter.  Enforces DOCUMENT_GENRES.md mechanically.

Created: 2026-07-12.  Updated: 2026-10-24.

No LLM, no network, no dependencies.  Runs in a few seconds.

Genres
------
immutable      prompt_logs/ (except its README)
public record  the released versions of the paper, paper/main_YYYY-MM-DD.tex
               and their PDFs
current state  every other document.  Material received from others and
               kept as received (the external review of 13 October 2026,
               its folder and its write-up; the PDFs in source_documents/)
               is never edited; checks 1 and 6 skip it, and --list labels
               it "verbatim".

Checks
------
1. GENRE       A current-state document describes what IS: no dated addenda,
               notes, correction or status blocks, no correction notices, no
               settlement or supersession narratives, no change-history
               footers.
2. IMMUTABLE   Files under prompt_logs/ that are already tracked by git must
               not be modified or deleted.  New entries are fine, and so is
               an append to an entry whose committed version is still open
               (it has no closing "Status:" line).
3. PUBLIC      A released version of the paper must not be modified or
               deleted.
4. DANGLING    Every repository path named in an index document must exist.
5. INDEX       Every program in python_project/src/ and gap_project/ must
               appear in the README that indexes it.
6. DATES       A current-state document that is changed (modified or new)
               carries the line "Created: YYYY-MM-DD.  Updated: YYYY-MM-DD.",
               and its Updated date is the date of the change (the file's
               modification date).  For the paper's working source the
               subtitle date plays this role.

Usage
-----
    python3 tools/lint_docs.py            # check
    python3 tools/lint_docs.py --list     # print the genre of every file

Exit code 0 if clean, 1 if any error was found.
"""
from __future__ import annotations

import datetime as dt
import fnmatch
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# ----------------------------------------------------------------------------
# Genre membership (see DOCUMENT_GENRES.md)
# ----------------------------------------------------------------------------

IMMUTABLE_DIRS = ["prompt_logs"]

PUBLIC_RECORD_GLOBS = [
    "paper/main_[0-9]*.tex",
    "paper/main_[0-9]*.pdf",
]

VERBATIM_GLOBS = [
    "paper/reviews/2026-10-13_claude_opus_5_5_review_v8/*",
    "paper/reviews/2026-10-13_claude_opus_5_5_review_v8.md",
    "source_documents/*.pdf",
]

# Documents in the sense of checks 1 and 6 (data files and outputs are not).
DOC_SUFFIXES = {".md", ".txt", ".py", ".g", ".gp", ".tex", ".sh"}

# The paper's working source carries its date in the subtitle (check 6).
PAPER_SOURCE = "paper/main.tex"

# Index documents: check 4 resolves every repository path they name.
INDEX_DOCS = [
    "MANIFESTO.md",
    "DOCUMENT_GENRES.md",
    "TRIAL_METHODOLOGY.md",
    "CURRENT_STATE.md",
    "README.md",
    "evidence_and_reasoning/terminology.md",
    "evidence_and_reasoning/editorial_standards.md",
    "evidence_and_reasoning/research_result.md",
    "evidence_and_reasoning/README.md",
    "evidence_and_reasoning/key_claims/README.md",
    "evidence_and_reasoning/references/README.md",
    "paper/README.md",
    "python_project/README.md",
    "python_project/src/README.md",
    "python_project/tests/README.md",
    "gap_project/README.md",
    "prompt_logs/README.md",
    "source_documents/README.md",
    "tools/README.md",
]

# research_statement.md is a special case: its "Corrections applied" section
# records how the research question was sharpened (DOCUMENT_GENRES.md).
GENRE_EXEMPT = {"evidence_and_reasoning/research_statement.md"}

# A rule book has to be able to name the thing it forbids.  These documents
# DEFINE the genre vocabulary, so they are exempt from check 1 (only).
SPEC_FILES = {"DOCUMENT_GENRES.md", "MANIFESTO.md", "tools/README.md", "tools/lint_docs.py"}

# ----------------------------------------------------------------------------
# 1. GENRE: patterns that must not appear in a current-state document
# ----------------------------------------------------------------------------

DATE = r"(?:\d{1,2}\s+\w+\s+20\d\d|20\d\d-\d\d-\d\d|\w+\s+\d{1,2},?\s+20\d\d)"

BANNED = [
    (rf"\(\s*[Cc]orrected\s+{DATE}", "dated correction notice"),
    (rf"\*\*\s*Correction\s*\(\s*{DATE}", "dated correction notice"),
    (r"[Cc]orrection of record", "dated correction notice"),
    (rf"\bsettled\s+(?:on\s+)?{DATE}", "settlement date"),
    (rf"\bestablished\s+{DATE}", "settlement date"),
    (rf"\bexcluded\s+(?:on\s+)?{DATE}", "settlement date"),
    (r"\bwas\s+OPEN\b", "settlement narrative"),
    (r"was open through", "settlement narrative"),
    (r"—\s*SETTLED\b|-- SETTLED\b", "settlement narrative"),
    (r"\bis now proved\b|\bnow settled\b", "settlement narrative"),
    (r"\bsupersede[sd]?\b", "supersession narrative"),
    (r"earlier version[s]? of this file", "supersession narrative"),
    (r"an earlier claim (?:was )?(?:corrected|refuted)", "supersession narrative"),
    (r"^Last updated:.*\(", "change-history footer"),
    (r"^Previous update", "change-history footer"),
    # A dated block opens its line (after list, heading, bold, quote or comment
    # decoration) and carries a date within a few words.
    (rf"^[ \t\-#*>(\"'%]*(?:(?:FORWARD\s+)?CORRECTIONS?|ADDENDUM|Note(?:\s+added)?|Updated|Status|Revision\s*\d*)\b.{{0,40}}?{DATE}",
     "dated addendum, note, correction or status block"),
    (r"^[ \t\-#*>(\"'%]*(?:DEPRECATED|SUPERSEDED)\b", "deprecation marker"),
]
BANNED = [(re.compile(p, re.M | re.I if "ADDENDUM" in p else re.M), why) for p, why in BANNED]

# A review compares versions of the paper by nature, so words like
# "superseded" there are about the paper, not about the review itself; in
# reviews only the dated-block, deprecation and footer patterns apply.
REVIEW_DIR = "paper/reviews/"
BLOCK_KINDS = {"dated addendum, note, correction or status block",
               "deprecation marker", "change-history footer"}

# The date stamp every current-state document carries (check 6); it is the
# one dated line a current-state document may have.
STAMP = re.compile(r"Created:\s*(20\d\d-\d\d-\d\d)\.?\s+Updated:\s*(20\d\d-\d\d-\d\d)\.?")

# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------


def git(*args: str) -> str:
    try:
        return subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def matches(rel: str, globs: list[str]) -> bool:
    return any(fnmatch.fnmatch(rel, g) or rel.startswith(g.rstrip("*")) and g.endswith("/*")
               for g in globs)


def genre_of(rel: str) -> str:
    if any(rel.startswith(d + "/") for d in IMMUTABLE_DIRS) and not rel.endswith("README.md"):
        return "immutable"
    if matches(rel, PUBLIC_RECORD_GLOBS):
        return "public-record"
    if matches(rel, VERBATIM_GLOBS):
        return "verbatim"
    return "current-state"


def all_files() -> list[str]:
    tracked = git("ls-files").splitlines()
    new = git("ls-files", "--others", "--exclude-standard").splitlines()
    return sorted({f for f in tracked + new if (ROOT / f).is_file()})


def changed_files() -> dict[str, str]:
    """Path -> porcelain status code, for modified, added, renamed or new files."""
    out: dict[str, str] = {}
    for line in git("status", "--porcelain", "--untracked-files=all").splitlines():
        code, path = line[:2], line[3:].strip()
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        out[path] = code
    return out


def is_doc(rel: str) -> bool:
    return Path(rel).suffix in DOC_SUFFIXES


errors: list[str] = []
warnings: list[str] = []


def err(msg: str) -> None:
    errors.append(msg)


def warn(msg: str) -> None:
    warnings.append(msg)


def read(rel: str) -> str:
    try:
        return (ROOT / rel).read_text(encoding="utf-8")
    except (UnicodeDecodeError, FileNotFoundError):
        return ""


# ----------------------------------------------------------------------------
# check 1: genre violations in current-state documents
# ----------------------------------------------------------------------------


def check_genre(files: list[str]) -> None:
    for rel in files:
        if genre_of(rel) != "current-state" or not is_doc(rel):
            continue
        if rel in SPEC_FILES or rel in GENRE_EXEMPT:
            continue
        text = read(rel)
        lines = text.splitlines()
        for rx, why in BANNED:
            if rel.startswith(REVIEW_DIR) and why not in BLOCK_KINDS:
                continue
            for m in rx.finditer(text):
                line = text[: m.start()].count("\n") + 1
                snippet = lines[line - 1].strip() if line - 1 < len(lines) else ""
                if STAMP.search(snippet):
                    continue
                err(f"{rel}:{line}: {why}: {snippet[:90]!r}\n"
                    f"    A current-state document states what IS. Delete the narrative; git carries it.")


# ----------------------------------------------------------------------------
# check 2: prompt_logs are immutable
# ----------------------------------------------------------------------------


CLOSED = re.compile(r"^Status:", re.M)


def is_open_append(path: str) -> bool:
    """True when the committed entry is still open and the change only appends."""
    committed = git("show", f"HEAD:{path}")
    if not committed or CLOSED.search(committed):
        return False
    current = read(path)
    return current.startswith(committed)


def check_immutable() -> None:
    for line in git("status", "--porcelain", "--", "prompt_logs").splitlines():
        code, path = line[:2], line[3:].strip()
        if path.endswith("README.md"):
            continue
        if "M" in code and not ("D" in code or "R" in code) and is_open_append(path):
            continue
        if "M" in code or "D" in code or "R" in code:
            err(f"{path}: prompt log modified or deleted. Prompt logs are IMMUTABLE.\n"
                f"    Corrections go in a LATER entry, never in an old one.")


# ----------------------------------------------------------------------------
# check 3: public records are never edited
# ----------------------------------------------------------------------------


def check_public(changed: dict[str, str]) -> None:
    for rel, code in changed.items():
        if code == "??" or code.startswith("A"):
            continue
        if genre_of(rel) == "public-record":
            err(f"{rel}: a released version of the paper was modified or deleted.\n"
                f"    A public record is never edited; a correction goes into the next version.")


# ----------------------------------------------------------------------------
# check 4: dangling repository paths in index documents
# ----------------------------------------------------------------------------

PLACEHOLDER = re.compile(r"NNN|YYYY|MM-DD|_to_v|/v[^0-9a-z]|\.\.\.")

PATH_RX = re.compile(
    r"[`(\[]((?:paper|python_project|gap_project|evidence_and_reasoning|prompt_logs|"
    r"source_documents|tools)/[A-Za-z0-9_./{}-]*[A-Za-z0-9_}])"
)


def check_dangling() -> None:
    for rel in INDEX_DOCS:
        p = ROOT / rel
        if not p.exists():
            continue
        for i, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1):
            for m in PATH_RX.finditer(line):
                cand = m.group(1)
                if "{" in cand or "*" in cand:      # brace expansion / glob
                    continue
                if PLACEHOLDER.search(cand):        # NNN_, YYYY-MM-DD, v*_to_v*
                    continue
                if cand.endswith("/"):
                    ok = (ROOT / cand).is_dir()
                else:
                    ok = (ROOT / cand).exists() or (ROOT / cand).is_dir()
                if not ok:
                    err(f"{rel}:{i}: dangling path {cand!r} (no such file or directory)")


# ----------------------------------------------------------------------------
# check 5: every program is indexed
# ----------------------------------------------------------------------------

INDEXES = [
    ("python_project/src/*.py", "python_project/src/README.md"),
    ("gap_project/**/*.g", "gap_project/README.md"),
]
INDEX_SKIP = {"__init__.py"}


def check_index() -> None:
    for pattern, index in INDEXES:
        idx = ROOT / index
        if not idx.exists():
            continue
        text = idx.read_text(encoding="utf-8")
        for f in sorted(ROOT.glob(pattern)):
            if f.name in INDEX_SKIP:
                continue
            if f.name not in text:
                err(f"{index}: {f.name} exists but is not indexed.\n"
                    f"    An index README lists every file in its folder.")


# ----------------------------------------------------------------------------
# check 6: a changed document carries its dates, and Updated is the day of the change
# ----------------------------------------------------------------------------

SUBTITLE = re.compile(r"\\subtitle\{(\d{1,2}\s+\w+\s+20\d\d)\s*\(v\d+\)\}")


def mtime_date(rel: str) -> dt.date:
    return dt.date.fromtimestamp((ROOT / rel).stat().st_mtime)


def check_dates(changed: dict[str, str]) -> None:
    for rel, code in sorted(changed.items()):
        if "D" in code or not (ROOT / rel).is_file():
            continue
        if genre_of(rel) != "current-state" or not is_doc(rel):
            continue
        day = mtime_date(rel)
        if rel == PAPER_SOURCE:
            m = SUBTITLE.search(read(rel))
            if not m:
                err(f"{rel}: no \\subtitle{{D Month YYYY (vN)}} date found.")
            elif dt.datetime.strptime(m.group(1), "%d %B %Y").date() != day:
                err(f"{rel}: subtitle date {m.group(1)!r} is not the day of the last change ({day}).\n"
                    f"    Set the paper's date in the same edit as the change.")
            continue
        m = STAMP.search(read(rel))
        if not m:
            err(f"{rel}: changed, but carries no 'Created: YYYY-MM-DD.  Updated: YYYY-MM-DD.' line.\n"
                f"    A current-state document gets its dates at its next edit.")
            continue
        if dt.date.fromisoformat(m.group(2)) != day:
            err(f"{rel}: Updated {m.group(2)} is not the day of the last change ({day}).\n"
                f"    Check the date (`date`) and set Updated in the same edit as the change.")


# ----------------------------------------------------------------------------


def main() -> int:
    files = all_files()
    if "--list" in sys.argv:
        for rel in files:
            print(f"{genre_of(rel):15s} {rel}")
        return 0

    changed = changed_files()
    check_genre(files)
    check_immutable()
    check_public(changed)
    check_dangling()
    check_index()
    check_dates(changed)

    for w in warnings:
        print(f"warning: {w}")
    for e in errors:
        print(f"ERROR: {e}")

    n, m = len(errors), len(warnings)
    if n:
        print(f"\n{n} error(s), {m} warning(s). See DOCUMENT_GENRES.md.")
        return 1
    print(f"documentation clean ({m} warning(s)).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
