#!/usr/bin/env python3
"""Documentation linter.  Enforces DOCUMENT_GENRES.md mechanically.

No LLM, no network, no dependencies.  Runs in about a second.

Checks
------
1. GENRE       Current-state documents must not narrate their own history:
               no dated correction notices, no settlement dates, no
               supersession narratives, no change-history footers.
2. IMMUTABLE   Files under prompt_logs/ that are already tracked by git must
               not be modified or deleted.  (New entries are fine.)
3. FROZEN      Frozen artifacts must not be rewritten in place: if one is
               modified relative to HEAD, the diff must ADD a dated addendum
               (or a deprecation marker), and must not merely rewrite prose.
4. DANGLING    Every repository path named in a current-state document must
               exist.
5. INDEX       Every program in python_project/src/ and gap_project/ must
               appear in the README that indexes it.

Usage
-----
    python3 tools/lint_docs.py            # check
    python3 tools/lint_docs.py --list     # print the genre of every file

Exit code 0 if clean, 1 if any error was found.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# ----------------------------------------------------------------------------
# Genre membership (see DOCUMENT_GENRES.md)
# ----------------------------------------------------------------------------

IMMUTABLE_DIRS = ["prompt_logs"]

FROZEN_GLOBS = [
    "CLOSING_NOTE.md",
    "paper/main_*.tex",
    "paper/companion.tex",
    "paper/automorphism_group_*.tex",
    "paper/hermes_mod2_structure_supplement.tex",
    "paper/update_brief_*.tex",
    "paper/kirmse_1924_exposition.tex",
    "paper/historical_appendix*.tex",
    "paper/Baez-Egan_*.tex",
    "paper/reviews/*.md",
    "paper/v*_summary.md",
    "paper/v*_changelog.md",
    "evidence_and_reasoning/key_claims/[0-9]*.txt",
    "evidence_and_reasoning/trial_*_results.md",
    "evidence_and_reasoning/[0-9]*_plan.md",
    "python_project/src/*.py",
    "python_project/tests/*.py",
    "gap_project/**/*.g",
]

CURRENT_STATE_FILES = [
    "MANIFESTO.md",
    "DOCUMENT_GENRES.md",
    "TRIAL_METHODOLOGY.md",
    "CURRENT_STATE.md",
    "README.md",
    "LICENSE.md",
    "LICENSE-CODE",
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

# research_statement.md is a special case: see DOCUMENT_GENRES.md.
EXEMPT = {"evidence_and_reasoning/research_statement.md"}

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
]
BANNED = [(re.compile(p, re.M), why) for p, why in BANNED]

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


def frozen_files() -> set[str]:
    out: set[str] = set()
    for pat in FROZEN_GLOBS:
        out |= {str(p.relative_to(ROOT)) for p in ROOT.glob(pat) if p.is_file()}
    return out


def genre_of(rel: str) -> str:
    if any(rel.startswith(d + "/") for d in IMMUTABLE_DIRS) and not rel.endswith("README.md"):
        return "immutable"
    if rel in EXEMPT:
        return "exempt"
    if rel in CURRENT_STATE_FILES:
        return "current-state"
    if rel in frozen_files():
        return "frozen"
    return "unclassified"


errors: list[str] = []
warnings: list[str] = []


def err(msg: str) -> None:
    errors.append(msg)


def warn(msg: str) -> None:
    warnings.append(msg)


# ----------------------------------------------------------------------------
# check 1: genre violations in current-state documents
# ----------------------------------------------------------------------------


# A rule book has to be able to name the thing it forbids.  These three
# documents DEFINE the genre vocabulary, so they are exempt from check 1
# (and from check 1 only).
SPEC_FILES = {"DOCUMENT_GENRES.md", "MANIFESTO.md", "tools/README.md"}


def check_genre() -> None:
    for rel in CURRENT_STATE_FILES:
        if rel in SPEC_FILES:
            continue
        p = ROOT / rel
        if not p.exists():
            continue
        text = p.read_text(encoding="utf-8")
        for rx, why in BANNED:
            for m in rx.finditer(text):
                line = text[: m.start()].count("\n") + 1
                snippet = text.splitlines()[line - 1].strip()[:90]
                err(f"{rel}:{line}: {why}: {snippet!r}\n"
                    f"    A current-state document states what IS. Delete the narrative; git carries it.")


# ----------------------------------------------------------------------------
# check 2: prompt_logs are immutable
# ----------------------------------------------------------------------------


def check_immutable() -> None:
    for line in git("status", "--porcelain", "--", "prompt_logs").splitlines():
        code, _, path = line[:2], line[2], line[3:].strip()
        if path.endswith("README.md"):
            continue
        if "M" in code or "D" in code or "R" in code:
            err(f"{path}: prompt log modified or deleted. Prompt logs are IMMUTABLE.\n"
                f"    Corrections go in a LATER entry, never in an old one.")


# ----------------------------------------------------------------------------
# check 3: frozen artifacts are not rewritten in place
# ----------------------------------------------------------------------------

ADDENDUM = re.compile(
    r"(FORWARD CORRECTION|ADDENDUM|CORRECTION|Note added|DEPRECATED|SUPERSEDED)",
    re.I,
)

# A version changelog freezes only "on release of the later version"
# (DOCUMENT_GENRES.md).  While it declares itself in progress it is a living
# document, edited in place like the current-state docs; the marker is removed
# when its version is frozen, at which point the frozen check re-engages.
IN_PROGRESS = re.compile(r"in progress, not frozen|not yet frozen", re.I)


def check_frozen() -> None:
    changed = [
        l[3:].strip()
        for l in git("status", "--porcelain").splitlines()
        if l[:2].strip() in {"M", "MM"}
    ]
    fz = frozen_files()
    for rel in changed:
        if rel not in fz:
            continue
        if rel.endswith(("_changelog.md", "_summary.md")) and (
            IN_PROGRESS.search((ROOT / rel).read_text(encoding="utf-8"))
            or IN_PROGRESS.search(git("show", f"HEAD:{rel}"))
        ):
            # In-progress changelog: a living document until its version
            # freezes.  The working-text match covers the development phase;
            # the HEAD-text match covers the one-time freeze transition (the
            # edit that drops the marker and stamps the freeze).  Once frozen
            # is committed, neither has the marker and the check re-engages.
            continue
        diff = git("diff", "HEAD", "--", rel)
        added = [l[1:] for l in diff.splitlines() if l.startswith("+") and not l.startswith("+++")]
        removed = [l[1:] for l in diff.splitlines() if l.startswith("-") and not l.startswith("---")]
        has_addendum = any(ADDENDUM.search(l) for l in added)
        if not has_addendum and removed:
            warn(f"{rel}: frozen artifact rewritten in place "
                 f"(+{len(added)}/-{len(removed)} lines, no dated addendum).\n"
                 f"    Frozen artifacts take a DATED ADDENDUM, or are marked deprecated/superseded.")


# ----------------------------------------------------------------------------
# check 4: dangling repository paths in current-state documents
# ----------------------------------------------------------------------------

PLACEHOLDER = re.compile(r"NNN|YYYY|MM-DD|_to_v|/v[^0-9a-z]|\.\.\.")

PATH_RX = re.compile(
    r"[`(\[]((?:paper|python_project|gap_project|evidence_and_reasoning|prompt_logs|"
    r"source_documents|tools)/[A-Za-z0-9_./{}-]*[A-Za-z0-9_}])"
)


def check_dangling() -> None:
    for rel in CURRENT_STATE_FILES:
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
    ("gap_project/*.g", "gap_project/README.md"),
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


def main() -> int:
    if "--list" in sys.argv:
        seen = set()
        for p in sorted(ROOT.rglob("*")):
            if not p.is_file() or ".git" in p.parts:
                continue
            rel = str(p.relative_to(ROOT))
            g = genre_of(rel)
            if g != "unclassified" and rel not in seen:
                print(f"{g:15s} {rel}")
                seen.add(rel)
        return 0

    check_genre()
    check_immutable()
    check_frozen()
    check_dangling()
    check_index()

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
