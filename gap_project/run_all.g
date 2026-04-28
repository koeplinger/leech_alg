#############################################################################
##
##  run_all.g — Driver script for the GAP/LOOPS verification suite.
##
##  Usage (from the repository root):
##      gap -q gap_project/run_all.g
##
##  Or, from any directory:
##      gap -q -c 'GAP_PROJECT_ROOT := "/abs/path/to/leech_alg/gap_project";' \
##              gap_project/run_all.g
##
##  Loads the LOOPS package, then runs each tests/*.g in sequence, tracking
##  pass/fail counts.  Prints a final summary and exits with code 0 if all
##  tests pass, non-zero otherwise.
##
##  Requires:  GAP 4.x and the LOOPS package (`apt install gap` plus the
##  LOOPS package via gap_packages, or install via GAP's PackageManager).
##
#############################################################################

# Determine GAP_PROJECT_ROOT: either set by the caller, or derived from this
# file's location.  For simplicity we require either:
#   (a) the working directory equals the repository root, in which case
#       "gap_project" is the relative path; OR
#   (b) GAP_PROJECT_ROOT is already set as an absolute path.

if not IsBound(GAP_PROJECT_ROOT) then
    GAP_PROJECT_ROOT := "gap_project";
fi;

LoadPackage("loops");

# Pre-bind harness counters so the syntax-checker for RunTestFile sees them.
PASSED := 0;
FAILED := 0;
FAIL_DESCS := [];

GLOBAL_PASS := 0;
GLOBAL_FAIL := 0;
GLOBAL_FAIL_DESCS := [];

RunTestFile := function(name)
    local file;
    file := Concatenation(GAP_PROJECT_ROOT, "/tests/", name);
    # Each test file calls ResetTestCounters() at the top, so after Read
    # the counters PASSED, FAILED, FAIL_DESCS hold the totals for that file.
    Read(file);
    GLOBAL_PASS := GLOBAL_PASS + PASSED;
    GLOBAL_FAIL := GLOBAL_FAIL + FAILED;
    Append(GLOBAL_FAIL_DESCS, FAIL_DESCS);
end;

# Load harness once so PASSED/FAILED/FAIL_DESCS exist.  Each test file calls
# ResetTestCounters() so counters are local to that file; we accumulate the
# delta into GLOBAL_*.
Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
ResetTestCounters();

RunTestFile("test_octonion.g");
RunTestFile("test_e8_wilson.g");
RunTestFile("test_sublattices.g");
RunTestFile("test_twist.g");
RunTestFile("test_lemmas.g");
RunTestFile("test_leech_wilson.g");
RunTestFile("test_triple_product.g");
RunTestFile("test_companion_examples.g");

Print("\n##############################################################\n");
Print("##  GAP/LOOPS verification suite — final summary\n");
Print("##############################################################\n");
Print("Total tests passed: ", GLOBAL_PASS, "\n");
Print("Total tests failed: ", GLOBAL_FAIL, "\n");
if GLOBAL_FAIL > 0 then
    Print("Failed tests:\n");
    Perform(GLOBAL_FAIL_DESCS, function(d) Print("  - ", d, "\n"); end);
    Print("##############################################################\n");
    QUIT_GAP(1);
fi;
Print("##############################################################\n");
QUIT_GAP(0);
