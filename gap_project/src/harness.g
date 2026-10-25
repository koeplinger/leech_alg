#############################################################################
##
##  harness.g — Tiny test harness used by the test scripts in tests/.
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
##
##  Usage in a test file (as in tests/*.g):
##      Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
##      LoadAllSrc(GAP_PROJECT_ROOT);
##      ResetTestCounters();
##      CHECK("description of the check", boolean_expression);
##      ...
##      PrintTestSummary();
##  PrintTestSummary() prints the totals of the test file.  run_all.g adds
##  up the totals of all test files, prints the sum, and exits with code 0
##  if every check passed and 1 otherwise.
##
#############################################################################

PASSED := 0;
FAILED := 0;
FAIL_DESCS := [];

CHECK := function(description, ok)
    if ok = true then
        PASSED := PASSED + 1;
    else
        FAILED := FAILED + 1;
        Add(FAIL_DESCS, description);
        Print("  FAIL: ", description, "\n");
    fi;
end;

ResetTestCounters := function()
    PASSED := 0;
    FAILED := 0;
    FAIL_DESCS := [];
end;

PrintTestSummary := function()
    Print("\n--------------------------------------------------------------\n");
    Print("Tests passed: ", PASSED, "    failed: ", FAILED, "\n");
    if FAILED > 0 then
        Print("Failed tests:\n");
        Perform(FAIL_DESCS, function(d) Print("  - ", d, "\n"); end);
    fi;
    Print("--------------------------------------------------------------\n");
end;

# Convenience: load all src/ modules in dependency order.
LoadAllSrc := function(root)
    Read(Concatenation(root, "/src/octonion.g"));
    Read(Concatenation(root, "/src/e8_wilson.g"));
    Read(Concatenation(root, "/src/sublattices.g"));
    Read(Concatenation(root, "/src/twist.g"));
    Read(Concatenation(root, "/src/leech_wilson.g"));
    Read(Concatenation(root, "/src/triple_product.g"));
end;
