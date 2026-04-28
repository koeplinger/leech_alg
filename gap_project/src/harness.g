#############################################################################
##
##  harness.g — Tiny test harness used by the test scripts in tests/.
##
##  Usage in a test file:
##      Read("src/harness.g");                     # initialises counters
##      TEST("description of the check", boolean_expression);
##      ...
##  At the end of run_all.g (or a single test file), the totals are printed
##  via PrintTestSummary(); a non-zero failure count causes a non-zero exit
##  (when running as a script, see run_all.g).
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
