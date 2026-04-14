# Key Claims

This folder contains one file per major claim, lemma, or reasoning step in the research. Together these files form the logical backbone of the investigation.

## Purpose

Breaking the argument into discrete, documented claims allows:
- Independent verification of each step
- Clear identification of which steps are established facts vs. conjectural
- A structured foundation for writing the formal paper

## Naming Convention

Files are named descriptively:

```
NNN_short_claim_name.txt
```

The three-digit prefix preserves logical ordering where a claim depends on earlier ones.

## Status Labels

Each claim file should carry one of the following status labels:

- `ESTABLISHED` — Follows directly from cited references with no novel reasoning required
- `DERIVED` — Follows from established claims by reasoning documented in this project; second-pass verified
- `CONJECTURE` — Stated as a hypothesis; not yet proven
- `UNCERTAIN` — Flagged for review; do not build further claims on this until resolved

## Contents

| File | Claim | Status |
|---|---|---|
| [001_wilson_dixon_multiplication_tables_compared.txt](001_wilson_dixon_multiplication_tables_compared.txt) | Wilson and Dixon use identical octonion multiplication tables (up to index relabelling) | ESTABLISHED |
| [002_wilson_e8_verified.txt](002_wilson_e8_verified.txt) | Wilson's 240-element octonionic construction is a valid E8 lattice | ESTABLISHED |
| [003_dixon_e8_verified.txt](003_dixon_e8_verified.txt) | Dixon's two 240-element constructions (Ξ^even, A^odd) are valid E8 minimal-vector shells | ESTABLISHED |
| [004_wilson_leech_verified.txt](004_wilson_leech_verified.txt) | Wilson's three families (720+11520+184320=196560) are the Leech lattice minimal vectors | ESTABLISHED |
| [005_dixon_leech_verified.txt](005_dixon_leech_verified.txt) | Dixon's three families (720+11520+184320=196560) form a different embedding of the Leech minimal shell; 17,232 vectors shared with Wilson's | ESTABLISHED |
| [006_double_check_wilson_dixon.txt](006_double_check_wilson_dixon.txt) | Double-check: Wilson fully established; Dixon structural verification only (even/self-dual gaps); order question explicitly OPEN | ESTABLISHED / OPEN |
| [007_triple_octonion_ruled_out.txt](007_triple_octonion_ruled_out.txt) | Triple-octonion algebra (O₁ ⊕ O₂ ⊕ O₃) comprehensively ruled out as order on Λ — 4 trials, all discrete/continuous variants fail on type3×type3 condition 3 | ESTABLISHED |
| [008_transposition_twist_order.txt](008_transposition_twist_order.txt) | Transposition-twisted triple octonion product (·_σ applied in all three blocks with Z₃ routing) is an order on Λ — symbolic proof via 4 lemmas, plus 12.5M+ random-pair verifications | ESTABLISHED |
