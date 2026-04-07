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
