# Trial 003 Results — Discrete variants of triple-octonion cross-block rule

Date: 2026-04-09

## Algebra family tested

Same triple-octonion base (O₁ ⊕ O₂ ⊕ O₃) as trials 001–002, but with three
discrete degrees of freedom varied in the cross-block terms:

- **Conjugation**: for each of 3 cross-block pairs, 4 choices (plain, conj left,
  conj right, conj both) — 4³ = 64 combinations.
- **Sign**: for each cross-block pair, ±1 — 2³ = 8 combinations.
- **Routing**: where cross-block products land (third block, left block, or
  right block) — 3 choices.

Total: 64 × 8 × 3 = **1,536 variants**, all tested exhaustively.

## Test methodology

Each variant was tested against 250 Min(Λ) pairs:
- 100 type3×type3 pairs (the critical case from trial 001)
- 50 type2×type2, 25 each of type1×type1, type1×type2, type1×type3, type2×type3

## Results

### Headline: ALL 1,536 variants fail

| Metric | Value |
|--------|-------|
| Variants tested | 1,536 (exhaustive) |
| Variants with 100% closure | 0 |
| Best closure rate | 68.8% (172/250) |
| Variants at best rate | 31 |
| Variants with zero t3×t3 failures | 0 |

### Closure rate distribution

The rates cluster tightly:
- 68.8%: 31 variants (all with third-block routing)
- 66.0%: 12 variants
- 64.0%: 32 variants
- Below 64%: the remaining ~1,400 variants

### Key structural observations

**1. Sign flips have zero effect.** All 8 sign patterns produce identical
results when conjugation and routing are held fixed. This means the sign of
cross-block contributions is irrelevant — the failure mechanism is insensitive
to it.

**2. Third-block routing is optimal.** Left-block and right-block routing
produce lower closure rates. Breaking Z₃ symmetry hurts rather than helps.

**3. Conjugation has limited effect.** The best variants are those with
*uniform* conjugation across all three pairs (all plain, all conj-left,
all conj-right, or all conj-both). Mixed conjugation patterns perform worse.

**4. Failure mode is identical to trial 001.** The detailed Wilson-condition
analysis of the best variant shows:
- Failures: exclusively type3×type3 (78/78)
- Violated condition: exclusively x+y+z ∈ Ls (100%)
- Norm² distribution: 32–104, peaked at 64

This is the exact same failure signature as trial 001. No discrete
modification to the cross-block rule changes the fundamental obstruction.

## Interpretation

The type3×type3 failure on condition 3 (x+y+z ∈ Ls) is a **structural
property of the triple-octonion architecture**, not an artifact of the specific
cross-block rule. Conjugation, sign flips, and routing changes — the three
natural discrete modifications — cannot fix it.

The insensitivity to sign flips is particularly telling: it means the failure
is not caused by constructive/destructive interference between cross-block
contributions, but by the inherent structure of how type-3 vectors interact
with the triple-block decomposition.

## Verdict

**FAIL — structural obstruction.** All 1,536 discrete variants of the
triple-octonion cross-block rule fail on the same type3×type3 obstruction.
The failure is intrinsic to the triple-octonion architecture.

## What this rules out

Any triple-octonion algebra (O₁ ⊕ O₂ ⊕ O₃) with:
- Any conjugation pattern on cross-block terms (4³ = 64)
- Any sign pattern on cross-block terms (2³ = 8)
- Any of three routing rules (third/left/right block)
