# Trial 003 Results — Discrete variants of triple-octonion cross-block rule

Created: 2026-04-09.  Updated: 2026-10-28.

Date: 2026-04-09

## Algebra family tested

Same triple-octonion base (O₁ ⊕ O₂ ⊕ O₃) as trials 001–002, but with three
discrete degrees of freedom varied in the cross-block terms:

- **Conjugation**: for each of 3 cross-block pairs, 4 choices (plain, conj left,
  conj right, conj both): 4³ = 64 combinations.
- **Sign**: for each cross-block pair, ±1: 2³ = 8 combinations.
- **Routing**: where cross-block products land (third block, left block, or
  right block): 3 choices.

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
- 68.8%: 31 variants (24 with third-block, 4 with left-block and 3 with
  right-block routing)
- 67.8% to 68.7%: 41 variants
- 65.6% to 66.0%: 24 variants
- 64.0% to 64.4%: 48 variants
- Below 64%: the remaining 1,392 variants

### Key structural observations

**1. Sign flips do not change the failure count.** With conjugation and
routing held fixed, all 8 sign patterns give the same number of failures. The
closure rates can still differ, because some sign patterns produce zero
products, which are not counted.

**2. Routing does not help.** Third-block routing reaches the best rate,
68.8%, in 24 variants; left-block and right-block routing reach it in 4 and
3 variants, and no routing does better.

**3. Conjugation has limited effect.** The best variants have *uniform*
conjugation across all three pairs: all plain, all conj-left or all
conj-right. All conj-both reaches at most 66.0%, and no mixed pattern
reaches 68.8%.

**4. Failure mode is identical to trial 001.** The detailed Wilson-condition
analysis of the best variant shows:
- Failures: exclusively type3×type3 (78/78)
- Violated condition: exclusively x+y+z ∈ Ls (100%)
- Norm² distribution: 32–104, peaked at 64

This is the exact same failure signature as trial 001. No discrete
modification to the cross-block rule changes the fundamental obstruction.

## Interpretation

The type3×type3 failure on condition 3 (x+y+z ∈ Ls) is a **property of the
standard octonion product in the three blocks**, not an artifact of the
specific cross-block rule. Conjugation, sign flips and routing changes, the
three natural discrete modifications, cannot fix it. For the base algebra
the reason is exact (key claim 007): the block sum of a product is the
product of the block sums, and Ls is not closed under the standard product
(paper Remark 4.5). With the σ-twisted product in all three blocks the
triple product closes (trial 007; paper Theorem 1.1).

The insensitivity to sign flips points the same way: the failure is not
caused by constructive or destructive interference between cross-block
contributions.

## Verdict

**FAIL.** All 1,536 discrete variants of the triple-octonion cross-block rule
fail, every one of them on type3×type3 pairs. The failure belongs to the
standard octonion product in the three blocks: with the σ-twisted product the
same architecture closes (trial 007).

## What this rules out

Any triple-octonion algebra (O₁ ⊕ O₂ ⊕ O₃) with:
- Any conjugation pattern on cross-block terms (4³ = 64)
- Any sign pattern on cross-block terms (2³ = 8)
- Any of three routing rules (third/left/right block)
