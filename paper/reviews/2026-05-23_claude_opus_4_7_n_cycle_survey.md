# Exact-arithmetic n-cycle closure survey

**Date:** 2026-05-23
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 153 — double-check the 4-cycle total failure noted in
`evidence_and_reasoning/trial_007_results.md` (Phase 3 survey, 0/5).
**Script:** `python_project/src/verify_all_cycles_exact.py`.

---

## The question

`evidence_and_reasoning/trial_007_results.md` Phase 3 reports, from a
small-sample (~5 per cycle type) survey:

| Cycle type | Tested | Worked |
|---|---|---|
| 4-cycles | 0 / 5 | 0 % |
| 5-cycles | 4 / 5 | 80 % |
| 7-cycles | 1 / 5 | 20 % |

The feedback flagged the 4-cycle total failure as a striking outlier.
This note confirms (or refutes) it by exhaustive exact-arithmetic
enumeration: for each n-cycle π in S_7, the π-twisted product is tested
on all 24 × 24 = 576 basis-pair products of a Z-basis of Λ; π closes on
Λ iff all 576 products lie in Λ.

## Method

For each n ∈ {3, 4, 5, 6, 7}, generate **every** n-cycle π in S_7 (acting
on indices {1,...,7}, fixing 0 = the octonion identity). Apply the
twisted-product closure test from `verify_consecutive_twists_exact.py`
(the same `closes(perm)` function used for the paper's footnote).

Counts of n-cycles in S_7 (k-cycles in S_n: C(n,k) · (k−1)!):
- 3-cycles: 35 · 2 = 70
- 4-cycles: 35 · 6 = 210
- 5-cycles: 21 · 24 = 504
- 6-cycles: 7 · 120 = 840
- 7-cycles: 1 · 720 = 720
- **Total: 2344**

## Result

```
cycle type      tested     close       pct
------------------------------------------
3-cycles            70        35    50.00%
4-cycles           210         0     0.00%
5-cycles           504       252    50.00%
6-cycles           840       210    25.00%
7-cycles           720       336    46.67%
------------------------------------------
total             2344       833    35.54%
```

**The 4-cycle total failure is real and exhaustive: 0 of 210 close.**
The trial_007 Phase 3 result (0/5) generalises to 0/210 — not a
sampling artifact.

## Striking aspects

1. **4-cycles are the unique pure-n-cycle type with 0 % closure.** Every
   other cycle length in S_7 has at least 25 % closure.
2. **3-cycles and 5-cycles both give exactly 50 %.** (35/70 and 252/504.)
3. **6-cycles give exactly 25 %** (210/840 = 1/4).
4. **7-cycles give exactly 7/15 ≈ 46.67 %** (336/720).

The 4-cycle result is the only one that is *categorical* (uniformly
fails) rather than fractional.

## Structural notes that may help an analytic argument

(Empirical only — not a proof.)

### (i) The square of a 4-cycle is a (2,2)-double on the same support

For π = (a b c d), π² = (a c)(b d), which is a (2,2)-double on
{a, b, c, d} fixing the same three indices {1,...,7} \ {a, b, c, d} that
π fixes. So a 4-cycle inherits the structure of its (2,2)-square.

### (ii) The Fano-line rule for (2,2)-doubles, lifted to 4-cycles

The (2,2)-double rule (from the consecutive-twists analysis): a
(2,2)-double whose fixed triple forms a Fano line never closes; the
other (2,2)-doubles split 42 of 84 closing.

There are C(7,4) = 35 four-element subsets of {1,...,7}, of which 7
have complement forming a Fano line. Each support carries
4!/4 = 6 distinct 4-cycles, so:
- 7 supports × 6 cycles = **42 four-cycles** whose square is a
  Fano-line-fixing (2,2)-double — these would inherit non-closure.
- 28 supports × 6 cycles = **168 four-cycles** whose square is a
  non-Fano-line (2,2)-double — half of those (2,2)-doubles *do*
  close on Λ.

If 4-cycle closure were only constrained by π²-closure, one would
expect roughly 84 of 168 of the non-Fano 4-cycles to close. The
exhaustive count finds **0**. So the 4-cycle has a *strictly
additional* obstruction beyond what π² inherits — closing the
(2,2)-square is necessary but not sufficient.

### (iii) Parity / order observations

| Cycle | Sign | Order | Closure |
|---|---|---|---|
| transposition | odd | 2 | 21 / 21 |
| 3-cycle | even | 3 | 35 / 70 |
| 4-cycle | odd | 4 | 0 / 210 |
| 5-cycle | even | 5 | 252 / 504 |
| 6-cycle | odd | 6 | 210 / 840 |
| 7-cycle | even | 7 | 336 / 720 |

Neither sign nor order alone predicts the closure pattern (e.g.,
transposition and 4-cycle are both odd; one closes universally and the
other never does).

### (iv) What a structural argument would need to show

For a 4-cycle π = (a b c d), the twisted product fails to close on Λ
iff there exist u, v ∈ Λ with u ⋆_π v ∉ Λ. Working through the
analog of paper §4:

- Lemma 4.1 analog: π(L) = L. Holds for any S_7 element. ✓
- Lemma 4.2 analog: L · L ⊆ L. Independent of π. ✓
- Lemma 4.3 analog (condition 2): L · π(Lsbar) ⊆ π(Lsbar). Depends on π.
- Lemma 4.4 analog (condition 3): π(Ls) · π(Ls) ⊆ π(Ls). Depends on π.

A structural argument for 4-cycle total failure would show that one of
the latter two inclusions *always* fails when π is a 4-cycle. The script
`verify_consecutive_twists_exact.py` already noted that the
condition-2 inclusion is *sufficient but not necessary*; the genuine
failure mode is the joint Λ-membership test, not either sublattice
inclusion alone. So the structural argument needs to operate on Λ
directly (or on a finer invariant that captures Λ).

The mod-2 picture of Appendix A (V = σ(Lsbar)/2L is a left ideal,
W = σ(Ls)/2L is a subalgebra of L/2L) is the natural place to look:
for which π is π(Ls)/2L still a subalgebra? For which π is L · π(Lsbar)
still in π(Lsbar) mod 2L? If 4-cycles uniformly violate one of these
mod-2 conditions, that would give the structural argument.

## Cross-references

- `evidence_and_reasoning/trial_007_results.md`: Phase 3 small-sample
  survey (the 0/5 figure that prompted this check).
- `python_project/src/verify_consecutive_twists_exact.py`: the
  `closes()` test reused here.
- `paper/reviews/2026-05-22_claude_opus_4_7_consecutive_twists.md`:
  earlier analysis covering products of two transpositions (3-cycles
  and (2,2)-doubles).
- `prompt_logs/153_4cycle_failure_check.txt`: prompt log.

## Status

**Empirical:** 4-cycle total failure (0/210) is verified by exhaustive
exact-arithmetic enumeration on a Z-basis of Λ. Not a sampling artifact.

**Structural argument:** open; user will look at the parallel between
3-cycles and 4-cycles. Starting points (i)–(iv) above identify where
the additional 4-cycle obstruction lives beyond the inherited
(2,2)-square structure.
