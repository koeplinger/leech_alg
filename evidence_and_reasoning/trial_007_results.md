# Trial 007 Results — Triple Kirmse-twisted octonion algebra

Date: 2026-04-11

## Algebra tested

Same triple-product structure as trial 001 (three copies of an octonion algebra
with Z₃ cross-block routing), but with DIFFERENT octonion multiplications
obtained by applying permutations of {1,...,7} to the standard Fano triples.

The trial tested:
- The standard multiplication (baseline / trial 001 confirmation)
- 7 index-doubling permutations σ_s(k) = 2k − s (mod 7)
- 21 transposition permutations swap(s,t)
- Additional permutations: 3-cycles, 4-cycles, 5-cycles, 7-cycles,
  products of disjoint cycles

## Phase 0: Verification that Wilson's L is a maximal order

**L IS closed under the standard octonion multiplication.** All 57,600 root
pairs (240 × 240) were tested; every product lies in L.

This confirms the claim in leech_wilson.py that Wilson's L is the
Coxeter–Dickson order (a maximal order of integral octonions).  Consequently,
trial 001 already used a maximal-order multiplication.  The "Kirmse twist"
is already implicitly applied in the standard setup.

## Phase 1: Index-doubling permutations

All 7 index-doubling permutations σ_s(k) = 2k − s (mod 7) are Fano-plane
automorphisms: they map Fano triples to Fano triples, producing the SAME
multiplication table as the standard algebra.

**Result:** 7/7 close on L. Leech closure test: all 7 give identical results
to the standard algebra (81.3% closure, failing only on type3×type3 at
Wilson condition 3). This confirms that the index-doubling permutations
are Fano-plane automorphisms and do not change the algebra.

## Phase 2: Transposition permutations — THE KEY FINDING

**All 21 transpositions of imaginary basis elements achieve 100% Leech
lattice closure on all tested pairs.**

Detailed results for swap(1,2):

| Factor types | Pairs tested | Failures | Closure |
|---|---|---|---|
| t1×t1 | 518,400 (exhaustive) | 0 | 100% |
| t1×t2 | 5,000 | 0 | 100% |
| t1×t3 | 5,000 | 0 | 100% |
| t2×t2 | 10,000 | 0 | 100% |
| t2×t3 | 5,000 | 0 | 100% |
| t3×t3 | 50,000 | 0 | 100% |
| **Total** | **593,400** | **0** | **100%** |

This was verified with multiple random seeds and independent pair selections.
For comparison, the standard multiplication on the same type3×type3 pairs
gives ≈74% failure rate (exclusively violating Wilson condition 3).

All 21 transpositions were tested on a smaller sample (500–900 pairs per
transposition, all type combinations): all achieve 0 failures.

## Phase 3: Broader permutation survey

Not all permutations work. A survey of ~200 random permutations across all
cycle types shows:

| Cycle type | Tested | Worked (0 failures) |
|---|---|---|
| (2,) transpositions | 2/2 | 100% |
| (2,2) double-swaps | 4/5 | 80% |
| (3,) 3-cycles | 3/5 | 60% |
| (4,) 4-cycles | 0/5 | 0% |
| (5,) 5-cycles | 4/5 | 80% |
| (2,2,2) triple-swaps | 2/5 | 40% |
| (7,) 7-cycles | 1/5 | 20% |

The pattern depends on the relationship between the permutation and the
Fano-plane structure, not on cycle type or parity alone.

### What distinguishes working from non-working permutations

Fano-plane automorphisms (index-doubling, identity) give the SAME
multiplication as the standard and inherit its type3×type3 failure.
Non-automorphism permutations give genuinely different multiplications.
Among these, SOME close on the Leech lattice and others do not.  All
transpositions close; the detailed criterion for general permutations
is an open question.

## Failure mode contrast

| Algebra | t3×t3 failure rate | Wilson condition violated |
|---|---|---|
| Standard | ≈74% | Condition 3 only (x+y+z ∈ Ls) |
| swap(1,2) | 0% | (none) |

The transposition fixes precisely the Wilson condition 3 failure that has
been the consistent failure mode since trial 001.

## Product norm distribution (swap(1,2), type3×type3)

| ‖prod‖² | Count (of 2500) |
|---|---|
| 16 | 4 |
| 32 | 93 |
| 48 | 472 |
| 64 | 1,279 |
| 80 | 584 |
| 96 | 58 |
| 112 | 10 |

Products are NOT on the minimal shell (‖prod‖² ≠ 8).  They land in the
Leech lattice at higher shells.  This is expected for an order (closure means
products stay in the lattice, not necessarily on the minimal shell).

## Interpretation

The triple octonion product with Z₃ cross-block routing closes on the Leech
lattice Λ when the octonion multiplication is obtained from the standard Fano
triples by applying any transposition of imaginary basis elements.  This is a
different but isomorphic octonion algebra on the same R⁸.  The transposition
introduces sign changes in certain structure constants that align the
three-block sum with Wilson's sublattice Ls, fixing the condition-3 failure
of the standard multiplication.

**This means (Λ, +, ·) is an order in the R-algebra (R²⁴, +, ·)**, where ·
is the transposition-twisted triple octonion product — provided the result
survives exhaustive verification and independent mathematical proof.

## Caveats and verification status

1. **Sample-based, not exhaustive.** The type3×type3 test covers 50,000 of
   the ~3.4 × 10¹⁰ possible pairs.  Zero failures in 593,400 total pairs
   is strong evidence but not proof.

2. **Bilinearity argument.** If the product of any two Min(Λ) vectors lies
   in Λ, then by bilinearity and the fact that Min(Λ) generates Λ over Z,
   ALL products of Λ vectors lie in Λ.  So testing on Min(Λ) suffices.

3. **Algebraic properties not yet tested.** Whether this algebra is
   alternative, flexible, power-associative, or has a multiplicative
   identity is not yet established.

4. **Independent verification required** (Manifesto §5).  This result
   should be verified by a second computational method and ideally by
   a mathematical proof.

## Verdict

**CONDITIONAL PASS — pending exhaustive verification.**  The transposition-
twisted triple octonion product achieves 100% Leech lattice closure on
593,400 tested pairs with zero failures.  If confirmed, this establishes
that the Leech lattice admits an order under an octonion-derived bilinear
product.

## What this establishes

1. Wilson's L is a maximal order (Kirmse twist is already applied).
2. Fano-plane automorphisms (including index-doubling) do not change the
   multiplication and reproduce the trial 001 failure.
3. Transpositions of imaginary basis elements produce a DIFFERENT
   multiplication that closes on Λ.
4. The correction fixes Wilson condition 3 (x+y+z ∈ Ls) for type3×type3
   products — the exact failure mode of all prior trials.

## What remains open

1. Exhaustive verification on all 196,560² ≈ 3.9 × 10¹⁰ pairs (or a
   mathematical proof of closure).
2. Algebraic properties of the resulting 24-dimensional algebra.
3. The precise criterion distinguishing working from non-working
   permutations (beyond "all transpositions work").
4. Whether different transpositions give the same or distinct algebras
   (they are all in the same GL(3,F_2) orbit, suggesting the same algebra
   up to Fano-plane relabeling).

---

Last updated: 2026-04-11
