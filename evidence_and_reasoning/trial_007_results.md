# Trial 007 Results — Triple Kirmse-twisted octonion algebra

Created: 2026-04-11.  Updated: 2026-10-28.

## Algebra tested

Same triple-product structure as trial 001 (three copies of an octonion algebra
with Z₃ cross-block routing), but with DIFFERENT octonion multiplications
obtained by applying permutations of {1,...,7} to the standard Fano triples.
The "Kirmse twist" of the title is a transposition of two imaginary basis
units, a linear involution of R⁸: the map by which Bruck repaired Kirmse's
J₁ (Coxeter 1946, §4: "transposing two of the i's").

The trial tested:
- The standard multiplication (baseline / trial 001 confirmation)
- 7 index-doubling permutations σ_s(k) = 2k − s (mod 7)
- 21 transposition permutations swap(s,t)
- Additional permutations: 3-cycles, 4-cycles, 5-cycles, 7-cycles,
  products of disjoint cycles

## Phase 0: Verification that Wilson's L is closed

**L IS closed under the standard octonion multiplication.** All 57,600 root
pairs (240 × 240) were tested; every product lies in L.

L = D₈⁺ at root norm 2 is a nonunital order (e₀ ∉ L) and, as a lattice, a
√2-scaled copy of a maximal order of integral octonions: L = (e₁ + e₂)𝒪 with
𝒪 := ½(e₁ + e₂)L, one of the seven maximal orders that contain the eight
basis units (paper Section 2.2).  The transposition σ fixes L (paper
Lemma 4.1) and acts on Wilson's sublattices: Ls̄ = 2J₁, twice the analogue
J₁ of Kirmse's set, σ(Ls̄) = 2𝒪 and σ(Ls) = p𝒪 with
p = e₂ − e₄ + e₅ − e₆ (paper (4.1); ledger claims 019 and 021).

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

The pair counts in this phase and in the norm distribution below were
produced by the harness trial_007_triple_octonion_swap.py, which is not in
the repository, so they cannot be rerun.  The trial-007 harness in the
repository, python_project/src/trial_007_kirmse_twist.py, tests 900 pairs
per transposition.  All 21 transpositions close on Λ exactly: the 576
basis-pair products of a Z-basis of Λ lie in Λ
(python_project/src/verify_consecutive_twists_exact.py).  For swap(1,2),
trial_007_scaled_test.py, trial_007_fast.py and trial_007_exhaust.py test
12,000,000 further random minimal-vector pairs, with 0 failures.

## Phase 3: Broader permutation survey

Not all permutations work. A survey of permutations of seven cycle types
gave the following; no program in the repository reproduces it, and the
exact counts over all of S₇ are those of the paper's Section 5.6
(`sec:census`; `python_project/src/verify_all_permutations_exact.py`):

| Cycle type | Worked (0 failures) / tested | Rate |
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

### What distinguishes working from nonworking permutations

Fano-plane automorphisms (index-doubling, identity) give the SAME
multiplication as the standard and inherit its type3×type3 failure.
Permutations that are not automorphisms give genuinely different multiplications.
Among these, SOME close on the Leech lattice and others do not.  All
transpositions close.  The criterion for general permutations is the
paper's Proposition 5.3 (`prop:criterion`): the twisted product
x ·_π y := π⁻¹(π(x) · π(y)) gives a triple product that closes on Λ exactly
when π(Ls̄) and π(Ls) are closed under the standard octonion product.  The
exact census over all of S₇ is the paper's Section 5.6 (`sec:census`).

## Failure mode contrast

| Algebra | t3×t3 failure rate | Wilson condition violated |
|---|---|---|
| Standard | ≈74% | Condition 3 only (x+y+z ∈ Ls) |
| swap(1,2) | 0% | (none) |

The transposition removes precisely the Wilson condition 3 failure of the
untwisted product (trials 001 and 003).

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
changes entries of the multiplication table in sign or target (swap(1,2)
changes 30 of the 64 entries; trial_007_explanation.py), so that Ls is
closed under the twisted product (paper Lemma 4.4 and Proposition 5.3), which
removes the condition-3 failure of the standard multiplication.

**This means (Λ, +, ·) is an order in the R-algebra (R²⁴, +, ·)**, where ·
is the transposition-twisted triple octonion product, the paper's
Z₃-symmetric triple-octonion product.  Closure on all of Λ is proved in the
paper (Theorem 1.1, `thm:main`; Section 4).

## Caveats and verification status

1. **Sample-based, not exhaustive.** The type3×type3 test covers 50,000 of
   the ~3.4 × 10¹⁰ possible pairs.  Zero failures in 593,400 total pairs
   is evidence, not proof; the proof is the paper's (Section 4).

2. **Bilinearity argument.** If the product of any two Min(Λ) vectors lies
   in Λ, then by bilinearity and the fact that Min(Λ) generates Λ over Z,
   ALL products of Λ vectors lie in Λ.  So testing on Min(Λ) suffices.

3. **Algebraic properties.** The algebra is not alternative, flexible or
   power-associative, and it has no multiplicative identity (paper Section
   5.1, `sec:identities`, and Remark 5.2, `rem:no-identity`).

4. **Independent verification** (Manifesto §5).  Closure is proved
   symbolically (paper Section 4) and checked in exact arithmetic by
   python_project/src/verify_structural_proof.py,
   python_project/src/symbolic_proof_checks.py and
   gap_project/tests/test_lemmas.g.

## Verdict

**PASS.**  The transposition-twisted triple octonion product achieves 100%
Leech lattice closure on 593,400 tested pairs with zero failures, and
closure on all of Λ is proved (paper Theorem 1.1): the Leech lattice admits
an order under an octonion-derived bilinear product.

## What this establishes

1. Wilson's L is closed under the standard octonion product: a nonunital
   order, as a lattice a √2-scaled copy of a maximal order (Phase 0).
2. Fano-plane automorphisms (including index-doubling) do not change the
   multiplication and reproduce the trial 001 failure.
3. Transpositions of imaginary basis elements produce a DIFFERENT
   multiplication that closes on Λ.
4. The twist removes the failure of Wilson condition 3 (x+y+z ∈ Ls) on
   type3×type3 products, the failure mode of the untwisted product (trials
   001 and 003).

## Questions raised by this trial

1. Closure on all 196,560² ≈ 3.9 × 10¹⁰ pairs: proved for all of Λ (paper
   Theorem 1.1, Section 4).
2. Algebraic properties of the resulting 24-dimensional algebra: paper
   Section 5 (`sec:properties`).  The automorphism group of the order is
   C₆ × S₃, of order 36, strictly inside Co₀ (paper Section 5.4,
   `sec:autgroup`).
3. The criterion distinguishing working from nonworking permutations:
   paper Proposition 5.3 (`prop:criterion`), with the census of paper
   Section 5.6 (`sec:census`).  A reason for the closure pattern beyond the
   2-cycles, in particular for the failure of every 4-cycle and for the
   exact one-half rates, is open (paper Section 7, `sec:conclusion`).
4. Different transpositions give isomorphic orders: the 21 transpositions
   form one orbit under conjugation by the order-21 group of coordinate
   permutations preserving the oriented Fano triples, each an octonion
   automorphism that fixes s and s̄ and preserves L (paper Section 4, first
   footnote).
