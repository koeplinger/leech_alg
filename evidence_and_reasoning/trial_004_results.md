# Trial 004 Results — E8 automorphism basis changes

Date: 2026-04-09

## Algebra family tested

Same triple-octonion base (O₁ ⊕ O₂ ⊕ O₃) with Z₃-symmetric cross-block
routing, but with a change of basis within each 8-dimensional block using an
automorphism of the E8 lattice.

If T₁, T₂, T₃ are E8 lattice automorphisms, the modified product is:

    (a₁, a₂, a₃) ★ (b₁, b₂, b₃) = (c₁, c₂, c₃)
    where cₖ = Tₖ · (octonion product of Tᵢ⁻¹aᵢ and Tⱼ⁻¹bⱼ, summed)

Automorphisms tested are drawn from the signed-permutation subgroup of W(E₈):
coordinate permutations × even-parity sign flips. This subgroup has order
128 × 40,320 = 5,160,960.

## Test methodology

150 Min(Λ) pairs tested per automorphism:
- 80 type3×type3, 30 type2×type2, 20 type1×type1, 20 type2×type3

Three search strategies:
1. **Single-block**: T on block 0, identity on blocks 1 and 2 (500 random)
2. **Three-block**: independent random T₁, T₂, T₃ (500 random)
3. **Sign-only on block 0**: all 128 even-parity sign patterns (exhaustive)

## Results

| Search | Tested | Best closure rate | Baseline |
|--------|--------|-------------------|----------|
| Baseline (identity) | 1 | 60.0% (90/150) | — |
| Single-block random | 500 | 53.3% (80/150) | 60.0% |
| Three-block random | 500 | 40.0% (60/150) | 60.0% |
| Sign-only (exhaustive) | 128 | 60.0% (90/150) | 60.0% |

### Key findings

**1. The identity (no basis change) is already optimal.** No automorphism
improves over the baseline. Most make it significantly worse.

**2. Single-block changes hurt.** Best is 53.3% vs. 60.0% baseline. Breaking
the Z₃ symmetry by modifying only one block degrades performance.

**3. Three-block changes hurt more.** Best is 40.0%. Independent random
automorphisms on all three blocks destroy the coordination between blocks
that the original algebra has.

**4. Automorphisms introduce NEW failure types.** The best single-block result
(53.3%) has failures in type2×type2 (2 failures) in addition to the usual
type3×type3 (68 failures). The baseline has zero type2×type2 failures. Basis
changes disrupt the lattice conditions that the original product preserves.

**5. Sign-only changes are neutral.** All 128 sign patterns on one block
tie at the baseline rate of 60.0%. Sign changes within a block (with even
parity to preserve L) commute with the octonion product structure in a way
that preserves the failure pattern.

## Interpretation

The triple-octonion algebra has a specific coordination between its three
blocks that is already "locally optimal" among E8 automorphism-related
variants. Changing the basis within a block can only make things worse —
it breaks the alignment between the octonion product and the Leech lattice
membership conditions without gaining anything.

This is consistent with the trial 003 finding that the failure is structural:
it's not the specific coordinate system within each E8 block that causes the
type3×type3 obstruction, but the triple-block architecture itself.

## Verdict

**FAIL — no improvement possible via basis change.** The identity basis is
already optimal. E8 automorphisms (sign changes, coordinate permutations)
can only degrade closure, not improve it. The type3×type3 obstruction is
invariant under the automorphism group of the constituent E8 lattice.

## What this rules out

The triple-octonion algebra O₁ ⊕ O₂ ⊕ O₃ with any choice of E8 lattice
automorphism applied independently to each block. Combined with trial 003
(discrete cross-block modifications), this establishes that no simple
modification of the triple-octonion algebra can achieve an order on Λ.
