# Trial 004 Results — E8 automorphism basis changes

Created: 2026-04-09.  Updated: 2026-10-28.

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

**1. No sampled automorphism improves on the identity.** Most make it
significantly worse.  The family has a closing member that the sampling
missed: the same transposition σ of two imaginary units in all three blocks
(T₁ = T₂ = T₃ = σ) closes on Λ (trial 007; paper Theorem 1.1).

**2. Single-block changes hurt.** Best is 53.3% vs. 60.0% baseline. Breaking
the Z₃ symmetry by modifying only one block degrades performance.

**3. Three-block changes hurt more.** Best is 40.0%. Independent random
automorphisms on all three blocks destroy the coordination between blocks
that the original algebra has.

**4. Automorphisms introduce NEW failure types.** The best single-block result
(53.3%) has failures in type2×type2 (2 failures) in addition to the usual
type3×type3 (68 failures). The baseline has zero type2×type2 failures. Basis
changes disrupt the lattice conditions that the original product preserves.

**5. Sign-only changes do not help.** The best of the 128 sign patterns on
one block ties with the baseline at 60.0% (90/150). The program prints only
the best rate; on the same 150 pairs, 16 of the 128 patterns reach 60.0% and
the other 112 give 50.0% to 55.3% (75 to 83 of 150).

## Interpretation

Random basis changes in single blocks, or independent ones in the three
blocks, break the alignment between the octonion product and the Leech
lattice membership conditions.  The coordinate system within the blocks does
matter: the same transposition σ in all three blocks moves Wilson's
sublattices Ls̄ and Ls onto lattices closed under the octonion product, and
the triple product closes (paper Section 4).

## Verdict

**FAIL.** None of the 1,128 sampled basis changes closes, and none improves
on the identity.  The obstruction is not invariant under the automorphisms of
E8: the unsampled change T₁ = T₂ = T₃ = σ closes (trial 007).

## What this rules out

The basis changes tested: 500 single-block and 500 independent three-block
signed permutations, and all 128 even sign changes on one block.  It does not
rule out the family: the same transposition σ in all three blocks closes on
Λ (trial 007; paper Theorem 1.1).
