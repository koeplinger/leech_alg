# Research Result

## Summary

The Leech lattice Λ admits an order under a bilinear product derived from
the octonion algebra.

### The product

Let O denote the real octonion algebra with standard Fano-plane
multiplication, and let σ: O → O be the linear map that swaps two
imaginary basis elements (any transposition of {e₁, …, e₇}).  Define a
new multiplication ·_σ on O by

    x ·_σ y  =  σ( σ(x) · σ(y) )        (†)

where · is the standard octonion product.  The map σ is an isomorphism
from (O, ·) to (O, ·_σ): σ(x · y) = σ(x) ·_σ σ(y).  The two
multiplication tables differ in 30 of 64 structure constants.

Now set R²⁴ = O₁ ⊕ O₂ ⊕ O₃, three copies of (O, ·_σ) — all using the
**same** swapped multiplication — with Z₃ cross-block routing:

    Same block:      O_α × O_α → O_α
    Cross block:     O_α × O_β → O_γ      where {α,β,γ} = {1,2,3}

All nine block-pair products use the single multiplication ·_σ.

### The claim

For every u, v ∈ Min(Λ) (the 196,560 vectors of squared norm 8), the
product u · v lies in Λ.

Since Min(Λ) generates Λ over Z, bilinearity extends closure to all of Λ:
for every u, v ∈ Λ, the product u · v lies in Λ.  Therefore (Λ, +, ·) is
an order in the 24-dimensional R-algebra (R²⁴, +, ·).

### What makes this simple

The entire construction has three ingredients:

1. **One octonion algebra**, modified by a single transposition of two
   imaginary basis elements.
2. **Three identical copies** of that algebra on R²⁴.
3. **Z₃ cross-block routing**: products of vectors from different blocks
   land in the third block.

No rescaling, no conjugation, no distinct algebras on different blocks.
The same product on all three copies.

### Verification status

| Test | Pairs | Failures |
|---|---|---|
| Trial 007 initial (swap(1,2)) | 593,400 | 0 |
| Scaled test (4M random, all types) | 4,000,000 | 0 |
| Scaled test (4M, multiprocessor) | 4,000,000 | 0 |
| All 21 transpositions (small sample) | ~15,000 | 0 |

Exhaustive verification of all 196,560² ≈ 3.86 × 10¹⁰ pairs is
computationally feasible (~2 hours with 16 CPU cores) and is planned.

All 21 transpositions of imaginary basis elements produce the same result.
Since all transpositions lie in a single orbit under the Fano-plane
automorphism group GL(3, F₂), the construction is essentially unique up
to basis relabelling.

### Symbolic proof

Closure is proved symbolically via four lemmas (exact arithmetic,
`symbolic_proof_checks.py`):

- **Lemma A**: σ(L) = L (coordinate permutation preserves D₈⁺).
- **Lemma B**: L·L ⊆ L (L is a maximal order — Coxeter 1946).
- **Lemma C**: L · σ(Ls̄) ⊆ σ(Ls̄) (64 basis products verified exactly).
- **Lemma D**: σ(Ls) · σ(Ls) ⊆ σ(Ls) (64 basis products verified exactly).

The standard product fails condition 3 because Ls·Ls ⊄ Ls.  The twist
maps the condition-3 sublattice from Ls to σ(Ls), where closure holds.
The accompanying observation σ(Ls) ≠ Ls (verified by an explicit
witness in the script) is recorded as a remark in Section 4 of the
paper rather than a lemma: it confirms the construction is
non-trivial, but is not used in the closure argument.

### Algebraic properties

The order (Λ, +, ·) is non-unital, non-commutative, non-associative,
not alternative, not flexible, and not power-associative.  The product
is not norm-multiplicative.

---

Last updated: 2026-04-14
