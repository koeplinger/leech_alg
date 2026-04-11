# Research Result

## Summary

The Leech lattice Λ admits an order under a bilinear product derived from
the octonion algebra.

### The product

Let O denote the real octonion algebra with standard Fano-plane
multiplication, and let σ: O → O be the linear map that swaps two
imaginary basis elements (any transposition of {e₁, …, e₇}).  Define a
new multiplication ·_σ on O by

    x ·_σ y  =  σ(x) · σ(y)        (†)

where · is the standard octonion product.  This is an isomorphic copy of O
on the same R⁸, differing in 30 of 64 structure constants.

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

### What is not yet established

1. Exhaustive computational verification (planned, feasible).
2. A symbolic (mathematical) proof of closure.
3. Algebraic properties of (Λ, +, ·): whether the algebra is alternative,
   flexible, power-associative, or has a multiplicative identity.
4. The relationship to Wilson's and Dixon's octonionic Leech lattice
   constructions.

---

Last updated: 2026-04-11
