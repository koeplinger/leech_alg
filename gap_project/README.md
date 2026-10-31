# GAP / LOOPS verification suite

Created: 2026-04-28.  Updated: 2026-10-31.

A GAP re-implementation of the verification tests that went into the paper
[`paper/main.tex`](../paper/main.tex) and the companion note
[`paper/companion.tex`](../paper/companion.tex), running in parallel to the
Python suite in [`python_project/`](../python_project/).  The two standalone
scripts described at the end of this file also supply the group-theoretic half
of the automorphism-group note
[`paper/automorphism_group_2026-07-12.tex`](../paper/automorphism_group_2026-07-12.tex).

The goal of this directory is independent re-derivation: the same
mathematical claims, executed in a different language, on a different
runtime (GAP, developed by a different community of mathematicians).  All arithmetic is
exact (rational numbers in GAP); no floating point is used.

The LOOPS package is used in `tests/test_octonion.g` as an independent check
that the multiplication encoded by our Fano triples really does form a
Moufang loop of order 16, and that this loop is isomorphic to the
LOOPS-library `MoufangLoop(16, 3)` (the standard octonion loop).

## Prerequisites

- **GAP 4.x.**  On Debian/Ubuntu:  `sudo apt install gap`.
- **The LOOPS package.**  Install it in GAP via
  `LoadPackage("PackageManager"); InstallPackage("loops");`, or unpack the
  release from https://gap-packages.github.io/loops/ into `~/.gap/pkg/`.
- **The Small Groups library (SmallGrp)** (usually bundled with GAP) for
  `StructureDescription`/`IdGroup` in `aut_lambda_star.g` and
  `octonion_stabilisers.g`.  Without it, `aut_lambda_star.g` prints
  `Size(G) = 36` and then aborts in a break loop at the structure
  identification, and `octonion_stabilisers.g` aborts at its first group,
  before printing any order.
- **The AtlasRep package** (usually autoloaded with GAP) for one line of
  `octonion_stabilisers.g`, which confirms Stab(Ls̄) ≅ U₃(3).2 = G₂(2) in
  Aut(O, ·_σ) against
  `AtlasGroup("U3(3).2")`.  Only that single cross-check needs it; the orders
  and structure descriptions do not.

To check that LOOPS is available:

```bash
echo 'LoadPackage("loops");; Print(IsPackageLoaded("loops"), "\n"); QUIT;' | gap -q
```

This suite was developed against **GAP 4.12.1** and **LOOPS 3.4.4**.
Full citation entries for both, in the project's reference registry, are
in
[`evidence_and_reasoning/references/computer_algebra_systems.txt`](../evidence_and_reasoning/references/computer_algebra_systems.txt).
For convenience, the recommended short citations are:

> [GAP] *GAP — Groups, Algorithms, and Programming, Version 4.12.1*,
> The GAP Group (2022), https://www.gap-system.org/.
>
> [LOOPS] G.P. Nagy and P. Vojtěchovský, *LOOPS, Computing with
> quasigroups and loops in GAP, Version 3.4.4* (2024), GAP package,
> https://gap-packages.github.io/loops/.

## Running

From the repository root:

```bash
gap -q gap_project/run_all.g
```

That runs all eight test files and prints a final pass/fail total.  Exit
code is zero iff all checks pass.

To run a single test file (useful while debugging):

```bash
gap -q -c 'GAP_PROJECT_ROOT := "/abs/path/to/leech_alg/gap_project";' \
        gap_project/tests/test_lemmas.g
```

## Layout

```
gap_project/
├── README.md              this file
├── run_all.g              driver: runs every tests/*.g
├── aut_lambda_star.g          standalone: Size / StructureDescription of
│                              Aut(Lambda, +, star)  (order 36, C6 x S3)
├── aut_lambda_star_gens.g     its input: 3 generators, 24x24 integer matrices
│                              in the Z-basis of Lambda, written by
│                              verify_aut_lambda_star.py
├── octonion_stabilisers.g     standalone: identifies Stab(L), Stab(Ls),
│                              Stab(Lsbar) inside the compact G_2 = Aut(O, ._sigma)
├── octonion_stabilisers_gens.g  its input: 8x8 rational generators, written
│                              by verify_aut_octonion_crosscheck.py
├── src/                   shared modules (loaded by every test)
│   ├── harness.g          tiny CHECK/PrintTestSummary harness, LoadAllSrc()
│   ├── octonion.g         Fano triples, OctMult, OctConjugate, OctNormSq,
│   │                      LOOPS Cayley-table builder
│   ├── e8_wilson.g        L = D_8^+, IsInL, root enumerators, L's Z-basis
│   ├── sublattices.g      Ls, Lsbar bases; IsIntegerCombination
│   ├── twist.g            sigma = (e_1 e_2), OctMultSigma, sigma-images of
│   │                      sublattice bases
│   ├── leech_wilson.g     Wilson's three conditions, J = +-{e_0,...,e_7},
│   │                      type-1 minimal-vector enumerator
│   └── triple_product.g   Z_3-routed triple product on R^24, untwisted and
│                          sigma-twisted forms
└── tests/                 test files (one CHECK per claim)
    ├── test_octonion.g            paper Section 2.1, plus LOOPS check
    ├── test_e8_wilson.g           paper Section 2.2
    ├── test_sublattices.g         paper Section 2.3
    ├── test_twist.g               paper Section 3 (Definitions 3.1, 3.2:
    │                              Transposition and σ-twisted product)
    ├── test_lemmas.g              paper Section 4: the four lemmas of the
    │                              proof of Theorem 1.1
    ├── test_leech_wilson.g        paper Definition 2.1 / Wilson Section 3
    ├── test_triple_product.g     paper Theorem 1.1 spot-check + Section 5
    │                              (no multiplicative identity)
    └── test_companion_examples.g  the companion's explicit numerics
```

## What each test file verifies

### `test_octonion.g`: 20 checks
Identity, imaginary squares, anticommutativity, noncommutativity,
nonassociativity, composition law on all 64 basis pairs, conjugate / norm
identities, left/right alternativity, flexibility, total antisymmetry of
the associator on all 512 basis triples (this certifies alternativity for
all of $\mathbb{R}^8$ by trilinearity), Moufang identity on 512 basis triples,
division-algebra inverses, power-associativity, Dixon's standard rule
$e_a e_{a+1} = e_{a+3}$.  LOOPS-package checks: the 16 signed basis units
form an order-16 nonassociative Moufang loop, isomorphic to the LOOPS
library's `MoufangLoop(16, 3)`.

### `test_e8_wilson.g`: 22 checks
Root counts (112, 128, 240), distinctness, coordinate structure, all roots
of norm 2, lattice membership predicate, single $e_i \notin L$,
$2e_i \in L$, integer Gram matrix, rank 8, unimodular Z-basis with $|\det|=1$,
$s\in L$ with norm 2, $\bar s \notin L$ with norm 2 (paper Section 2.3),
even-lattice property, $L\cdot L \subseteq L$ on all
$240^2 = 57{,}600$ root pairs (Coxeter / Lemma 4.2), and $r\cdot s, r\cdot \bar s\in L$
for all 240 roots ($Ls, L\bar s\subseteq L$).

### `test_sublattices.g`: 16 checks
Ranks of $Ls$, $L\bar s$ bases; their Gram determinants $= 256 = 16^2$;
$\mathrm{IsInLs}$, $\mathrm{IsInLsBar}$ accept all 240 root products; both
reject every E8 root; $2L \subseteq Ls \cap L\bar s$; $Ls + L\bar s$ of rank 8;
Lemma 4.1 ($\sigma(L)=L$); Remark 4.5 ($\sigma(Ls) \ne Ls$, $Ls\cdot Ls
\not\subseteq Ls$); the explicit companion witness vector
$(-1,0,1,0,0,1,0,1)\in\sigma(Ls)\setminus Ls$.

### `test_twist.g`: 8 checks
$\sigma$ is an involution and an isometry; the algebra-isomorphism
identity $\sigma(xy) = \sigma(x)\cdot_\sigma\sigma(y)$ from
Definition 3.2 of the paper; composition law for $\cdot_\sigma$; $\cdot_\sigma$ disagrees with $\cdot$ on at least one basis
pair; the twist resolves the $Ls$ obstruction
($Ls\cdot Ls\not\subseteq Ls$ but $Ls\cdot_\sigma Ls\subseteq Ls$);
companion Section 10 ($L\cdot_\sigma L\subseteq L$, so the twist does not
undo Coxeter).

### `test_lemmas.g`: 7 checks
The four lemmas of the symbolic proof (paper Section 4):
$\sigma(L)=L$; $L\cdot L\subseteq L$;
$L\cdot\sigma(L\bar s)\subseteq\sigma(L\bar s)$;
$\sigma(Ls)\cdot\sigma(Ls)\subseteq\sigma(Ls)$.
Lemma 4.1 is checked on the eight basis vectors of $L$, Lemmas 4.2 to 4.4
on the 8×8 = 64 products of basis vectors of the two factors, for each
lemma, all exactly.  Plus Remark 4.5 (two nontriviality witnesses) and
the bonus untwisted condition 2 ($L\cdot L\bar s\subseteq L\bar s$).
The paper cites this file in Appendix C (paragraph "Symbolic proof and
numeric validation").

### `test_leech_wilson.g`: 14 checks
$|J|=16$, all unit norm; $L\bar s, Ls\subseteq L$ on roots;
$\mathrm{IsInLsBar}, \mathrm{IsInLs}$ reject all roots; the 720 type-1
minimal vectors of $\Lambda$ (full enumeration) all have ambient norm 8 and
satisfy Wilson's three conditions; samples of type-2 and type-3 minimal
vectors built from generators $\lambda, j, k$ all satisfy the three
conditions; nonmembers `(root, 0, 0)` and `(root, root, 0)` are correctly
rejected.

### `test_triple_product.g`: 6 checks
Theorem 1.1 spot-check: $u\star_\sigma v\in\Lambda$ for every type-1 $u$
paired with six right operands from the three families (three of type 1,
two of type 2, one of type 3; 720 × 6 = 4,320 products), plus
type-2 × {type-2, type-3} and type-3 × type-3 sample grids.  The standard
(untwisted) triple product fails on at least one pair of the type-3 sample
(companion Example 8.1).
Two checks illustrate the paper's Remark 5.2 (Section 5.1), that
$(\mathbb{R}^{24}, \star_\sigma)$ has no multiplicative identity:
$(e_0,e_0,e_0)\star_\sigma v$ has block 1 equal to $x'+y'+z'$, not $x'$, and
no $(a,b,c)$ with components in $J$ is a left identity.  The remark proves
the general statement by a block-matrix argument; on the Python side,
`python_project/src/verify_identity_witnesses.py` (item (G)) cross-checks it
on the linear equations for a left and a right identity, and
`python_project/src/verify_idempotent_classification.py` gives a second
route (an identity would be an idempotent; there are exactly eight, paper
Section 5.2, and none of them acts as one).

### `test_companion_examples.g`: 17 checks
The explicit hand-checkable numerics in the companion:
$e_1\cdot e_2 = e_4$, $e_1\cdot e_3 = e_7$, $e_2\cdot e_7 = -e_6$;
$s\cdot s = -\tfrac32 e_0 - \tfrac12(e_1+\cdots+e_7) \in L$;
the companion's Example 4.2 pair
$a = (-1,0,0,0,1,0,1,1)$, $b = (-1,1,0,0,0,1,0,1)$ in $Ls$ with
$a\cdot b = (0,-2,2,1,-2,-1,-1,-1) \in L\setminus Ls$;
$e_1\cdot_\sigma e_2 = -e_4$ (companion Section 5);
the companion's "Summary of concrete data" three-point list.

## Cross-references with the Python suite

| GAP file                          | Python equivalent (mostly)                                  |
|-----------------------------------|-------------------------------------------------------------|
| `src/octonion.g`                  | `python_project/src/octonions.py`                           |
| `src/e8_wilson.g`                 | `python_project/src/e8_wilson.py`                           |
| `src/sublattices.g`               | parts of `python_project/src/leech_wilson.py`               |
| `src/twist.g`                     | parts of `python_project/src/trial_007_kirmse_twist.py`     |
| `src/leech_wilson.g`              | `python_project/src/leech_wilson.py`                        |
| `src/triple_product.g`            | parts of `python_project/src/trial_001_triple_octonion.py`  |
| `tests/test_octonion.g`           | `python_project/tests/test_octonions.py`                    |
| `tests/test_e8_wilson.g`          | `python_project/tests/test_e8_wilson.py`                    |
| `tests/test_leech_wilson.g`       | `python_project/tests/test_leech_wilson.py`                 |
| `tests/test_sublattices.g`        | (companion Example 4.2 + `leech_wilson.py` invariants)      |
| `tests/test_lemmas.g`             | `python_project/src/symbolic_proof_checks.py` (Lemmas 4.1, 4.3, 4.4) and `verify_structural_proof.py` (Lemmas 4.2 to 4.4, item (E)) |
| `tests/test_twist.g`              | (paper Section 3)                                           |
| `tests/test_triple_product.g`     | parts of `python_project/src/trial_007_kirmse_twist.py`     |
| `tests/test_companion_examples.g` | (the companion note's examples, recomputed here)                    |
| `aut_lambda_star.g`               | `python_project/src/verify_aut_lambda_star.py`              |
| `octonion_stabilisers.g`          | `python_project/src/verify_aut_octonion_crosscheck.py`      |

## Standalone group-identification scripts

Two scripts sit outside the `run_all.g` suite because they consume matrices
produced by the Python side (they are the group-theoretic half of the
automorphism-group calculation, and GAP's Schreier–Sims machinery is the
independent second opinion on the orders):

| Script | Input | What it reports |
|---|---|---|
| `aut_lambda_star.g` | `aut_lambda_star_gens.g` (written by `python_project/src/verify_aut_lambda_star.py`) | `Size` and `StructureDescription` of Aut(Λ, +, ⋆) from three 24×24 integer generators: **order 36, C₆ × S₃**.  Also the center (C₆), the derived subgroup (C₃), and the fact that −I₂₄ is **not** in the group (so the containment in Co₀ is strict).  The paper cites this script in Section 5.4 |
| `octonion_stabilisers.g` | `octonion_stabilisers_gens.g` (written by `python_project/src/verify_aut_octonion_crosscheck.py`) | identification of the three lattice stabilizers in Aut(O, ·_σ), the automorphism group of the σ-twisted product: Stab(L) = 2³·L₃(2) of order 1344, Stab(Ls) = 2³:(7:3) of order 168, Stab(Ls̄) = U₃(3).2 = G₂(2) of order 12096 (not computed here: under the standard product the orders are 1344, 21 and 1344; ledger claim 010a) |

Run them from inside `gap_project/`, here starting from the repository root;
reading standard input from `/dev/null` ends GAP when the script is done:

```bash
cd gap_project
gap -q -b aut_lambda_star.g < /dev/null
gap -q -b octonion_stabilisers.g < /dev/null
```

Reading the output of `aut_lambda_star.g`: it prints **two** lists of orders,
and they are not the same list.

- `orders of all 36 elements` is the element-order distribution, printed
  `Collected`: 1 (×1), 2 (×7), 3 (×8), 6 (×20).  This is the one to quote for
  the group.
- `orders of conjugacy class representatives` is a list of 18 numbers, one per
  conjugacy class.  It is *not* the multiset of element orders, and must not be
  quoted as such.

## Notes on conventions

- Coordinates are 1-indexed (GAP convention).  Position 1 holds the $e_0$
  coefficient, position 2 holds $e_1$, …, position 8 holds $e_7$.
- The transposition $\sigma$ is $(e_1\;e_2)$ throughout (the canonical
  choice: the 21 transpositions of imaginary basis units form a single orbit
  under conjugation by the order-21 group of coordinate permutations
  preserving the oriented Fano triples, so the proof for $(e_1\;e_2)$ carries
  over to each of them, as the footnote to $\sigma = (1\;2)$ in paper
  Section 4 states; all 21 close on $\Lambda$, paper Appendix E).
- What $\sigma$ is, and is not.  It is a transposition of two imaginary basis
  units, extended linearly and fixing $e_0$.  It is **not** an octonion
  automorphism.  It carries the octonion product to the isomorphic twisted
  product $x \cdot_\sigma y := \sigma(\sigma(x) \cdot \sigma(y))$, so
  $\cdot_\sigma$ is again an octonion product and not a new algebra; closure on
  Wilson's $\Lambda$ turns on the twist.  $\sigma$ fixes $L$ setwise
  ($\sigma(L) = L$, `test_sublattices.g`) but **moves** the sublattices $Ls$
  and $L\bar s$, which is the whole point ($\sigma(Ls) \ne Ls$ in
  `test_sublattices.g` and `test_lemmas.g`; $\sigma(L\bar s) \ne L\bar s$
  on the Python side, in `python_project/src/symbolic_proof_checks.py`).
- All lattice-membership tests are exact: $L$ by its coordinate criterion
  or over its $\mathbb{Z}$-basis, and $Ls$, $L\bar s$ and their
  $\sigma$-images over their $\mathbb{Z}$-bases, in exact rational
  arithmetic (no floating-point arithmetic anywhere).
