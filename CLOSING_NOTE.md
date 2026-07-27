# Closing note

**27 July 2026.** This effort is closed. What follows records the result, the
reasons for closing, and what this repository is from here on.

## The result

The Leech lattice Λ admits an order: it is closed under a Z₃-symmetric
triple-octonion product on R²⁴, the same octonion product in all three blocks
with Z₃-symmetric cross-block routing. In Robert A. Wilson's sublattice
description of Λ, closure turns on a transposition σ of two imaginary basis
units, which leaves the E₈ lattice invariant and is an octonion-algebra
isomorphism, but moves Wilson's sublattices Ls̄ and Ls.

The result and its verification are published:

> J. Köplinger, *An order on the Leech lattice from a Z₃-symmetric
> triple-octonion product*,
> [doi:10.13140/RG.2.2.22093.19686](https://doi.org/10.13140/RG.2.2.22093.19686)
> (ResearchGate).

The corresponding source is frozen in this repository as
[`paper/main_2026-07-19.tex`](paper/main_2026-07-19.tex) (v7, 24 pages), with
the revision chain in [`paper/`](paper/) and the review record in
[`paper/reviews/`](paper/reviews/).

## Why the effort is closed

1. **Journal submission is not the productive path right now.** Journals appear
   to be swamped with AI-generated submissions, and it is hard to differentiate
   a structured approach from that volume. The path forward would be to
   solidify and understand the approach itself, and to advertise that instead.

2. **The prompt log's value is front-loaded.** The discovery was made in
   [prompt 025](prompt_logs/025_triple_kirmse_twisted_octonions.txt)
   (11 April 2026), where the construction was proposed: three copies of the
   octonions, each twisted on the same imaginary direction, combined with Z₃
   cross-block routing. The next few dozen prompts are the interesting ones,
   tracing how the result was checked and ultimately owned, and which parts
   were not. Logging every incremental update across the last 150 or so
   prompts has yielded diminishing value.

3. **The mathematics continues by other means.** Analysis of the structure is
   ongoing, with thanks for the feedback and interest received in the work
   since. That work now proceeds by largely conventional methods, led by human
   understanding, with AI as a supporting tool.

## What this repository is now

A record, not a work in progress. It holds the published result and the full
trail behind it: the prompt logs, the trial results and key claims, the Python
and GAP verification suites, the reference registry, and the paper with its
frozen versions, changelogs, and reviews. The document genres in
[`DOCUMENT_GENRES.md`](DOCUMENT_GENRES.md) continue to describe how each
artifact may be changed, and [`tools/lint_docs.py`](tools/lint_docs.py)
continues to enforce them mechanically.

The mathematical questions the paper leaves open are stated in its Sections 7
and 8, and summarized in [`CURRENT_STATE.md`](CURRENT_STATE.md). They remain
open; they are simply no longer pursued here.

— Recorded by Claude (Anthropic), at the direction of Jens Köplinger,
2026-07-27.
