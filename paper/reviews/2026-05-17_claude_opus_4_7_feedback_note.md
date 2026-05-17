# Feedback note — is $\Lambda/2\Lambda \cong \mathbb{F}_2[\mathbb{O}]^3$?

**Date:** 2026-05-17
**Reviewer:** Claude Opus 4.7 (Anthropic), at the direction of Jens Köplinger
**Prompt:** 097. Records the assessment of one item from the computer-aided
review at `github.com/koeplinger/leech_alg` (issue 1).

## The feedback item

> *"Is $\Lambda/2\Lambda \cong \mathbb{F}_2[\mathbb{O}]^3$? If yes, this
> connects the order to the octonion $\mathbb{F}_2$-algebra, which has known
> connections to the binary Golay code and the sporadic groups."*

A natural question to raise, and checking it was worthwhile. The verdict,
however, is **not actionable** — no change to the manuscript follows.

## Assessment

**As an algebra isomorphism, the claim is false.** The product $\star$
descends to a well-defined $\mathbb{F}_2$-bilinear product $\bar\star$ on the
24-dimensional space $\Lambda/2\Lambda$. Computing the full structure
constants (Wilson's $\Lambda$ built explicitly, mod-2 reduced) shows
$(\Lambda/2\Lambda,\bar\star)$ is **non-unital**: it has neither a left nor a
right identity, and every product lies in an 8-dimensional subspace
($\dim(A\cdot A)=8$). The octonion algebra over $\mathbb{F}_2$ is unital, so
$\mathbb{F}_2[\mathbb{O}]^3$ is unital with full product image. Unitality
alone settles it — the two are not isomorphic. As *vector spaces* both are
$\mathbb{F}_2^{24}$, but that isomorphism is trivial and carries no algebraic
content.

**The Golay/sporadic backdrop is classical.** The link between the Leech
lattice mod 2, the binary Golay code, and the Conway group is textbook
(Conway; Conway–Sloane, *SPLAG*). The Turyn-style construction of $\Lambda$
from three copies of $E_8$ — parallel to the construction of the Golay code
from three Hamming codes — is Conway–Sloane (1982), already cited in the
manuscript's Related Work. The order is octonionic by construction, so a
mod-2 reduction introduces no new bridge. Nothing here is a new discovery.

## One observation worth keeping

The reduction did surface a genuine, if modest, fact: $\bar\star$ on
$\Lambda/2\Lambda$ is **alternative**, even though $\star$ itself is *not*
alternative over $\mathbb{R}$ (Section 5 of the manuscript). A non-alternative
bilinear product whose mod-2 reduction is alternative is a small curiosity,
and an interesting one against the backdrop of the structure theory of
alternative algebras developed by A. A. Albert and others. It is an
observation, not a result — its meaning, if any, is unexplored — and at most
it merits a one-line future-work remark. (A specific Albert reference can be
added to `evidence_and_reasoning/references/` if this thread is pursued.)

## Acknowledgment

Thanks to the computer-aided review for prompting the question; ruling it out
cleanly sharpened the picture of what the mod-2 reduction does and does not
give.
