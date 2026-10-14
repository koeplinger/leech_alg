"""f8_ideal_classes.py -- ideal type of V in every closing polarization.

Added after the document check of 13 October 2026, to back two statements
of the review with the pipeline itself:
  (1) among closing twists, the analogue of the paper's Lemma 4.3
      (V_g = g(Lsbar)/2L a left ideal of L/2L) holds exactly in the
      index-2^16 class;
  (2) on most index-1 pairs V is neither a left nor a right ideal, so the
      paper's four-lemma template (or its mirror image, which needs a right
      ideal) does not apply there and the closure criterion is needed.

It reuses the definitions of f5_polarizations.py (imported as a module; the
heavy pair computation there runs only under __main__).  Two facts that
f5's main run verifies on all 17,280 ordered pairs are used here:
closure <=> both members are subalgebras, and index = h(W) * h(V)^2.
Pairs are ordered (i, j) = (V, W), as in f5.  Ideals are with respect to
the standard product on L/2L (mod-2 structure constants of L.L in L).
"""
from collections import Counter
import f5_polarizations as F

out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


closing = [(i, j) for (i, j) in F.pairs if F.sub[i] and F.sub[j]]
say("closing ordered pairs (both members subalgebras):", len(closing))


def kind(i):
    l, r = F.left[i], F.right[i]
    return "two-sided" if (l and r) else "left only" if l else "right only" if r else "neither"


table = Counter()
for i, j in closing:
    idx = F.H[j] * F.H[i] ** 2
    table[(idx, kind(i))] += 1

say("index of the pair, ideal type of V (first member), number of pairs:")
for (idx, k), n in sorted(table.items()):
    say("   index %-6d  V %-10s  %5d" % (idx, k, n))

left_iff_top = all((F.left[i]) == (F.H[j] * F.H[i] ** 2 == 2 ** 16) for i, j in closing)
say("V is a left ideal  <=>  index = 2^16, on every closing pair:", left_iff_top)
idx1 = [(i, j) for i, j in closing if F.H[j] * F.H[i] ** 2 == 1]
say("index-1 pairs: %d; V neither left nor right ideal on %d; V right (not left) ideal on %d"
    % (len(idx1), sum(1 for i, j in idx1 if kind(i) == "neither"),
       sum(1 for i, j in idx1 if kind(i) == "right only")))

open("f8_ideal_classes.out.txt", "w").write("\n".join(out) + "\n")
