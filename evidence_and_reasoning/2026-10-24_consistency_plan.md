# 2026-10-24 — Plan: consistency before the author's fine-comb read

Created: 2026-10-24.  Updated: 2026-10-24.

**Status: proposed; decisions taken.**  (s) Option 2: the keyword
"computer-assisted proof" is dropped and MSC 68V05 stays (P1, P3).  (t)
Option 1: the paragraph after the census table loses the three table facts it
restates (P2, P4).  Every edit of §5 was drafted read-only, checked against
its source, and applied to scratch copies that build or run as §4 states; the
plan was reviewed before hand-over by five readers, each finding checked by a
refuting verifier.  Drafted by Claude Opus 5.5 at the direction
of Jens Köplinger.

## 0. Goal and scope

The paper (`paper/main.tex`, the v8 draft) and the supporting documents with
a known error are brought into agreement with what the record now knows, so
that the author's fine-comb read finds a stable text.  The restructure of the
paper, starting with Section 5.6, follows that read, with its own plan.

- **In the paper:** only what the two decisions of §1 choose.  Nothing else in
  `paper/main.tex` is known to disagree with the record.
- **Outside the paper:** the automorphism note, which the paper links in
  Section 5.4 (§2); the Baez–Egan note, the companion note, the
  historical-appendix draft, the minimum-norm argument in `leech_wilson.py`
  with ledger claim 004 and a test docstring, eight GAP files and two index
  lines (§3).  None of these is linked from the paper except
  `gap_project/tests/test_lemmas.g`, which Appendix C cites; its edit Gm1
  changes header comments only.
- **Each edited document** gets its `Created`/`Updated` line and loses its
  history wording in the same edit (`DOCUMENT_GENRES.md`).

Out of scope: new results; the paper's prose beyond §1; documents without a
known error, which get their dates and cleanup at their next edit; the style
leftovers that the drafting noted in the edited notes (spaced dashes,
hyphenated non- words, British spellings), which are not errors.

## 1. Decisions for the author

### (s) The keyword "computer-assisted proof" and the MSC code 68V05

- The keywords end with "computer-assisted proof"; the secondary MSC codes
  are 17D05 and 68V05.  MSC 2020 names 68V05 "Computer assisted proofs of
  proofs-by-exhaustion type".
- The proof of Theorem 1.1 (Section 4) is symbolic.
- Results that rest on complete computation remain: the automorphism group
  (Section 5.4: "The derivation, the complete enumerations, and an independent
  verification are recorded separately in a companion note"), the census over
  $S_7$ and the span indices of Section 5.6, the classifications "by symbolic
  classification" of Sections 5.2 and 5.5, the span index of $\varphi$
  (Section 6) and the self-product norms of all minimal vectors (Appendix D).
  The footnotes cite their scripts, Section 5.4 refers to the companion note
  for the independent verification, and Appendix C states the verification
  protocol.

Options:
1. Keep both.
2. Drop the keyword and keep 68V05.  The keywords then name the subject of the
   main theorem, and the code still classifies the complete enumerations of
   Sections 5.4 and 5.6.
3. Drop both.

*Decided:* option 2, as recommended (edits P1 and P3).

### (t) The paragraph after the census table (Section 5.6)

The paragraph restates three facts that the table shows: the 4-cycles are
"the only one failing without exception", "the $(5,2)$ type has the lowest
nonzero rate", and three types "close at exactly half".  The rest of it says
what the table does not: the 2-cycles close by Theorem 1.1, the identity fails
by Proposition 5.3, what the column "Close" counts, the description of the
closing 3-cycles, and what the paper does not give.

Options:
1. Drop the three restated facts and put the sentence on the column "Close"
   before the 3-cycle description (edits P2 and P4).  The paragraph keeps all
   that the table does not show.  Section 7's third question still names the
   failure of every 4-cycle and the one-half rates, which the table shows.
2. Leave the paragraph to the restructure of Section 5.6.

*Decided:* option 1, as recommended: it is small, it removes the one known
redundancy in the section before the read, and the restructure can start
from it.

## 2. The automorphism note (linked from Section 5.4)

`paper/automorphism_group_2026-07-12.tex` still describes the paper as it was
at v6.  Edits A1 to A12:
- date lines (A1, A2);
- the paper's status: "(v6)" and "the author of the main paper will decide
  how, and how much of, this material is worked into it" become "whose
  Section 5.4 states the main results of this note and cites it" (A3); "The
  main paper leaves open the identification" becomes "This note determines"
  (A4); the subsection that grades the v6 Conclusion becomes "Relation to the
  main paper", which maps the statements of Section 5.4 to the note's results
  (A7);
- two to-dos about a docstring that already gives the noncircular argument
  (A5, A8);
- a pointer that has been wrong since the note was written: Route 1 is the
  note's own derivation, and `verify_star_algebra_structure.py` checks its
  input, the Fourier normal form (A6);
- the runtime of `verify_aut_lambda_star.py`, measured at about 5 minutes on
  one core (300 s and 274 s): the note's "About 5 minutes" stands and names
  the machine (A9), and the script's docstring, which said "well over ten
  minutes", is corrected (A12);
- the empty GAP path (A10) and the wording on the earlier probe (A11).

The PDF is rebuilt; the scratch build has 11 pages, no undefined references
and no overfull boxes.

## 3. Other supporting documents with known errors

Of these, only `gap_project/tests/test_lemmas.g` is linked from the paper
(Appendix C); its edit Gm1 changes header comments only.

**3.1 The Baez–Egan note** (`paper/Baez-Egan_Leech_Verification_addendum_2026-05-24.tex`;
file name kept, PDF rebuilt).  The known errors: "24-dimensional $L/2L$"
(it is 8-dimensional, B17); the polarization reading of $L/2L$ called new in
"our paper", while the paper's Appendix A credits it to Lepowsky and Meurman
(§2) and to Silverman, and claims as new only the product $\bar\mu$ (B12,
B14); "reproduced here" for Egan's proof, which the note gives in one line
(B13).  Found while drafting: $L\bar s$ and $Ls$ are images of $L$ under
right multiplication, not left (B6); the key lemma is used in showing that
$L_L/\sqrt2$, not $L_L$, is self-dual and without vectors of squared norm 2
(B7); $\sigma$
fixes $2L$ because $\sigma(L) = L$, not because it permutes coordinates of an
even lattice (B10); Lemma 4.3 gives a left ideal, and two-sidedness is
Remark 4.5 (B11).  The rest brings the note to the current state: date lines,
its description as a note rather than a dated addendum (B2 to B4), a pointer
to the Verification and Commentary (B5), pointers to v8 (B9, B15, B16), the
references (B18), US spelling (B8), and its line in
`paper/README.md` (B19).

**3.2 The companion note** (`paper/companion.tex`, PDF rebuilt).  The known
error: "Theorem 3.1 of the main paper" is Theorem 1.1 (C9).  Found while
drafting: the note's squared norm is what the paper calls the norm, and its
$\star$ is the untwisted product while the paper's is the twisted one; a
notation sentence says both (C4); its $\Lambda$ is the
paper's $\sqrt2$-scaled lattice (C8); the paper proves $L \cdot L \subseteq L$
as Lemma 4.2 (C5), Remark 4.5 holds the witness pair (C6), and the script
names the check Lemma 4.4, not "Lemma D" (C7); the paper's early title (C3);
date lines (C1, C2).

**3.3 The historical-appendix draft** (`paper/historical_appendix.tex`, PDF
rebuilt).  The known error: Mahler's ideal generator is an integer multiple of
an integral Cayley number of norm 1, 2 or 4 (H5).  Found while drafting:
Mahler's covering-radius theorem was checked by exact-arithmetic sampling,
not re-derived (H4); a comment about v4 (H3); date lines (H1, H2).

**3.4 The minimum-norm argument and Wilson's note on Dixon.**  The known
error: the sketch in `python_project/src/leech_wilson.py` that $\Lambda$ has
no vectors of norm below 8 jumps from $N_{\mathrm{std}} < 8$ to a zero
coordinate, which fails for three nonzero coordinates; Wilson's argument
(Section 4, p. 2189) replaces it, which first shows $N_{\mathrm{std}} \in 4\mathbb Z$
(W1).  Found while drafting: ledger claim 004 repeats the faulty step (K3) and
calls Wilson's bracketed note on Dixon a footnote about $s$ versus $\bar s$,
while the note (Section 5, p. 2189) concerns the multiplication by $j$ in the
first coordinate (K2; the module likewise, W2 and W3); the test docstring
writes the type-3 vectors with $\bar s$ where the module and Wilson use $s$
(T1 to T3); claim 004 and the test docstrings read a rank check as proof that
$L\bar s + Ls = L$, which is Wilson's statement (Section 2), and say so now
(K4, T4 to T6).  Claim 004 gets its date line (K1) and keeps its status,
ESTABLISHED.

**3.5 The GAP files** (comments and CHECK labels only; no code logic).  The
known stale pointers: eq. (1) is (2.1) (Go2, Gt2); Definition 2.4 is 2.1
(Gl1); Proposition 3.2 and Remark 3.5 are not in the paper (Gw1, Gw3, Gw4);
companion Example 4.5 is 4.2, Example 4.6 does not exist, Sections 7 and 8
are 5 or 10, "Sections 2-8" are 3 to 5 and 11, and the "final box" is
Section 11 (Gw2, Gw5, Gw6, Gs2 to Gs5, Gc1 to Gc10); `test_lemmas.g` names
where the paper cites it (Gm1); the usage example of `harness.g` calls
`TEST(...)`, where the harness defines `CHECK(...)` (Gh1).  Found while
drafting: `test_twist.g` called $\cdot_\sigma$ "not merely an
algebra-isomorphic relabelling" of the standard product, while the paper says
$\sigma$ is an algebra isomorphism (Gw1); `test_sublattices.g` labels a rank
check "$Ls + L\bar s = L$" (G-S6a to G-S6c; the README row, G-R).  Date lines
(Go1, Gh1, Gt1, Gw1, Gl1, Gs1, Gc1, Gm1).

## 4. Execution, checks and records

On the author's word, in this order:

1. §5.2 to §5.7.  The `Updated` dates in §5 read 2026-10-24; at execution
   each, and the existing `Updated` line of every other file edited in §5, is
   set to the day of the edit, checked with `date`.
2. §5.1: P1 and P3 for (s); P2 and P4 for (t).
   The paper's subtitle date is set to the day of the edit, together with its
   mentions in `paper/v7_to_v8_changelog.md` (the v8 line of its header and
   the subtitle under "Title block and preamble") and in `CURRENT_STATE.md`
   (the paper's date); `main.pdf` is built (22 pages expected), and the
   changed lines are read in the PDF.
3. The PDFs of the automorphism note, the Baez–Egan note, the companion note
   and the historical-appendix draft are rebuilt.
4. Checks: the GAP suite (110 passed; the scratch run with §5.7 applied gave
   110); `python_project/tests/test_leech_wilson.py` (31 passed in scratch
   with W1 to W3 and T1 to T6); every edited Python file parses; the linter
   is clean.
5. Records: this plan's status; its lines in `evidence_and_reasoning/README.md`
   and `CURRENT_STATE.md`; the prompt log.  The ledger's statuses do not
   change.
6. Before hand-over: an adversarial review of the applied edits, each finding
   checked by a refuting verifier; the confirmed findings applied.

Then the author's fine-comb read of the paper, and after it the plan for the
restructure, starting with Section 5.6.

## 5. Exact edits

Each edit replaces the old text, which occurs exactly once in its file, by the
new text; the edits of a file apply in the order given.

### 5.1 The paper

`paper/main.tex`

**P1.**  Decision (s), option 2: drop the keyword "computer-assisted proof"; MSC 68V05 stays.

Old:
```
\noindent\textbf{Keywords:} Leech lattice; integral octonions; $E_8$
lattice; nonassociative algebra; nonunital order; computer-assisted
proof.
```

New:
```
\noindent\textbf{Keywords:} Leech lattice; integral octonions; $E_8$
lattice; nonassociative algebra; nonunital order.
```

**P2.**  Decision (t), option 1: the paragraph after the census table without the three restated table facts; the sentence on the column "Close" moves before the 3-cycle description.

Old:
```
The $2$-cycles, the transpositions of Definition~\ref{def:twist},
close by Theorem~\ref{thm:main}.  By Proposition~\ref{prop:criterion},
the identity permutation does not close, since $Ls$ is not closed
under the standard octonion product (Remark~\ref{rem:Ls-not-closed}).
Among the cycle types in the table,
the $4$-cycles are the only one failing without exception; the $(5,2)$ type has the lowest
nonzero rate.  The $3$-cycles, $5$-cycles, and $(3,2,2)$ types
close at exactly half.  For the $3$-cycles the closing half admits a
description: for each three-element subset of imaginary indices, one
of the two cyclic orientations closes.  Each permutation~$\pi$
preserves~$L$, as $\sigma$ does (Lemma~\ref{lem:sigma-L}), so by the
same proposition the column ``Close'' counts the permutations~$\pi$
for which $\pi(L\bar s)$ and~$\pi(Ls)$ are closed under the standard
octonion product.  No description in terms of~$\pi$ like that of the
$3$-cycles is given here for the other types that close in part, nor a
reason for the failure of every $4$-cycle.
```

New:
```
The $2$-cycles, the transpositions of Definition~\ref{def:twist},
close by Theorem~\ref{thm:main}.  By Proposition~\ref{prop:criterion},
the identity permutation does not close, since $Ls$ is not closed
under the standard octonion product (Remark~\ref{rem:Ls-not-closed}).
Each permutation~$\pi$ preserves~$L$, as $\sigma$ does
(Lemma~\ref{lem:sigma-L}), so by the same proposition the column
``Close'' counts the permutations~$\pi$ for which $\pi(L\bar s)$
and~$\pi(Ls)$ are closed under the standard octonion product.  For the
$3$-cycles the closing half admits a description: for each
three-element subset of imaginary indices, one of the two cyclic
orientations closes.  No description in terms of~$\pi$ like that of
the $3$-cycles is given here for the other types that close in part,
nor a reason for the failure of every $4$-cycle.
```

`paper/v7_to_v8_changelog.md`

**P3.**  With P1 (decision (s), option 2): the changelog entry "Keywords and MSC".

Old:
```
Unchanged.  Keywords: "Leech lattice; integral octonions; $E_8$ lattice; nonassociative algebra; nonunital order; computer-assisted proof."  MSC 2020: Primary 11H06, 17A75; Secondary 17D05, 68V05.
```

New:
```
- **Keywords.**  "computer-assisted proof" is removed: the proof of Theorem 1.1 is symbolic (Section 4).  The keywords are "Leech lattice; integral octonions; $E_8$ lattice; nonassociative algebra; nonunital order."
- **MSC 2020.**  Unchanged: Primary 11H06, 17A75; Secondary 17D05, 68V05 ("Computer assisted proofs of proofs-by-exhaustion type"), for the complete enumerations of Sections 5.4 and 5.6.
```

**P4.**  With P2 (decision (t), option 1): the changelog bullet on the paragraph after the census table.

Old:
```
  - After the table: the 2-cycles close by Theorem 1.1; by Proposition 5.3 the identity permutation does not close, since $Ls$ is not closed (Remark 4.5), and the column "Close" counts the $\pi$ for which $\pi(L\bar s)$ and $\pi(Ls)$ are closed under the standard octonion product.  "The 4-cycles are the only type failing without exception" → "Among the cycle types in the table, the 4-cycles are the only one failing without exception"; "No such description is known for the other cycle types, and none is known for the $0\%$ of the 4-cycles." → "No description in terms of $\pi$ like that of the 3-cycles is given here for the other types that close in part, nor a reason for the failure of every 4-cycle."  The rest of the paragraph keeps v7's wording.
```

New:
```
  - After the table: the 2-cycles close by Theorem 1.1; by Proposition 5.3 the identity permutation does not close, since $Ls$ is not closed (Remark 4.5); the column "Close" counts the $\pi$ for which $\pi(L\bar s)$ and $\pi(Ls)$ are closed under the standard octonion product; v7's description of the closing 3-cycles follows.  v7's restatements of what the table shows are removed: "The 4-cycles are the only type failing without exception; the $(5,2)$ type has the lowest nonzero rate.  The 3-cycles, 5-cycles, and $(3,2,2)$ types close at exactly half."  "No such description is known for the other cycle types, and none is known for the $0\%$ of the 4-cycles." → "No description in terms of $\pi$ like that of the 3-cycles is given here for the other types that close in part, nor a reason for the failure of every 4-cycle."
```

### 5.2 The automorphism note and its script

`paper/automorphism_group_2026-07-12.tex`

**A1.**  Every edit needs the date stamp (DOCUMENT_GENRES.md), and in LaTeX it goes on a first comment line, not inside \date, where it would print in the title block.

Old:
```
\documentclass[11pt]{amsart}
```

New:
```
% Created: 2026-07-12.  Updated: 2026-10-24.
\documentclass[11pt]{amsart}
```

**A2.**  Shows the update in the visible date, in the same form as the Hermes supplement.

Old:
```
\date{12 July 2026}
```

New:
```
\date{12 July 2026; updated 24 October 2026}
```

**A3.**  Stale status: '(v6)', and the sentence about how much of this material will be worked into the paper. The v8 draft's Section 5.4 now states the results and cites this note.

Old:
```
from a $\ZZ_3$-symmetric triple-octonion product} (v6). It is not part of
that paper and is not in its referee scope; the author of the main paper
will decide how, and how much of, this material is worked into it.
```

New:
```
from a $\ZZ_3$-symmetric triple-octonion product}, whose Section~5.4
states the main results of this note and cites it. It is not part of that
paper and is not in its referee scope.
```

**A4.**  Stale: the paper no longer leaves the group open.

Old:
```
The main paper leaves open the identification of $\Aut(\LL, +, \star)$,
```

New:
```
This note determines $\Aut(\LL, +, \star)$,
```

**A5.**  Stale to-do: the docstring now gives the noncircular route. The remark itself stays as part of the note's verification record.

Old:
```
is what the second pass verified. The conclusion is unaffected. The
docstring of \texttt{python\_project/src/verify\_star\_algebra\_structure.py}
still carries the circular wording and should be corrected.
```

New:
```
is what the second pass verified. The conclusion is unaffected.
```

**A6.**  Wrong pointer. It has been wrong since the note was written, so it is not a staleness fix, and it can be dropped on its own. The script does not contain the Route 1 derivation. It checks Route 1's input and refers back to this note for the classification of Aut(S).

Old:
```
Route~2 was carried out in the second, independent pass and does not share a
step with Route~1. Route~1 is the derivation recorded in
\texttt{verify\_star\_algebra\_structure.py}. The two agree.
```

New:
```
Route~2 was carried out in the second, independent pass and does not share a
step with Route~1. Route~1 is the derivation given above; its input, the
Fourier normal form, is checked exactly in
\texttt{verify\_star\_algebra\_structure.py}. The two agree.
```

**A7.**  History and stale status: the subsection quotes and grades the v6 Conclusion ('should be updated', 'was missed', 'now a theorem'). The new text maps the paper's Section 5.4 statements to Theorems 13 and 15 and to $A_6$, and states the stabilizer fact directly. The last bullet is worded so that it sets without an overfull line. The label sec:v6 becomes sec:main-paper; grep finds nothing that references it.

Old:
```
\subsection{Relation to what the main paper (v6) asserts}
\label{sec:v6}

The v6 Conclusion lists as an open problem: ``The automorphism group of
$(\LL, +, \star)$, equivalently the stabilizer of the product tensor
inside~$\Co_0$, and its relationship to the Conway group. Empirically,
$S_3 \subseteq \Aut(\LL, +, \star) \subsetneq \Co_0$: the block-permutation
$S_3$ acts by automorphisms, and the inclusion in $\Co_0$ is strict because
$-I_{24} \in \Co_0$ does not preserve $\star$. Full identification of
$\Aut(\LL, +, \star)$ is open.'' Against the present calculation:

\begin{itemize}[leftmargin=*]
\item \textbf{``Full identification is open'' is now closed.}
  $\Aut(\LL,+,\star) \cong C_6 \times S_3$, of order~$36$, with the three
  explicit generators of \S\ref{sec:lattice}. The open-problem list should be
  updated.
\item \textbf{$S_3 \subseteq \Aut(\LL,+,\star)$ is correct}, and the $S_3$ is
  exactly the block-permutation subgroup. It is a proper subgroup of index~$6$:
  the blockwise $C_6$ was missed.
\item \textbf{The containment $\Aut(\LL,+,\star) \subseteq \Co_0$ is correct
  and is now a theorem}, not an empirical observation. But the paper's phrasing
  presents it as if it were free: the phrase ``equivalently the stabilizer of
  the product tensor inside~$\Co_0$'' already presupposes that every
  $\star$-automorphism of $\LL$ is an isometry. That presupposition is true, by
  Theorem~\ref{thm:co0}, and the only proof we have of it goes through the
  compactness of $\Aut(\OO,\cdot_s)$. If the sentence stays, it needs that
  justification attached.
\item \textbf{The strictness $\subsetneq \Co_0$ is correct}, and the reason
  given ($-I_{24}$ is not multiplicative) is correct. It is a considerable
  understatement: $36$ against $|\Co_0| = 8{,}315{,}553{,}613{,}086{,}720{,}000$.
\end{itemize}
```

New:
```
\subsection{Relation to the main paper}
\label{sec:main-paper}

Section~5.4 of the main paper states the order, structure and generators of
$\Aut(\LL,+,\star)$ and its strict containment in~$\Co_0$, and cites this note
for their derivation.

\begin{itemize}[leftmargin=*]
\item $\Aut(\LL,+,\star) \cong C_6 \times S_3$, of order~$36$
  (Theorem~\ref{thm:main}), with the three explicit generators of
  \S\ref{sec:lattice}. The signed permutation that the paper gives as the
  generator of the~$C_6$ is~$A_6$.
\item The block permutations form a subgroup~$S_3$ of index~$6$; the
  blockwise part is the cyclic group~$C_6$.
\item The containment $\Aut(\LL,+,\star) \subseteq \Co_0$ is
  Theorem~\ref{thm:co0}; the paper's derivation of it, from the ambient
  algebra, follows the proof of that theorem. Consequently
  $\Aut(\LL,+,\star)$ is the stabilizer of the product~$\star$ in~$\Co_0$.
\item The inclusion in~$\Co_0$ is strict because $-I_{24}$ is not
  multiplicative; $|\Aut(\LL,+,\star)| = 36$ against
  $|\Co_0| = 8{,}315{,}553{,}613{,}086{,}720{,}000$.
\end{itemize}
```

**A8.**  The same stale to-do in the note's item C1, plus 'now'. The item itself stays.

Old:
```
(Remark~\ref{rem:circular}). The conclusion survives; the proof is now via
  simplicity of the two ideals and a dimension count. The docstring of
  \texttt{verify\_star\_algebra\_structure.py} still carries the circular
  wording.
```

New:
```
(Remark~\ref{rem:circular}). The conclusion survives; the proof is via
  simplicity of the two ideals and a dimension count.
```

**A9.**  States what the runtimes were measured on, matching the docstring of verify_aut_lambda_star.py as edited by A12. 'The project's reference machine' is defined nowhere in the repository.

Old:
```
are on the project's reference machine. Python scripts live in
```

New:
```
are wall-clock times on one core of an Intel Core i7-12650H; the Python
scripts run single-threaded. Python scripts live in
```

**A10.**  The path was empty.

Old:
```
GAP scripts in \texttt{}; both
```

New:
```
GAP scripts in \texttt{gap\_project/}; both
```

**A11.**  History wording ('earlier', 'recorded', 'subsumed'). The new text says what the probe does.

Old:
```
The earlier probe \texttt{python\_project/src/probe\_aut\_lambda\_star.py}
recorded the $S_3$ of block permutations and the exclusion of $-I_{24}$; those
observations are subsumed here.
```

New:
```
The probe \texttt{python\_project/src/probe\_aut\_lambda\_star.py} tests
individual candidate maps against $\LL$ and~$\star$ on all $576$ basis pairs;
it finds the $S_3$ of block permutations and excludes~$-I_{24}$, both
contained in the results above.
```

`python_project/src/verify_aut_lambda_star.py`

**A12.**  The measured runtime is about 5 minutes, not well over ten. The new wording matches Section 8 of the note, and it drops a spaced dash aside.

Old:
```
Runtime: long -- well over ten minutes single-threaded (the complete
stabilizer enumerations dominate); the script prints its own runtime at
the end.
```

New:
```
Runtime: about 5 minutes on one core of an Intel Core i7-12650H (the
script runs single-threaded; the complete stabilizer enumerations take
about 4 of the 5 minutes).  It prints its own runtime at the end.
```

### 5.3 The Baez–Egan note

`paper/Baez-Egan_Leech_Verification_addendum_2026-05-24.tex`

**B1.**  Date line (current-state rule, DOCUMENT_GENRES.md).

Old:
```
\documentclass[11pt]{amsart}
```

New:
```
% Created: 2026-05-24.  Updated: 2026-10-24.
\documentclass[11pt]{amsart}
```

**B2.**  Self-description: a note on the L/2L parallel that accompanies the 16 April report, not a dated addendum (file name kept).

Old:
```
\title[Addendum (2026-05-24): the $L/2L$ parallel]%
{Addendum to\\
\emph{Verification and Commentary: The Baez--Egan Embedding of the
Leech Lattice in the Exceptional Jordan Algebra}\\
\large The $L/2L$ parallel and how our work builds on it}
```

New:
```
\title[The $L/2L$ parallel]%
{The $L/2L$ parallel and how our work builds on it\\
\large A note accompanying \emph{Verification and Commentary: The
Baez--Egan Embedding of the Leech Lattice in the Exceptional Jordan
Algebra}}
```

**B3.**  Visible date: created; updated.

Old:
```
\date{24 May 2026}
```

New:
```
\date{24 May 2026; updated 24 October 2026}
```

**B4.**  Abstract: current state, pointing to the current paper. Removes 'dated addendum', the 'Re-reading ... we notice' narration, and the v4/Section A.1 pointers.

Old:
```
This is a dated addendum to the original 16 April 2026 Verification
and Commentary on the Baez--Egan embedding
(\texttt{paper/Baez-Egan\_Leech\_Verification.pdf}, 8\,pp).  Re-reading
that work in light of our subsequent paper \emph{An order on the Leech
lattice from a $\ZZ_3$-symmetric triple-octonion product} (v4,
24 May 2026), we notice an explicit lattice-algebraic parallel
between Baez/Egan's ``key lemma'' (Part 9, $L_1 \cap L_2 = 2 L_0$)
and Wilson's identity ($L\bar s \cap Ls = 2 L$), and we identify
exactly how Section~A.1 of our paper (\emph{the mod-2-quotient
reading}) builds on top of that shared identity.  The purpose of
this addendum is to credit the shared lattice-algebraic foundation
correctly and to delineate what is new in our use of it.
```

New:
```
This note accompanies the Verification and Commentary of 16 April
2026 on the Baez--Egan embedding
(\texttt{paper/Baez-Egan\_Leech\_Verification.pdf}, 8\,pp).  It sets
out an explicit lattice-algebraic parallel between Baez/Egan's ``key
lemma'' (Part 9, $L_1 \cap L_2 = 2 L_0$) and Wilson's identity
($L\bar s \cap Ls = 2 L$), and how Appendix~A of our paper
\emph{An order on the Leech lattice from a $\ZZ_3$-symmetric
triple-octonion product} (\texttt{paper/main.tex}, v8), ``The closure
criterion modulo~$2$'', builds on that shared identity.  It credits
the shared lattice-algebraic foundation and the prior work on $L/2L$,
and states what our paper adds.
```

**B5.**  Pointer to the parent report without the addendum's 'original PDF' wording.

Old:
```
Baez/Egan (Part 9, Section~4 of original PDF): &
```

New:
```
Baez/Egan (Part 9; \S4 of the Verification and Commentary): &
```

**B6.**  Found while drafting: Ls-bar and Ls are images of L under RIGHT multiplication, not left; also US spelling 'realized'.

Old:
```
In Wilson's setup, $L = D_8^+$ is the octonion-realised $E_8$; the
sublattices~$L\bar s$ and~$Ls$ are the images of left-multiplication
by $\bar s$ and~$s$ respectively
```

New:
```
In Wilson's setup, $L = D_8^+$ is the octonion-realized $E_8$; the
sublattices~$L\bar s$ and~$Ls$ are the images of~$L$ under right
multiplication by $\bar s$ and~$s$ respectively
```

**B7.**  Found while drafting: wrong scale. The key lemma is used to show that L_L/sqrt2 (not L_L) is self-dual and has no vectors of squared norm 2. L_L has Gram determinant 2^24.

Old:
```
is used in showing both (i)~$L_L$ is even and self-dual, and (ii)~$L_L$
has no vector of squared norm~$<4$.
```

New:
```
is used in showing that $L_L/\sqrt{2}$ is (i)~self-dual and
(ii)~without vectors of squared norm~$2$.
```

**B8.**  Cosmetic, optional: US spelling (editorial_standards.md, standard 4).

Old:
```
and analyse the sublattice
```

New:
```
and analyze the sublattice
```

**B9.**  Pointer to the current Appendix A. It opens with the descent and uses the identity in its next sentence, so it 'uses the identity as the foundation' rather than 'starts from' it.

Old:
```
\section{What our Section A.1 builds on the same identity}
% ===================================================================

Section~A.1 of our paper (\emph{The mod-2-quotient reading of
Tables~A.2 and~A.3}, v4) starts from the same identity, in Wilson's
form:
\begin{equation}\label{eq:key}
  L\bar s \cap Ls \;=\; 2 L, \tag{$\dagger$}
\end{equation}
and uses it as the foundation for an explicitly algebraic reading.
```

New:
```
\section{What our Appendix A builds on the same identity}
% ===================================================================

Appendix~A of our paper (``The closure criterion modulo~$2$'',
\texttt{paper/main.tex}, v8) uses the same identity, in Wilson's
form,
\begin{equation}\label{eq:key}
  L\bar s \cap Ls \;=\; 2 L, \tag{$\dagger$}
\end{equation}
as the foundation for an explicitly algebraic reading.
```

**B10.**  Found while drafting: invalid reason. sigma fixes 2L because sigma(L) = L (Lemma 4.1); being a coordinate permutation of an even integer lattice does not imply this.

Old:
```
  applying~$\sigma$ (which fixes $2 L$, since it is a coordinate
  permutation of an even integer lattice),
```

New:
```
  applying~$\sigma$ (which fixes $2 L$, since $\sigma(L) = L$ by
  Lemma~4.1 in our paper),
```

**B11.**  Precision of a pointer into the current paper: Lemma 4.3 gives a left ideal, and two-sidedness comes from Remark 4.5.

Old:
```
  \(V\) is a (two-sided) ideal of $\bar L$;
```

New:
```
  \(V\) is a left ideal of $\bar L$, two-sided since
  $\sigma(L\bar s)$ is a two-sided ideal of~$L$ (Remark~4.5);
```

**B12.**  Stops calling the polarization reading new. It is prior work for Wilson's pair, and sigma, an isometry of L, carries that pair to (V,W). What the paper adds is the product mu-bar, under which V and W are a two-sided ideal and a subalgebra.

Old:
```
Step~(1) uses identity~\eqref{eq:key} directly.  Steps~(2)--(5) are
new in our paper: they descend the octonion product to an
$\FF_2$-bilinear product on $L / 2 L$, read the sublattices as
ideal/subalgebra in that descended algebra, and recognise the whole
picture as a polarization of an $\FF_2$-orthogonal space.  Baez/Egan
in Parts~9--10 do not take this $\FF_2$-algebra step.
```

New:
```
Step~(1) uses identity~\eqref{eq:key} directly.  Step~(5) reads
$\bar L = V \oplus W$ as a polarization of~$(\bar L, q)$, an ordered
pair of complementary Lagrangians.  For Wilson's pair
$(L\bar s/2L, Ls/2L)$ that reading is prior work (Section~4), and
$\sigma$, an isometry of~$L$, maps that pair to~$(V, W)$.  What our
paper adds is the product~$\bar\mu$ of step~(2), under which $V$
and~$W$ are a two-sided ideal and a subalgebra (step~(4)).
Baez/Egan in Parts~9--10 do not take this $\FF_2$-algebra step.
```

**B13.**  The note gives Egan's proof in one line (Section 1); the proof itself is in section 4 of the report.

Old:
```
  Egan's duality-based proof of this is in the original PDF~\S4 of
  the 16~April 2026 Verification and Commentary, and is reproduced
  here for the reader's convenience.
```

New:
```
  Egan's duality-based proof of this is in~\S4 of the Verification
  and Commentary of 16~April 2026; Section~1 above summarizes it in
  one line.
```

**B14.**  Replaces the 'not in Baez/Egan and new in our paper' list with the credits of the current Appendix A: Lepowsky-Meurman section 2 (Tits per their introduction), Silverman (part 12 and comments), Calderbank. What is new is stated once, in B12.

Old:
```
What is not in Baez/Egan and is new in our paper:

\begin{itemize}
\item The descent of the octonion product modulo $2 L$ to an
  $\FF_2$-bilinear product on $L / 2L$.
\item The identification of $V$ and $W$ as ideal and subalgebra of
  this descended algebra.
\item The Witt-decomposition (polarization) reading of $L / 2L$ into
  complementary Lagrangians under the natural plus-type quadratic
  form.
\end{itemize}
```

New:
```
The polarization reading of~$\bar L$ is prior work.  Lepowsky and
Meurman (1982, \S2) identify the lattices between $2L$ and~$L$ that
are isometric to $\sqrt{2}\,L$ with the preimages of the Lagrangians
of~$(\bar L, q)$, and construct a lattice isometric to
$\sqrt{2}\,\LL$ from any polarization; their introduction notes that
work of Tits contains what amounts to the same construction.  Egan's
construction is of this form.  In Baez's series, Tim Silverman
arrived at the identification of these sublattices with the
Lagrangians independently (Part~12).  Silverman's comments of 16 and
19~February 2016 on Part~12 read Egan's construction modulo~$2$ and
find that Wilson's pair $(L\bar s/2L, Ls/2L)$ is a polarization.
David Calderbank related Egan's construction to that of Lepowsky and
Meurman (comment of 2~April 2023 on Part~9).
```

**B15.**  Pointer to the current paper. Removes the 'v4' and 'made explicit in this addendum' narration and the spaced ' -- ' aside.

Old:
```
The main paper (v4, Section~A.1) cites Wilson~2009 for the
sublattice identities and notes Baez/Egan in~\S6 as part of the
related-work discussion.  The cross-credit recorded here -- that
our Section~A.1 builds on the same lattice-algebraic key lemma
that Egan's Part~9 stripped to its geometric core -- is implicit
in the main paper's references but is made explicit in this
addendum.
```

New:
```
The main paper (\texttt{paper/main.tex}, v8) states the sublattice
identities in Wilson's setting (Section~2.3) and Egan's construction
in Section~6, for sublattices $L_1, L_2 \subset L$ isometric to
$\sqrt{2}\,L$ with $L_1 \cap L_2 = 2L$; Wilson's lattice is the case
$(L_1, L_2) = (L\bar s, Ls)$.  The paragraph ``Relation to prior
work'' of its Appendix~A gives the credits above.  This note sets out
the parallel between Egan's key lemma and Wilson's
identity~\eqref{eq:key}, on which Appendix~A builds.
```

**B16.**  Pointer to the current Appendix A; US spelling 'generalized'.

Old:
```
$\FF_2$-algebra picture of Section~A.1?  In their setup the
generalised lattice
```

New:
```
$\FF_2$-algebra picture of Appendix~A?  In their setup the
generalized lattice
```

**B17.**  L/2L is 8-dimensional over F_2, not 24.

Old:
```
within $24$-dimensional $L / 2L$
```

New:
```
within the $8$-dimensional $L / 2L$
```

**B18.**  References: names the report without 'Original', points to the current paper, and adds the sources the new credits use (Silverman, Calderbank, Lepowsky-Meurman).

Old:
```
\item Original Verification PDF:
  \texttt{paper/Baez-Egan\_Leech\_Verification.pdf}, 16 April 2026
  (Claude, Anthropic).
\item Our main paper:
  \emph{An order on the Leech lattice from a $\ZZ_3$-symmetric
  triple-octonion product}, v4 (24 May 2026), Section~A.1
  ``The mod-2-quotient reading of Tables A.2 and A.3''.
\item Wilson, R.~A., \emph{Octonions and the Leech lattice},
  J.~Algebra \textbf{322} (2009), 2186--2190.
\item Baez, J.~C.\ (with G.~Egan and others),
  \emph{Integral octonions} (parts 9 and 10),
  $n$-Category Caf\'e, November--December 2014.
```

New:
```
\item Verification and Commentary:
  \texttt{paper/Baez-Egan\_Leech\_Verification.pdf}, 16 April 2026
  (Claude, Anthropic).
\item Our main paper:
  \emph{An order on the Leech lattice from a $\ZZ_3$-symmetric
  triple-octonion product}, v8 (\texttt{paper/main.tex}), Appendix~A
  ``The closure criterion modulo~$2$'' and Section~6 ``Related work''.
\item Wilson, R.~A., \emph{Octonions and the Leech lattice},
  J.~Algebra \textbf{322} (2009), 2186--2190.
\item Baez, J.~C.\ (with G.~Egan, T.~Silverman and others),
  \emph{Integral octonions} (parts 9, 10 and~12),
  $n$-Category Caf\'e, November--December 2014 (parts 9 and~10) and
  31~January 2016 (part~12, a guest post by T.~Silverman); comments
  by T.~Silverman of 16 and 19~February 2016 on part~12, and by
  D.~Calderbank of 2~April 2023 on part~9.
\item Lepowsky, J., and Meurman, A., \emph{An $E_8$-approach to the
  Leech lattice and the Conway group}, J.~Algebra \textbf{77} (1982),
  484--504.
```

`paper/README.md`

**B19.**  Optional: the index line still calls the note an addendum after B2 and B4.

Old:
```
and the addendum comparing their L/2L treatment with ours
```

New:
```
and a note on the L/2L parallel and how our paper builds on it
```

### 5.4 The companion note

`paper/companion.tex`

**C1.**  Adds the date line that DOCUMENT_GENRES.md requires when a current-state document is edited, as the first-line comment used by the LaTeX notes.

Old:
```
\documentclass[11pt]{amsart}
```

New:
```
% Created: 2026-04-14.  Updated: 2026-10-24.
\documentclass[11pt]{amsart}
```

**C2.**  Puts the visible date in the repository form '<created>; updated 24 October 2026' and removes the history clutter 'revised 23 May 2026 (v4 companion)'.

Old:
```
\date{16 April 2026; revised 23 May 2026 (v4 companion)}
```

New:
```
\date{14 April 2026; updated 24 October 2026}
```

**C3.**  Removes the dated revision note from the abstract and replaces the main paper's early title with its current title.

Old:
```
(drafted by
Claude Opus~4.6 in max-thinking mode, revised 23 May 2026 by
Claude Opus~4.7, Anthropic) at the direction of Jens K\"oplinger,
the author of the accompanying paper ``An order on the Leech lattice
from transposition-twisted octonions.''
```

New:
```
(drafted by
Claude Opus~4.6 in max-thinking mode, revised by
Claude Opus~4.7, Anthropic) at the direction of Jens K\"oplinger,
the author of the accompanying paper ``An order on the Leech lattice
from a $\ZZ_3$-symmetric triple-octonion product.''
```

**C4.**  The pointer 'follows the conventions of the main paper' is false for the symbol ⋆. In the note, ⋆ is the untwisted product and ⋆_σ the twisted one, while the main paper's ⋆ is the twisted product. Without this note, Section 11's sentence 'neither says that Λ is closed under ⋆, which is false' reads as a denial of the paper's Theorem 1.1. The note's squared norm ‖v‖² is also what the main paper calls the norm (the squared length), so the sentence names both differences.

Old:
```
This note follows the conventions of the main paper.  We work with
```

New:
```
This note follows the conventions of the main paper, with two
differences of notation: the squared norm~$\|v\|^2$ here is what the
main paper calls the norm, and $\star$ here is the untwisted triple
product and $\star_\sigma$ the $\sigma$-twisted one, which the main
paper writes $\star_{\mathrm{id}}$ and~$\star$.  We work with
```

**C5.**  Stale: the main paper no longer establishes L·L ⊆ L by computing 64 products. It proves it as Lemma 4.2 from the closure of the maximal order and records the basis check as validation.

Old:
```
This is a classical result of Coxeter~(1946),
verified in the main paper by computing all $8 \times 8 = 64$
products of a $\ZZ$-basis of~$L$ in exact arithmetic and checking
that each lies again in~$L$.
```

New:
```
This is, at the scale where roots have squared norm~$2$, a classical
result of Coxeter~(1946): the maximal orders of integral octonions
are closed under multiplication.  The main paper proves it as
Lemma~4.2, from the closure of one maximal order, and records an
independent check of all $8 \times 8 = 64$ products of a
$\ZZ$-basis of~$L$ in exact arithmetic.
```

**C6.**  Updates the pointer: the main paper states this pair as the witness in Remark 4.5. 'The main paper's exact-arithmetic verification' no longer describes anything in the paper.

Old:
```
Both $a$ and $b$ are in $Ls$ (they appear as basis vectors of~$Ls$
produced in the main paper's exact-arithmetic verification).
```

New:
```
Both $a$ and $b$ are in $Ls$ (they are the elements $u_1$
and~$u_2$ of Remark~4.5 of the main paper).
```

**C7.**  Stale: the script no longer has a 'Lemma D' (it labels this check 'LEMMA 4.4', the paper's number). The inclusion is the paper's Lemma 4.4, which the paper proves from the closure of the maximal order, not from this basis check.

Old:
```
Using the basis generated in the main paper's verification, all
$8 \times 8 = 64$ products of basis vectors of~$\sigma(Ls)$, under
the standard octonion product~$\cdot$, lie again in~$\sigma(Ls)$.
The verification is in
\href{../python_project/src/symbolic_proof_checks.py}{\texttt{symbolic\_proof\_checks.py}}
as Lemma~D.

Hence $\sigma(Ls) \cdot \sigma(Ls) \subseteq \sigma(Ls)$.
```

New:
```
All $8 \times 8 = 64$ products of these basis vectors, under the
standard octonion product~$\cdot$, lie again in~$\sigma(Ls)$.  The
verification, in exact arithmetic, is in
\href{../python_project/src/symbolic_proof_checks.py}{\texttt{symbolic\_proof\_checks.py}}.

Hence $\sigma(Ls) \cdot \sigma(Ls) \subseteq \sigma(Ls)$.  This is
Lemma~4.4 of the main paper, which proves it from the closure of one
maximal order of integral octonions.
```

**C8.**  The note used Λ both for the unit-scale Leech lattice (minimal norm 4) and for Wilson's lattice in L^3 (display eq:leech), which is the √2-scaled copy. The main paper's convention is the scaled copy, with minimal norm 8.

Old:
```
The Leech lattice~$\LL$ is the unique even unimodular lattice of
rank~$24$ with no roots (i.e.\ no vectors of squared norm~$2$).  Its
minimal vectors have squared norm~$4$, and there are exactly
$196{,}560$ of them.
```

New:
```
The Leech lattice is the unique even unimodular lattice of
rank~$24$ with no roots (i.e.\ no vectors of squared norm~$2$).  Its
minimal vectors have squared norm~$4$, and there are exactly
$196{,}560$ of them.  As in the main paper, $\LL$ denotes its
$\sqrt{2}$-scaled copy, at the scale of Wilson's construction below,
whose $196{,}560$ minimal vectors have squared norm~$8$.
```

**C9.**  Corrects the stale theorem number (the known correction).

Old:
```
This is Theorem~3.1 of the main paper.
```

New:
```
This is Theorem~1.1 of the main paper.
```

### 5.5 The historical appendix draft

`paper/historical_appendix.tex`

**H1.**  Adds the date line that DOCUMENT_GENRES.md requires when a current-state document is edited.

Old:
```
\documentclass[11pt]{amsart}
```

New:
```
% Created: 2026-05-23.  Updated: 2026-10-24.
\documentclass[11pt]{amsart}
```

**H2.**  Puts the visible date in the repository form '<created>; updated 24 October 2026'.

Old:
```
\date{23 May 2026}
```

New:
```
\date{23 May 2026; updated 24 October 2026}
```

**H3.**  Removes history clutter: the comment narrates a past plan for v4. The replacement points to the paper's current Appendix B.

Old:
```
% Body intended for v4 as an appendix.  Target length in the
% main paper: <= 2 pages (target 1 page).
```

New:
```
% The main paper carries this note, in its own wording, as
% Appendix B (label app:history).
```

**H4.**  Accuracy: 're-derived in exact arithmetic' overstates the Mahler check. The covering-radius theorem was checked by exact-arithmetic sampling, as the paper's Appendix B states.

Old:
```
verification notes in
\texttt{paper/reviews/}).
```

New:
```
verification notes in
\texttt{paper/reviews/}), for Mahler's covering-radius theorem by
exact-arithmetic sampling.
```

**H5.**  Corrects Mahler's Theorem 2 (the known correction): the generator is an integer multiple δG* of an integral Cayley number G* of norm 1, 2 or 4.

Old:
```
and every left or right ideal is principal, generated by an integral
Cayley number of norm $1$, $2$, or $4$.
```

New:
```
and every left or right ideal is principal, generated by an integer
multiple of an integral Cayley number of norm $1$, $2$, or $4$.
```

### 5.6 The minimum-norm argument and Wilson's note on Dixon

`python_project/src/leech_wilson.py`

**W1.**  The step "at least one of x, y, z must be 0" does not follow from N_std < 8. Three nonzero coordinates can give N_std = 6, and 61,440 root triples, such as (λ, λ, λ), satisfy conditions 1 and 2. The two-nonzero case also never reached a conclusion. The replacement is Wilson's own argument: evenness from an identity that uses condition 3, then the zero-coordinate cases.

Old:
```
Suppose (x, y, z) ∈ Λ has Wilson-norm N < 4 (i.e., N_std < 8).
Then at least one of x, y, z must be 0 (since each nonzero coordinate from L
has N ≥ 2).  If, say, z = 0, then conditions 2 and 3 give x+y ∈ Ls̄ and
x+y ∈ Ls, so x+y ∈ Ls̄ ∩ Ls = 2L.  Also y+z = y ∈ Ls̄, and x+z = x ∈ Ls̄,
so both x,y ∈ Ls̄ = 2B ⊂ L.  Then y ∈ Ls̄ and x = (x+y) − y ∈ 2L + 2B = 2B.
In the only-one-nonzero case (y=z=0), condition 3 gives x ∈ Ls, and
condition 2 gives x ∈ Ls̄, so x ∈ Ls ∩ Ls̄ = 2L, hence N(x) ≥ N_std(2L) = 8.
Contradiction.  Therefore Λ has no vectors of Wilson-norm < 4 (N_std < 8).
```

New:
```
For (x, y, z) ∈ Λ, expanding the norms gives
    N_std(x, y, z) = N(x+y) + N(x+z) + N(y+z) − N(x+y+z).
By conditions 2 and 3 the terms on the right are norms of elements of Ls̄
and Ls, and N(ls̄) = N(ls) = 2N(l) for l ∈ L; as L is even with minimum
norm 2, these lie in 4Z and are at least 4 when l ≠ 0.  So N_std ∈ 4Z: Λ is
even under Wilson's norm N_Wilson = N_std/2.  A nonzero vector of Λ with
N_std < 8 therefore has N_std = 4, so one coordinate is 0 (three nonzero
ones would give N_std ≥ 6); by symmetry, say z = 0.  Then x = x+z and
y = y+z lie in Ls̄, so exactly one of them is nonzero (two would give
N_std ≥ 8), say x.  Then x = x+y+z ∈ Ls as well, so x ∈ Ls̄ ∩ Ls = 2L and
N(x) ≥ 8, a contradiction.  Hence every nonzero vector of Λ has N_std ≥ 8
(N_Wilson ≥ 4).
```

`evidence_and_reasoning/key_claims/004_wilson_leech_verified.txt`

**K1.**  Date line (current-state rule); the file has a bare Date line.

Old:
```
Date: 2026-04-07
```

New:
```
Created: 2026-04-07.  Updated: 2026-10-24.
```

**K2.**  Wilson's bracketed note on Dixon (Section 5, p. 2189) concerns the multiplication by j in the first coordinate, not the s versus s-bar distinction, and it is not a footnote; quoted from the published text.

Old:
```
Wilson's membership conditions.  This distinction is verified computationally
and is also explicitly noted in Wilson's paper:
  "Dixon's formula is equivalent to (λs, ±λk, ±(λj)k), that is, omitting the
   (necessary) multiplication by j in the first coordinate."
  — [Wilson2009 p. 2189, footnote on Dixon's construction]
```

New:
```
Wilson's membership conditions; this is verified computationally.  Wilson's
bracketed note on Dixon's conjectured formula [Wilson2009 Section 5, p. 2189]
concerns a different point, the multiplication by j in the first coordinate:
  "in effect, his formula is equivalent to (λs, ±λk, ±(λj)k), that is,
   omitting the (necessary) multiplication by j in the first coordinate."
```

**K3.**  The old argument jumps from N_std < 8 to a zero coordinate, which fails for three nonzero coordinates (N_std = 6 is possible for conditions 1 and 2); Wilson's own argument (Section 4, p. 2189) first shows N_std in 4Z.

Old:
```
No norm-2 vectors (proven analytically):
  If (x, y, z) ∈ Λ has N_std < 8, at least one coordinate is 0.
  If z = 0: conditions (2) and (3) give x+y ∈ Ls̄ and x+y ∈ Ls, so
    x+y ∈ Ls̄ ∩ Ls = 2L.  Also x ∈ Ls̄ and y ∈ Ls̄.
  If y = z = 0: condition (3) gives x ∈ Ls, and condition (2) gives x ∈ Ls̄,
    so x ∈ 2L, hence N(x) ≥ 8.  Contradiction.
```

New:
```
No norm-2 vectors (Wilson's argument, [Wilson2009 Section 4, p. 2189]):
  For (x, y, z) ∈ Λ, N_std(x, y, z) = N(x+y) + N(x+z) + N(y+z) − N(x+y+z).
  The terms on the right are norms of elements of Ls̄ and Ls, and
  N(l s̄) = N(l s) = 2N(l) for l ∈ L, which lies in 4Z and is at least 4 for
  l ≠ 0; so N_std ∈ 4Z.  A nonzero vector with N_std < 8 therefore has
  N_std = 4, so one coordinate is 0 (three nonzero ones give N_std ≥ 6), say
  z = 0.  Then x and y lie in Ls̄, so only one of them is nonzero (two would
  give N_std ≥ 8), say x; then x = x+y+z ∈ Ls as well, so x ∈ Ls̄ ∩ Ls = 2L
  and N(x) ≥ 8, a contradiction.
```

**K4.**  A rank check does not prove Ls̄ + Ls = L (rank 8 and index 16 would hold even if Ls̄ = Ls); the identity is Wilson's (Section 2, p. 2187), as G-S6a to G-S6c say for the GAP test.

Old:
```
  Ls̄ + Ls = L   (combined generators span L; verified by rank check)
```

New:
```
  Ls̄ + Ls = L   (Wilson Section 2; the combined generators have rank 8)
```

`python_project/src/leech_wilson.py`

**W2.**  Same misattribution as K2: Wilson's note concerns the multiplication by j; the s versus s-bar distinction is this project's computational check (key claim 004).

Old:
```
  Note: type-2 uses s̄; type-3 uses s.  This is the distinction Wilson proves is
  necessary: see [Wilson2009 p. 2189] footnote on Dixon's related construction.
```

New:
```
  Note: type-2 uses s̄; type-3 uses s.  With s̄ in the first coordinate the
  type-3 triples have norm 8 but fail Wilson's conditions (key claim 004).
```

**W3.**  Quote Wilson's bracketed note exactly and place it (Section 5, p. 2189).

Old:
```
    Wilson notes: "Dixon's formula is equivalent to (λs, ±λk, ±(λj)k),
    omitting the (necessary) multiplication by j in the first coordinate."
    [Wilson2009, p. 2189 footnote]
```

New:
```
    Wilson's bracketed note on Dixon's conjectured formula ([Wilson2009]
    Section 5, p. 2189) concerns the multiplication by j in the first
    coordinate: "in effect, his formula is equivalent to (λs, ±λk, ±(λj)k),
    that is, omitting the (necessary) multiplication by j in the first
    coordinate."
```

`python_project/tests/test_leech_wilson.py`

**T1.**  The type-3 vectors use s, as the module builds them (and Wilson, p. 2188).

Old:
```
 Type-3 minimal vectors  —  ((λs̄)j, λk, (λj)k) and 2 cyclic perms:
```

New:
```
 Type-3 minimal vectors  —  ((λs)j, λk, (λj)k) and 2 cyclic perms:
```

**T2.**  As T1.

Old:
```
    Type-3: ((λs̄)j, λk, (λj)k) and 2 cyclic permutations, λ ∈ L, j,k ∈ J.
```

New:
```
    Type-3: ((λs)j, λk, (λj)k) and 2 cyclic permutations, λ ∈ L, j,k ∈ J.
```

**T3.**  As T1 (the value is unchanged, since N(s) = N(s̄) = 2).

Old:
```
    Norm:  N_std = N(λ)N(s̄)N(j) + N(λ)N(k) + N(λ)N(j)N(k) = 4 + 2 + 2 = 8.
```

New:
```
    Norm:  N_std = N(λ)N(s)N(j) + N(λ)N(k) + N(λ)N(j)N(k) = 4 + 2 + 2 = 8.
```

**T4.**  As K4.

Old:
```
15.  Combined Z-span of Ls̄ and Ls vectors has rank 8  (Ls̄ + Ls spans L)
```

New:
```
15.  Combined Z-span of Ls̄ and Ls vectors has rank 8
```

**T5.**  As K4.

Old:
```
    The Ls̄ + Ls = L direction is tested via rank: combining all Ls̄ and Ls
    generators, the matrix rank is 8 (= rank of L).
```

New:
```
    The rank check shows only that the combined Ls̄ and Ls generators have
    rank 8 (= rank of L); Ls̄ + Ls = L is Wilson's statement (Section 2).
```

**T6.**  As K4: rank 8 and the index bound do not imply Ls̄ + Ls = L.

Old:
```
        This confirms Ls̄ + Ls ⊇ L in the vector-space sense; together with
        the index bound this implies Ls̄ + Ls = L.
```

New:
```
        This shows only that the combined generators have rank 8;
        Ls̄ + Ls = L is Wilson's statement (Section 2).
```

### 5.7 The GAP files

`gap_project/tests/test_sublattices.g`

**G-S6a.**  Item 6 checks the rank of the joint span only; Ls + Lsbar = L is stated in the paper (Section 2.3).

Old:
```
##    6.  Ls + Lsbar generates L (joint Z-span has rank 8).
```

New:
```
##    6.  The joint Z-span of Ls and Lsbar has rank 8 (Ls + Lsbar = L is
##        paper Section 2.3).
```

**G-S6b.**  As G-S6a.

Old:
```
# 6. Ls + Lsbar generates L (joint span has rank 8).
```

New:
```
# 6. The joint span of Ls and Lsbar has rank 8 (Ls + Lsbar = L: paper Section 2.3).
```

**G-S6c.**  As G-S6a: the CHECK label names what is checked.

Old:
```
CHECK("Z-span of LS_BASIS u LS_BAR_BASIS has rank 8  (Ls + Lsbar = L)",
```

New:
```
CHECK("Z-span of LS_BASIS u LS_BAR_BASIS has rank 8",
```

`gap_project/README.md`

**G-R.**  As G-S6a.

Old:
```
reject every E8 root; $2L \subseteq Ls \cap L\bar s$; $Ls + L\bar s = L$;
```

New:
```
reject every E8 root; $2L \subseteq Ls \cap L\bar s$; $Ls + L\bar s$ of rank 8;
```

`gap_project/src/octonion.g`

**Go1.**  Adds the date line DOCUMENT_GENRES.md requires on an edited document. Created is the file's first commit.

Old:
```
##  octonion.g — Octonion algebra over Q^8 with the standard Fano-plane
##  multiplication used in the paper.
```

New:
```
##  octonion.g — Octonion algebra over Q^8 with the standard Fano-plane
##  multiplication used in the paper.
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
```

**Go2.**  Stale equation number: the Fano triples are paper eq. (2.1).

Old:
```
##  Standard Fano triples (paper eq. (1)):
```

New:
```
##  Standard Fano triples (paper eq. (2.1), `eq:fano`):
```

`gap_project/src/harness.g`

**Gh1.**  The usage example calls TEST(...), but the harness defines CHECK(...). The example and the paragraph after it now show what the test files and run_all.g actually do: run_all.g prints the suite totals and sets the exit code, not PrintTestSummary. Adds the date line.

Old:
```
##  harness.g — Tiny test harness used by the test scripts in tests/.
##
##  Usage in a test file:
##      Read("src/harness.g");                     # initialises counters
##      TEST("description of the check", boolean_expression);
##      ...
##  At the end of run_all.g (or a single test file), the totals are printed
##  via PrintTestSummary(); a non-zero failure count causes a non-zero exit
##  (when running as a script, see run_all.g).
```

New:
```
##  harness.g — Tiny test harness used by the test scripts in tests/.
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
##
##  Usage in a test file (as in tests/*.g):
##      Read(Concatenation(GAP_PROJECT_ROOT, "/src/harness.g"));
##      LoadAllSrc(GAP_PROJECT_ROOT);
##      ResetTestCounters();
##      CHECK("description of the check", boolean_expression);
##      ...
##      PrintTestSummary();
##  PrintTestSummary() prints the totals of the test file.  run_all.g adds
##  up the totals of all test files, prints the sum, and exits with code 0
##  if every check passed and 1 otherwise.
```

`gap_project/tests/test_octonion.g`

**Gt1.**  Adds the date line to the edited document.

Old:
```
##  test_octonion.g — Octonion algebra basics (paper Section 2.1).
##
```

New:
```
##  test_octonion.g — Octonion algebra basics (paper Section 2.1).
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
##
```

**Gt2.**  Stale equation number in the CHECK label: the Fano triples are paper eq. (2.1).

Old:
```
for a = 1..7  (paper eq. (1))", ok);
```

New:
```
for a = 1..7  (paper eq. (2.1))", ok);
```

`gap_project/tests/test_twist.g`

**Gw1.**  The paper has no Proposition 3.2 or Remark 3.5. Section 3 has Definitions 3.1 to 3.3; the identity sigma(xy) = sigma(x) ._sigma sigma(y) is stated after Definition 3.2; the equivalence of the transpositions is the first footnote of Section 4. Item 5 said *_sigma is 'not merely an algebra-isomorphic relabelling', but the paper, and item 3 of this test, make sigma an algebra isomorphism from (O, .) to (O, ._sigma). The title uses the term of Definition 3.2, 'sigma-twisted octonion product'. Adds the date line.

Old:
```
##  test_twist.g — The transposition-twisted octonion algebra (paper
##  Definition 3.1, Proposition 3.2, Remark 3.5).
##
##  Verifies:
##    1.  ApplySigma is an involution: sigma o sigma = id  on R^8.
##    2.  ApplySigma is an isometry: it preserves the Euclidean norm.
##    3.  Proposition 3.2 (paper):  sigma( x * y )  =  sigma(x) *_sigma sigma(y).
##    4.  *_sigma satisfies the composition law N(x *_sigma y) = N(x) N(y) on
##        all 64 basis pairs.
##    5.  *_sigma differs from * on at least one basis pair (so it is a
##        DIFFERENT bilinear product on R^8, not merely an algebra-isomorphic
##        relabelling of the same product).
```

New:
```
##  test_twist.g — The sigma-twisted octonion product (paper Definitions
##  3.1 and 3.2, `def:twist` and `def:twisted-product`), for
##  sigma = (1 2).  All transposition twists are equivalent up to
##  relabeling (paper Section 4, `sec:proof`, footnote to sigma = (1 2)).
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
##
##  Verifies:
##    1.  ApplySigma is an involution: sigma o sigma = id  on R^8.
##    2.  ApplySigma is an isometry: it preserves the Euclidean norm.
##    3.  sigma( x * y )  =  sigma(x) *_sigma sigma(y): sigma is an algebra
##        isomorphism from (O, *) to (O, *_sigma), as the paper states after
##        Definition 3.2.
##    4.  *_sigma satisfies the composition law N(x *_sigma y) = N(x) N(y) on
##        all 64 basis pairs.
##    5.  *_sigma differs from * on at least one basis pair: it is a
##        different bilinear product on R^8, isomorphic to * via sigma
##        (item 3).
```

**Gw2.**  Stale companion section: L ._sigma L in L is the Proposition of companion Section 10, 'Could the sigma-twist undo Coxeter's fix?'.

Old:
```
##    7.  Companion Section 8: the twisted product preserves L:
```

New:
```
##    7.  Companion Section 10: the twisted product preserves L:
```

**Gw3.**  The paper has no Proposition 3.2; the identity is stated after Definition 3.2.

Old:
```
# 3. Proposition 3.2:  sigma(x*y) = sigma(x) *_sigma sigma(y).
```

New:
```
# 3. sigma(x*y) = sigma(x) *_sigma sigma(y)  (paper, after Definition 3.2).
```

**Gw4.**  The CHECK label cited a Proposition 3.2 that does not exist. It now names the statement; the paper pointer is in the comment above the check and in the header.

Old:
```
CHECK("Proposition 3.2:  sigma(x*y) = sigma(x) *_sigma sigma(y)  on all 64 basis pairs",
```

New:
```
CHECK("sigma is an algebra isomorphism:  sigma(x*y) = sigma(x) *_sigma sigma(y)  on all 64 basis pairs",
```

**Gw5.**  Stale companion section; the content is in Section 10. The comment now states what the companion's Proposition says, instead of naming a failure regime after Kirmse (fair-attribution rule).

Old:
```
# 7. L *_sigma L  subset L  (companion Section 8: Kirmse's regime not reached).
```

New:
```
# 7. L *_sigma L  subset L  (companion Section 10: the twist keeps L closed).
```

**Gw6.**  Stale companion section in the CHECK label.

Old:
```
CHECK("Companion Section 8:  L *_sigma L  subset L   (twist preserves L)",
```

New:
```
CHECK("Companion Section 10:  L *_sigma L  subset L   (twist preserves L)",
```

`gap_project/tests/test_leech_wilson.g`

**Gl1.**  Stale definition number: Wilson's Leech lattice is paper Definition 2.1. Adds the date line.

Old:
```
##  test_leech_wilson.g — Wilson's Leech lattice Lambda inside L^3
##  (paper Definition 2.4, Wilson [Wilson2009] Section 3).
```

New:
```
##  test_leech_wilson.g — Wilson's Leech lattice Lambda inside L^3
##  (paper Definition 2.1, `def:wilson-leech`; Wilson [Wilson2009] Section 3).
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
```

`gap_project/tests/test_sublattices.g`

**Gs1.**  Adds the date line to the edited document.

Old:
```
##  test_sublattices.g — The sublattices Ls, Lsbar, sigma(Ls), sigma(Lsbar).
##
```

New:
```
##  test_sublattices.g — The sublattices Ls, Lsbar, sigma(Ls), sigma(Lsbar).
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
##
```

**Gs2.**  Stale companion section. Companion Section 5 ('What is twisted and what is invariant') states sigma(Ls) != Ls with the witness (-1,0,1,0,0,1,0,1). Section 7 is Wilson's Leech lattice and does not state it.

Old:
```
##    8.  sigma(Ls) ne Ls  (paper Remark 4.5; companion Section 7).
```

New:
```
##    8.  sigma(Ls) ne Ls  (paper Remark 4.5; companion Section 5).
```

**Gs3.**  Stale companion example number: 'Standard product fails on Ls' is companion Example 4.2.

Old:
```
##    9.  Ls * Ls  not subset Ls  (companion Example 4.5; main paper Remark 4.5).
```

New:
```
##    9.  Ls * Ls  not subset Ls  (companion Example 4.2; main paper Remark 4.5).
```

**Gs4.**  Stale companion example number.

Old:
```
# 9. Ls * Ls  not subset Ls  (companion Example 4.5).
```

New:
```
# 9. Ls * Ls  not subset Ls  (companion Example 4.2).
```

**Gs5.**  Stale companion example number in the CHECK label.

Old:
```
CHECK("Companion Example 4.5:  Ls * Ls  NOT subset Ls   (explicit witness pair found)",
```

New:
```
CHECK("Companion Example 4.2:  Ls * Ls  NOT subset Ls   (explicit witness pair found)",
```

`gap_project/tests/test_companion_examples.g`

**Gc1.**  The checks come from companion Sections 3, 4, 5 and 11, not 2 to 8. Adds the date line.

Old:
```
##  to give the reader a fully independent re-derivation.
##
##  Verifies (companion Sections 2-8):
```

New:
```
##  to give the reader a fully independent re-derivation.
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
##
##  Verifies (companion Sections 3 to 5 and 11):
```

**Gc2.**  Stale companion example number.

Old:
```
##    5.  Companion Example 4.5:
```

New:
```
##    5.  Companion Example 4.2 (`ex:Ls-fails`):
```

**Gc3.**  There is no companion Example 4.6. The product e_1 ._sigma e_2 = -e_4 is worked out in the 'Concrete' part of Section 5, which has no numbered example.

Old:
```
##    6.  Companion Example 4.6 / Section 5: under  sigma = (1 2),
```

New:
```
##    6.  Companion Section 5 (`sec:swap`): under  sigma = (1 2),
```

**Gc4.**  Stale companion section. sigma(L) = L, sigma(Ls) != Ls and the witness are stated in Section 5 ('What is twisted and what is invariant').

Old:
```
##    7.  Companion Section 8: sigma fixes L  (sigma(L) = L) but moves Ls
```

New:
```
##    7.  Companion Section 5: sigma fixes L  (sigma(L) = L) but moves Ls
```

**Gc5.**  The summary is companion Section 11, not a box.

Old:
```
##    8.  Companion final box ("Summary of concrete data"):
```

New:
```
##    8.  Companion Section 11 ("Summary of concrete data"):
```

**Gc6.**  Stale companion example number.

Old:
```
# 5. Companion Example 4.5: explicit a, b in Ls with a*b not in Ls.
```

New:
```
# 5. Companion Example 4.2: explicit a, b in Ls with a*b not in Ls.
```

**Gc7.**  Stale companion example number in five CHECK labels. Only the labels change; the code lines between them stay as they are.

Old:
```
CHECK("Companion Example 4.5: a in Ls",
      IsInLs(a));
CHECK("Companion Example 4.5: b in Ls",
      IsInLs(b));
CHECK("Companion Example 4.5: a * b = (0, -2, 2, 1, -2, -1, -1, -1)",
      OctMult(a, b) = ab_expected);
CHECK("Companion Example 4.5: a * b is in L",
      IsInL(OctMult(a, b)));
CHECK("Companion Example 4.5: a * b is NOT in Ls   (so Ls * Ls not subset Ls)",
```

New:
```
CHECK("Companion Example 4.2: a in Ls",
      IsInLs(a));
CHECK("Companion Example 4.2: b in Ls",
      IsInLs(b));
CHECK("Companion Example 4.2: a * b = (0, -2, 2, 1, -2, -1, -1, -1)",
      OctMult(a, b) = ab_expected);
CHECK("Companion Example 4.2: a * b is in L",
      IsInL(OctMult(a, b)));
CHECK("Companion Example 4.2: a * b is NOT in Ls   (so Ls * Ls not subset Ls)",
```

**Gc8.**  Stale companion section. 'Lemma A' names no lemma in either the paper or the companion; sigma(L) = L is paper Lemma 4.1.

Old:
```
# 7. Companion Section 8: sigma fixes L (Lemma A) but moves Ls (sigma(Ls) ne Ls).
CHECK("Companion Section 8: sigma(L) = L   (sigma(b) is a Z-comb of L_BASIS for all 8 basis b)",
```

New:
```
# 7. Companion Section 5: sigma fixes L (paper Lemma 4.1) but moves Ls (sigma(Ls) ne Ls).
CHECK("Companion Section 5: sigma(L) = L   (sigma(b) is a Z-comb of L_BASIS for all 8 basis b)",
```

**Gc9.**  Stale companion section in two CHECK labels.

Old:
```
CHECK("Companion Section 8 witness  (-1, 0, 1, 0, 0, 1, 0, 1)  in sigma(Ls)",
      IsInSigmaLs(witness));
CHECK("Companion Section 8 witness  (-1, 0, 1, 0, 0, 1, 0, 1)  NOT in Ls",
```

New:
```
CHECK("Companion Section 5 witness  (-1, 0, 1, 0, 0, 1, 0, 1)  in sigma(Ls)",
      IsInSigmaLs(witness));
CHECK("Companion Section 5 witness  (-1, 0, 1, 0, 0, 1, 0, 1)  NOT in Ls",
```

**Gc10.**  The summary is companion Section 11, not a box.

Old:
```
# 8. Companion final box: three concrete data points.
```

New:
```
# 8. Companion Section 11: three concrete data points.
```

`gap_project/tests/test_lemmas.g`

**Gm1.**  A program states where the paper cites it (DOCUMENT_GENRES.md, 'Kinds of current-state document'), and this is one of the two GAP files the paper cites. Adds the date line. In the title line inside this block, 'non-triviality' closes up to 'nontriviality' (prefix rule).

Old:
```
##  paper Section 4 (and the accompanying non-triviality remark).
##
##  Mirrors python_project/src/symbolic_proof_checks.py.
```

New:
```
##  paper Section 4 (and the accompanying nontriviality remark).
##
##  Created: 2026-04-28.  Updated: 2026-10-24.
##
##  The paper cites this file in Appendix C (`sec:methodology`), paragraph
##  "Symbolic proof and numeric validation", with
##  python_project/src/symbolic_proof_checks.py and
##  python_project/src/verify_structural_proof.py: the four lemmas
##  verified on Z-bases in a Python and a GAP implementation.
##
##  Mirrors python_project/src/symbolic_proof_checks.py.
```
