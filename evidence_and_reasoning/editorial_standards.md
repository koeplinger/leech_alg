# Editorial standards

*The five principles below govern every prose
change to every *durable artifact* in this repository: the main
paper and any other documents under `paper/` (update briefs, companion,
historical appendix, exposition, verification and review notes);
governance and methodology files (`MANIFESTO.md`, this file,
`evidence_and_reasoning/terminology.md`, `README.md`); plan documents and
research statements; and the prose parts of verification scripts and
their accompanying notes (docstrings, headers, comment blocks of any
length). Ephemeral working notes (e.g. /tmp scratch) are out of scope.
Prompt logs are durable but forward-evolving (MANIFESTO §12): the
standards apply to newly written logs but are not retroactively enforced
on past ones.*

The intended reader of any of these artifacts is an expert in the field,
or a future investigator (human or AI) who must reconstruct what was
done and why. Their time is to be respected: every sentence should earn
its place, and the logic should be followable on a first reading.

## The five standards

1. **Lead programmatically, then develop.** State the finding or claim
   first; build the apparatus that supports it afterwards. A paragraph
   announces what it is about before it elaborates.

2. **Introduce one thing at a time, in dependency order.** Never use a
   term before it is defined. Each sentence should bring in at most one
   new object, and only objects the reader already has. Do not break the
   thread to visit a side topic.

3. **Cut every redundancy.** Say each thing once, in its proper place.
   If a fact is stated where it belongs, do not restate it elsewhere.

4. **One consistent vocabulary.** Fixed conventions, used the same way
   throughout. For example: *division algebra* = no zero divisors;
   *third-/fourth-/power-associativity* kept distinct; `--` (en-dash)
   never `---`.  Dash-as-punctuation (spaced ` -- ` asides) is avoided
   altogether: split the sentence in two, use a
   colon or parentheses, or move the aside to a footnote if it is out
   of the main flow.  Unspaced en-dashes remain for name joins
   (Baez--Egan) and numeric ranges (1923--1946).
   Spelling is US English: *-ize*/*-ization*
   (characterize, polarization, stabilizer), *center*, *labeled*.
   Quoted titles and verbatim quotations keep their source spelling.

5. **Precision over good-sounding vagueness.** Define terms; do not
   over-claim; state plainly what is open. A sentence that sounds
   substantial but asserts nothing checkable is to be cut or made
   concrete.

## Corollaries

- Keep at most one or two forward references outside the introduction of
  any document.
- An introduction (where one exists) should read as a compressed version
  of the document.

## How to apply

Check every prose change against all five. A practical test for
standard 2: read the paragraph as a first-time reader and, at each
sentence, ask: *do I have every term used here, and is this sentence
still on the topic the paragraph announced?* If the answer is no, the
logical order is wrong, or a side topic has intruded.

A practical test for standard 5: for each sentence, ask what would
change for the reader if it were deleted. A sentence whose deletion
costs the reader nothing checkable is vague; cut it or make it
concrete.
