# TCR input for exercise 10 (NetTCR)

`tcrs_example.csv` is the missing input that NetTCR needs and that a protein
sequence can never provide.

## Why this file has to exist

Every other exercise in this module starts from one thing, the fhaC sequence.
NetTCR is the exception. It scores a **pair**, a peptide against a T cell
receptor, so it needs a second input that is not in the FASTA:

| Input | Comes from |
|---|---|
| peptide | the fhaC sequence, tiled into 9-mers |
| CDR3 alpha and CDR3 beta | TCR sequencing of actual T cells |

There is no way to derive a TCR from an antigen. If nobody has sequenced T
cells against your protein, you do not have this input, and that is a real
experimental gap rather than a gap in the pipeline.

## What is in the example file

Twenty real, published CDR3 pairs from the NetTCR-2.0 training data, two for
each of the ten best covered peptides. Nothing is invented.

| Column | Meaning |
|---|---|
| `name` | a label for the row |
| `CDR3a` | CDR3 of the alpha chain |
| `CDR3b` | CDR3 of the beta chain |
| `measured_against` | the peptide this TCR was actually measured against |
| `source` | where the sequence comes from |

The `measured_against` column is there on purpose. It is the reminder that
none of these TCRs has anything to do with fhaC. They are in the file so you
can see the format and run the tooling end to end, not because pairing them
with fhaC peptides would mean anything.

## Using your own TCRs

Replace the file with your own, keeping at least `name`, `CDR3a` and `CDR3b`:

```csv
name,CDR3a,CDR3b
donor1_clone3,AVRPLLDGTYIPT,ASSYSGTGSYEQY
```

Then:

```bash
python solutions/ex10_nettcr_prepare.py data/tcrs_example.csv 20
```

That crosses your TCRs with the top PRIME candidates and writes the input file
NetTCR expects. Run NetTCR itself per its own documentation.

## Read exercise 13 before you trust any of this

Exercise 13 shows, with numbers, that NetTCR-2.0 only works on the 17 peptides
in its training data. On fhaC peptides its score moves with the CDR3 and barely
moves with the peptide, which means it is not answering the question being
asked. For a new antigen the honest output is "no evidence", which is exactly
what the prioritisation exercise records.

To fill that column for real you would need to generate the data: sort T cells
with peptide MHC tetramers built from your own candidates, sequence their
receptors, and only then is there anything for NetTCR to score.
