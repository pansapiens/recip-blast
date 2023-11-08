# Recip_BLAST

A quick and dirty reciprocal BLAST to find orthologues between two (bacterial?) strains, made by guiding ChatGPT4.

You are probably better off using something like [Roary](https://github.com/sanger-pathogens/Roary) or [Panaroo](https://github.com/gtonkinhill/panaroo).

> NOT WELL TESTED, BUT APPEARS TO WORK ?

*Requires:* NCBI BLAST+, `makeblastdb`, `blastp` and/or `blastn`

*Input:* Two FASTA files of protein or nucleotide sequences.

*Output:*
 
  - A table `reciprocal_best_hits.tsv` with pairs of orthologues common between the two strains
  - a bunch of intermediate file from running BLAST

# Quickstart

Example:
```bash
# sudo apt install ncbi-blast+

cd example

../recip_blast.py --threads=4 --identity=90 GCF_000827025.1.faa GCF_017821535.1.faa
```

Output files go to the directory `output` by default. 
Look at the file `reciprocal_best_hits.tsv` for pairs of orthologues common between the two strains.

# Tests

You can run some simple unit tests like:
```bash
python tests.py
```

_These should be expanded so we can be more confident of the result_
