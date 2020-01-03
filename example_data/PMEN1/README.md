# PMEN1 dataset
The PMEN1 dataset (12 Mbyte multi-FASTA alignment file) which was used in the paper was downloaded from:
```
ftp://ftp.sanger.ac.uk/pub/project/pathogens/gubbins/PMEN1.aln.gz
```
and was run through gubbins with all of the default parameters:
```
gunzip PMEN1.aln.gz
run_gubbins.py PMEN1.aln
```
This directory contains the output files generated.
It required 60 Mbytes of memory and took 50 seconds to run.
