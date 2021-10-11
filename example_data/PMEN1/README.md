# PMEN1 dataset
The PMEN1 dataset (12 Mbyte multi-FASTA alignment file) which was used in the paper was downloaded from:
```
ftp://ftp.sanger.ac.uk/pub/project/pathogens/gubbins/PMEN1.aln.gz
```
and was run through gubbins with all of the default parameters:
```
gunzip PMEN1.aln.gz
run_gubbins.py --prefix PMEN1 --first-tree-builder rapidnj --first-model JC --tree-builder raxmlng --model GTR PMEN1.aln
```
This directory contains the output files generated. The analysis took ~20 seconds to run.
