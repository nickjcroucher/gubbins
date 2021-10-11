# ST239 dataset
The ST239 dataset (7 Mbyte multi-FASTA alignment file) which was used in the paper was downloaded from:
```
ftp://ftp.sanger.ac.uk/pub/project/pathogens/gubbins/ST239.aln.gz
```
and was run through gubbins with all of the default parameters:
```
gunzip ST239.aln.gz
run_gubbins.py --prefix ST239 --first-tree-builder rapidnj --first-model JC --tree-builder raxmlng --model GTR ST239.aln
```
This directory contains the output files generated. The analysis took ~30 seconds to run.
