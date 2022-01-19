# ST239 dataset
The ST239 dataset (11.3 Mb multi-FASTA alignment file) which was used in the paper was downloaded from:
```
https://figshare.com/ndownloader/files/33472622
```
and was run through gubbins with all of the default parameters:
```
gunzip ST239.aln.gz
run_gubbins.py --prefix ST239 --first-tree-builder rapidnj --first-model JC --tree-builder raxmlng --model GTR ST239.aln
```
This directory contains the output files generated. The analysis took ~30 seconds to run. When rendered in Phandango, the output resembles:

![ST239 analysis](ST239_output.png)
