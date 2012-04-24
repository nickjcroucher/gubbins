test:
	gcc -c -g *.c  -I .
	gcc -o run_all_tests *.o ../src/seqUtil.o ../src/Newickform.o ../src/alignment_file.o ../src/block_tab_file.o ../src/gff_file.o ../src/parse_vcf.o ../src/snp_searching.c ../src/vcf.o ../src/branch_sequences.o ../src/gubbins.o ../src/parse_phylip.o ../src/phylib_of_snp_sites.o ../src/snp_sites.o ../src/tree_statistics.o ../src/fasta_of_snp_sites.o -lm -lz -L /usr/local/lib/libcheck.dylib
	./run_all_tests
