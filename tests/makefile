test:
	gcc -c check_snp_sites.c  -I .
	gcc -o check_snp_sites check_snp_sites.o ../seqUtil.o ../Newickform.o ../alignment_file.o ../block_tab_file.o ../parse_vcf.o ../snp_searching.c ../vcf.o ../branch_sequences.o ../gubbins.o ../parse_phylip.o ../phylib_of_snp_sites.o ../fasta_of_snp_sites.o ../snp_sites.o -lm -lz -L /usr/local/lib/libcheck.dylib
	./check_snp_sites

