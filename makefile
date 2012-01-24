gubbins: seqUtil.o Newickform.o alignment_file.o main.o parse_vcf.o snp_searching.c vcf.o branch_sequences.o gubbins.o parse_phylip.o phylib_of_snp_sites.o snp_sites.o block_tab_file.o fasta_of_snp_sites.o
	gcc -o gubbins seqUtil.o Newickform.o alignment_file.o block_tab_file.o main.o parse_vcf.o snp_searching.c vcf.o branch_sequences.o gubbins.o parse_phylip.o phylib_of_snp_sites.o snp_sites.o fasta_of_snp_sites.o -lm -lz

seqUtil.o: seqUtil.c
	gcc -c seqUtil.c
Newickform.o: Newickform.c
	gcc -c Newickform.c

alignment_file.o: alignment_file.c
	gcc -c alignment_file.c

block_tab_file.o: block_tab_file.c
	gcc -c block_tab_file.c

branch_sequences.o: branch_sequences.c
	gcc -c branch_sequences.c

fasta_of_snp_sites.o: 	fasta_of_snp_sites.c
	gcc -c fasta_of_snp_sites.c

gubbins.o: gubbins.c
	gcc -c gubbins.c

main.o: main.c
	gcc -c main.c

parse_phylip.o: parse_phylip.c
	gcc -c parse_phylip.c

parse_vcf.o: 	parse_vcf.c
	gcc -c parse_vcf.c

phylib_of_snp_sites.o: 	phylib_of_snp_sites.c
	gcc -c phylib_of_snp_sites.c

snp_searching.o: snp_searching.c
	gcc -c snp_searching.c

snp_sites.o: snp_sites.c
	gcc -c snp_sites.c

vcf.o: vcf.c
	gcc -c vcf.c

clean:
	-rm *.o

test:
	cd tests && make
