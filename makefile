gubbins: seqUtil.o Newickform.o alignment_file.o main.o parse_vcf.o snp_searching.o vcf.o branch_sequences.o gff_file.o gubbins.o parse_phylip.o phylib_of_snp_sites.o snp_sites.o block_tab_file.o fasta_of_snp_sites.o tree_statistics.o
	gcc -o gubbins seqUtil.o Newickform.o alignment_file.o block_tab_file.o main.o parse_vcf.o gff_file.o snp_searching.c vcf.o branch_sequences.o gubbins.o parse_phylip.o tree_statistics.o phylib_of_snp_sites.o snp_sites.o fasta_of_snp_sites.o -lm -lz

seqUtil.o: seqUtil.c
	gcc -c -g seqUtil.c
Newickform.o: Newickform.c
	gcc -c -g Newickform.c

alignment_file.o: alignment_file.c
	gcc -c -g alignment_file.c

block_tab_file.o: block_tab_file.c
	gcc -c -g block_tab_file.c

branch_sequences.o: branch_sequences.c
	gcc -c -g branch_sequences.c

fasta_of_snp_sites.o: 	fasta_of_snp_sites.c
	gcc -c -g fasta_of_snp_sites.c

gff_file.o: gff_file.c
	gcc -c -g gff_file.c

gubbins.o: gubbins.c
	gcc -c -g gubbins.c

main.o: main.c
	gcc -c -g main.c

parse_phylip.o: parse_phylip.c
	gcc -c -g parse_phylip.c

parse_vcf.o: 	parse_vcf.c
	gcc -c -g parse_vcf.c

phylib_of_snp_sites.o: 	phylib_of_snp_sites.c
	gcc -c -g phylib_of_snp_sites.c

snp_searching.o: snp_searching.c
	gcc -c -g snp_searching.c

snp_sites.o: snp_sites.c
	gcc -c -g snp_sites.c

tree_statistics.o: tree_statistics.c
	gcc -c -g tree_statistics.c

vcf.o: vcf.c
	gcc -c -g vcf.c

clean:
	-rm *.o

test: clean gubbins
	cd tests && make
