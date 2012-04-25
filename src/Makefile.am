gubbins_no_main= Newickform.c alignment_file.c block_tab_file.c branch_sequences.c fasta_of_snp_sites.c gff_file.c gubbins.c parse_phylip.c parse_vcf.c phylib_of_snp_sites.c seqUtil.c snp_searching.c snp_sites.c tree_statistics.c vcf.c

bin_PROGRAMS = gubbins
gubbins_SOURCES = $(gubbins_no_main) main.c
dist_data_DATA = Newickform.h alignment_file.h block_tab_file.h branch_sequences.h fasta_of_snp_sites.h gff_file.h gubbins.h kseq.h parse_phylip.h parse_vcf.h phylib_of_snp_sites.h seqUtil.h snp_searching.h snp_sites.h tree_statistics.h vcf.h
LDFLAGS=-lz -lm

TESTS = $(check_PROGRAMS)
check_PROGRAMS = run_all_tests
run_all_tests_SOURCES = $(gubbins_no_main) \
	../tests/check_branch_sequences.c \
	../tests/check_calculate_ancestor_sequence.c \
	../tests/check_gubbins.c \
	../tests/check_parse_phylip.c \
	../tests/check_snp_searching.c \
	../tests/check_snp_sites.c \
	../tests/check_vcf_parsing.c \
	../tests/helper_methods.c \
	../tests/run_all_tests.c
run_all_tests_CFLAGS = -I../tests
run_all_tests_LDADD = -lm -lz -lcheck # there's a better way of lookign for libcheck