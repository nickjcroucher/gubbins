
# ship these headers in the "make dist" target
dist_data_DATA = Newickform.h alignment_file.h block_tab_file.h branch_sequences.h fasta_of_snp_sites.h gff_file.h gubbins.h kseq.h parse_phylip.h parse_vcf.h phylib_of_snp_sites.h seqUtil.h snp_searching.h snp_sites.h tree_statistics.h vcf.h

# gubbins is our top level progra
bin_PROGRAMS = gubbins
gubbins_SOURCES = main.c
gubbins_LDADD=-lz -lm libgubbins.la

# libgubbins.so is our library
lib_LTLIBRARIES = libgubbins.la libPyGubbins.la
libgubbins_la_SOURCES = Newickform.c alignment_file.c block_tab_file.c branch_sequences.c fasta_of_snp_sites.c gff_file.c gubbins.c parse_phylip.c parse_vcf.c phylib_of_snp_sites.c seqUtil.c snp_searching.c snp_sites.c tree_statistics.c vcf.c
libgubbins_la_LDFLAGS= -version-info 0:1

libPyGubbins_la_SOURCES = PyGubbins.cpp
libPyGubbins_la_LDFLAGS = -version-info 0:1
libPyGubbins_la_CPPFLAGS = -I/usr/include/python2.7 # yeah, there are better ways to check (and to particularly not hardcode the py version)
libPyGubbins_la_LIBADD = -lpython2.7 -lboost_python -lz -lm libgubbins.la


# "make check" target
TESTS = $(check_PROGRAMS)
check_PROGRAMS = run_all_tests
run_all_tests_SOURCES = \
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
run_all_tests_LDADD = -lm -lz -lcheck libgubbins.la