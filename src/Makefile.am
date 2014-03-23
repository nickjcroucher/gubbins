
# ship these headers in the "make dist" target
dist_data_DATA = Newickform.h alignment_file.h block_tab_file.h string_cat.h branch_sequences.h fasta_of_snp_sites.h gff_file.h gubbins.h kseq.h parse_phylip.h parse_vcf.h phylip_of_snp_sites.h seqUtil.h snp_searching.h snp_sites.h tree_statistics.h tree_scaling.h vcf.h

# gubbins is our top level progra
bin_PROGRAMS = gubbins
gubbins_SOURCES = main.c
gubbins_LDADD=libgubbins.la -lz -lm

# libgubbins.so is our library
lib_LTLIBRARIES = libgubbins.la
libgubbins_la_SOURCES = Newickform.c alignment_file.c block_tab_file.c branch_sequences.c string_cat.c fasta_of_snp_sites.c gff_file.c gubbins.c parse_phylip.c parse_vcf.c phylip_of_snp_sites.c seqUtil.c snp_searching.c snp_sites.c tree_statistics.c tree_scaling.c vcf.c
libgubbins_la_LDFLAGS= -version-info 0:1
libgubbins_la_LIBADD = -lm

#libPyGubbins_la_SOURCES = PyGubbins.cpp
#libPyGubbins_la_LDFLAGS = -version-info 0:1
#libPyGubbins_la_CPPFLAGS = $(PYTHON_CPPFLAGS)
#libPyGubbins_la_LIBADD = $(PYTHON_LDFLAGS) -lboost_python libgubbins.la -lz -lm


# "make check" target
TESTS = $(check_PROGRAMS)
check_PROGRAMS = run_all_tests
run_all_tests_SOURCES = \
	../tests/check_branch_sequences.c \
	../tests/check_gubbins.c \
	../tests/check_parse_phylip.c \
	../tests/check_snp_searching.c \
	../tests/check_snp_sites.c \
	../tests/check_vcf_parsing.c \
	../tests/helper_methods.c \
	../tests/run_all_tests.c
run_all_tests_CFLAGS = -I../tests
run_all_tests_LDADD = -lcheck libgubbins.la -lz -lm
