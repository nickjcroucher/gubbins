test:
	#export DYLD_LIBRARY_PATH=${CONDA_PREFIX}/lib/
	gcc -c -g *.c  -I . -I ../src ${CFLAGS}
	gcc -o run_all_tests *.o ../src/seqUtil.o ../src/string_cat.o ../src/Newickform.o ../src/alignment_file.o ../src/block_tab_file.o ../src/gff_file.o ../src/parse_vcf.o ../src/snp_searching.o ../src/vcf.o ../src/branch_sequences.o ../src/gubbins.o ../src/parse_phylip.o ../src/phylip_of_snp_sites.o ../src/snp_sites.o ../src/tree_statistics.o ../src/fasta_of_snp_sites.o ../src/tree_scaling.o -lm -lz ${LDFLAGS} --coverage -lcheck -I . -fprofile-arcs -ftest-coverage
	./run_all_tests
	gcov ./run_all_tests
	lcov -c  -d . -o gubbins_c_cov.info
