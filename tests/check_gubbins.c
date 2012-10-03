#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_parse_phylip.h"
#include "helper_methods.h"
#include "gubbins.h"

START_TEST (check_gubbins_no_recombinations)
{
	remove("../tests/data/no_recombinations.tre");
	cp("../tests/data/no_recombinations.tre", "../tests/data/no_recombinations.original.tre");
	run_gubbins("../tests/data/no_recombinations.aln.vcf", "../tests/data/no_recombinations.tre","../tests/data/no_recombinations.aln.snp_sites.aln",3);
	fail_unless(file_exists("../tests/data/no_recombinations.tre.tab") == 1);
	fail_unless(file_exists("../tests/data/no_recombinations.tre.vcf") == 1);
	fail_unless(file_exists("../tests/data/no_recombinations.tre.phylip") == 1);
	fail_unless(file_exists("../tests/data/no_recombinations.tre.snp_sites.aln") == 1);
	fail_unless(file_exists("../tests/data/no_recombinations.tre.stats") == 1);
	fail_unless(file_exists("../tests/data/no_recombinations.tre.gff") == 1);
	
  fail_unless(number_of_recombinations_in_file("../tests/data/no_recombinations.tre.tab") == 0);
  fail_unless(compare_files("../tests/data/no_recombinations.tre","../tests/data/no_recombinations.expected.tre") == 1);

  remove("../tests/data/no_recombinations.tre");
	remove("../tests/data/no_recombinations.tre.tab");
	remove("../tests/data/no_recombinations.tre.vcf");
	remove("../tests/data/no_recombinations.tre.phylip");
	remove("../tests/data/no_recombinations.tre.stats");
	remove("../tests/data/no_recombinations.tre.gff");
	remove("../tests/data/no_recombinations.tre.snp_sites.aln");
}
END_TEST

START_TEST (check_gubbins_one_recombination)
{
	remove("../tests/data/one_recombination.tre");
	cp("../tests/data/one_recombination.tre", "../tests/data/one_recombination.original.tre");
	run_gubbins("../tests/data/one_recombination.aln.vcf", "../tests/data/one_recombination.tre","../tests/data/one_recombination.aln.snp_sites.aln",3);
	fail_unless(file_exists("../tests/data/one_recombination.tre.tab") == 1);
	fail_unless(file_exists("../tests/data/one_recombination.tre.vcf") == 1);
	fail_unless(file_exists("../tests/data/one_recombination.tre.phylip") == 1);
	fail_unless(file_exists("../tests/data/one_recombination.tre.stats") == 1);
	fail_unless(file_exists("../tests/data/one_recombination.tre.gff") == 1);
	fail_unless(file_exists("../tests/data/one_recombination.tre.snp_sites.aln") == 1);
	
	fail_unless(number_of_recombinations_in_file("../tests/data/one_recombination.tre.tab") == 1);
	fail_unless(compare_files("../tests/data/one_recombination.tre.vcf","../tests/data/one_recombination.expected.vcf") == 1);
	fail_unless(compare_files("../tests/data/one_recombination.tre.stats","../tests/data/one_recombination.expected.stats") == 1);
  
  remove("../tests/data/one_recombination.tre");
	remove("../tests/data/one_recombination.tre.tab");
	remove("../tests/data/one_recombination.tre.vcf");
	remove("../tests/data/one_recombination.tre.phylip");
	remove("../tests/data/one_recombination.tre.stats");
	remove("../tests/data/one_recombination.tre.gff");
	remove("../tests/data/one_recombination.tre.snp_sites.aln");
}
END_TEST

START_TEST (check_gubbins_multiple_recombinations)
{
	remove("../tests/data/multiple_recombinations.tre");
	cp("../tests/data/multiple_recombinations.tre", "../tests/data/multiple_recombinations.original.tre");
	run_gubbins("../tests/data/multiple_recombinations.aln.vcf", "../tests/data/multiple_recombinations.tre","../tests/data/multiple_recombinations.aln.snp_sites.aln",3);
	fail_unless(file_exists("../tests/data/multiple_recombinations.tre.tab") == 1);
	fail_unless(file_exists("../tests/data/multiple_recombinations.tre.vcf") == 1);
	fail_unless(file_exists("../tests/data/multiple_recombinations.tre.phylip") == 1);
	fail_unless(file_exists("../tests/data/multiple_recombinations.tre.stats") == 1);
	fail_unless(file_exists("../tests/data/multiple_recombinations.tre.gff") == 1);
	fail_unless(file_exists("../tests/data/multiple_recombinations.tre.snp_sites.aln") == 1);

	fail_unless(number_of_recombinations_in_file("../tests/data/multiple_recombinations.tre.tab") == 3);
  fail_unless(compare_files("../tests/data/multiple_recombinations.tre","../tests/data/multiple_recombinations.expected.tre") == 1);

	remove("../tests/data/multiple_recombinations.tre");
	remove("../tests/data/multiple_recombinations.tre.tab");
	remove("../tests/data/multiple_recombinations.tre.vcf");
	remove("../tests/data/multiple_recombinations.tre.phylip");
	remove("../tests/data/multiple_recombinations.tre.stats");
	remove("../tests/data/multiple_recombinations.tre.gff");
	remove("../tests/data/multiple_recombinations.tre.snp_sites.aln");
}
END_TEST



Suite * run_gubbins_suite(void)
{
  Suite *s = suite_create ("Checking the gubbins functionality");
  TCase *tc_gubbins = tcase_create ("check_gubbins_produces_files");
	tcase_add_test (tc_gubbins, check_gubbins_no_recombinations);
	tcase_add_test (tc_gubbins, check_gubbins_one_recombination);
	tcase_add_test (tc_gubbins, check_gubbins_multiple_recombinations);
  suite_add_tcase (s, tc_gubbins);
  return s;
}



