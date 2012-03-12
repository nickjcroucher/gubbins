#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_parse_phylip.h"
#include "helper_methods.h"
#include "../gubbins.h"

START_TEST (check_gubbins_no_recombinations)
{
	run_gubbins("data/no_recombinations.aln.vcf", "data/no_recombinations.tre", "data/no_recombinations.aln.phylip","data/no_recombinations.aln.snp_sites.aln");
	fail_unless(file_exists("data/no_recombinations.tre.tab") == 1);
	fail_unless(file_exists("data/no_recombinations.tre.vcf") == 1);
	fail_unless(file_exists("data/no_recombinations.tre.phylip") == 1);
	fail_unless(file_exists("data/no_recombinations.tre.snp_sites.aln") == 1);
	
  fail_unless(number_of_recombinations_in_file("data/no_recombinations.tre.tab") == 0);
  
	remove("data/no_recombinations.tre.tab");
	remove("data/no_recombinations.tre.vcf");
	remove("data/no_recombinations.tre.phylip");
	remove("data/no_recombinations.tre.snp_sites.aln");
}
END_TEST

START_TEST (check_gubbins_one_recombination)
{
	run_gubbins("data/one_recombination.aln.vcf", "data/one_recombination.tre", "data/one_recombination.aln.phylip","data/one_recombination.aln.snp_sites.aln");
	fail_unless(file_exists("data/one_recombination.tre.tab") == 1);
	fail_unless(file_exists("data/one_recombination.tre.vcf") == 1);
	fail_unless(file_exists("data/one_recombination.tre.phylip") == 1);
	fail_unless(file_exists("data/one_recombination.tre.snp_sites.aln") == 1);
	
	fail_unless(number_of_recombinations_in_file("data/one_recombination.tre.tab") == 1);

	remove("data/one_recombination.tre.tab");
	remove("data/one_recombination.tre.vcf");
	remove("data/one_recombination.tre.phylip");
	remove("data/one_recombination.tre.snp_sites.aln");
}
END_TEST

START_TEST (check_gubbins_multiple_recombinations)
{
	run_gubbins("data/multiple_recombinations.aln.vcf", "data/multiple_recombinations.tre", "data/multiple_recombinations.aln.phylip","data/multiple_recombinations.aln.snp_sites.aln");
	fail_unless(file_exists("data/multiple_recombinations.tre.tab") == 1);
	fail_unless(file_exists("data/multiple_recombinations.tre.vcf") == 1);
	fail_unless(file_exists("data/multiple_recombinations.tre.phylip") == 1);
	fail_unless(file_exists("data/multiple_recombinations.tre.snp_sites.aln") == 1);

	fail_unless(number_of_recombinations_in_file("data/multiple_recombinations.tre.tab") == 4);

	remove("data/multiple_recombinations.tre.tab");
	remove("data/multiple_recombinations.tre.vcf");
	remove("data/multiple_recombinations.tre.phylip");
	remove("data/multiple_recombinations.tre.snp_sites.aln");
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



