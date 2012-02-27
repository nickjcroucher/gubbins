#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_parse_phylip.h"
#include "helper_methods.h"
#include "../gubbins.h"

START_TEST (check_gubbins_produces_files)
{
	run_gubbins("data/alignment_file_one_line_per_sequence.aln.vcf", "data/alignment_file_one_line_per_sequence.tre", "data/alignment_file_one_line_per_sequence.aln.phylip","data/alignment_file_one_line_per_sequence.aln.snp_sites.aln");
	fail_unless(file_exists("data/alignment_file_one_line_per_sequence.tre.tab") == 1);
	fail_unless(file_exists("data/alignment_file_one_line_per_sequence.tre.vcf") == 1);
	fail_unless(file_exists("data/alignment_file_one_line_per_sequence.tre.phylip") == 1);
	fail_unless(file_exists("data/alignment_file_one_line_per_sequence.tre.snp_sites.aln") == 1);
	
	remove("data/alignment_file_one_line_per_sequence.tre.tab");
	remove("data/alignment_file_one_line_per_sequence.tre.vcf");
	remove("data/alignment_file_one_line_per_sequence.tre.phylip");
	remove("data/alignment_file_one_line_per_sequence.tre.snp_sites.aln");
}
END_TEST

Suite * run_gubbins_suite(void)
{
  Suite *s = suite_create ("Checking the gubbins functionality");
  TCase *tc_gubbins = tcase_create ("check_gubbins_produces_files");
  tcase_add_test (tc_gubbins, check_gubbins_produces_files);
  suite_add_tcase (s, tc_gubbins);
  return s;
}

