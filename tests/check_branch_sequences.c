
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_branch_sequences.h"
#include "helper_methods.h"

#include "../branch_sequences.h"


START_TEST (check_exclude_recombination_windows)
{
	int * snp_site_coords;
	int window_start_coordinate = 3;
	int window_end_coordinate = 10;
	int number_of_branch_snps = 5;
	
	snp_site_coords = (int *) malloc((number_of_branch_snps)*sizeof(int));
	snp_site_coords[0] = 1;
	snp_site_coords[1] = 3;
	snp_site_coords[2] = 5;
	snp_site_coords[3] = 10;
	snp_site_coords[4] = 12;
	
	fail_unless(exclude_snp_sites_in_block(window_start_coordinate, window_end_coordinate, snp_site_coords,number_of_branch_snps) == 3);
	fail_unless(snp_site_coords[0] == 1);
	fail_unless(snp_site_coords[1] == 10);
	fail_unless(snp_site_coords[2] == 12);
}
END_TEST

Suite * check_branch_sequences_suite (void)
{
  Suite *s = suite_create ("checking branch sequences");

  TCase *tc_branch_sequences = tcase_create ("excluding_recombinations");
  tcase_add_test (tc_branch_sequences, check_exclude_recombination_windows );
  suite_add_tcase (s, tc_branch_sequences);

  return s;
}
