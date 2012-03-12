
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

//merge_adjacent_blocks(int ** block_coordinates, int number_of_blocks)
START_TEST (check_merge_adjacent_blocks_not_adjacent)
{
	int * blocks_not_adjacent[2]; 
	blocks_not_adjacent[0]   = (int *)  malloc((3)*sizeof(int)); 
	blocks_not_adjacent[1]   = (int *)  malloc((3)*sizeof(int)); 
	blocks_not_adjacent[0][0] = 10; 
	blocks_not_adjacent[1][0] = 20; 
	blocks_not_adjacent[0][1] = 1000; 
	blocks_not_adjacent[1][1] = 1200; 
	
	fail_unless(merge_adjacent_blocks(blocks_not_adjacent,2) == 2);
	fail_unless(blocks_not_adjacent[0][0] == 10  );
	fail_unless(blocks_not_adjacent[1][0] == 20  );
	fail_unless(blocks_not_adjacent[0][1] == 1000); 
	fail_unless(blocks_not_adjacent[1][1] == 1200);
}
END_TEST

START_TEST (check_merge_adjacent_blocks_beside_each_other)
{

	int * blocks_beside_each_other[2]; 
	blocks_beside_each_other[0]   = (int *)  malloc((3)*sizeof(int)); 
	blocks_beside_each_other[1]   = (int *)  malloc((3)*sizeof(int)); 
	blocks_beside_each_other[0][0] = 10; 
	blocks_beside_each_other[1][0] = 20; 
	blocks_beside_each_other[0][1] = 20; 
	blocks_beside_each_other[1][1] = 30;
	
	fail_unless(merge_adjacent_blocks(blocks_beside_each_other,2) == 1);
	fail_unless(blocks_beside_each_other[0][0] == 10  );
	fail_unless(blocks_beside_each_other[1][0] == 30  );
	fail_unless(blocks_beside_each_other[0][1] == 0); 
	fail_unless(blocks_beside_each_other[1][1] == 0);
}
END_TEST

START_TEST (check_merge_adjacent_blocks_near_each_other)
{	
	int * blocks_near_each_other[2]; 
	blocks_near_each_other[0]   = (int *)  malloc((3)*sizeof(int)); 
	blocks_near_each_other[1]   = (int *)  malloc((3)*sizeof(int)); 
	blocks_near_each_other[0][0] = 10; 
	blocks_near_each_other[1][0] = 20; 
	blocks_near_each_other[0][1] = 21; 
	blocks_near_each_other[1][1] = 30;
	
	fail_unless(merge_adjacent_blocks(blocks_near_each_other,2) == 1);
	fail_unless(blocks_near_each_other[0][0] == 10  );
	fail_unless(blocks_near_each_other[1][0] == 30  );
	fail_unless(blocks_near_each_other[0][1] == 0); 
	fail_unless(blocks_near_each_other[1][1] == 0);
}
END_TEST

START_TEST (check_merge_adjacent_blocks_overlapping)
{
	int * blocks_overlapping[2]; 
	blocks_overlapping[0]   = (int *)  malloc((3)*sizeof(int)); 
	blocks_overlapping[1]   = (int *)  malloc((3)*sizeof(int)); 
	blocks_overlapping[0][0] = 10; 
	blocks_overlapping[1][0] = 20; 
	blocks_overlapping[0][1] = 19; 
	blocks_overlapping[1][1] = 30;
	
	fail_unless(merge_adjacent_blocks(blocks_overlapping,2) == 1);
	fail_unless(blocks_overlapping[0][0] == 10  );
	fail_unless(blocks_overlapping[1][0] == 30  );
	fail_unless(blocks_overlapping[0][1] == 0); 
	fail_unless(blocks_overlapping[1][1] == 0);

	
}
END_TEST
		

Suite * check_branch_sequences_suite (void)
{
  Suite *s = suite_create ("checking branch sequences");

  TCase *tc_branch_sequences = tcase_create ("excluding_recombinations");
  tcase_add_test (tc_branch_sequences, check_exclude_recombination_windows );
	tcase_add_test (tc_branch_sequences, check_merge_adjacent_blocks_not_adjacent);
	tcase_add_test (tc_branch_sequences, check_merge_adjacent_blocks_beside_each_other);
	tcase_add_test (tc_branch_sequences, check_merge_adjacent_blocks_near_each_other);
	tcase_add_test (tc_branch_sequences, check_merge_adjacent_blocks_overlapping);
  suite_add_tcase (s, tc_branch_sequences);

  return s;
}
