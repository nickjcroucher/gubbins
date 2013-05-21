
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_branch_sequences.h"
#include "helper_methods.h"

#include "branch_sequences.h"






START_TEST (check_exclude_snp_sites_in_block)
{
	int number_of_branch_snps = 8;
	int * snp_sites;  
	snp_sites    = (int *)  malloc((number_of_branch_snps+1)*sizeof(int)); 
	snp_sites[0] = 1;
	snp_sites[1] = 3;
	snp_sites[2] = 5;
	snp_sites[3] = 6;
	snp_sites[4] = 7;
	snp_sites[5] = 8;
	snp_sites[6] = 10;
	snp_sites[7] = 11;
	
	fail_unless(exclude_snp_sites_in_block(0,2,  snp_sites, number_of_branch_snps)   == 7);
	fail_unless(exclude_snp_sites_in_block(5,7,  snp_sites, number_of_branch_snps-1) == 4);
	fail_unless(exclude_snp_sites_in_block(8,11, snp_sites, number_of_branch_snps-4) == 1);
	fail_unless(exclude_snp_sites_in_block(3,3,  snp_sites, number_of_branch_snps-7) == 0);
}
END_TEST

START_TEST (check_copy_and_concat_2d_integer_arrays)
{
	int ** block_coords;  
	block_coords  = (int **) malloc(2*sizeof(int*));
	block_coords[0] = (int*) malloc((3)*sizeof(int ));
	block_coords[1] = (int*) malloc((3)*sizeof(int ));
	
	int ** block_coords_2;  
	block_coords_2  = (int **) malloc(2*sizeof(int*));
	block_coords_2[0] = (int*) malloc((2)*sizeof(int ));
	block_coords_2[1] = (int*) malloc((2)*sizeof(int ));
	
	int ** block_coords_out;  
	block_coords_out  = (int **) malloc(2*sizeof(int*));
	block_coords_out[0] = (int*) malloc((5)*sizeof(int ));
	block_coords_out[1] = (int*) malloc((5)*sizeof(int ));
	
	block_coords[0][0] = 5;
	block_coords[1][0] = 10;
	block_coords[0][1] = 100;
	block_coords[1][1] = 110;
	block_coords[0][2] = 7;
	block_coords[1][2] = 15;
	
	block_coords_2[0][0] = 200;
	block_coords_2[1][0] = 204;
	block_coords_2[0][1] = 2;
	block_coords_2[1][1] = 8;
	
	int output_size = 0;
	output_size = copy_and_concat_2d_integer_arrays(block_coords, 3, block_coords_2, 2, block_coords_out) ;
	fail_unless(output_size == 5); 
	fail_unless(block_coords_out[0][0] == 5);
	fail_unless(block_coords_out[1][0] == 10);
	fail_unless(block_coords_out[0][2] == 7);
	fail_unless(block_coords_out[1][2] == 15);
	fail_unless(block_coords_out[0][4] == 2);
	fail_unless(block_coords_out[1][4] == 8);
	
}
END_TEST

int test_bases_in_recombinations(int block_size)
{
	int ** block_coords;  
	block_coords  = (int **) malloc(2*sizeof(int*));
	block_coords[0] = (int*) malloc((4)*sizeof(int ));
	block_coords[1] = (int*) malloc((4)*sizeof(int ));
	block_coords[0][0] = 5;
	block_coords[1][0] = 10;
	block_coords[0][1] = 100;
	block_coords[1][1] = 110;
	block_coords[0][2] = 15;
	block_coords[1][2] = 20;
	block_coords[0][3] = 7;
	block_coords[1][3] = 15;
	char * child_sequence      = "AAAAAAAAATAA";
	int snp_locations[12] = {1,2,3,5,7,10,11,15,20,30,100,110};
	calculate_number_of_bases_in_recombations_excluding_gaps(block_coords, block_size, child_sequence, snp_locations,12);
}

int test_bases_in_recombinations_with_gaps(int block_size)
{
	int ** block_coords;  
	block_coords  = (int **) malloc(2*sizeof(int*));
	block_coords[0] = (int*) malloc((4)*sizeof(int ));
	block_coords[1] = (int*) malloc((4)*sizeof(int ));
	block_coords[0][0] = 5;
	block_coords[1][0] = 10;
	block_coords[0][1] = 100;
	block_coords[1][1] = 110;
	block_coords[0][2] = 15;
	block_coords[1][2] = 20;
	block_coords[0][3] = 7;
	block_coords[1][3] = 15;
	char * child_sequence      =  "--A---AAAAAAAAAAAAAT";
	int snp_locations[16] = {1,4,5,6,7,8,9,10,11,15,20,30,40,50,100,110};
	calculate_number_of_bases_in_recombations_excluding_gaps(block_coords, block_size, child_sequence, snp_locations,16);
}


START_TEST (check_calculate_number_of_bases_in_recombations)
{
	fail_unless(test_bases_in_recombinations(4) == 25);
	fail_unless(test_bases_in_recombinations(3) == 20);
	fail_unless(test_bases_in_recombinations(2) == 15);
	fail_unless(test_bases_in_recombinations(1) == 5);
	
	fail_unless(test_bases_in_recombinations_with_gaps(4) == 22);
	fail_unless(test_bases_in_recombinations_with_gaps(3) == 17);
	fail_unless(test_bases_in_recombinations_with_gaps(2) == 12);
	fail_unless(test_bases_in_recombinations_with_gaps(1) == 2);
}
END_TEST




Suite * check_branch_sequences_suite (void)
{
  Suite *s = suite_create ("checking branch sequences");

  TCase *tc_branch_sequences = tcase_create ("excluding_recombinations");
	tcase_add_test (tc_branch_sequences, check_exclude_snp_sites_in_block);
	tcase_add_test (tc_branch_sequences, check_copy_and_concat_2d_integer_arrays);
	tcase_add_test (tc_branch_sequences, check_calculate_number_of_bases_in_recombations);
  suite_add_tcase (s, tc_branch_sequences);

  return s;
}
