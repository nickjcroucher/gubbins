#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_snp_searching.h"
#include "helper_methods.h"

#include "snp_searching.h"


START_TEST (check_find_starting_index_windows)
{
	int coords_empty[0] = {};
	int coords_one[1] = {1};
	int coords_odd[3] = {1,3,5};
	int coords_even[4] = {1,3,5,7};
	
	fail_unless( find_starting_index(0, coords_empty, 0, 0) == 0);
  fail_unless( find_starting_index(0, coords_one, 0, 1) == 0);
	fail_unless( find_starting_index(1, coords_one, 0, 1) == 0);
	
  fail_unless( find_starting_index(1, coords_odd, 0, 3) == 0);
  fail_unless( find_starting_index(3, coords_odd, 0, 3) == 1);
  fail_unless( find_starting_index(5, coords_odd, 0, 3) == 2);
  
  fail_unless( find_starting_index(0, coords_odd, 0, 3) == 0);
  fail_unless( find_starting_index(2, coords_odd, 0, 3) == 0);
  fail_unless( find_starting_index(4, coords_odd, 0, 3) == 1);
  
  fail_unless( find_starting_index(1, coords_even, 0, 4) == 0);
  fail_unless( find_starting_index(3, coords_even, 0, 4) == 1);
  fail_unless( find_starting_index(5, coords_even, 0, 4) == 2);
  fail_unless( find_starting_index(7, coords_even, 0, 4) == 3);
  fail_unless( find_starting_index(9, coords_even, 0, 4) == 3);
  
  fail_unless( find_starting_index(0, coords_even, 0, 4) == 0);
  fail_unless( find_starting_index(2, coords_even, 0, 4) == 0);
  fail_unless( find_starting_index(4, coords_even, 0, 4) == 1);
  fail_unless( find_starting_index(6, coords_even, 0, 4) == 2);
  fail_unless( find_starting_index(8, coords_even, 0, 4) == 3);
	
}
END_TEST

START_TEST (check_advance_window_start_to_next_snp)
{
  int coords_empty[0] = {};
  int coords_one[1]   = {1};
  int coords_odd[3]   = {1,3,5};
  int coords_even[4]  = {1,3,5,7};
	char *child_sequence_without_gaps = "ACGT";
	char *child_sequence_with_gaps = "-AG-";
	
	// Without gaps
	fail_unless( advance_window_start_to_next_snp(0, coords_empty,child_sequence_without_gaps,  0) == 0);
  fail_unless( advance_window_start_to_next_snp(0, coords_one,  child_sequence_without_gaps,  1) == 1);
	fail_unless( advance_window_start_to_next_snp(1, coords_one,  child_sequence_without_gaps,  1) == 1);
  fail_unless( advance_window_start_to_next_snp(1, coords_odd,  child_sequence_without_gaps,  3) == 1);
  fail_unless( advance_window_start_to_next_snp(3, coords_odd,  child_sequence_without_gaps,  3) == 3);
  fail_unless( advance_window_start_to_next_snp(5, coords_odd,  child_sequence_without_gaps,  3) == 5);
  fail_unless( advance_window_start_to_next_snp(0, coords_odd,  child_sequence_without_gaps,  3) == 1);
  fail_unless( advance_window_start_to_next_snp(2, coords_odd,  child_sequence_without_gaps,  3) == 3);
  fail_unless( advance_window_start_to_next_snp(4, coords_odd,  child_sequence_without_gaps,  3) == 5);
  fail_unless( advance_window_start_to_next_snp(1, coords_even, child_sequence_without_gaps,  4) == 1);
  fail_unless( advance_window_start_to_next_snp(3, coords_even, child_sequence_without_gaps,  4) == 3);
  fail_unless( advance_window_start_to_next_snp(5, coords_even, child_sequence_without_gaps,  4) == 5);
  fail_unless( advance_window_start_to_next_snp(7, coords_even, child_sequence_without_gaps,  4) == 7);
  fail_unless( advance_window_start_to_next_snp(9, coords_even, child_sequence_without_gaps,  4) == 9);
  fail_unless( advance_window_start_to_next_snp(0, coords_even, child_sequence_without_gaps,  4) == 1);
  fail_unless( advance_window_start_to_next_snp(2, coords_even, child_sequence_without_gaps,  4) == 3);
  fail_unless( advance_window_start_to_next_snp(4, coords_even, child_sequence_without_gaps,  4) == 5);
  fail_unless( advance_window_start_to_next_snp(6, coords_even, child_sequence_without_gaps,  4) == 7);
  fail_unless( advance_window_start_to_next_snp(8, coords_even, child_sequence_without_gaps,  4) == 8);

  // With gaps
	fail_unless( advance_window_start_to_next_snp(0, coords_empty,child_sequence_with_gaps,  0) == 0);
  fail_unless( advance_window_start_to_next_snp(0, coords_one,  child_sequence_with_gaps,  1) == 0);
	fail_unless( advance_window_start_to_next_snp(1, coords_one,  child_sequence_with_gaps,  1) == 1);
  fail_unless( advance_window_start_to_next_snp(1, coords_odd,  child_sequence_with_gaps,  3) == 3);
  fail_unless( advance_window_start_to_next_snp(3, coords_odd,  child_sequence_with_gaps,  3) == 3);
  fail_unless( advance_window_start_to_next_snp(5, coords_odd,  child_sequence_with_gaps,  3) == 5);
  fail_unless( advance_window_start_to_next_snp(0, coords_odd,  child_sequence_with_gaps,  3) == 3);
  fail_unless( advance_window_start_to_next_snp(2, coords_odd,  child_sequence_with_gaps,  3) == 3);
  fail_unless( advance_window_start_to_next_snp(4, coords_odd,  child_sequence_with_gaps,  3) == 5);
  fail_unless( advance_window_start_to_next_snp(1, coords_even, child_sequence_with_gaps,  4) == 3);
  fail_unless( advance_window_start_to_next_snp(3, coords_even, child_sequence_with_gaps,  4) == 3);
  fail_unless( advance_window_start_to_next_snp(5, coords_even, child_sequence_with_gaps,  4) == 5);
  fail_unless( advance_window_start_to_next_snp(7, coords_even, child_sequence_with_gaps,  4) == 7);
  fail_unless( advance_window_start_to_next_snp(9, coords_even, child_sequence_with_gaps,  4) == 9);
  fail_unless( advance_window_start_to_next_snp(0, coords_even, child_sequence_with_gaps,  4) == 3);
  fail_unless( advance_window_start_to_next_snp(2, coords_even, child_sequence_with_gaps,  4) == 3);
  fail_unless( advance_window_start_to_next_snp(4, coords_even, child_sequence_with_gaps,  4) == 5);
  fail_unless( advance_window_start_to_next_snp(6, coords_even, child_sequence_with_gaps,  4) == 6);
  fail_unless( advance_window_start_to_next_snp(8, coords_even, child_sequence_with_gaps,  4) == 8);
}
END_TEST

START_TEST (check_rewind_window_end_to_last_snp)
{
  int coords_empty[0] = {};
  int coords_one[1]   = {1};
  int coords_odd[3]   = {1,3,5};
  int coords_even[4]  = {1,3,5,7};
	char *child_sequence_without_gaps = "ACGT";
	char *child_sequence_with_gaps = "-AG-";
	
	// Without gaps
	fail_unless( rewind_window_end_to_last_snp(0, coords_empty,child_sequence_without_gaps,  0) == 0);
  fail_unless( rewind_window_end_to_last_snp(0, coords_one,  child_sequence_without_gaps,  1) == 0);
	fail_unless( rewind_window_end_to_last_snp(1, coords_one,  child_sequence_without_gaps,  1) == 1);
  fail_unless( rewind_window_end_to_last_snp(1, coords_odd,  child_sequence_without_gaps,  3) == 1);
  fail_unless( rewind_window_end_to_last_snp(3, coords_odd,  child_sequence_without_gaps,  3) == 3);
  fail_unless( rewind_window_end_to_last_snp(5, coords_odd,  child_sequence_without_gaps,  3) == 5);
  fail_unless( rewind_window_end_to_last_snp(0, coords_odd,  child_sequence_without_gaps,  3) == 0);
  fail_unless( rewind_window_end_to_last_snp(2, coords_odd,  child_sequence_without_gaps,  3) == 1);
  fail_unless( rewind_window_end_to_last_snp(4, coords_odd,  child_sequence_without_gaps,  3) == 3);
  fail_unless( rewind_window_end_to_last_snp(1, coords_even, child_sequence_without_gaps,  4) == 1);
  fail_unless( rewind_window_end_to_last_snp(3, coords_even, child_sequence_without_gaps,  4) == 3);
  fail_unless( rewind_window_end_to_last_snp(5, coords_even, child_sequence_without_gaps,  4) == 5);
  fail_unless( rewind_window_end_to_last_snp(7, coords_even, child_sequence_without_gaps,  4) == 7);
  fail_unless( rewind_window_end_to_last_snp(9, coords_even, child_sequence_without_gaps,  4) == 7);
  fail_unless( rewind_window_end_to_last_snp(0, coords_even, child_sequence_without_gaps,  4) == 0);
  fail_unless( rewind_window_end_to_last_snp(2, coords_even, child_sequence_without_gaps,  4) == 1);
  fail_unless( rewind_window_end_to_last_snp(4, coords_even, child_sequence_without_gaps,  4) == 3);
  fail_unless( rewind_window_end_to_last_snp(6, coords_even, child_sequence_without_gaps,  4) == 5);
  fail_unless( rewind_window_end_to_last_snp(8, coords_even, child_sequence_without_gaps,  4) == 7);


	fail_unless( rewind_window_end_to_last_snp(0, coords_empty,child_sequence_with_gaps,  0) == 0);
  fail_unless( rewind_window_end_to_last_snp(0, coords_one,  child_sequence_with_gaps,  1) == 0);
	fail_unless( rewind_window_end_to_last_snp(1, coords_one,  child_sequence_with_gaps,  1) == 1);
  fail_unless( rewind_window_end_to_last_snp(1, coords_odd,  child_sequence_with_gaps,  3) == 1);
  fail_unless( rewind_window_end_to_last_snp(3, coords_odd,  child_sequence_with_gaps,  3) == 3);
  fail_unless( rewind_window_end_to_last_snp(5, coords_odd,  child_sequence_with_gaps,  3) == 5);
  fail_unless( rewind_window_end_to_last_snp(0, coords_odd,  child_sequence_with_gaps,  3) == 0);
  fail_unless( rewind_window_end_to_last_snp(2, coords_odd,  child_sequence_with_gaps,  3) == 2);
  fail_unless( rewind_window_end_to_last_snp(4, coords_odd,  child_sequence_with_gaps,  3) == 3);
  fail_unless( rewind_window_end_to_last_snp(1, coords_even, child_sequence_with_gaps,  4) == 1);
  fail_unless( rewind_window_end_to_last_snp(3, coords_even, child_sequence_with_gaps,  4) == 3);
  fail_unless( rewind_window_end_to_last_snp(5, coords_even, child_sequence_with_gaps,  4) == 5);
  fail_unless( rewind_window_end_to_last_snp(7, coords_even, child_sequence_with_gaps,  4) == 5);
  fail_unless( rewind_window_end_to_last_snp(9, coords_even, child_sequence_with_gaps,  4) == 5);
  fail_unless( rewind_window_end_to_last_snp(0, coords_even, child_sequence_with_gaps,  4) == 0);
  fail_unless( rewind_window_end_to_last_snp(2, coords_even, child_sequence_with_gaps,  4) == 2);
  fail_unless( rewind_window_end_to_last_snp(4, coords_even, child_sequence_with_gaps,  4) == 3);
  fail_unless( rewind_window_end_to_last_snp(6, coords_even, child_sequence_with_gaps,  4) == 5);
  fail_unless( rewind_window_end_to_last_snp(8, coords_even, child_sequence_with_gaps,  4) == 5);
}
END_TEST

START_TEST (check_get_window_end_coordinates_excluding_gaps)
{
  int coords_empty[0] = {};
  int coords_one[1]   = {1};
  int coords_odd[3]   = {1,3,5};
  int coords_even[8]  = {1,3,5,7,11,13,17,19};
	char *child_sequence_without_gaps = "ACGTACGT";
	char *child_sequence_with_gaps = "-AC-GT-A";
	
	//int get_window_end_coordinates_excluding_gaps(1, 3, coords_even, char * child_sequence, int number_of_snps)

	fail_unless( get_window_end_coordinates_excluding_gaps(0, 3, coords_empty,child_sequence_without_gaps,  0) == 3);
  fail_unless( get_window_end_coordinates_excluding_gaps(0, 3, coords_one,  child_sequence_without_gaps,  1) == 3);
	fail_unless( get_window_end_coordinates_excluding_gaps(1, 3, coords_one,  child_sequence_without_gaps,  1) == 4);
  fail_unless( get_window_end_coordinates_excluding_gaps(1, 3, coords_odd,  child_sequence_without_gaps,  3) == 4);
  fail_unless( get_window_end_coordinates_excluding_gaps(3, 3, coords_odd,  child_sequence_without_gaps,  3) == 6);
  fail_unless( get_window_end_coordinates_excluding_gaps(5, 3, coords_odd,  child_sequence_without_gaps,  3) == 8);
  fail_unless( get_window_end_coordinates_excluding_gaps(0, 3, coords_odd,  child_sequence_without_gaps,  3) == 3);
  fail_unless( get_window_end_coordinates_excluding_gaps(2, 3, coords_odd,  child_sequence_without_gaps,  3) == 5);
  fail_unless( get_window_end_coordinates_excluding_gaps(4, 3, coords_odd,  child_sequence_without_gaps,  3) == 7);
  fail_unless( get_window_end_coordinates_excluding_gaps(1, 3, coords_even, child_sequence_without_gaps,  8) == 4);
  fail_unless( get_window_end_coordinates_excluding_gaps(3, 3, coords_even, child_sequence_without_gaps,  8) == 6);
  fail_unless( get_window_end_coordinates_excluding_gaps(5, 3, coords_even, child_sequence_without_gaps,  8) == 8);
  fail_unless( get_window_end_coordinates_excluding_gaps(7, 3, coords_even, child_sequence_without_gaps,  8) == 10);
  fail_unless( get_window_end_coordinates_excluding_gaps(9, 3, coords_even, child_sequence_without_gaps,  8) == 12);
  fail_unless( get_window_end_coordinates_excluding_gaps(0, 3, coords_even, child_sequence_without_gaps,  8) == 3);
  fail_unless( get_window_end_coordinates_excluding_gaps(2, 3, coords_even, child_sequence_without_gaps,  8) == 5);
  fail_unless( get_window_end_coordinates_excluding_gaps(4, 3, coords_even, child_sequence_without_gaps,  8) == 7);
  fail_unless( get_window_end_coordinates_excluding_gaps(6, 3, coords_even, child_sequence_without_gaps,  8) == 9);
  fail_unless( get_window_end_coordinates_excluding_gaps(8, 3, coords_even, child_sequence_without_gaps,  8) == 11);
           
	fail_unless( get_window_end_coordinates_excluding_gaps(0, 3, coords_empty,child_sequence_with_gaps,  0) == 3);
  fail_unless( get_window_end_coordinates_excluding_gaps(0, 3, coords_one,  child_sequence_with_gaps,  1) == 4);
	fail_unless( get_window_end_coordinates_excluding_gaps(1, 3, coords_one,  child_sequence_with_gaps,  1) == 5);
  fail_unless( get_window_end_coordinates_excluding_gaps(1, 3, coords_odd,  child_sequence_with_gaps,  3) == 5);
  fail_unless( get_window_end_coordinates_excluding_gaps(3, 3, coords_odd,  child_sequence_with_gaps,  3) == 6);
  fail_unless( get_window_end_coordinates_excluding_gaps(5, 3, coords_odd,  child_sequence_with_gaps,  3) == 8);
  fail_unless( get_window_end_coordinates_excluding_gaps(0, 3, coords_odd,  child_sequence_with_gaps,  3) == 4);
  fail_unless( get_window_end_coordinates_excluding_gaps(2, 3, coords_odd,  child_sequence_with_gaps,  3) == 5);
  fail_unless( get_window_end_coordinates_excluding_gaps(4, 3, coords_odd,  child_sequence_with_gaps,  3) == 7);
  fail_unless( get_window_end_coordinates_excluding_gaps(1, 3, coords_even, child_sequence_with_gaps,  8) == 5);
  fail_unless( get_window_end_coordinates_excluding_gaps(3, 3, coords_even, child_sequence_with_gaps,  8) == 6);
  fail_unless( get_window_end_coordinates_excluding_gaps(5, 3, coords_even, child_sequence_with_gaps,  8) == 9);
  fail_unless( get_window_end_coordinates_excluding_gaps(7, 3, coords_even, child_sequence_with_gaps,  8) == 11);
  fail_unless( get_window_end_coordinates_excluding_gaps(9, 3, coords_even, child_sequence_with_gaps,  8) == 12);
  fail_unless( get_window_end_coordinates_excluding_gaps(0, 3, coords_even, child_sequence_with_gaps,  8) == 4);
  fail_unless( get_window_end_coordinates_excluding_gaps(2, 3, coords_even, child_sequence_with_gaps,  8) == 5);
  fail_unless( get_window_end_coordinates_excluding_gaps(4, 3, coords_even, child_sequence_with_gaps,  8) == 7);
  fail_unless( get_window_end_coordinates_excluding_gaps(6, 3, coords_even, child_sequence_with_gaps,  8) == 10);
  fail_unless( get_window_end_coordinates_excluding_gaps(8, 3, coords_even, child_sequence_with_gaps,  8) == 11);
	
}
END_TEST

START_TEST (check_find_number_of_snps_in_block)
{
  int coords_empty[0] = {};
  int coords_even[8]  = {1,3,5,7,11,13,17,19};
	char *child_sequence = "AAAAA-AAAAAAAAAAAAA";
	
	fail_unless( find_number_of_snps_in_block(1,3,  coords_empty, child_sequence, 8) == 0);
	fail_unless( find_number_of_snps_in_block(2,2,  coords_even, child_sequence, 8) == 0);
	fail_unless( find_number_of_snps_in_block(1,3,  coords_even, child_sequence, 8) == 1);
  fail_unless( find_number_of_snps_in_block(1,4,  coords_even, child_sequence, 8) == 2);
  fail_unless( find_number_of_snps_in_block(1,19, coords_even, child_sequence, 8) == 6);
  fail_unless( find_number_of_snps_in_block(0,20, coords_even, child_sequence, 8) == 7);
	
}
END_TEST

	//int calculate_number_of_snps_excluding_gaps(char * ancestor_sequence, char * child_sequence, int child_sequence_size, int * branch_snp_coords, int * snp_locations)
START_TEST (check_calculate_number_of_snps_excluding_gaps)
{
  char * reference_sequence           = "TTTTTTTTTT";
  char * ancestor_sequence_no_snps    = "AAAAAAAAAA";
	char * ancestor_sequence_one_snp    = "AAAAA.AAAA";
	char * ancestor_sequence_many_snps  = ".AAA.AAA..";
	char * ancestor_sequence_gaps       = ".AAA..AAA.";
	
	char * child_sequence_no_snps       = "AAAAAAAAAA";
	char * child_sequence_one_snp       = "AAAAACAAAA";
  char * child_sequence_many_snps     = "CAAACAAACC";
  char * child_sequence_gaps          = "-AAAC----C";
  char * child_ancestor_sequence      = ".AAAC....C";
  char * ancestor_sequence_with_dots  = ".C...C...C";

 	int child_sequence_size = 10;
	int * branch_snp_coords;
	char * branch_snp_sequence ;

	int snp_locations[10] = {0,4,8,20,30,40,50,88,90,100};
	
	branch_snp_coords   = (int *)  malloc((child_sequence_size+1)*sizeof(int)); 
	branch_snp_sequence = (char *) malloc((child_sequence_size+1)*sizeof(char));
	fail_unless( calculate_number_of_snps_excluding_gaps(ancestor_sequence_no_snps, child_sequence_no_snps, child_sequence_size, branch_snp_coords, snp_locations,reference_sequence,branch_snp_sequence ) == 0 );
	fail_unless(strcmp(branch_snp_sequence, "") == 0);
	
  branch_snp_coords   = (int *)  malloc((child_sequence_size+1)*sizeof(int)); 
  branch_snp_sequence = (char *) malloc((child_sequence_size+1)*sizeof(char));
  fail_unless( calculate_number_of_snps_excluding_gaps(ancestor_sequence_one_snp, child_sequence_one_snp, child_sequence_size, branch_snp_coords, snp_locations,reference_sequence,branch_snp_sequence ) == 1 );
  fail_unless(strcmp(branch_snp_sequence, "C") == 0);
	fail_unless(branch_snp_coords[0] == 40);

  branch_snp_coords   = (int *)  malloc((child_sequence_size+1)*sizeof(int)); 
  branch_snp_sequence = (char *) malloc((child_sequence_size+1)*sizeof(char));
  fail_unless( calculate_number_of_snps_excluding_gaps(ancestor_sequence_many_snps, child_sequence_many_snps, child_sequence_size, branch_snp_coords, snp_locations,reference_sequence,branch_snp_sequence ) == 4 );
	fail_unless(strcmp(branch_snp_sequence, "CCCC") == 0);
	fail_unless(branch_snp_coords[0] == 0);
	fail_unless(branch_snp_coords[3] == 100);

  branch_snp_coords   = (int *)  malloc((child_sequence_size+1)*sizeof(int)); 
  branch_snp_sequence = (char *) malloc((child_sequence_size+1)*sizeof(char));
  fail_unless( calculate_number_of_snps_excluding_gaps(ancestor_sequence_gaps, child_sequence_gaps, child_sequence_size, branch_snp_coords, snp_locations,reference_sequence,branch_snp_sequence ) == 2 );
	fail_unless(strcmp(branch_snp_sequence, "CC") == 0);
	fail_unless(branch_snp_coords[0] == 30);
	fail_unless(branch_snp_coords[1] == 100);

  branch_snp_coords   = (int *)  malloc((child_sequence_size+1)*sizeof(int)); 
  branch_snp_sequence = (char *) malloc((child_sequence_size+1)*sizeof(char));
  fail_unless( calculate_number_of_snps_excluding_gaps(ancestor_sequence_with_dots, child_ancestor_sequence, child_sequence_size, branch_snp_coords, snp_locations,reference_sequence,branch_snp_sequence ) == 3);
	fail_unless(strcmp(branch_snp_sequence, "AAC") == 0);
	fail_unless(branch_snp_coords[0] == 8);
	fail_unless(branch_snp_coords[1] == 20);
}
END_TEST

Suite * check_snp_searching_suite (void)
{
  Suite *s = suite_create ("snp_searching");

  TCase *tc_snp_searching = tcase_create ("snp_searching");
  tcase_add_test (tc_snp_searching, check_find_starting_index_windows );
	tcase_add_test (tc_snp_searching, check_advance_window_start_to_next_snp);
	tcase_add_test (tc_snp_searching, check_rewind_window_end_to_last_snp);
	tcase_add_test (tc_snp_searching, check_get_window_end_coordinates_excluding_gaps);
	tcase_add_test (tc_snp_searching, check_find_number_of_snps_in_block);
	tcase_add_test (tc_snp_searching, check_calculate_number_of_snps_excluding_gaps);
  suite_add_tcase (s, tc_snp_searching);

  return s;
}
