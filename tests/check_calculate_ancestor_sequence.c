#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_calculate_ancestor_sequence.h"
#include "helper_methods.h"
#include "gubbins.h"

START_TEST (check_calculate_ancestor_sequence)
{
  char ancestor_sequence[6];
  
  // one snp
  char *child_sequences[] ={"AAAGA","AAAGA","ACAGA"};
  fail_unless( strcmp(calculate_ancestor_sequence(ancestor_sequence, child_sequences,5, 3), "A.AGA") == 0);  
  
  // all are snps
  char *child_sequences_all_snps[] ={"AAAAA","CCCCC","-----"};
  fail_unless( strcmp(calculate_ancestor_sequence(ancestor_sequence, child_sequences_all_snps,5, 3), ".....") == 0);  
  
  // no snps
  char *child_sequences_no_snps[] ={"AAAAA","AAAAA", "AAAAA"};
  fail_unless( strcmp(calculate_ancestor_sequence(ancestor_sequence, child_sequences_no_snps,5, 3), "AAAAA") == 0);

  // child not a leaf node and has snps accounted for
  char *child_sequences_not_leaf_and_no_snps[] ={"AAAAA","A..AA", "AAAAA"};
  fail_unless( strcmp(calculate_ancestor_sequence(ancestor_sequence, child_sequences_not_leaf_and_no_snps,5, 3), "AAAAA") == 0);

  // child not a leaf node and has some snps
  char *child_sequences_not_leaf_and_has_snps[] ={"AAAAA","A..AA", "AACAA"};
  fail_unless( strcmp(calculate_ancestor_sequence(ancestor_sequence, child_sequences_not_leaf_and_has_snps,5, 3), "AA.AA") == 0);

  // children not a leaf nodes 
  char *child_sequence_with_no_snps[] ={".....","AACAA"};
  fail_unless( strcmp(calculate_ancestor_sequence(ancestor_sequence, child_sequence_with_no_snps,5, 2), "AACAA") == 0);

  // children not a leaf nodes 
  char *child_sequence_where_first_seq_has_dots[] ={".....","AA.AA","AAGAA","AATAA"};
  fail_unless( strcmp(calculate_ancestor_sequence(ancestor_sequence, child_sequence_where_first_seq_has_dots,5, 4), "AA.AA") == 0);
}
END_TEST

START_TEST (check_find_first_real_base)
{
  char *child_sequence_with_dots[] ={"..A.",".CC.","GGG."};
	fail_unless( find_first_real_base(0,3,child_sequence_with_dots) == 'G');
	fail_unless( find_first_real_base(1,3,child_sequence_with_dots) == 'C');
	fail_unless( find_first_real_base(2,3,child_sequence_with_dots) == 'A');
	fail_unless( find_first_real_base(3,3,child_sequence_with_dots) == '.');
	
  char *child_sequence_with_ns[] ={"NNAN","NCCN","GGGN"};
	fail_unless( find_first_real_base(0,3,child_sequence_with_ns) == 'G');
	fail_unless( find_first_real_base(1,3,child_sequence_with_ns) == 'C');
	fail_unless( find_first_real_base(2,3,child_sequence_with_ns) == 'A');
	fail_unless( find_first_real_base(3,3,child_sequence_with_ns) == 'N');
	
	char *child_sequence_with_gaps[] ={"--A-","-CC-","GGG-"};
	fail_unless( find_first_real_base(0,3,child_sequence_with_gaps) == 'G');
	fail_unless( find_first_real_base(1,3,child_sequence_with_gaps) == 'C');
	fail_unless( find_first_real_base(2,3,child_sequence_with_gaps) == 'A');
	fail_unless( find_first_real_base(3,3,child_sequence_with_gaps) == '-');
}
END_TEST


Suite * calculate_ancestor_sequence_suite(void)
{
   Suite *s = suite_create ("Calculate an ancestor sequence");
   TCase *tc_phylip = tcase_create ("check_calculate_ancestor_sequence");
   tcase_add_test (tc_phylip, check_calculate_ancestor_sequence);
	 tcase_add_test (tc_phylip, check_find_first_real_base);
   suite_add_tcase (s, tc_phylip);
   return s;
}


