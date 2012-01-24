#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_calculate_ancestor_sequence.h"
#include "helper_methods.h"
#include "../gubbins.h"

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
}
END_TEST

Suite * calculate_ancestor_sequence_suite(void)
{
   Suite *s = suite_create ("Calculate an ancestor sequence");
   TCase *tc_phylip = tcase_create ("check_calculate_ancestor_sequence");
   tcase_add_test (tc_phylip, check_calculate_ancestor_sequence);
   suite_add_tcase (s, tc_phylip);
   return s;
}