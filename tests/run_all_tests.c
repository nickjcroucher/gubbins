#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_parse_phylip.h"
#include "check_snp_sites.h"
#include "check_vcf_parsing.h"
#include "check_branch_sequences.h"
#include "check_gubbins.h"
#include "check_snp_searching.h"

int main (void)
{
  int number_failed;
  Suite *s = snp_sites_suite ();
  SRunner *sr = srunner_create (s);
  srunner_add_suite (sr, parse_phylip_suite());
  srunner_add_suite (sr, parse_vcf_suite());
  srunner_add_suite (sr, check_branch_sequences_suite());
  srunner_add_suite (sr, run_gubbins_suite());
  srunner_add_suite (sr, check_snp_searching_suite());

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
