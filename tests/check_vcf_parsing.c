#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_vcf_parsing.h"
#include "helper_methods.h"
#include "../parse_vcf.h"


START_TEST (check_parsing_of_vcf_files)
{
  // Todo insert tests
}
END_TEST


Suite * parse_vcf_suite(void)
{
  Suite *s = suite_create ("Parsing a vcf file");
  TCase *tc_parse_vcf = tcase_create ("check_parsing_of_vcf_files");
  tcase_add_test (tc_parse_vcf, check_parsing_of_vcf_files);
  suite_add_tcase (s, tc_parse_vcf);
  return s;
}
