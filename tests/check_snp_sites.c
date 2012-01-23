
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "../snp_sites.h"
	

START_TEST (valid_alignment_with_one_line_per_sequence)
{
  generate_snp_sites("data/alignment_file_one_line_per_sequence.aln");
  fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.vcf", "alignment_file_one_line_per_sequence.aln.vcf" ) == 1, 'Invalid VCF file for 1 line per seq' );
  fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.phylip", "alignment_file_one_line_per_sequence.aln.phylip" ) == 1, 'Invalid Phylip file for 1 line per seq' );
  fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_one_line_per_sequence.aln.snp_sites.aln" ) == 1 ,'Invalid ALN file for 1 line per seq');
  remove("alignment_file_one_line_per_sequence.aln.vcf");
  remove("alignment_file_one_line_per_sequence.aln.phylip");
  remove("alignment_file_one_line_per_sequence.aln.snp_sites.aln");
}
END_TEST

Suite * snp_sites_suite (void)
{
  Suite *s = suite_create ("Snp_sites");

  /* Core test case */
  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, valid_alignment_with_one_line_per_sequence);
  suite_add_tcase (s, tc_core);

  return s;
}

int compare_files(char expected_output_filename[],char actual_output_filename[] )
{
  FILE *expected_output_fh;
  FILE *actual_output_fh;
  
  char    *expected_buffer;
  char    *actual_buffer;
  long    numbytes;
  
  expected_output_fh = fopen(expected_output_filename, "r");
  actual_output_fh = fopen(actual_output_filename, "r");
  
  fseek(expected_output_fh, 0L, SEEK_END);
  numbytes = ftell(expected_output_fh);
  fseek(expected_output_fh, 0L, SEEK_SET);	
  expected_buffer = (char*)calloc(numbytes, sizeof(char));	
  fread(expected_buffer, sizeof(char), numbytes, expected_output_fh);
  fclose(expected_output_fh);
  
  fseek(actual_output_fh, 0L, SEEK_END);
  numbytes = ftell(actual_output_fh);
  fseek(actual_output_fh, 0L, SEEK_SET);	
  actual_buffer = (char*)calloc(numbytes, sizeof(char));	
  fread(actual_buffer, sizeof(char), numbytes, actual_output_fh);
  fclose(actual_output_fh);
  
  if(strcmp(expected_buffer,actual_buffer) == 0)
  { 
    free(expected_buffer);
    free(actual_buffer);
    return 1;
  }

  free(expected_buffer);
  free(actual_buffer);
  
  return 0;
}


int main (void)
{
  int number_failed;
  Suite *s = snp_sites_suite ();
  SRunner *sr = srunner_create (s);
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

