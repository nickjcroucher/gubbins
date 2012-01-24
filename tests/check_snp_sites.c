
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "../snp_sites.h"
#include "../alignment_file.h"
#include "../parse_phylip.h"
	

START_TEST (valid_alignment_with_one_line_per_sequence)
{
  generate_snp_sites("data/alignment_file_one_line_per_sequence.aln");
  fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.vcf", "alignment_file_one_line_per_sequence.aln.vcf" ) == 1, "Invalid VCF file for 1 line per seq" );
  fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.phylip", "alignment_file_one_line_per_sequence.aln.phylip" ) == 1, "Invalid Phylip file for 1 line per seq" );
  fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_one_line_per_sequence.aln.snp_sites.aln" ) == 1 , "Invalid ALN file for 1 line per seq");
  remove("alignment_file_one_line_per_sequence.aln.vcf");
  remove("alignment_file_one_line_per_sequence.aln.phylip");
  remove("alignment_file_one_line_per_sequence.aln.snp_sites.aln");
}
END_TEST


START_TEST (valid_alignment_with_one_line_per_sequence_gzipped)
{
  generate_snp_sites("data/alignment_file_one_line_per_sequence.aln.gz");
  fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.vcf", "alignment_file_one_line_per_sequence.aln.gz.vcf" ) == 1, "Invalid VCF file for 1 line per seq" );
  fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.phylip", "alignment_file_one_line_per_sequence.aln.gz.phylip" ) == 1, "Invalid Phylip file for 1 line per seq" );
  fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_one_line_per_sequence.aln.gz.snp_sites.aln" ) == 1 , "Invalid ALN file for 1 line per seq");
  remove("alignment_file_one_line_per_sequence.aln.gz.vcf");
  remove("alignment_file_one_line_per_sequence.aln.gz.phylip");
  remove("alignment_file_one_line_per_sequence.aln.gz.snp_sites.aln");
}
END_TEST

START_TEST (valid_alignment_with_multiple_lines_per_sequence)
{
    generate_snp_sites("data/alignment_file_multiple_lines_per_sequence.aln");
    fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.vcf", "alignment_file_multiple_lines_per_sequence.aln.vcf" ) == 1, "Invalid VCF file for multiple lines per seq" );
    fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.phylip", "alignment_file_multiple_lines_per_sequence.aln.phylip" ) == 1, "Invalid Phylip file for multiple lines per seq" );
    fail_unless( compare_files("data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_multiple_lines_per_sequence.aln.snp_sites.aln" ) == 1 ,"Invalid ALN file for multiple lines per seq");
    remove("alignment_file_multiple_lines_per_sequence.aln.vcf");
    remove("alignment_file_multiple_lines_per_sequence.aln.phylip");
    remove("alignment_file_multiple_lines_per_sequence.aln.snp_sites.aln");
}
END_TEST


START_TEST (valid_genome_length)
{
  fail_unless( genome_length("data/alignment_file_one_line_per_sequence.aln") == 2000 );
}
END_TEST

START_TEST (valid_genome_length_with_multiple_lines_per_sequence)
{
  fail_unless( genome_length("data/alignment_file_multiple_lines_per_sequence.aln") == 2000 );
}
END_TEST

START_TEST (valid_number_of_sequences_in_file)
{
  fail_unless( number_of_sequences_in_file("data/alignment_file_one_line_per_sequence.aln") == 109 );
}
END_TEST

START_TEST (valid_number_of_sequences_in_file_with_multiple_lines_per_sequence)
{
  fail_unless( number_of_sequences_in_file("data/alignment_file_multiple_lines_per_sequence.aln") == 109 );
}
END_TEST

START_TEST (valid_initial_reference_sequence)
{
  char actual_reference_sequence[2001];
  char *expected_reference_sequence = "-------------------------CTATATAGAGATCTTTTTATTAGATCTACTATTAAGGAGCAGGATCTTTGTGGATAAGTGAAAAATGATCAACAAGATCATGCGATTCAGAAGGATCAGATCGTGTGATCAACCACTGATCTGTTCAAGGATTAGCTGGGATCAAAAACCTATGTTATACACAGCCACCTTGGGATCTAAAACTTGTTATATGGATAACTATAGGAAGATCACCGGATAATCGTATAGTTATCCACATGAGATTTGATTGAAAAAGCATCAATCAATTTTTTCACTACCGTTAAATTTATCCACAATCCAAAAAAAAGAGCGGCATTAAGCCGCTCTGCATGGAATAGGTCATTATTTAGAAGCGATTGATGACGCGTTTGAGCCAAGCTTCAGCGGCATCTTCAGGCACTGGGTGCTCTTGTACATCGATGGTAAAGCAGTTGGCCAGAGGTTTAGCACCAATATCCCCCAGCAGCTGATAGGCATGTTTACCTGCCGCGCAGAAAGTATCGTAGCTTGAATCACCAATCGCGACCACGGCATAACGTAGTGCAGAGGTATTCGGTGGTGTATTCTGCAGAGCCTGAATAAAGGGCTGGATATTATCCGGGTACTCACCAGCCCCGTGGGTTGAGGTGATGATCAGCCAAGTCCCTTTAGCAGGGATCTCACTCATGTTGGGCTGGTTATGAATTTTGGTGTCAAAGCCTTGTTCTTGCAGTAAATCACTCAGGTGGTCACCCACATATTCCGCACCGCCTAGGGTGCTGCCAGTAATGATATGAATCATAGCGTTACTCTATTTCCCAATACAGAATGATGAAAAAATGCGGCCAAGCAGATCATCGGAGCTGAACTCGCCCGTAATTTCGTTAAGGTGTTGCTGGGCTATACGCAGCTCTTCGGCGAGGATTTCTCCGGCCATATAGCCTTCAAGTTGTTGCTGGCCAATCGCTAAGTGCTCTGCGGCTCGCTCTAGGGCATCGAGATGACGGCGGCGTGCCATAAAGCCACCTTCCTGATTGCCTGAAAAACCCATGCACTCTTTGAGGTGCTGACGCAAGGCATCGACCCCTTGGCCTGTTTTGGCTGATAGGCGGATCAAGGTGGGTTGATTAACATGGCAGATCCCAAGGGGCTCACCAGTTTGATCGGCTTTATTACGGATCACAGTGATCCCAATATTCTCTGGCAGTTTGTCAACAAAATCAGGCCAGATGTCCTGTGGATCGGTGGCCTCTGTGGTGGTGCCATCGACCATAAACAGTACGCGATCGGCTTGGCGGATCTCTTCCCATGCGCGCTCAATACCAATTTTTTCTACCGCATCAGAAGCGTCTCGTAGTCCCGCAGTATCGATGATGTGCAGCGGCATCCCATCAATATGGATATGCTCACGCAGAACATCACGGGTGGTACCGGCAATGTCGGTAACGATGGCAGACTCTTTACCTGAAAGCGCATTGAGTAGGCTCGATTTACCCGCATTAGGACGCCCAGCAATCACCACCTTCATCCCTTCGCGCATAATGGCGCCTTGGTTGGCTTCACGGCGCACTGCGGCAAGATTATCTATGATGGTTTGCAGATCAGCGGAAACCTTACCATCGGCCAGAAAATCGATCTCTTCTTCTGGGAAATCAATTGCGGCTTCAACATAGATGCGCAGGTGAATCAGCGATTCCACCAAGGTATGGATGCGTTTAGAAAACTCGCCTTGCAGTGATTGCAGCGCGGATTTCGCGGCTTGCTCAGAGCTGGCATCAATCAGGTCTGCGATGGCTTCCGCTTGGGTTAAATCCATCTTGTCATTGAGGAAAGCGCGTTCTGAGAATTCACCGGGACGGGCTGGGCGCACTCCTTTAATCTGCAAAATACGGCGGATCAGCATATCCATGACGACCGGGCCACCGTGACCTTGCAGCTCAAGCACATCTTCACCGGTAAATGAATGAGGATTGGGGAAAAACAGCGCAATGCCTTG";
  build_reference_sequence(actual_reference_sequence, "data/alignment_file_multiple_lines_per_sequence.aln") ;
  fail_unless( strcmp(actual_reference_sequence,expected_reference_sequence) == 0 );
}
END_TEST  

START_TEST (number_of_snps_detected)
{
  char actual_reference_sequence[2001];
  build_reference_sequence(actual_reference_sequence, "data/alignment_file_multiple_lines_per_sequence.aln") ;
  fail_unless(  detect_snps(actual_reference_sequence, "data/alignment_file_multiple_lines_per_sequence.aln", 2000) == 5);
}
END_TEST

START_TEST (number_of_snps_detected_small)
{
  char actual_reference_sequence[9];
  build_reference_sequence(actual_reference_sequence, "data/small_alignment.aln");
  fail_unless(  detect_snps(actual_reference_sequence, "data/small_alignment.aln", 8) == 1);
}
END_TEST


START_TEST (sample_names_from_alignment_file)
{
  char *expected_sequence_names[] ={"reference_sequence","comparison_sequence","another_comparison_sequence"};
  char* sequence_names[3];
  int i = 0;
	sequence_names[3-1] = '\0';
  for(i = 0; i < 3; i++)
	{
		sequence_names[i] = malloc(30*sizeof(char));
	}
  get_sample_names_for_header("data/small_alignment.aln",sequence_names, 3);
  
  for(i =0; i< 3; i++)
  {
    fail_unless( strcmp(expected_sequence_names[i], sequence_names[i]) ==0 );
  }
}
END_TEST



START_TEST (phylip_read_in_small_file)
{
  load_sequences_from_phylib_file("data/small_phylip_file.phylip");
  
  fail_unless( number_of_samples_from_parse_phylip() == 3);
  fail_unless( find_sequence_index_from_sample_name("2956_6_1") == 0);
  fail_unless( find_sequence_index_from_sample_name("2956_6_2") == 1);
  fail_unless( find_sequence_index_from_sample_name("2956_6_3") == 2);
  
  char *reference_bases = "*ACG*";
  char *filtered_bases_for_snps[3];

  filter_sequence_bases_and_rotate(reference_bases, filtered_bases_for_snps, 3);
  fail_unless( strcmp(filtered_bases_for_snps[0], "AAT") == 0 );
  fail_unless( strcmp(filtered_bases_for_snps[1], "CGT") == 0 );
  fail_unless( strcmp(filtered_bases_for_snps[2], "GGT") == 0 );
  
  char *sample_names[3];
  get_sample_names_from_parse_phylip(sample_names);
  fail_unless( strcmp(sample_names[0],"2956_6_1") == 0 );
  fail_unless( strcmp(sample_names[1],"2956_6_2") == 0 );
  fail_unless( strcmp(sample_names[2],"2956_6_3") == 0 );
  
  fail_unless( does_column_contain_snps(0, 'A') == 0);
  fail_unless( does_column_contain_snps(1, 'A') == 1);
  fail_unless( does_column_contain_snps(2, 'A') == 1);
  // bad reference base
  fail_unless( does_column_contain_snps(0, 'X') == 1);
  
  char sequence_bases[10];
  get_sequence_for_sample_name(sequence_bases, "2956_6_2");
  fail_unless( strcmp(sequence_bases, "AAGGC") == 0);
  
  update_sequence_base('X', 1, 4);
  get_sequence_for_sample_name(sequence_bases, "2956_6_2");
  fail_unless( strcmp(sequence_bases, "AAGGX") == 0);
  
}
END_TEST


Suite * snp_sites_suite (void)
{
  Suite *s = suite_create ("Creating_SNP_Sites");

  TCase *tc_alignment_file = tcase_create ("alignment_file");
  tcase_add_test (tc_alignment_file, valid_genome_length);
  tcase_add_test (tc_alignment_file, valid_genome_length_with_multiple_lines_per_sequence);
  tcase_add_test (tc_alignment_file, valid_number_of_sequences_in_file);
  tcase_add_test (tc_alignment_file, valid_number_of_sequences_in_file_with_multiple_lines_per_sequence);
  tcase_add_test (tc_alignment_file, valid_initial_reference_sequence);
  tcase_add_test (tc_alignment_file, number_of_snps_detected_small);
  tcase_add_test (tc_alignment_file, number_of_snps_detected);
  tcase_add_test (tc_alignment_file, sample_names_from_alignment_file);
  suite_add_tcase (s, tc_alignment_file);
  

  TCase *tc_snp_sites = tcase_create ("snp_sites");
  tcase_add_test (tc_snp_sites, valid_alignment_with_one_line_per_sequence);
  tcase_add_test (tc_snp_sites, valid_alignment_with_multiple_lines_per_sequence);
  tcase_add_test (tc_snp_sites, valid_alignment_with_one_line_per_sequence_gzipped);
  suite_add_tcase (s, tc_snp_sites);

  TCase *tc_phylip = tcase_create ("phylip_files");
  tcase_add_test (tc_phylip, phylip_read_in_small_file);
  suite_add_tcase (s, tc_phylip);

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

