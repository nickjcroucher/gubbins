
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_snp_sites.h"
#include "helper_methods.h"

#include "snp_sites.h"
#include "alignment_file.h"
	
	
START_TEST (valid_alignment_with_large_numbers_of_snps)
{
  generate_snp_sites("../tests/data/alignment_file_with_large_number_of_snps.aln",0,"");
  fail_unless( compare_files("../tests/data/alignment_file_with_large_number_of_snps.aln.vcf", "alignment_file_with_large_number_of_snps.aln.vcf" ) == 1, "Invalid VCF file for 1 line per seq" );
  fail_unless( compare_files("../tests/data/alignment_file_with_large_number_of_snps.aln.phylip", "alignment_file_with_large_number_of_snps.aln.phylip" ) == 1, "Invalid Phylip file for 1 line per seq" );
  fail_unless( compare_files("../tests/data/alignment_file_with_large_number_of_snps.aln.snp_sites.aln","alignment_file_with_large_number_of_snps.aln.snp_sites.aln" ) == 1 , "Invalid ALN file for 1 line per seq");
  remove("alignment_file_with_large_number_of_snps.aln.vcf");
  remove("alignment_file_with_large_number_of_snps.aln.phylip");
  remove("alignment_file_with_large_number_of_snps.aln.snp_sites.aln");
}
END_TEST
	

START_TEST (valid_alignment_with_one_line_per_sequence)
{
  generate_snp_sites("../tests/data/alignment_file_one_line_per_sequence.aln",0,"");
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.phylip", "alignment_file_one_line_per_sequence.aln.phylip" ) == 1, "Invalid Phylip file for 1 line per seq" );
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.vcf", "alignment_file_one_line_per_sequence.aln.vcf" ) == 1, "Invalid VCF file for 1 line per seq" );
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_one_line_per_sequence.aln.snp_sites.aln" ) == 1 , "Invalid ALN file for 1 line per seq");
  remove("alignment_file_one_line_per_sequence.aln.vcf");
  remove("alignment_file_one_line_per_sequence.aln.phylip");
  remove("alignment_file_one_line_per_sequence.aln.snp_sites.aln");
}
END_TEST


START_TEST (valid_alignment_with_n_as_gap)
{
	  generate_snp_sites("../tests/data/alignment_file_with_n.aln",0,"");
	  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.vcf", "alignment_file_with_n.aln.vcf" ) == 1, "Invalid VCF file for 1 line per seq" );
	  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.phylip", "alignment_file_with_n.aln.phylip" ) == 1, "Invalid Phylip file for 1 line per seq" );
	  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_with_n.aln.snp_sites.aln" ) == 1 , "Invalid ALN file for 1 line per seq");
	  remove("alignment_file_with_n.aln.vcf");
	  remove("alignment_file_with_n.aln.phylip");
	  remove("alignment_file_with_n.aln.snp_sites.aln");
}
END_TEST



START_TEST (valid_alignment_with_one_line_per_sequence_gzipped)
{
  generate_snp_sites("../tests/data/alignment_file_one_line_per_sequence.aln.gz",0,"");
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.vcf", "alignment_file_one_line_per_sequence.aln.gz.vcf" ) == 1, "Invalid VCF file for 1 line per seq" );
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.phylip", "alignment_file_one_line_per_sequence.aln.gz.phylip" ) == 1, "Invalid Phylip file for 1 line per seq" );
  fail_unless( compare_files("../tests/data/alignment_file_one_line_per_sequence.aln.snp_sites.aln","alignment_file_one_line_per_sequence.aln.gz.snp_sites.aln" ) == 1 , "Invalid ALN file for 1 line per seq");
  remove("alignment_file_one_line_per_sequence.aln.gz.vcf");
  remove("alignment_file_one_line_per_sequence.aln.gz.phylip");
  remove("alignment_file_one_line_per_sequence.aln.gz.snp_sites.aln");
}
END_TEST

START_TEST (valid_alignment_with_multiple_lines_per_sequence)
{
    generate_snp_sites("../tests/data/alignment_file_multiple_lines_per_sequence.aln",0,"");
    fail_unless( compare_files("../tests/data/alignment_file_multiple_lines_per_sequence.aln.vcf", "alignment_file_multiple_lines_per_sequence.aln.vcf" ) == 1, "Invalid VCF file for multiple lines per seq" );
    fail_unless( compare_files("../tests/data/alignment_file_multiple_lines_per_sequence.aln.phylip", "alignment_file_multiple_lines_per_sequence.aln.phylip" ) == 1, "Invalid Phylip file for multiple lines per seq" );
    fail_unless( compare_files("../tests/data/alignment_file_multiple_lines_per_sequence.aln.snp_sites.aln","alignment_file_multiple_lines_per_sequence.aln.snp_sites.aln" ) == 1 ,"Invalid ALN file for multiple lines per seq");
    remove("alignment_file_multiple_lines_per_sequence.aln.vcf");
    remove("alignment_file_multiple_lines_per_sequence.aln.phylip");
    remove("alignment_file_multiple_lines_per_sequence.aln.snp_sites.aln");
}
END_TEST

START_TEST (two_sequences)
{
    generate_snp_sites("../tests/data/two_sequences.aln",0,"");
    fail_unless( compare_files("../tests/data/two_sequences.aln.vcf", "two_sequences.aln.vcf" ) == 1, "Invalid VCF file for multiple lines per seq" );
    fail_unless( compare_files("../tests/data/two_sequences.aln.phylip", "two_sequences.aln.phylip" ) == 1, "Invalid Phylip file for multiple lines per seq" );
    fail_unless( compare_files("../tests/data/two_sequences.aln.snp_sites.aln","two_sequences.aln.snp_sites.aln" ) == 1 ,"Invalid ALN file for multiple lines per seq");
    remove("two_sequences.aln.vcf");
    remove("two_sequences.aln.phylip");
    remove("two_sequences.aln.snp_sites.aln");
}
END_TEST


START_TEST (valid_genome_length)
{
  fail_unless( genome_length("../tests/data/alignment_file_one_line_per_sequence.aln") == 2000 );
}
END_TEST

START_TEST (valid_genome_length_with_multiple_lines_per_sequence)
{
  fail_unless( genome_length("../tests/data/alignment_file_multiple_lines_per_sequence.aln") == 2000 );
}
END_TEST

START_TEST (valid_number_of_sequences_in_file)
{
  fail_unless( number_of_sequences_in_file("../tests/data/alignment_file_one_line_per_sequence.aln") == 109 );
}
END_TEST

START_TEST (valid_number_of_sequences_in_file_with_multiple_lines_per_sequence)
{
  fail_unless( number_of_sequences_in_file("../tests/data/alignment_file_multiple_lines_per_sequence.aln") == 109 );
}
END_TEST

START_TEST (valid_initial_reference_sequence)
{
  char actual_reference_sequence[2001];
  char *expected_reference_sequence = "-------------------------CTATATAGAGATCTTTTTATTAGATCTACTATTAAGGAGCAGGATCTTTGTGGATAAGTGAAAAATGATCAACAAGATCATGCGATTCAGAAGGATCAGATCGTGTGATCAACCACTGATCTGTTCAAGGATTAGCTGGGATCAAAAACCTATGTTATACACAGCCACCTTGGGATCTAAAACTTGTTATATGGATAACTATAGGAAGATCACCGGATAATCGTATAGTTATCCACATGAGATTTGATTGAAAAAGCATCAATCAATTTTTTCACTACCGTTAAATTTATCCACAATCCAAAAAAAAGAGCGGCATTAAGCCGCTCTGCATGGAATAGGTCATTATTTAGAAGCGATTGATGACGCGTTTGAGCCAAGCTTCAGCGGCATCTTCAGGCACTGGGTGCTCTTGTACATCGATGGTAAAGCAGTTGGCCAGAGGTTTAGCACCAATATCCCCCAGCAGCTGATAGGCATGTTTACCTGCCGCGCAGAAAGTATCGTAGCTTGAATCACCAATCGCGACCACGGCATAACGTAGTGCAGAGGTATTCGGTGGTGTATTCTGCAGAGCCTGAATAAAGGGCTGGATATTATCCGGGTACTCACCAGCCCCGTGGGTTGAGGTGATGATCAGCCAAGTCCCTTTAGCAGGGATCTCACTCATGTTGGGCTGGTTATGAATTTTGGTGTCAAAGCCTTGTTCTTGCAGTAAATCACTCAGGTGGTCACCCACATATTCCGCACCGCCTAGGGTGCTGCCAGTAATGATATGAATCATAGCGTTACTCTATTTCCCAATACAGAATGATGAAAAAATGCGGCCAAGCAGATCATCGGAGCTGAACTCGCCCGTAATTTCGTTAAGGTGTTGCTGGGCTATACGCAGCTCTTCGGCGAGGATTTCTCCGGCCATATAGCCTTCAAGTTGTTGCTGGCCAATCGCTAAGTGCTCTGCGGCTCGCTCTAGGGCATCGAGATGACGGCGGCGTGCCATAAAGCCACCTTCCTGATTGCCTGAAAAACCCATGCACTCTTTGAGGTGCTGACGCAAGGCATCGACCCCTTGGCCTGTTTTGGCTGATAGGCGGATCAAGGTGGGTTGATTAACATGGCAGATCCCAAGGGGCTCACCAGTTTGATCGGCTTTATTACGGATCACAGTGATCCCAATATTCTCTGGCAGTTTGTCAACAAAATCAGGCCAGATGTCCTGTGGATCGGTGGCCTCTGTGGTGGTGCCATCGACCATAAACAGTACGCGATCGGCTTGGCGGATCTCTTCCCATGCGCGCTCAATACCAATTTTTTCTACCGCATCAGAAGCGTCTCGTAGTCCCGCAGTATCGATGATGTGCAGCGGCATCCCATCAATATGGATATGCTCACGCAGAACATCACGGGTGGTACCGGCAATGTCGGTAACGATGGCAGACTCTTTACCTGAAAGCGCATTGAGTAGGCTCGATTTACCCGCATTAGGACGCCCAGCAATCACCACCTTCATCCCTTCGCGCATAATGGCGCCTTGGTTGGCTTCACGGCGCACTGCGGCAAGATTATCTATGATGGTTTGCAGATCAGCGGAAACCTTACCATCGGCCAGAAAATCGATCTCTTCTTCTGGGAAATCAATTGCGGCTTCAACATAGATGCGCAGGTGAATCAGCGATTCCACCAAGGTATGGATGCGTTTAGAAAACTCGCCTTGCAGTGATTGCAGCGCGGATTTCGCGGCTTGCTCAGAGCTGGCATCAATCAGGTCTGCGATGGCTTCCGCTTGGGTTAAATCCATCTTGTCATTGAGGAAAGCGCGTTCTGAGAATTCACCGGGACGGGCTGGGCGCACTCCTTTAATCTGCAAAATACGGCGGATCAGCATATCCATGACGACCGGGCCACCGTGACCTTGCAGCTCAAGCACATCTTCACCGGTAAATGAATGAGGATTGGGGAAAAACAGCGCAATGCCTTG";
  build_reference_sequence(actual_reference_sequence, "../tests/data/alignment_file_multiple_lines_per_sequence.aln") ;
  fail_unless( strcmp(actual_reference_sequence,expected_reference_sequence) == 0 );
}
END_TEST  

START_TEST (number_of_snps_detected)
{
  char actual_reference_sequence[2001];
  build_reference_sequence(actual_reference_sequence, "../tests/data/alignment_file_multiple_lines_per_sequence.aln") ;
  fail_unless(  detect_snps(actual_reference_sequence, "../tests/data/alignment_file_multiple_lines_per_sequence.aln", 2000,1) == 5);
}
END_TEST

START_TEST (number_of_snps_detected_include_gaps)
{
  char actual_reference_sequence[2001];
  build_reference_sequence(actual_reference_sequence, "../tests/data/alignment_file_multiple_lines_per_sequence.aln") ;
  fail_unless(  detect_snps(actual_reference_sequence, "../tests/data/alignment_file_multiple_lines_per_sequence.aln", 2000,0) == 1975);
}
END_TEST

START_TEST (number_of_snps_detected_small)
{
  char actual_reference_sequence[9];
  build_reference_sequence(actual_reference_sequence, "../tests/data/small_alignment.aln");
  fail_unless(  detect_snps(actual_reference_sequence, "../tests/data/small_alignment.aln", 8,0) == 1);
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
  get_sample_names_for_header("../tests/data/small_alignment.aln",sequence_names, 3);
  
  for(i =0; i< 3; i++)
  {
    fail_unless( strcmp(expected_sequence_names[i], sequence_names[i]) ==0 );
  }
}
END_TEST

START_TEST (check_strip_directory_from_filename_without_directory)
{
	char *input_filename_without_directory =  "my_file_name.aln";
	char output_filename[30];
	strip_directory_from_filename(input_filename_without_directory, output_filename);
	fail_unless( strcmp(input_filename_without_directory, output_filename) ==0 );
}
END_TEST

START_TEST (check_strip_directory_from_filename_with_directory)
{
	char *input_filename_without_directory =  "/some/directory/name/my_file_name.aln";
	char output_filename[30];
	strip_directory_from_filename(input_filename_without_directory, output_filename);
	fail_unless( strcmp("my_file_name.aln", output_filename) ==0 );
}
END_TEST



Suite * snp_sites_suite (void)
{
  Suite *s = suite_create ("Creating_SNP_Sites");

  TCase *tc_alignment_file = tcase_create ("alignment_file");
	tcase_add_test (tc_alignment_file, valid_alignment_with_large_numbers_of_snps);
  tcase_add_test (tc_alignment_file, valid_genome_length);
  tcase_add_test (tc_alignment_file, valid_genome_length_with_multiple_lines_per_sequence);
  tcase_add_test (tc_alignment_file, valid_number_of_sequences_in_file);
  tcase_add_test (tc_alignment_file, valid_number_of_sequences_in_file_with_multiple_lines_per_sequence);
  tcase_add_test (tc_alignment_file, valid_initial_reference_sequence);
  tcase_add_test (tc_alignment_file, number_of_snps_detected_small);
  tcase_add_test (tc_alignment_file, number_of_snps_detected);
	tcase_add_test (tc_alignment_file, number_of_snps_detected_include_gaps);
  tcase_add_test (tc_alignment_file, sample_names_from_alignment_file);
  tcase_add_test (tc_alignment_file, check_strip_directory_from_filename_without_directory);
	tcase_add_test (tc_alignment_file, check_strip_directory_from_filename_with_directory);
  suite_add_tcase (s, tc_alignment_file);
  

  TCase *tc_snp_sites = tcase_create ("snp_sites");
  tcase_add_test (tc_snp_sites, valid_alignment_with_one_line_per_sequence);
  tcase_add_test (tc_snp_sites, valid_alignment_with_multiple_lines_per_sequence);
  tcase_add_test (tc_snp_sites, valid_alignment_with_one_line_per_sequence_gzipped);
	tcase_add_test (tc_snp_sites, valid_alignment_with_n_as_gap);
	tcase_add_test (tc_snp_sites, two_sequences);
  suite_add_tcase (s, tc_snp_sites);

  return s;
}


