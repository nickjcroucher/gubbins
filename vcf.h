#ifndef _VCF_H_
#define _VCF_H_

void output_vcf_header( FILE * vcf_file_pointer, char ** sequence_names, int number_of_samples);
void create_vcf_file(char filename[],  FILE * alignment_file_pointer, int snp_locations[], int length_of_genome, int number_of_snps);
void output_vcf_snps(FILE * vcf_file_pointer, char ** bases_for_snps, int * snp_locations, int number_of_snps, int number_of_samples);
void output_vcf_row(FILE * vcf_file_pointer, char * bases_for_snp, int snp_location, int number_of_samples);
void output_vcf_row_samples_bases(FILE * vcf_file_pointer, char reference_base, char * bases_for_snp, int number_of_samples);
void alternative_bases(char reference_base, char * bases_for_snp, char alt_bases[], int number_of_samples);
int check_if_char_in_string(char search_string[], char target_char, int search_string_length);


#endif