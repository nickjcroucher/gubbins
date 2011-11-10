#ifndef _PARSE_VCF_H_
#define _PARSE_VCF_H_

void get_sequence_from_column_in_vcf(FILE * vcf_file_pointer, int snp_locations[], char * sequence_bases, int number_of_snps, int column_number);
int get_number_of_snps(FILE * vcf_file_pointer);
int get_number_of_samples(FILE * vcf_file_pointer);
char *split_string_and_return_specific_index(char * result, char *input_string, int token_index);
int get_number_of_columns(char * column_header);
int get_number_of_columns_from_file(FILE * vcf_file_pointer);
void get_column_names(FILE * vcf_file_pointer, char ** column_names, int number_of_columns);
int column_number_for_column_name(char ** column_names, char * column_name, int number_of_columns);

#endif