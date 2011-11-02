#ifndef _VCF_H_
#define _VCF_H_

void output_vcf_header( FILE * vcf_file_pointer);
void create_vcf_file(char filename[],  FILE * alignment_file_pointer, int snp_locations[], int length_of_genome, int number_of_snps);
void output_vcf_snps(FILE * vcf_file_pointer, char ** bases_for_snps, int snp_locations[]);
void output_vcf_row(FILE * vcf_file_pointer, char * bases_for_snp, int snp_location);
void output_vcf_row_samples_bases(FILE * vcf_file_pointer, char reference_base, char * bases_for_snp);

#endif