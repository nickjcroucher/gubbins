
#ifndef _SNP_SITES_H_
#define _SNP_SITES_H_

int build_reference_sequence(char reference_sequence[], FILE * alignment_file_pointer);
int detect_snps(char reference_sequence[], FILE * alignment_file_pointer, int length_of_genome);
void build_snp_locations(int snp_locations[], char reference_sequence[]);
void get_bases_for_each_snp(FILE * alignment_file_pointer, int snp_locations[], char ** bases_for_snps, int length_of_genome, int number_of_snps);
int generate_snp_sites(char filename[]);

#endif