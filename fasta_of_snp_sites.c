#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "fasta_of_snp_sites.h"

void create_fasta_of_snp_sites(char filename[], int number_of_snps, char ** bases_for_snps, char ** sequence_names, int number_of_samples)
{
	FILE *fasta_file_pointer;
	int sample_counter;
	int snp_counter; 
	
	fasta_file_pointer = fopen(strcat(filename,".snps.aln"), "w");
	
	for(sample_counter=0; sample_counter< number_of_samples; sample_counter++)
	{
		fprintf( fasta_file_pointer, ">%s\n", sequence_names[sample_counter]);
		for(snp_counter=0; snp_counter< number_of_snps; snp_counter++)
		{
			fprintf( fasta_file_pointer, "%c", bases_for_snps[snp_counter][sample_counter]);
		}
		fprintf( fasta_file_pointer, "\n");
	}
}
