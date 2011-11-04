#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "phylib_of_snp_sites.h"


void create_phylib_of_snp_sites(char filename[], int number_of_snps, char ** bases_for_snps, char ** sequence_names, int number_of_samples)
{
	FILE *fasta_file_pointer;
	int sample_counter;
	int snp_counter;
	char * base_filename;
	
	base_filename = (char *) malloc(256*sizeof(char));
	strcpy(base_filename, filename);
	
	fasta_file_pointer = fopen(strcat(base_filename,".phylip"), "w");
	
	fprintf( fasta_file_pointer, "%d %d\n", number_of_samples, number_of_snps);
	
	for(sample_counter=0; sample_counter< number_of_samples; sample_counter++)
	{
		// sequence_name can be more than 10 (relaxed phylib format) and contain [\w\s]
		//TODO check for illegal characters [^\w\s]
		fprintf( fasta_file_pointer, "%s\t", sequence_names[sample_counter]);
		
		for(snp_counter=0; snp_counter< number_of_snps; snp_counter++)
		{
			fprintf( fasta_file_pointer, "%c", bases_for_snps[snp_counter][sample_counter]);
		}
		fprintf( fasta_file_pointer, "\n");
	}
}

