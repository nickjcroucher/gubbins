/*
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "vcf.h"
#include "alignment_file.h"
#include "snp_sites.h"


int build_reference_sequence(char reference_sequence[], FILE * alignment_file_pointer)
{
	int i;
	
	rewind(alignment_file_pointer);
    advance_to_sequence(alignment_file_pointer);
	
    read_line(reference_sequence, alignment_file_pointer);
	
	for(i = 0; reference_sequence[ i ]; i++)
	{
		reference_sequence[i] = toupper(reference_sequence[i]);
	}
	
	return 1;
}


int detect_snps(char reference_sequence[], FILE * alignment_file_pointer, int length_of_genome)
{
	char * comparison_sequence;
	int i;
	int number_of_snps = 0;
	
	do{
		comparison_sequence = (char *) malloc(length_of_genome*sizeof(char));
		advance_to_sequence(alignment_file_pointer);
		read_line(comparison_sequence, alignment_file_pointer);
		
		if(comparison_sequence[0] == '\0')
		{
			break;
		}
		
		// Set the reference base to * if 
		for(i = 0; reference_sequence[ i ]; i++)
		{
			if(reference_sequence[i] != '*' && comparison_sequence[i] != '-' && reference_sequence[i] != toupper(comparison_sequence[i]))
			{
				reference_sequence[i] = '*';
				number_of_snps++;
			}
		}
	}while(comparison_sequence[0] != '\0');
	
	free(comparison_sequence);
	return number_of_snps;
}


void build_snp_locations(int snp_locations[], char reference_sequence[])
{
	int i;
	int snp_counter = 0;
	
	for(i = 0; reference_sequence[i]; i++)
    {
		if(reference_sequence[i] == '*')
		{
			snp_locations[snp_counter] = i;
			snp_counter++;
		}
	}
}


void get_bases_for_each_snp(FILE * alignment_file_pointer, int snp_locations[], char ** bases_for_snps, int length_of_genome)
{
	int i;
	int sequence_number = 0;
	char * comparison_sequence;
	rewind(alignment_file_pointer);
	
	do{
		comparison_sequence = (char *) malloc(length_of_genome*sizeof(char));
		
		advance_to_sequence(alignment_file_pointer);
		read_line(comparison_sequence, alignment_file_pointer);
		
		if(comparison_sequence[0] == '\0')
		{
			break;
		}
		
		for(i = 0; snp_locations[i]; i++)
		{
			bases_for_snps[i][sequence_number] = toupper(((char *) comparison_sequence)[snp_locations[i]]);
		}
		
		sequence_number++;
	}while(comparison_sequence[0] != '\0');
	
	free(comparison_sequence);
}


int generate_snp_sites(char filename[])
{
	FILE *alignment_file_pointer;
	int length_of_genome;
	char * reference_sequence;
	int number_of_snps;
	int * snp_locations;
	
	alignment_file_pointer=fopen(filename, "r");
	
	if(validate_alignment_file(alignment_file_pointer) == 0)
	{
		return 0;
	}
	
	length_of_genome = genome_length(alignment_file_pointer);
	reference_sequence = (char *) malloc(length_of_genome*sizeof(char));
	
	build_reference_sequence(reference_sequence,alignment_file_pointer);
	number_of_snps = detect_snps(reference_sequence, alignment_file_pointer, length_of_genome);
	
	snp_locations = (int *) malloc(number_of_snps*sizeof(int));
	build_snp_locations(snp_locations, reference_sequence);
	free(reference_sequence);
	 
	//printf("Number of SNPs: %d\n", number_of_snps);
	create_vcf_file(filename, alignment_file_pointer, snp_locations,length_of_genome, number_of_snps);
	
	free(snp_locations);
	return 1;
}

