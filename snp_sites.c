/*
 *  Wellcome Trust Sanger Institute
 *  Copyright (C) 2011  Wellcome Trust Sanger Institute
 *  
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "vcf.h"
#include "alignment_file.h"
#include "snp_sites.h"
#include "fasta_of_snp_sites.h"
#include "phylib_of_snp_sites.h"


int build_reference_sequence(char reference_sequence[], FILE * alignment_file_pointer)
{
	int i;
	
	rewind(alignment_file_pointer);
    advance_to_sequence(alignment_file_pointer);
	
    read_line(reference_sequence, alignment_file_pointer);
	
	for(i = 0; reference_sequence[i]; i++)
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
		for(i = 0; i < length_of_genome; i++)
		{
			// If there is an indel in the reference sequence, replace with the first proper base you find
			if(reference_sequence[i] == '-' && comparison_sequence[i] != '-' )
			{
				reference_sequence[i] = toupper(comparison_sequence[i]);
			}
			
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


void get_bases_for_each_snp(FILE * alignment_file_pointer, int snp_locations[], char ** bases_for_snps, int length_of_genome, int number_of_snps)
{
	int i;
	int sequence_number = 0;
	char * comparison_sequence;
	rewind(alignment_file_pointer);
	
	// initialise the strings in the array
	for(i = 0; i < number_of_snps; i++)
	{
		strcpy(bases_for_snps[i], "");
	}

	
	do{
		comparison_sequence = (char *) malloc(length_of_genome*sizeof(char));
		
		advance_to_sequence(alignment_file_pointer);
		read_line(comparison_sequence, alignment_file_pointer);
		
		if(comparison_sequence[0] == '\0')
		{
			break;
		}
		
		for(i = 0; i< number_of_snps; i++)
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
	int number_of_samples;
	int i;
	
	
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
	 
	rewind(alignment_file_pointer);
	number_of_samples = count_lines_in_file(alignment_file_pointer)/2;
	
	// Find out the names of the sequences
	char* sequence_names[number_of_samples];
	sequence_names[number_of_samples-1] = '\0';
	for(i = 0; i < number_of_samples; i++)
	{
		sequence_names[i] = malloc(21*sizeof(char));
		strcpy(sequence_names[i],"");
	}
	
	get_sample_names_for_header(alignment_file_pointer, sequence_names, number_of_samples);
	
	char* bases_for_snps[number_of_snps];
	
	for(i = 0; i < number_of_snps; i++)
	{
		bases_for_snps[i] = malloc(number_of_samples*sizeof(char));
	}
	
	get_bases_for_each_snp(alignment_file_pointer, snp_locations, bases_for_snps, length_of_genome, number_of_snps);
	
	create_vcf_file(filename, snp_locations, number_of_snps, bases_for_snps, sequence_names, number_of_samples);
	
	create_fasta_of_snp_sites(filename, number_of_snps, bases_for_snps, sequence_names, number_of_samples);
	
	create_phylib_of_snp_sites(filename, number_of_snps, bases_for_snps, sequence_names, number_of_samples);
	
	
	free(snp_locations);
	return 1;
}

