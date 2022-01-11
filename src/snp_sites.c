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
#include "phylip_of_snp_sites.h"
#include "parse_phylip.h"
#include "string_cat.h"
#include "fasta_of_snp_sites.h"


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


int generate_snp_sites(char filename[],  int exclude_gaps, char suffix[])
{
	int length_of_genome;
	char * reference_sequence;
	int number_of_snps;
	int * snp_locations;
	int number_of_samples;
	int i;
	
	length_of_genome = genome_length(filename);
	reference_sequence = (char *) calloc((length_of_genome+1),sizeof(char));
	
	build_reference_sequence(reference_sequence,filename);
	number_of_snps = detect_snps(reference_sequence, filename, length_of_genome, exclude_gaps);
	
	snp_locations = (int *) calloc((number_of_snps+1),sizeof(int));
	build_snp_locations(snp_locations, reference_sequence);
	free(reference_sequence);
	
	number_of_samples = number_of_sequences_in_file(filename);
	
	// Find out the names of the sequences
	char* sequence_names[number_of_samples];
	sequence_names[number_of_samples-1] = "\0";
	for(i = 0; i < number_of_samples; i++)
	{
		sequence_names[i] = calloc(MAX_SAMPLE_NAME_SIZE,sizeof(char));
	}
	
	get_sample_names_for_header(filename, sequence_names, number_of_samples);
	
	int internal_nodes[number_of_samples];
	int a = 0;
	for(a =0; a < number_of_samples; a++)
	{
		internal_nodes[a] = 0;
	}
	
	char** bases_for_snps = malloc(number_of_snps * sizeof(char *));
	
	for(i = 0; i < number_of_snps; i++)
	{
		bases_for_snps[i] = calloc((number_of_samples+1),sizeof(char));
	}
	
	get_bases_for_each_snp(filename, snp_locations, bases_for_snps, length_of_genome, number_of_snps);
	
    char filename_without_directory[MAX_FILENAME_SIZE];
    strip_directory_from_filename(filename, filename_without_directory);
	
	concat_strings_created_with_malloc(filename_without_directory,suffix);
	
	create_vcf_file(filename_without_directory, snp_locations, number_of_snps, bases_for_snps, sequence_names, number_of_samples,internal_nodes,1,length_of_genome);
	create_phylip_of_snp_sites(filename_without_directory, number_of_snps, bases_for_snps, sequence_names, number_of_samples,internal_nodes);
	create_fasta_of_snp_sites(filename_without_directory, number_of_snps, bases_for_snps, sequence_names, number_of_samples,internal_nodes);

	free(snp_locations);
	return 1;
}

// Inefficient
void strip_directory_from_filename(char * input_filename, char * output_filename)
{
  int i;
  int end_index = 0;
  int last_forward_slash_index = -1;
  for(i = 0; i< MAX_FILENAME_SIZE; i++)
  {
    if(input_filename[i] == '/')
    {
      last_forward_slash_index = i;
    }
    
    if(input_filename[i] == '\0' || input_filename[i] == '\n')
    {
      end_index = i;
      break;
    }
  }
  
  int current_index = 0;
  for(i = last_forward_slash_index+1; i< end_index; i++)
  {
    output_filename[current_index] = input_filename[i];
    current_index++;
  }
  output_filename[current_index] = '\0';
}

// return new number of snps
int refilter_existing_snps(char * reference_bases, int number_of_snps, int * snp_locations, int * filtered_snp_locations, int internal_nodes[])
{
	// go through each snp column and check to see if there is still variation
	int i;
	int number_of_filtered_snps = number_of_snps;
	for(i = 0; i < number_of_snps; i++)
	{
		if( does_column_contain_snps(i, reference_bases[i]) == 0)
		{
			snp_locations[i] = -1;
			reference_bases[i] = '*';
			
			number_of_filtered_snps--;
		}
	}
	
	remove_filtered_snp_locations(filtered_snp_locations, snp_locations, number_of_snps);
	return number_of_filtered_snps;
}

void remove_filtered_snp_locations(int * filtered_snp_locations, int * snp_locations, int number_of_snps)
{
	int i;
	int filtered_counter=0;
	for(i = 0; i< number_of_snps; i++)
	{
		if(snp_locations[i] != -1)
		{
			filtered_snp_locations[filtered_counter] = snp_locations[i];
			filtered_counter++;
		}
	}
}













