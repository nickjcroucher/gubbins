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
#include "gubbins.h"
#include "parse_vcf.h"
#include "parse_phylip.h"
#include "alignment_file.h"

#include "seqUtil.h"
#include "Newickform.h"



// get reference sequence from VCF, and store snp locations
// given a sample name extract the sequences from the vcf
// compare two sequences to get pseudo sequnece and fill in with difference from reference sequence

void run_gubbins(char vcf_filename[], char tree_filename[], char phylip_filename[],char multi_fasta_filename[])
{
	load_sequences_from_phylib_file(phylip_filename);
	extract_sequences(vcf_filename, tree_filename, multi_fasta_filename);
}


int get_length_of_genome_from_alignment_file(char multi_fasta_filename[])
{
	FILE *alignment_file_pointer;
	alignment_file_pointer=fopen(multi_fasta_filename, "r");
	return genome_length(alignment_file_pointer);	
}

void extract_sequences(char vcf_filename[], char tree_filename[],char multi_fasta_filename[])
{
	FILE *vcf_file_pointer;
	vcf_file_pointer=fopen(vcf_filename, "r");
	
	char * reference_bases;
	
	int number_of_snps;
	int number_of_columns;
	int i;
	int reference_column_number;
	int length_of_original_genome;
	length_of_original_genome = get_length_of_genome_from_alignment_file(multi_fasta_filename);
	
	number_of_columns = get_number_of_columns_from_file(vcf_file_pointer);
	char* column_names[number_of_columns];
	for(i = 0; i < number_of_columns; i++)
	{
		column_names[i] = malloc(MAX_SAMPLE_NAME_SIZE*sizeof(char));
	}
	get_column_names(vcf_file_pointer, column_names, number_of_columns);
	
	number_of_snps  = get_number_of_snps(vcf_file_pointer);
	reference_bases = (char *) malloc(number_of_snps*sizeof(char));
	
	int snp_locations[number_of_snps];
	
	get_integers_from_column_in_vcf(vcf_file_pointer, snp_locations, number_of_snps, column_number_for_column_name(column_names, "POS", number_of_columns));
	
	// get reference sequence from VCF
	reference_column_number = column_number_for_column_name(column_names, "REF", number_of_columns);
	get_sequence_from_column_in_vcf(vcf_file_pointer, reference_bases, number_of_snps, reference_column_number);
	
	build_newick_tree(tree_filename, vcf_file_pointer,snp_locations, number_of_snps, column_names, number_of_columns, reference_bases,length_of_original_genome);
}

// If there are snps between the child sequences, fill in with the reference sequence
char *calculate_ancestor_sequence(char * ancestor_sequence, char ** child_sequences, int sequence_length, int number_of_child_sequences)
{
	int base_position;
	int sequence_number;
	int found_snp;
	strcpy(ancestor_sequence,"");
		   
		   
	for(base_position = 0; base_position < sequence_length; base_position++)
	{
		if(child_sequences[0][base_position] == '\0')
		{
			ancestor_sequence[base_position] == '\0';
			break;
		}
		
		found_snp = 0;
		
		for(sequence_number = 1; sequence_number < number_of_child_sequences; sequence_number++)	
		{
			if(child_sequences[0][base_position] != child_sequences[sequence_number][base_position])
			{
				found_snp = 1;
				break;
			}
		}
		if(found_snp == 1)
		{
			ancestor_sequence[base_position] = '.';
		}
		else
		{
			ancestor_sequence[base_position] = child_sequences[0][base_position];
		}
		
	}
	
	if(ancestor_sequence[sequence_length] != '\0')
	{
		ancestor_sequence[sequence_length] = '\0';
	}
	
	
	return ancestor_sequence;
}








