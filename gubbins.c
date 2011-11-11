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

#include "seqUtil.h"
#include "Newickform.h"



// get reference sequence from VCF, and store snp locations
// given a sample name extract the sequences from the vcf
// compare two sequences to get pseudo sequnece and fill in with difference from reference sequence

void run_gubbins(char vcf_filename[], char tree_filename[])
{
	//extract_sequences(vcf_filename);
	build_newick_tree(tree_filename);

	

}


void extract_sequences(char vcf_filename[])
{
	FILE *vcf_file_pointer;
	vcf_file_pointer=fopen(vcf_filename, "r");
	
	char * reference_bases;
	int * snp_locations;
	int number_of_snps;
	int number_of_columns;
	int i;
	int reference_column_number;
	
	number_of_columns = get_number_of_columns_from_file(vcf_file_pointer);
	char* column_names[number_of_columns];
	for(i = 0; i < number_of_columns; i++)
	{
		column_names[i] = malloc(100*sizeof(char));
	}
	get_column_names(vcf_file_pointer, column_names, number_of_columns);
	
	number_of_snps  = get_number_of_snps(vcf_file_pointer);
	snp_locations   = (int *)  malloc(number_of_snps*sizeof(int));
	reference_bases = (char *) malloc(number_of_snps*sizeof(char));
	// get reference sequence from VCF
	reference_column_number = column_number_for_column_name(column_names, "REF", number_of_columns);
	get_sequence_from_column_in_vcf(vcf_file_pointer, snp_locations, reference_bases, number_of_snps, reference_column_number);
	
	
	char * seq1;
	char * seq2;
	char * seq3;
	char * ancestor_sequence;
	seq1= (char *) malloc(number_of_snps*sizeof(char));
	seq2= (char *) malloc(number_of_snps*sizeof(char));
	seq3= (char *) malloc(number_of_snps*sizeof(char));
	ancestor_sequence = (char *) malloc(number_of_snps*sizeof(char));
	char * child_sequences[3];
	child_sequences[0] = seq1;
	child_sequences[1] = seq2;
	child_sequences[2] = seq3;
	get_sequence_from_column_in_vcf(vcf_file_pointer, snp_locations, seq1, number_of_snps,  column_number_for_column_name(column_names, "4232_3_2", number_of_columns));
	get_sequence_from_column_in_vcf(vcf_file_pointer, snp_locations, seq2, number_of_snps,  column_number_for_column_name(column_names, "4882_8_8", number_of_columns));
	get_sequence_from_column_in_vcf(vcf_file_pointer, snp_locations, seq3, number_of_snps,  column_number_for_column_name(column_names, "4882_8_11", number_of_columns));
	
	calculate_ancestor_sequence(ancestor_sequence, reference_bases,child_sequences, number_of_snps, 3);
}

// If there are snps between the child sequences, fill in with the reference sequence
void calculate_ancestor_sequence(char * ancestor_sequence, char * reference_sequence, char ** child_sequences, int sequence_length, int number_of_child_sequences)
{
	int base_position;
	int sequence_number;
	int found_snp;
	
	for(base_position = 0; base_position < sequence_length; base_position++)
	{
		found_snp = 0; 
		if(child_sequences[0][base_position] != '.')
		{
		
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
				ancestor_sequence[base_position] = reference_sequence[base_position];
			}
			else
			{
				ancestor_sequence[base_position] = child_sequences[0][base_position];
			}
		}
		else
		{
			ancestor_sequence[base_position] = reference_sequence[base_position];
		}
		
	}
}





















