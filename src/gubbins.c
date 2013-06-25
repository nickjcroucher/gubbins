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
#include "snp_sites.h"
#include "vcf.h"
#include "phylip_of_snp_sites.h"
#include "tree_scaling.h"
#include "seqUtil.h"
#include "Newickform.h"



// get reference sequence from VCF, and store snp locations
// given a sample name extract the sequences from the vcf
// compare two sequences to get pseudo sequnece and fill in with difference from reference sequence

void run_gubbins(char vcf_filename[], char tree_filename[],char multi_fasta_filename[], int min_snps, char original_multi_fasta_filename[])
{
	load_sequences_from_multifasta_file(multi_fasta_filename);
	extract_sequences(vcf_filename, tree_filename, multi_fasta_filename,min_snps,original_multi_fasta_filename);
	create_tree_statistics_file(tree_filename,get_sample_statistics(),number_of_samples_from_parse_phylip());
	freeup_memory();

}


void extract_sequences(char vcf_filename[], char tree_filename[],char multi_fasta_filename[],int min_snps, char original_multi_fasta_filename[])
{
	FILE *vcf_file_pointer;
	vcf_file_pointer=fopen(vcf_filename, "r");
	
	newick_node* root_node;
	int number_of_snps;
	int number_of_columns;
	int i;
	int length_of_original_genome;
	length_of_original_genome = genome_length(original_multi_fasta_filename);	
	
	number_of_columns = get_number_of_columns_from_file(vcf_file_pointer);
	char* column_names[number_of_columns];
	for(i = 0; i < number_of_columns; i++)
	{
		column_names[i] = malloc(MAX_SAMPLE_NAME_SIZE*sizeof(char));
	}
	get_column_names(vcf_file_pointer, column_names, number_of_columns);
	
	number_of_snps  = number_of_snps_in_phylip();
	
	int snp_locations[number_of_snps];
	
	get_integers_from_column_in_vcf(vcf_file_pointer, snp_locations, number_of_snps, column_number_for_column_name(column_names, "POS", number_of_columns));

	root_node = build_newick_tree(tree_filename, vcf_file_pointer,snp_locations, number_of_snps, column_names, number_of_columns, length_of_original_genome,min_snps);

	int filtered_snp_locations[number_of_snps];
	int number_of_filtered_snps;
	int number_of_samples = number_of_samples_from_parse_phylip();

	char * sample_names[number_of_samples];
  get_sample_names_from_parse_phylip(sample_names);

  char * reference_sequence_bases;
  reference_sequence_bases = (char *) malloc((number_of_snps+1)*sizeof(char));

	get_sequence_for_sample_name(reference_sequence_bases, sample_names[0]);
	int internal_nodes[number_of_samples];
	int a = 0;
	for(a =0; a < number_of_samples; a++)
	{
		internal_nodes[a] = get_internal_node(a);
	}

	number_of_filtered_snps = refilter_existing_snps(reference_sequence_bases, number_of_snps, snp_locations, filtered_snp_locations,internal_nodes);
	char * filtered_bases_for_snps[number_of_filtered_snps];

	filter_sequence_bases_and_rotate(reference_sequence_bases, filtered_bases_for_snps, number_of_filtered_snps);
	free(reference_sequence_bases);
	create_phylip_of_snp_sites(tree_filename, number_of_filtered_snps, filtered_bases_for_snps, sample_names, number_of_samples,internal_nodes);
	create_vcf_file(tree_filename, filtered_snp_locations, number_of_filtered_snps, filtered_bases_for_snps, sample_names, number_of_samples,internal_nodes,0);
	create_fasta_of_snp_sites(tree_filename, number_of_filtered_snps, filtered_bases_for_snps, sample_names, number_of_samples,internal_nodes);
	
	// Create an new tree with updated distances
	scale_branch_distances(root_node, number_of_filtered_snps);

	FILE *output_tree_pointer;
	output_tree_pointer=fopen(tree_filename, "w");
	print_tree(root_node,output_tree_pointer);
	fprintf(output_tree_pointer,";");
	fflush(output_tree_pointer);
	fclose(output_tree_pointer);
	cleanup_node_memory(root_node);

}


char find_first_real_base(int base_position,  int number_of_child_sequences, char ** child_sequences)
{
	int i = 0; 
	for(i =0; i< number_of_child_sequences; i++)
	{
	  if(child_sequences[i][base_position] != 'N' && child_sequences[i][base_position] != '-' && child_sequences[i][base_position] != '.')
		{
			return child_sequences[i][base_position];
		}
	}
	return child_sequences[0][base_position];
}








