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
#include <math.h>
#include "seqUtil.h"
#include "Newickform.h"
#include "branch_sequences.h"
#include "gubbins.h"
#include "parse_vcf.h"
#include "parse_phylip.h"
#include "snp_searching.h"
#include "block_tab_file.h"
#include "gff_file.h"

int node_counter = 0;


// Order is not preserved.
int copy_and_concat_integer_arrays(int * array_1, int array_1_size, int * array_2, int array_2_size, int * output_array)
{
	int array_1_counter=0;
	int array_2_counter=0;
	
	for(array_1_counter = 0; array_1_counter< array_1_size; array_1_counter++)
	{
		output_array[array_1_counter] = array_1[array_1_counter];
	}
	
	for(array_2_counter = 0; array_2_counter < array_2_size; array_2_counter++)
	{
		output_array[array_2_counter+array_1_size] = array_2[array_2_counter];
	}
	return array_1_size+array_2_size;
}

int copy_and_concat_2d_integer_arrays(int ** array_1, int array_1_size, int ** array_2, int array_2_size, int ** output_array)
{
	int array_1_counter=0;
	int array_2_counter=0;
	
	for(array_1_counter = 0; array_1_counter< array_1_size; array_1_counter++)
	{
		output_array[0][array_1_counter] = array_1[0][array_1_counter];
		output_array[1][array_1_counter] = array_1[1][array_1_counter];
	}
	
	for(array_2_counter = 0; array_2_counter < array_2_size; array_2_counter++)
	{
		output_array[0][array_2_counter+array_1_size] = array_2[0][array_2_counter];
		output_array[1][array_2_counter+array_1_size] = array_2[1][array_2_counter];
	}
	return array_1_size+array_2_size;
}


// Go through the tree and build up the recombinations list from root to branch. Print out each sample name and a list of recombinations
void fill_in_recombinations_with_gaps(newick_node *root, int * parent_recombinations, int parent_num_recombinations, int current_total_snps,int num_blocks, int ** current_block_coordinates,int length_of_original_genome,int * snp_locations )
{
	newick_child *child;
	int * current_recombinations;
	int num_current_recombinations = 0 ;
	
	current_recombinations = (int *) malloc((root->num_recombinations+1+parent_num_recombinations)*sizeof(int));
	num_current_recombinations = copy_and_concat_integer_arrays(root->recombinations, root->num_recombinations,parent_recombinations, parent_num_recombinations, current_recombinations);
	
 	// overwrite the bases of snps with N's
 	int i;
 	int sequence_index;
 	sequence_index = find_sequence_index_from_sample_name(root->taxon);
 	
 	set_number_of_recombinations_for_sample(root->taxon,root->num_recombinations);
 	set_number_of_snps_for_sample(root->taxon,root->number_of_snps);
 	set_number_of_blocks_for_sample(root->taxon, root->number_of_blocks);
 	
	char * child_sequence = (char *) malloc((length_of_original_genome +1)*sizeof(char));
	strcpy(child_sequence,"");
	get_sequence_for_sample_name(child_sequence, root->taxon);

  // This should take the merged block coordinates?
 	set_number_of_bases_in_recombinations(root->taxon, calculate_number_of_bases_in_recombations_excluding_gaps(root->block_coordinates, root->number_of_blocks, child_sequence, snp_locations,current_total_snps));
	free(child_sequence); 	

 	for(i = 0; i < num_current_recombinations; i++)
 	{
 		int snp_index;
 		snp_index = current_recombinations[i];
 		
 		update_sequence_base('N', sequence_index, snp_index);
 	}

	if (root->childNum > 0)
	{
		child = root->child;
		set_internal_node(1,sequence_index);

		while (child != NULL)
		{
			// recursion
			int ** merged_block_coordinates;
			merged_block_coordinates = (int **) malloc(3*sizeof(int *));
			merged_block_coordinates[0] = (int*) malloc((num_blocks + root->number_of_blocks+1)*sizeof(int ));
			merged_block_coordinates[1] = (int*) malloc((num_blocks + root->number_of_blocks+1)*sizeof(int ));
			copy_and_concat_2d_integer_arrays(current_block_coordinates,num_blocks,root->block_coordinates, root->number_of_blocks,merged_block_coordinates );
			fill_in_recombinations_with_gaps(child->node, current_recombinations, num_current_recombinations,(current_total_snps + root->number_of_snps),(num_blocks + root->number_of_blocks),merged_block_coordinates,length_of_original_genome, snp_locations );
			child = child->next;
			
			free(merged_block_coordinates[0]);
			free(merged_block_coordinates[1]);
			free(merged_block_coordinates);
		}
	}
	else
	{
	set_internal_node(0,sequence_index);	
	}
	free(current_recombinations);
}

int calculate_number_of_bases_in_recombations_excluding_gaps(int ** block_coordinates, int num_blocks,char * child_sequence, int * snp_locations,int length_of_original_genome)
{
	int total_bases = 0;
	int current_block = 1;
	int start_block = 0;

	for(start_block = 0; start_block < num_blocks; start_block++)
	{
		if(block_coordinates[0][start_block] == -1 || block_coordinates[1][start_block] == -1)
		{
			continue;
		}
		
		for(current_block = 0 ; current_block < num_blocks; current_block++)
		{ 
			if(current_block == start_block)
			{
				continue;	
			}
			
			if(block_coordinates[0][current_block] == -1 || block_coordinates[1][current_block] == -1)
			{
				continue;
			}
			
			
			int found_overlap = 0;
		  if(block_coordinates[0][start_block] >=  block_coordinates[0][current_block] && block_coordinates[0][start_block] <= block_coordinates[1][current_block] )
		  {
				block_coordinates[0][start_block] = block_coordinates[0][current_block];
				found_overlap = 1;
			}
			
			if(block_coordinates[1][start_block] >=  block_coordinates[0][current_block]  && block_coordinates[1][start_block] <= block_coordinates[1][current_block])
		  {
				block_coordinates[1][start_block] = block_coordinates[1][current_block];
				found_overlap = 1;
			}
			
			if(found_overlap == 1)
			{
				block_coordinates[0][current_block] = -1;
				block_coordinates[1][current_block] = -1;
			}
		}	
		
	}
	for(start_block = 0; start_block < num_blocks; start_block++)
	{
		if(block_coordinates[0][start_block] == -1 || block_coordinates[1][start_block] == -1)
		{
			continue;
		}
		
		
	  total_bases += calculate_block_size_without_gaps(child_sequence,  snp_locations, block_coordinates[0][start_block], block_coordinates[1][start_block], length_of_original_genome);
  }
	
	return total_bases;
}

void carry_unambiguous_gaps_up_tree(newick_node *root)
{
	if(root->childNum > 0)
	{
		newick_child *child;
		int parent_sequence_index =  find_sequence_index_from_sample_name(root->taxon);
		
		child = root->child;
		int child_sequence_indices[root->childNum];
		int child_counter = 0;
		while (child != NULL)
		{
			child_sequence_indices[child_counter] = find_sequence_index_from_sample_name(child->node->taxon);
			carry_unambiguous_gaps_up_tree(child->node);
			child = child->next;
			child_counter++;
		}
		
		// compare the parent sequence to the each child sequence and update the gaps
		fill_in_unambiguous_gaps_in_parent_from_children(parent_sequence_index, child_sequence_indices,child_counter);
		fill_in_unambiguous_bases_in_parent_from_children_where_parent_has_a_gap(parent_sequence_index, child_sequence_indices,child_counter);
	}
}

char *generate_branch_sequences(newick_node *root, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns, char * leaf_sequence, int length_of_original_genome, FILE * block_file_pointer, FILE * gff_file_pointer,int min_snps,FILE * branch_snps_file_pointer)
{
	newick_child *child;
	int child_counter = 0;
	int current_branch =0;
	int branch_genome_size = 0;
	int number_of_branch_snps=0;
  root->current_node_id = ++node_counter;
	
	if (root->childNum == 0)
	{
		leaf_sequence = (char *) malloc((number_of_snps +1)*sizeof(char));
		strcpy(leaf_sequence, "");
		get_sequence_for_sample_name(leaf_sequence, root->taxon);
		
		root->taxon_names = (char *) malloc(MAX_SAMPLE_NAME_SIZE*sizeof(char));
    strcpy(root->taxon_names, root->taxon);

    // Save some statistics about the sequence
		branch_genome_size = calculate_size_of_genome_without_gaps(leaf_sequence, 0,number_of_snps, length_of_original_genome);
		set_genome_length_without_gaps_for_sample(root->taxon,branch_genome_size);
		int number_of_gaps = length_of_original_genome-branch_genome_size;
		
		return leaf_sequence;
	}
	else
	{
		child = root->child;
		char * child_sequences[root->childNum];
		newick_node * child_nodes[root->childNum];
		root->taxon_names = (char *) malloc(MAX_SAMPLE_NAME_SIZE*number_of_columns*sizeof(char));
    strcpy(root->taxon_names, "");

		// generate pointers for each child seuqn
		

		
		while (child != NULL)
		{
			// recursion
			child_sequences[child_counter] = generate_branch_sequences(child->node, vcf_file_pointer, snp_locations, number_of_snps, column_names, number_of_columns,  child_sequences[child_counter],length_of_original_genome, block_file_pointer,gff_file_pointer,min_snps,branch_snps_file_pointer);
			child_nodes[child_counter] = child->node;
      strcat(root->taxon_names, " ");
      strcat(root->taxon_names, child_nodes[child_counter]->taxon_names);
			
			child_counter++;
			child = child->next;
		}
		
		// For all bases update the parent sequence with N if all child sequences.
		
		
		leaf_sequence = (char *) malloc((number_of_snps +1)*sizeof(char));
		strcpy(leaf_sequence, "");
		// All child sequneces should be available use them to find the ancestor sequence
		get_sequence_for_sample_name(leaf_sequence, root->taxon);
		
		branch_genome_size = calculate_size_of_genome_without_gaps(leaf_sequence, 0,number_of_snps, length_of_original_genome);
		set_genome_length_without_gaps_for_sample(root->taxon,branch_genome_size);
		
		int * branches_snp_sites[root->childNum];
		
		for(current_branch = 0 ; current_branch< (root->childNum); current_branch++)
		{
			branches_snp_sites[current_branch] = (int *) malloc((number_of_snps +1)*sizeof(int));
			char * branch_snp_sequence;
			char * branch_snp_ancestor_sequence;
			branch_snp_sequence = (char *) malloc((number_of_snps +1)*sizeof(char));
			branch_snp_ancestor_sequence = (char *) malloc((number_of_snps +1)*sizeof(char));
			
			branch_genome_size = calculate_size_of_genome_without_gaps(child_sequences[current_branch], 0,number_of_snps, length_of_original_genome);
			number_of_branch_snps = calculate_number_of_snps_excluding_gaps(leaf_sequence, child_sequences[current_branch], number_of_snps, branches_snp_sites[current_branch], snp_locations,branch_snp_sequence,branch_snp_ancestor_sequence);
			
			print_branch_snp_details(branch_snps_file_pointer, child_nodes[current_branch]->taxon,root->taxon, branches_snp_sites[current_branch], number_of_branch_snps, branch_snp_sequence, branch_snp_ancestor_sequence,child_nodes[current_branch]->taxon_names);
			
			get_likelihood_for_windows(child_sequences[current_branch], number_of_snps, branches_snp_sites[current_branch], branch_genome_size, number_of_branch_snps,snp_locations, child_nodes[current_branch], block_file_pointer, root, branch_snp_sequence,gff_file_pointer,min_snps,length_of_original_genome);
			free(branch_snp_sequence);
			free(branch_snp_ancestor_sequence);
			free(child_sequences[current_branch]);
			free(branches_snp_sites[current_branch]);
		}
		
		return leaf_sequence;
	}
}


// Windows need to be of a fixed size
// calculate window size
// starting at coord of first snp, count number of snps which fall into window
// if region is blank, move on
int calculate_window_size(int branch_genome_size, int number_of_branch_snps)
{
	int window_size = 0;
	if(number_of_branch_snps == 0)
	{
		return MIN_WINDOW_SIZE;
	}
	
	window_size = (int) ((branch_genome_size*1.0)/(number_of_branch_snps*1.0/WINDOW_SNP_MODE_TARGET));
	
	if(window_size < MIN_WINDOW_SIZE)
	{
		return MIN_WINDOW_SIZE;
	}
	else if(window_size > MAX_WINDOW_SIZE)
	{
		return 	MAX_WINDOW_SIZE;
	}

	return window_size;
}


void get_likelihood_for_windows(char * child_sequence, int length_of_sequence, int * snp_site_coords, int branch_genome_size, int number_of_branch_snps, int * snp_locations, newick_node * current_node, FILE * block_file_pointer, newick_node *root, char * branch_snp_sequence, FILE * gff_file_pointer, int min_snps, int length_of_original_genome)
{
	int i = 0;
	int window_size = 0;
	int window_start_coordinate = 0;
	int window_end_coordinate = 0;
	int number_of_snps_in_block = 0;
	int block_genome_size_without_gaps = 0;
	double branch_snp_density = 0.0;
	double block_snp_density = 0.0;
	int number_of_blocks = 0 ;
	int original_branch_genome_size = branch_genome_size;

	// place to store coordinates of recombinations snps
	current_node->recombinations = (int *) malloc((number_of_branch_snps+1)*sizeof(int));
	
	int number_of_windows = (branch_genome_size/MIN_WINDOW_SIZE) + 1;
	int * block_coordinates[4];
	block_coordinates[0] = (int *) malloc((number_of_windows+1)*sizeof(int));
	block_coordinates[1] = (int *) malloc((number_of_windows+1)*sizeof(int));
	block_coordinates[2] = (int *) malloc((number_of_windows+1)*sizeof(int));
	block_coordinates[3] = (int *) malloc((number_of_windows+1)*sizeof(int));
	
	double * block_likelihoods;	
	block_likelihoods = (double *) malloc((number_of_windows+1)*sizeof(double));

	while(number_of_branch_snps > min_snps)
	{
		if(number_of_branch_snps <= min_snps)
		{
			free(block_coordinates[0]) ;
			free(block_coordinates[1]) ;
			free(block_coordinates[2]) ;
			free(block_coordinates[3]) ;
			free(block_likelihoods);
			return;
		}
		branch_snp_density = snp_density(branch_genome_size, number_of_branch_snps);
		
		window_size = calculate_window_size(branch_genome_size, number_of_branch_snps);

    int cutoff = calculate_cutoff(branch_genome_size, window_size, number_of_branch_snps);

		number_of_blocks = get_blocks(block_coordinates, length_of_original_genome, snp_site_coords, number_of_branch_snps,  window_size, cutoff);


		for(i = 0; i < number_of_blocks; i++)
		{
			number_of_snps_in_block = find_number_of_snps_in_block(block_coordinates[0][i], block_coordinates[1][i], snp_site_coords, branch_snp_sequence, number_of_branch_snps);
			block_genome_size_without_gaps = calculate_block_size_without_gaps(child_sequence, snp_locations, block_coordinates[0][i], block_coordinates[1][i], length_of_sequence);

			// minimum number of snps to be statistically significant in block
			if(number_of_snps_in_block <= min_snps)
			{
				block_coordinates[0][i] = -1;
				block_coordinates[1][i] = -1;
				continue;
			}
			
			block_snp_density = snp_density(block_genome_size_without_gaps, number_of_snps_in_block);
			// region with low number of snps so skip over
			if(block_snp_density <= branch_snp_density)
			{
				block_coordinates[0][i] = -1;
				block_coordinates[1][i] = -1;
				continue;	
			}
			
			block_likelihoods[i] = get_block_likelihood(branch_genome_size, number_of_branch_snps, block_genome_size_without_gaps, number_of_snps_in_block);
			block_coordinates[2][i] = (int) block_likelihoods[i];
			block_coordinates[3][i] = block_genome_size_without_gaps;
		}

		move_blocks_inwards_while_likelihood_improves(number_of_blocks,block_coordinates, min_snps, snp_site_coords, number_of_branch_snps, branch_snp_sequence, snp_locations, branch_genome_size, child_sequence, length_of_sequence,block_likelihoods,cutoff);

		int * candidate_blocks[4];
		candidate_blocks[0] = (int *) malloc((number_of_blocks+1)*sizeof(int));
		candidate_blocks[1] = (int *) malloc((number_of_blocks+1)*sizeof(int));
		candidate_blocks[2] = (int *) malloc((number_of_blocks+1)*sizeof(int));
		candidate_blocks[3] = (int *) malloc((number_of_blocks+1)*sizeof(int));
		
		double * candidate_block_likelihoods;
		candidate_block_likelihoods = (double *) malloc((number_of_blocks+1)*sizeof(double));
		
		int number_of_candidate_blocks = 0;

		for(i = 0 ; i < number_of_blocks; i++)
		{
			if(block_coordinates[0][i] == -1 || block_coordinates[1][i] == -1)
			{
				continue;
			}
			int current_start = block_coordinates[0][i];
			int current_end = block_coordinates[1][i];
		  int block_snp_count = find_number_of_snps_in_block(current_start, current_end, snp_site_coords, branch_snp_sequence, number_of_branch_snps);
		  int block_genome_size_without_gaps = block_coordinates[3][i];

			if(p_value_test(branch_genome_size, block_genome_size_without_gaps, number_of_branch_snps, block_snp_count, min_snps) == 1)
			{
				
				candidate_blocks[0][number_of_candidate_blocks] = block_coordinates[0][i];
				candidate_blocks[1][number_of_candidate_blocks] = block_coordinates[1][i];
				// TODO use a float in a struct here, should be okay for the moment but assumes that there will be a clear integer difference between best and second best
				candidate_blocks[2][number_of_candidate_blocks] = block_coordinates[2][i];
				candidate_blocks[3][number_of_candidate_blocks] = block_genome_size_without_gaps;
				
				candidate_block_likelihoods[number_of_candidate_blocks] = block_likelihoods[i];
				number_of_candidate_blocks++;
			}
		}
		if(number_of_candidate_blocks == 0 )
		{
			free(block_coordinates[0]) ;
			free(block_coordinates[1]) ;
			free(block_coordinates[2]) ;
			free(block_coordinates[3]) ;
			free(block_likelihoods);
			free(candidate_blocks[0]);
		  free(candidate_blocks[1]);
		  free(candidate_blocks[2]);
		  free(candidate_blocks[3]);
			free(candidate_block_likelihoods);
			
			int new_recombination_size = (current_node->num_recombinations+1)*sizeof(int);
			if(new_recombination_size > 1024)
			{
			  current_node->recombinations = (int *) realloc(current_node->recombinations, new_recombination_size);
		  }
			return;	
		}
		number_of_branch_snps = flag_smallest_log_likelihood_recombinations(candidate_blocks, number_of_candidate_blocks, number_of_branch_snps, snp_site_coords, current_node->recombinations, current_node->num_recombinations,current_node, block_file_pointer, root, snp_locations, length_of_sequence,gff_file_pointer,candidate_block_likelihoods );
		branch_genome_size = original_branch_genome_size  - current_node->total_bases_removed_excluding_gaps;
	  free(candidate_blocks[0]);
	  free(candidate_blocks[1]);
	  free(candidate_blocks[2]);
	  free(candidate_blocks[3]);
		free(candidate_block_likelihoods);
	
	}
	free(block_coordinates[0]) ;
	free(block_coordinates[1]) ;
	free(block_coordinates[2]) ;
	free(block_coordinates[3]) ;
	free(block_likelihoods);
	int new_recombination_size = (current_node->num_recombinations+1)*sizeof(int);
	if(new_recombination_size > 1024)
	{
	  current_node->recombinations = (int *) realloc(current_node->recombinations, new_recombination_size);
  }
}


int get_blocks(int ** block_coordinates, int genome_size,int * snp_site_coords,int number_of_branch_snps, int window_size, int cutoff)
{
	// Set up the window counter with 1 value per base in the branch
 	int * window_count;
	window_count = (int *) malloc((genome_size+1)*sizeof(int));
	int i;
	for(i =0; i< genome_size; i++)
	{
		window_count[i] = 0;
	}
	
	
	// create the pileup of the snps and their sphere of influence
	int snp_counter = 0;
	for(snp_counter = 0; snp_counter < number_of_branch_snps; snp_counter++)
	{
		// Lower bound of the window around a snp
		int snp_sliding_window_counter = snp_site_coords[snp_counter]-(window_size/2);
		if(snp_sliding_window_counter < 0)
		{
			snp_sliding_window_counter = 0;
		}
		
		// Upper bound of the window around a snp
		int max_snp_sliding_window_counter = snp_site_coords[snp_counter]+(window_size/2);
		if(max_snp_sliding_window_counter>genome_size)
		{
			max_snp_sliding_window_counter = genome_size;
		}
		
		int j = 0;
		for(j = snp_sliding_window_counter; j < max_snp_sliding_window_counter; j++)
		{
			window_count[j] += 1;	
		}
	}
	
	int number_of_blocks = 0;
	int in_block = 0;
	int block_lower_bound = 0;
	// Scan across the pileup and record where blocks are above the cutoff
	for(i = 0; i < genome_size; i++)
	{
		// Just entered the start of a block
		if(window_count[i] > cutoff && in_block == 0)
		{
			block_lower_bound = i;
			in_block = 1;
    }

		// Just left a block
		if(window_count[i] <= cutoff && in_block == 1)
		{
			block_coordinates[0][number_of_blocks] = block_lower_bound;
			block_coordinates[1][number_of_blocks] = i-1;
			number_of_blocks++;
			in_block = 0;
		}
		
	}
	free(window_count);
	return number_of_blocks;

}



void move_blocks_inwards_while_likelihood_improves(int number_of_blocks,int ** block_coordinates, int min_snps, int * snp_site_coords,  int number_of_branch_snps,char * branch_snp_sequence, int * snp_locations, int branch_genome_size,char * child_sequence, int length_of_sequence, double * block_likelihoods, int cutoff_value)
{
	int i;
	
	int previous_start;
	int previous_end;
	
	for(i = 0 ; i < number_of_blocks; i++)
	{
		int current_start = block_coordinates[0][i];
		int current_end = block_coordinates[1][i];
		int start_index = find_starting_index( current_start, snp_site_coords,0, number_of_branch_snps);
    int end_index   = find_starting_index( current_end, snp_locations, start_index, number_of_branch_snps);
		block_coordinates[0][i] = advance_window_start_to_next_snp_with_start_end_index(current_start, snp_site_coords, branch_snp_sequence, number_of_branch_snps,start_index,end_index);
		block_coordinates[1][i] = rewind_window_end_to_last_snp_with_start_end_index(current_end, snp_site_coords, branch_snp_sequence, number_of_branch_snps,start_index,end_index);
		
		if( i == 0)
		{
			previous_start = block_coordinates[0][i];
			previous_end = block_coordinates[1][i];
		}
		else if(previous_start == block_coordinates[0][i] && previous_end == block_coordinates[1][i] && i > 0)
		{
			block_coordinates[0][i] = -1;
			block_coordinates[1][i] = -1;
		}
		else
		{
			previous_start = block_coordinates[0][i];
			previous_end = block_coordinates[1][i];
		}
	}
	
	
	for(i = 0 ; i < number_of_blocks; i++)
	{
		if(block_coordinates[0][i] == -1 || block_coordinates[1][i] == -1)
		{
			continue;
		}
		
		int current_start = block_coordinates[0][i];
		int current_end = block_coordinates[1][i];
		double current_block_likelihood = block_coordinates[2][i]*1.0;
		int block_genome_size_without_gaps = block_coordinates[3][i];
		int block_snp_count;
		

	  int next_start_position = current_start;
		int start_index = find_starting_index( current_start, snp_site_coords,0, number_of_branch_snps);
    int end_index =  find_starting_index( current_end, snp_locations, start_index, number_of_branch_snps);

		block_snp_count = find_number_of_snps_in_block_with_start_end_index(current_start, current_end, snp_site_coords, branch_snp_sequence, number_of_branch_snps,start_index,end_index);
    
		if(current_block_likelihood < 0 || block_genome_size_without_gaps == -1)
		{
			
			block_genome_size_without_gaps = calculate_block_size_without_gaps(child_sequence, snp_locations, current_start, current_end, length_of_sequence);
			current_block_likelihood = get_block_likelihood(branch_genome_size, number_of_branch_snps, block_genome_size_without_gaps, block_snp_count);
			block_coordinates[2][i] = current_block_likelihood;
			block_coordinates[3][i] = block_genome_size_without_gaps;
		}


		// Move left inwards while the likelihood gets better
		while(current_start < current_end && block_snp_count >= min_snps && block_snp_count >= cutoff_value && block_genome_size_without_gaps > MIN_WINDOW_SIZE)
		{
			  next_start_position++;
			  next_start_position = advance_window_start_to_next_snp_with_start_end_index(next_start_position, snp_site_coords, branch_snp_sequence, number_of_branch_snps,start_index,end_index);
			
				if(next_start_position == current_start)
				{
					break;
				}
			
				int previous_block_snp_count = block_snp_count;
				int previous_block_genome_size_without_gaps = block_genome_size_without_gaps;
			  block_snp_count = find_number_of_snps_in_block_with_start_end_index(next_start_position, current_end, snp_site_coords, branch_snp_sequence, number_of_branch_snps,start_index,end_index);
			  block_genome_size_without_gaps = calculate_block_size_without_gaps_with_start_end_index(child_sequence, snp_locations, next_start_position, current_end, length_of_sequence,start_index,end_index);
			  
			  double next_block_likelihood = get_block_likelihood(branch_genome_size, number_of_branch_snps, block_genome_size_without_gaps, block_snp_count);
				
			if(next_block_likelihood <= current_block_likelihood)
			  {
			  	current_block_likelihood = next_block_likelihood;
					current_start = next_start_position;
					start_index++;
			  }
			  else
			  {
					block_snp_count = previous_block_snp_count;
					block_genome_size_without_gaps = previous_block_genome_size_without_gaps;
					break;
			  }
		}
		
		int next_end_position = current_end;
		
		// Move Right inwards while the likelihood gets better
		while(current_start < current_end && block_snp_count >= min_snps && block_snp_count >= cutoff_value && block_genome_size_without_gaps > MIN_WINDOW_SIZE)
		{
			  next_end_position--;
				int next_end_position_proposed = rewind_window_end_to_last_snp_with_start_end_index(next_end_position, snp_site_coords, branch_snp_sequence, number_of_branch_snps,start_index,end_index);
				
				if(next_end_position_proposed > next_end_position)
				{
					next_end_position = next_end_position_proposed;
				}
				if(next_end_position == current_end)
				{
					break;
				}
				
			  block_snp_count = find_number_of_snps_in_block_with_start_end_index(current_start, current_end, snp_site_coords, branch_snp_sequence, number_of_branch_snps,start_index,end_index);
			  block_genome_size_without_gaps = calculate_block_size_without_gaps_with_start_end_index(child_sequence, snp_locations, current_start, current_end, length_of_sequence,start_index,end_index);
			  
			  double next_block_likelihood = get_block_likelihood(branch_genome_size, number_of_branch_snps, block_genome_size_without_gaps, block_snp_count);
			  if(next_block_likelihood <= current_block_likelihood)
			  {
			  	current_block_likelihood = next_block_likelihood;
					current_end = next_end_position;
					end_index--;
			  }
			  else
			  {
					break;
			  }
			
		}
		
		  block_coordinates[0][i] = current_start;
		  block_coordinates[1][i] = current_end;
	  	block_coordinates[2][i] = (int) current_block_likelihood;
		  block_coordinates[3][i] = block_genome_size_without_gaps;
		
		  block_likelihoods[i] = current_block_likelihood;
	}	
}



int exclude_snp_sites_in_block(int window_start_coordinate, int window_end_coordinate, int * snp_site_coords, int number_of_branch_snps)
{
	int i;
	// inefficient
	int updated_snp_site_coords[number_of_branch_snps];
	int number_of_branch_snps_excluding_block = 0 ;
	
	for(i = 0 ; i< number_of_branch_snps; i++)
	{
		if(snp_site_coords[i]>= window_start_coordinate && snp_site_coords[i] <= window_end_coordinate)
		{
		}
		else
		{
			updated_snp_site_coords[number_of_branch_snps_excluding_block] = snp_site_coords[i];
			number_of_branch_snps_excluding_block++;
		}
	}
	
	for(i = 0; i < number_of_branch_snps_excluding_block; i++)
	{
	  snp_site_coords[i] = updated_snp_site_coords[i];
	}
	
	snp_site_coords = realloc(snp_site_coords, (number_of_branch_snps_excluding_block+1)*sizeof(int));
	return number_of_branch_snps_excluding_block;
}

int flag_smallest_log_likelihood_recombinations(int ** candidate_blocks, int number_of_candidate_blocks, int number_of_branch_snps, int * snp_site_coords, int * recombinations, int number_of_recombinations,newick_node * current_node, FILE * block_file_pointer, newick_node *root,int * snp_locations, int total_num_snps, FILE * gff_file_pointer, double * block_likelihooods)
{
	int number_of_branch_snps_excluding_block = number_of_branch_snps;
	if(number_of_candidate_blocks > 0)
	{
		int smallest_index = 0;
    int number_of_recombinations_in_window = 0;
		smallest_index = get_smallest_log_likelihood(block_likelihooods, number_of_candidate_blocks);
		number_of_recombinations_in_window = flag_recombinations_in_window(candidate_blocks[0][smallest_index], candidate_blocks[1][smallest_index],number_of_branch_snps, snp_site_coords, recombinations, number_of_recombinations,snp_locations,total_num_snps);	
    number_of_recombinations += number_of_recombinations_in_window;
		number_of_branch_snps_excluding_block = exclude_snp_sites_in_block(candidate_blocks[0][smallest_index],candidate_blocks[1][smallest_index], snp_site_coords,number_of_branch_snps);
		
		current_node->num_recombinations = number_of_recombinations;

		print_block_details(block_file_pointer, candidate_blocks[0][smallest_index], candidate_blocks[1][smallest_index],  number_of_recombinations_in_window, current_node->taxon,  root->taxon, current_node->taxon_names, current_node->childNum);
		print_gff_line(gff_file_pointer, candidate_blocks[0][smallest_index], candidate_blocks[1][smallest_index],  number_of_recombinations_in_window, current_node->taxon,  root->taxon, current_node->taxon_names);
		current_node->number_of_blocks = current_node->number_of_blocks + 1;
		
		current_node->total_bases_removed_excluding_gaps = current_node->total_bases_removed_excluding_gaps  + candidate_blocks[3][smallest_index];

		current_node->block_coordinates[0] = realloc((int *)current_node->block_coordinates[0], ((int)current_node->number_of_blocks +1)*sizeof(int));
		current_node->block_coordinates[1] = realloc((int *)current_node->block_coordinates[1], ((int)current_node->number_of_blocks +1)*sizeof(int));
		
		current_node->block_coordinates[0][current_node->number_of_blocks -1] = candidate_blocks[0][smallest_index];
		current_node->block_coordinates[1][current_node->number_of_blocks -1] = candidate_blocks[1][smallest_index];
	}
	current_node->number_of_snps = number_of_branch_snps_excluding_block;
	
	return number_of_branch_snps_excluding_block;
}

// candidate blocks contains, start coordinate, end_coordinate and log likelihood
int get_smallest_log_likelihood(double * candidate_blocks, int number_of_candidate_blocks)
{
	int i;
	int smallest_index = 0 ; 
	
	for(i=0; i< number_of_candidate_blocks; i++)
	{
		if(candidate_blocks[i] < candidate_blocks[smallest_index] && candidate_blocks[i] > 0)
		{
		   smallest_index = i;
		}
	}
	return smallest_index;
}


double snp_density(int length_of_sequence, int number_of_snps)
{
	return number_of_snps*1.0/length_of_sequence;
}




double calculate_threshold(int branch_genome_size, int window_size)
{
	return 1-(RANDOMNESS_DAMPNER/((branch_genome_size*1.0)/((window_size*1.0)/WINDOW_SNP_MODE_TARGET)));
}

int calculate_cutoff(int branch_genome_size, int window_size, int num_branch_snps)
{
	double threshold = 0.0;
	int cutoff = 0;
	double pvalue = 0.0;
	double part1, part2, part3 = 0.0;
	
	threshold = calculate_threshold(branch_genome_size, window_size);
	
	while( pvalue <= threshold)
	{
		part1 = reduce_factorial(window_size,cutoff)-reduce_factorial(cutoff,cutoff);
		part2 = log10((num_branch_snps*1.0)/branch_genome_size)*cutoff;
		part3 = log10(1.0-((num_branch_snps*1.0)/branch_genome_size))*(window_size-cutoff);
		pvalue = pvalue + pow(10,(part1 + part2 + part3));
		cutoff++;
	}
	cutoff--;

	
	return cutoff;
}

int p_value_test(int branch_genome_size, int window_size, int num_branch_snps, int block_snp_count, int min_snps)
{
	double threshold = 0.0;
	int cutoff = 0;
	double pvalue = 0.0;
	double part1, part2, part3 = 0.0;
	
	if( block_snp_count < min_snps)
	{
		return 0;	
	}
	
	threshold = 0.05/branch_genome_size;
	
	while( cutoff < block_snp_count)
	{
		part1 = reduce_factorial(window_size,cutoff)-reduce_factorial(cutoff,cutoff);
		part2 = log10((num_branch_snps*1.0)/branch_genome_size)*cutoff;
		part3 = log10(1.0-((num_branch_snps*1.0)/branch_genome_size))*(window_size-cutoff);
		pvalue = pvalue + pow(10,(part1 + part2 + part3));
		cutoff++;
	}
	
	// There should be rounding here to 10 decimal places
	pvalue = 1.0-pvalue;
	
	if(pvalue < threshold)
	{
		// the block could contain a recombination
		return 1;
	}
	
	return 0;
}


double reduce_factorial(int l, int i)
{
	double factorial;
	int x;
	
	factorial = log10(1.0);
	
	for(x = l-(i-1); x < l+1 ; x++)
	{
		factorial = factorial + log10(x);
	}
	return factorial;
}


// N = branch_genome_size
// C = number_of_branch_snps
// n = block_genome_size_without_gaps
// c = number_of_block_snps
double get_block_likelihood(int branch_genome_size, int number_of_branch_snps, int block_genome_size_without_gaps, int number_of_block_snps)
{
	double part1, part2, part3, part4;
	
	if(block_genome_size_without_gaps == 0)
	{
		return 0.0;
	}
	if(number_of_block_snps == 0)
	{
		return 0.0;
	}
	
	part1 = log10(number_of_block_snps*1.0/block_genome_size_without_gaps*1.0)*number_of_block_snps;
	
	if((block_genome_size_without_gaps-number_of_block_snps) == 0)
	{
		part2 = 0.0;	
	}
	else
	{
		part2 = log10( ((block_genome_size_without_gaps-number_of_block_snps)*1.0)/block_genome_size_without_gaps*1.0 )*(block_genome_size_without_gaps-number_of_block_snps);
	}
	
	if((number_of_branch_snps-number_of_block_snps) == 0)
	{
		part3 = 0.0;
	}
	else
	{
		part3 = log10(((number_of_branch_snps-number_of_block_snps)*1.0)/((branch_genome_size-block_genome_size_without_gaps)*1.0))*(number_of_branch_snps-number_of_block_snps);
	}
			
	if(((branch_genome_size-block_genome_size_without_gaps)-(number_of_branch_snps-number_of_block_snps))==0)
	{
	  part4 = 0.0;
	}
	else
	{
		part4=log10(((branch_genome_size-block_genome_size_without_gaps)-(number_of_branch_snps-number_of_block_snps)*1.0 )/((branch_genome_size-block_genome_size_without_gaps)*1.0)) * ((branch_genome_size-block_genome_size_without_gaps)-(number_of_branch_snps-number_of_block_snps));
	}
	
	return (part1+part2+part3+part4)*-1;
}


