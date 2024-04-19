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
#include <stdint.h>
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
#include "string_cat.h"

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


// Loop over all recombination blocks and return a list of all snp indices that fall within those blocks.
int get_list_of_snp_indices_which_fall_in_downstream_recombinations(int ** current_block_coordinates,int num_blocks, int * snp_locations,int number_of_snps, int * snps_in_recombinations)
{
	int num_snps_in_recombinations =0;
	int i = 0;
  
  // loop over each block
	for(i = 0; i<num_blocks; i++ )
	{
		int current_index = 0;
    // convert the starting coordinates of block to the nearest SNP index
		current_index = find_starting_index(current_block_coordinates[0][i],snp_locations,0, number_of_snps);
		
    //make sure that the index begins at start of block
		int beginning_j = current_index;
    for(beginning_j = current_index; snp_locations[beginning_j] < current_block_coordinates[0][i];beginning_j++)
    {
    }
    
    int j;
    // starting at the begining index of block, count all the snps until the end of the bock.
		for(j = beginning_j; (j < number_of_snps && snp_locations[j] <= current_block_coordinates[1][i]); j++)
		{
				snps_in_recombinations[num_snps_in_recombinations] = j;
				num_snps_in_recombinations++;
		}
	}

  // may contain duplications
	return num_snps_in_recombinations;
}


// Go through the tree and build up the recombinations list from root to branch. Print out each sample name and a list of recombinations
void fill_in_recombinations_with_gaps(newick_node *root, int * parent_recombinations, int parent_num_recombinations, int current_total_snps,int num_blocks, int ** current_block_coordinates,int length_of_original_genome,int * snp_locations, int number_of_snps )
{
	newick_child *child;
	int * current_recombinations;
	int num_current_recombinations = 0 ;
	char * child_sequence = (char *) calloc((length_of_original_genome +1),sizeof(char));
	
	current_recombinations = (int *) calloc((root->num_recombinations+1+parent_num_recombinations),sizeof(int));
	num_current_recombinations = copy_and_concat_integer_arrays(root->recombinations,
                                                                root->num_recombinations,
                                                                parent_recombinations,
                                                                parent_num_recombinations,
                                                                current_recombinations);
	
 	// overwrite the bases of snps with N's
 	int i;
 	int sequence_index;
 	sequence_index = find_sequence_index_from_sample_name(root->taxon);
 	
 	set_number_of_recombinations_for_sample(root->taxon,root->num_recombinations);
 	set_number_of_snps_for_sample(root->taxon,root->number_of_snps);
	
	get_sequence_for_sample_name(child_sequence, root->taxon);
	int genome_length_excluding_blocks_and_gaps = calculate_genome_length_excluding_blocks_and_gaps(child_sequence,
                                                                                                    length_of_original_genome,
                                                                                                    current_block_coordinates,
                                                                                                    num_blocks);
	
	set_genome_length_excluding_blocks_and_gaps_for_sample(root->taxon,
                                                           genome_length_excluding_blocks_and_gaps);
	
	int ** merged_block_coordinates;
	merged_block_coordinates = (int **) calloc(3,sizeof(int *));
	merged_block_coordinates[0] = (int*) calloc((num_blocks + root->number_of_blocks+1),sizeof(int ));
	merged_block_coordinates[1] = (int*) calloc((num_blocks + root->number_of_blocks+1),sizeof(int ));
	copy_and_concat_2d_integer_arrays(current_block_coordinates,
                                      num_blocks,
                                      root->block_coordinates,
                                      root->number_of_blocks,
                                      merged_block_coordinates
                                      );
	
	set_number_of_blocks_for_sample(root->taxon, root->number_of_blocks);
    set_number_of_branch_bases_in_recombinations(root->taxon,
                                                 calculate_number_of_bases_in_recombations_excluding_gaps(merged_block_coordinates,
                                                                                                          root->number_of_blocks,
                                                                                                          child_sequence,
                                                                                                          snp_locations,
                                                                                                          current_total_snps)
                                                 );
 	set_number_of_bases_in_recombinations(root->taxon,
                                          calculate_number_of_bases_in_recombations_excluding_gaps(merged_block_coordinates,
                                                                                                   (num_blocks + root->number_of_blocks),
                                                                                                   child_sequence,
                                                                                                   snp_locations,
                                                                                                   current_total_snps)
                                          );
	free(child_sequence); 	

 	for(i = 0; i < num_current_recombinations; i++)
 	{
 		update_sequence_base('N', sequence_index, current_recombinations[i]);
 	}

    // TODO: The stats for the number of snps in recombinations will need to be updated.
	int * snps_in_recombinations = (int *) calloc((number_of_snps +1),sizeof(int));
	int num_snps_in_recombinations = get_list_of_snp_indices_which_fall_in_downstream_recombinations(merged_block_coordinates,
                                                                                                     (num_blocks + root->number_of_blocks),
                                                                                                     snp_locations,
                                                                                                     number_of_snps,
                                                                                                     snps_in_recombinations);
 	for(i = 0; i < num_snps_in_recombinations; i++)
 	{
 		update_sequence_base('N', sequence_index, snps_in_recombinations[i]);
 	}
	free(snps_in_recombinations); 	

	if (root->childNum > 0)
	{
		child = root->child;
		set_internal_node(1,sequence_index);

		while (child != NULL)
		{
			fill_in_recombinations_with_gaps(child->node,
                                             current_recombinations,
                                             num_current_recombinations,
                                             (current_total_snps + root->number_of_snps),
                                             (num_blocks + root->number_of_blocks),
                                             merged_block_coordinates,
                                             length_of_original_genome,
                                             snp_locations,
                                             number_of_snps
                                             );
			child = child->next;

		}
	}
	else
	{
        set_internal_node(0,sequence_index);
	}
	free(current_recombinations);
	free(merged_block_coordinates[0]);
	free(merged_block_coordinates[1]);
	free(merged_block_coordinates);
}

int calculate_number_of_bases_in_recombations_excluding_gaps(int ** block_coordinates, int num_blocks,char * child_sequence, int * snp_locations,int current_total_snps)
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
		
		
	  total_bases += calculate_block_size_without_gaps(child_sequence,
                                                       snp_locations,
                                                       block_coordinates[0][start_block],
                                                       block_coordinates[1][start_block],
                                                       current_total_snps);
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
		//fill_in_unambiguous_bases_in_parent_from_children_where_parent_has_a_gap(parent_sequence_index, child_sequence_indices,child_counter);
	}
}

void generate_branch_sequences(newick_node *node, char ** node_sequences, char ** node_names, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns, int length_of_original_genome, int num_stored_nodes, FILE * block_file_pointer, FILE * gff_file_pointer,int min_snps,FILE * branch_snps_file_pointer, int window_min, int window_max, float uncorrected_p_value, float trimming_ratio, int extensive_search_flag)
{
	newick_child *child;
	int child_counter = 0;
	int current_branch = 0;
	int branch_genome_size = 0;
	int number_of_branch_snps = 0;
	
  char * node_sequence = (char *) calloc((number_of_snps +1),sizeof(char));
  
	if (node->childNum == 0)
	{
    
		get_sequence_for_sample_name(node_sequence, node->taxon);
		
    node->taxon_names = (char *) calloc(MAX_SAMPLE_NAME_SIZE,sizeof(char));
		memcpy(node->taxon_names, node->taxon, size_of_string(node->taxon)+1);

    // Save some statistics about the sequence
		branch_genome_size = calculate_size_of_genome_without_gaps(node_sequence, 0,number_of_snps, length_of_original_genome);
		set_genome_length_without_gaps_for_sample(node->taxon,branch_genome_size);
		
	}
	else
	{
		child = node->child;
		char * child_sequences[node->childNum];
		newick_node * child_nodes[node->childNum];
    node->taxon_names = (char *) calloc(MAX_SAMPLE_NAME_SIZE*number_of_columns,sizeof(char));

		// generate pointers for each child sequence

		while (child != NULL)
		{
			// Retrieve child sequences from store
      for (int seq_store_index = 0; seq_store_index  < num_stored_nodes; ++seq_store_index)
      {
        if (node_names[seq_store_index] == *child->node->taxon) {
          child_sequences[child_counter] = &node_sequences[seq_store_index];
          break;
        }
      }
			child_nodes[child_counter] = child->node;
			
      // Remove from store as cannot be children of any other nodes
      // TO DO
      
			char delimiter_string[3] = {" "};
			concat_strings_created_with_malloc(node->taxon_names, delimiter_string);
			concat_strings_created_with_malloc(node->taxon_names, child_nodes[child_counter]->taxon_names);
			
			child_counter++;
			child = child->next;
		}
		
    // Get sequence reconstructed at internal node
		node_sequence = (char *) calloc((number_of_snps +1),sizeof(char));
		get_sequence_for_sample_name(node_sequence, node->taxon);
		branch_genome_size = calculate_size_of_genome_without_gaps(node_sequence, 0,number_of_snps, length_of_original_genome);
		set_genome_length_without_gaps_for_sample(node->taxon,branch_genome_size);
		
    // Identify recombinations on descendant branches
		for(current_branch = 0 ; current_branch< (node->childNum); current_branch++)
		{
			int * branches_snp_sites;
			branches_snp_sites = (int *) calloc((number_of_snps +1),sizeof(int));
			char * branch_snp_sequence;
			char * branch_snp_ancestor_sequence;
			branch_snp_sequence = (char *) calloc((number_of_snps +1),sizeof(char));
			branch_snp_ancestor_sequence = (char *) calloc((number_of_snps +1),sizeof(char));
			
			branch_genome_size = calculate_size_of_genome_without_gaps(child_sequences[current_branch], 0,number_of_snps, length_of_original_genome);
			number_of_branch_snps = calculate_number_of_snps_excluding_gaps(node_sequence, child_sequences[current_branch], number_of_snps, branches_snp_sites, snp_locations,branch_snp_sequence,branch_snp_ancestor_sequence);
			
			
			child_nodes[current_branch]->number_of_snps = number_of_branch_snps;
			print_branch_snp_details(branch_snps_file_pointer, child_nodes[current_branch]->taxon,node->taxon, branches_snp_sites, number_of_branch_snps, branch_snp_sequence, branch_snp_ancestor_sequence,child_nodes[current_branch]->taxon_names);
			
			get_likelihood_for_windows(child_sequences[current_branch],
                                       number_of_snps,
                                       branches_snp_sites,
                                       branch_genome_size,
                                       number_of_branch_snps,
                                       snp_locations,
                                       child_nodes[current_branch],
                                       block_file_pointer,
                                       node,
                                       branch_snp_sequence,
                                       gff_file_pointer,
                                       min_snps,
                                       length_of_original_genome,
                                       node_sequence,
                                       window_min,
                                       window_max,
                                       uncorrected_p_value,
                                       trimming_ratio,
                                       extensive_search_flag);
			free(branch_snp_sequence);
			free(branch_snp_ancestor_sequence);
			free(branches_snp_sites);
			
		}

	}
  
  // Store node sequence
  for (int seq_store_index = 0; seq_store_index  < num_stored_nodes; ++seq_store_index)
  {
    if (node_names[seq_store_index] == ' ') {
      node_names[seq_store_index]  = node->taxon;
      node_sequences[seq_store_index] = node_sequence;
      break;
    }
  }
  
}


// Windows need to be of a fixed size
// calculate window size
// starting at coord of first snp, count number of snps which fall into window
// if region is blank, move on
int calculate_window_size(int branch_genome_size, int number_of_branch_snps,int window_min, int window_max, int min_snps, int window_factor)
{
	int window_size = 0;
	if(number_of_branch_snps == 0)
	{
		return window_min;
	}
	
	window_size = (int) ((branch_genome_size*1.0)/(window_factor*number_of_branch_snps*1.0/(min_snps - 1)));
	
	if(window_size < window_min)
	{
		return window_min;
	}
	else if(window_size > window_max)
	{
		return 	window_max;
	}

	return window_size;
}


void get_likelihood_for_windows(char * child_sequence, int current_num_snps, int * snp_site_coords, int branch_genome_size, int number_of_branch_snps, int * snp_locations, newick_node * current_node, FILE * block_file_pointer, newick_node *root, char * branch_snp_sequence, FILE * gff_file_pointer, int min_snps, int length_of_original_genome, char * original_sequence,int window_min, int window_max, float uncorrected_p_value, float trimming_ratio, int extensive_search_flag)
{
    
    // define variables
    int i = 0;
    int window_size = window_max;
    int number_of_snps_in_block = 0;
    int block_genome_size_without_gaps = 0;
    double branch_snp_density = 0.0;
    double block_snp_density = 0.0;
    int number_of_blocks = 0 ;
    int original_branch_genome_size = branch_genome_size;

    // place to store coordinates of recombinations snps
    current_node->recombinations = (int *) calloc((number_of_branch_snps+1),sizeof(int));

    // place to store candidate recombination information
    int number_of_windows = (branch_genome_size/window_min) + 1;
    int * block_coordinates[4];
    block_coordinates[0] = (int *) calloc((number_of_windows+1),sizeof(int));
    block_coordinates[1] = (int *) calloc((number_of_windows+1),sizeof(int));
    block_coordinates[2] = (int *) calloc((number_of_windows+1),sizeof(int));
    block_coordinates[3] = (int *) calloc((number_of_windows+1),sizeof(int));

    // place to store candidate block likelihoods
    double * block_likelihoods;
    block_likelihoods = (double *) calloc((number_of_windows+1),sizeof(double));

    // iterate over SNPs
    // Keep searching while there is the possibility of detecting a small recombination containing
    // the minimum number of SNPs
    int window_factor = 1;
    int cutoff = min_snps - 1;
    int previous_cutoff = number_of_branch_snps + 1;
    while(number_of_branch_snps >= min_snps && window_size > window_min)
    {
        
        // return SNP density as double
        branch_snp_density = snp_density(branch_genome_size, number_of_branch_snps);

        // return sensible window size
        window_size = calculate_window_size(branch_genome_size,
                                            number_of_branch_snps,
                                            window_min,
                                            window_max,
                                            min_snps,
                                            window_factor);

        // return a cutoff number of SNPs in a window
        // for extensive search, this is every window containing min SNPs count
        // otherwise focus only on windows likely to exceed the statistical threshold
        // for detecting recombination (faster)
        if (extensive_search_flag == 0)
        {
            cutoff = calculate_cutoff(branch_genome_size,
                                          window_size,
                                          number_of_branch_snps,
                                          min_snps,
                                          uncorrected_p_value);
        }
        
        // If returned cutoff == 0, then exit the program and error
        if (cutoff == 0)
        {
            fprintf(stderr,
                    "Cannot identify recombinations on at least one branch with window size %d; try increasing this value\n",
                    window_size);
            exit(EXIT_FAILURE);
        }

        // Test if reduction in window size has reduced the cutoff
        int number_of_candidate_blocks = 0;
        if (cutoff < previous_cutoff)
        {
            // populate block coordinate data structure by identifying windows containing
            // a greater number of SNPs than the threshold and trimming them based on SNP
            // positions
            number_of_blocks = get_blocks(block_coordinates,
                                          length_of_original_genome,
                                          snp_site_coords,
                                          number_of_branch_snps,
                                          window_size,
                                          cutoff,
                                          child_sequence,
                                          snp_locations,
                                          current_num_snps);
            
            // iterate over blocks
            for(i = 0; i < number_of_blocks; i++)
            {
                // get number of SNPs in block
                number_of_snps_in_block = find_number_of_snps_in_block(block_coordinates[0][i],
                                                                       block_coordinates[1][i],
                                                                       snp_site_coords,
                                                                       branch_snp_sequence,
                                                                       number_of_branch_snps);
                
                // get number of bases in block
                block_genome_size_without_gaps = calculate_block_size_without_gaps(child_sequence,
                                                                                   snp_locations,
                                                                                   block_coordinates[0][i],
                                                                                   block_coordinates[1][i],
                                                                                   current_num_snps);

                // minimum number of snps to be statistically significant in block
                if(number_of_snps_in_block < min_snps)
                {
                    block_coordinates[0][i] = -1;
                    block_coordinates[1][i] = -1;
                    continue;
                }

                // calculate SNP density of block
                block_snp_density = snp_density(block_genome_size_without_gaps, number_of_snps_in_block);
                
                // region with low number of snps so skip over
                if(block_snp_density <= branch_snp_density)
                {
                    block_coordinates[0][i] = -1;
                    block_coordinates[1][i] = -1;
                    continue;
                }

                // calculate block log likelihood under null model
                block_likelihoods[i] = get_block_likelihood(branch_genome_size,
                                                            number_of_branch_snps,
                                                            block_genome_size_without_gaps,
                                                            number_of_snps_in_block);
                block_coordinates[2][i] = (int) block_likelihoods[i]; // casts double log likelihood to int
                block_coordinates[3][i] = block_genome_size_without_gaps;
                
            }

            // trim the edges of candidate recombinations
            move_blocks_inwards_while_likelihood_improves(number_of_blocks,
                                                          block_coordinates,
                                                          min_snps,
                                                          snp_site_coords,
                                                          number_of_branch_snps,
                                                          branch_snp_sequence,
                                                          snp_locations,
                                                          branch_genome_size,
                                                          child_sequence,
                                                          current_num_snps,
                                                          block_likelihoods,
                                                          cutoff,
                                                          trimming_ratio);

            int * candidate_blocks[4];
            double * candidate_block_likelihoods;
            
            candidate_blocks[0] = (int *) calloc((number_of_blocks+1),sizeof(int));
            candidate_blocks[1] = (int *) calloc((number_of_blocks+1),sizeof(int));
            candidate_blocks[2] = (int *) calloc((number_of_blocks+1),sizeof(int));
            candidate_blocks[3] = (int *) calloc((number_of_blocks+1),sizeof(int));

            candidate_block_likelihoods = (double *) calloc((number_of_blocks+1),sizeof(double));

            for(i = 0 ; i < number_of_blocks; i++)
            {
                if(block_coordinates[0][i] == -1 || block_coordinates[1][i] == -1)
                {
                    continue;
                }
                int current_start = block_coordinates[0][i];
                int current_end = block_coordinates[1][i];
                int block_snp_count = find_number_of_snps_in_block(current_start,
                                                                   current_end,
                                                                   snp_site_coords,
                                                                   branch_snp_sequence,
                                                                   number_of_branch_snps);
                int block_genome_size_without_gaps = block_coordinates[3][i];

                if(p_value_test(branch_genome_size, block_genome_size_without_gaps, number_of_branch_snps, block_snp_count, min_snps, uncorrected_p_value) == 1)
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
            
            if(number_of_candidate_blocks > 0)
            {
                // remove recombination with smallest log likelihood and
                // correspondingly reduce the number of branch SNPs
                number_of_branch_snps = flag_smallest_log_likelihood_recombinations(candidate_blocks,
                                                                                    number_of_candidate_blocks,
                                                                                    number_of_branch_snps,
                                                                                    snp_site_coords,
                                                                                    current_node->recombinations,
                                                                                    current_node->num_recombinations,
                                                                                    current_node,
                                                                                    block_file_pointer,
                                                                                    root,
                                                                                    snp_locations,
                                                                                    current_num_snps,
                                                                                    gff_file_pointer,
                                                                                    candidate_block_likelihoods);
                
                branch_genome_size = original_branch_genome_size  - current_node->total_bases_removed_excluding_gaps;

            }
            
            free(candidate_blocks[0]);
            free(candidate_blocks[1]);
            free(candidate_blocks[2]);
            free(candidate_blocks[3]);
            free(candidate_block_likelihoods);

        }
        
        if(number_of_candidate_blocks == 0)
        {
            window_factor = window_factor * 2;
            previous_cutoff = cutoff;
        }
        
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

int extend_upper_part_of_window(int starting_coord, int initial_max_coord, int genome_size, int8_t * gaps_in_original_genome_space)
{
		int max_snp_sliding_window_counter = initial_max_coord;
		int upper_offset = 0;
		int j = 0; 
		for(j = starting_coord; (j < initial_max_coord + upper_offset) && (j < genome_size); j++)
		{
		  if(gaps_in_original_genome_space[j]  == 1)
			{
				upper_offset++;
			}
		}

		max_snp_sliding_window_counter = initial_max_coord + upper_offset;
		return max_snp_sliding_window_counter;
}

int extend_lower_part_of_window(int starting_coord, int initial_min_coord, int genome_size, int8_t * gaps_in_original_genome_space)
{
		int lower_offset = 0;
		int snp_sliding_window_counter = initial_min_coord;
		int j = 0; 
		for(j = starting_coord; (j >= 0) && j > (initial_min_coord - lower_offset) && (initial_min_coord - lower_offset >= 0); j--)
		{
			if(gaps_in_original_genome_space[j]  == 1)
			{
				lower_offset++;
			}
		}
		snp_sliding_window_counter = snp_sliding_window_counter - lower_offset;
		return snp_sliding_window_counter;
}


int get_blocks(int ** block_coordinates, int genome_size,int * snp_site_coords,int number_of_branch_snps, int window_size, int cutoff, char * original_sequence, int * snp_locations, int number_of_snps)
{
    // Set up the window counter with 1 value per base in the branch
    int * window_count;
    window_count = (int *) calloc((genome_size+1),sizeof(int));

    // Integer array with location of gaps
    int8_t * gaps_in_original_genome_space;
    gaps_in_original_genome_space = (int8_t *) calloc((genome_size+1),sizeof(int8_t));
    int x = 0;
    for(x=0; x< number_of_snps; x++)
    {
        if((original_sequence[x] == 'N' || original_sequence[x] == '-' ) && snp_locations[x] != 0)
        {
            gaps_in_original_genome_space[snp_locations[x]-1] = 1;
        }
    }

    // create the pileup of the snps and their sphere of influence
    int snp_counter = 0;
    int min_postion = genome_size;
    int max_position = 0;
    for(snp_counter = 0; snp_counter < number_of_branch_snps; snp_counter++)
    {
        int j = 0;
        // Lower bound of the window around a snp
        int snp_sliding_window_counter = snp_site_coords[snp_counter]-(window_size/2);

        snp_sliding_window_counter = extend_lower_part_of_window(snp_site_coords[snp_counter] - 1 ,
                                                                 snp_sliding_window_counter,
                                                                 genome_size,
                                                                 gaps_in_original_genome_space);

        if(snp_sliding_window_counter < 0)
        {
            snp_sliding_window_counter = 0;
        }

        if(snp_sliding_window_counter < min_postion)
        {
            min_postion = snp_sliding_window_counter;
        }
      
        // Upper bound of the window around a snp
        int max_snp_sliding_window_counter = snp_site_coords[snp_counter]+(window_size/2);
        max_snp_sliding_window_counter = extend_upper_part_of_window(snp_site_coords[snp_counter] + 1,
                                                                     max_snp_sliding_window_counter,
                                                                     genome_size,
                                                                     gaps_in_original_genome_space);

        if(max_snp_sliding_window_counter>genome_size)
        {
            max_snp_sliding_window_counter = genome_size;
        }
      
        if(max_snp_sliding_window_counter > max_position)
        {
          max_position= max_snp_sliding_window_counter;
        }

        for(j = snp_sliding_window_counter; j < max_snp_sliding_window_counter; j++)
        {
            window_count[j] += 1;
        }
        
    }

    int number_of_blocks = 0;
    int in_block = 0;
    int block_lower_bound = 0;
    // Scan across the pileup and record where blocks are above the cutoff
    int i;
    for(i = min_postion; i <= max_position; i++)
    {
        // Just entered the start of a block
        if(window_count[i] > cutoff && in_block == 0)
        {
            block_lower_bound = i;
            in_block = 1;
        }

        // Reached end of genome
        if (in_block == 1)
        {
          if(i == genome_size)
          {
              block_coordinates[0][number_of_blocks] = block_lower_bound;
              block_coordinates[1][number_of_blocks] = i;
              number_of_blocks++;
              in_block = 0;
          }
          // Just left a block
          else if(window_count[i] <= cutoff)
          {
              block_coordinates[0][number_of_blocks] = block_lower_bound;
              block_coordinates[1][number_of_blocks] = i-1;
              number_of_blocks++;
              in_block = 0;
          }
        }

    }

    // Move blocks inwards to next SNP
    for(i = 0; i < number_of_blocks; i++)
    {
        for(snp_counter = 0; snp_counter < number_of_branch_snps; snp_counter++)
        {
            if(snp_site_coords[snp_counter] >= block_coordinates[0][i] )
            {
                block_coordinates[0][i] = snp_site_coords[snp_counter];
                break;
            }
        }

        for(snp_counter = number_of_branch_snps-1; snp_counter >= 0 ; snp_counter--)
        {
            if(snp_site_coords[snp_counter] <= block_coordinates[1][i] )
            {
                block_coordinates[1][i] = snp_site_coords[snp_counter];
                break;
            }
        }
    }

    free(gaps_in_original_genome_space);
    free(window_count);
    return number_of_blocks;

}



void move_blocks_inwards_while_likelihood_improves(int number_of_blocks,int ** block_coordinates, int min_snps, int * snp_site_coords,  int number_of_branch_snps,char * branch_snp_sequence, int * snp_locations, int branch_genome_size,char * child_sequence, int current_num_snps, double * block_likelihoods, int cutoff_value, float trimming_ratio)
{
	int i;
	
	int previous_start = -1;
	int previous_end = -1;
	
	for(i = 0 ; i < number_of_blocks; i++)
	{

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
		double current_block_likelihood = 0;
		int block_genome_size_without_gaps = block_coordinates[3][i];
		int block_snp_count;
		
        int next_start_position = current_start;
		int start_index = find_starting_index( current_start, snp_site_coords,0, number_of_branch_snps);
        int end_index =  find_starting_index( current_end, snp_site_coords, start_index, number_of_branch_snps);

		block_snp_count = find_number_of_snps_in_block_with_start_end_index(current_start, current_end, snp_site_coords, branch_snp_sequence, number_of_branch_snps,start_index,end_index);
    
		if(block_genome_size_without_gaps == -1)
		{
			
			block_genome_size_without_gaps = calculate_block_size_without_gaps(child_sequence, snp_locations, current_start, current_end, current_num_snps);
			block_coordinates[2][i] = current_block_likelihood;
			block_coordinates[3][i] = block_genome_size_without_gaps;
		}
		current_block_likelihood = get_block_likelihood(branch_genome_size, number_of_branch_snps, block_genome_size_without_gaps, block_snp_count);


		// Move left inwards while the likelihood gets better
        // to avoid the current asymmetry in which the two edges are treated, we
        // should trim SNPs from each edge in an alternating order
		while(current_start < current_end && block_snp_count >= min_snps)
		{
            next_start_position++;
            next_start_position = advance_window_start_to_next_snp(next_start_position, snp_site_coords, branch_snp_sequence, number_of_branch_snps);
			
            if(next_start_position >= current_end)
            {
                break;
            }
            if(next_start_position <= current_start)
            {
                break;
            }
			
            int previous_block_snp_count = block_snp_count;
            int previous_block_genome_size_without_gaps = block_genome_size_without_gaps;
            block_snp_count = find_number_of_snps_in_block_with_start_end_index(next_start_position, current_end, snp_site_coords, branch_snp_sequence, number_of_branch_snps,start_index,end_index);
            block_genome_size_without_gaps = calculate_block_size_without_gaps(child_sequence, snp_locations, next_start_position, current_end, current_num_snps);
			  
            double next_block_likelihood = get_block_likelihood(branch_genome_size, number_of_branch_snps, block_genome_size_without_gaps, block_snp_count);
				
            if(next_block_likelihood*trimming_ratio <= current_block_likelihood && block_snp_count >= min_snps)
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
		while(current_start < current_end && block_snp_count >= min_snps)
		{
            next_end_position--;
            next_end_position = rewind_window_end_to_last_snp(next_end_position, snp_site_coords, branch_snp_sequence, number_of_branch_snps);
				
            if(next_end_position <= current_start )
            {
                break;
            }
            if(next_end_position >= current_end)
            {
                break;
            }
				
            int previous_block_snp_count = block_snp_count;
            int previous_block_genome_size_without_gaps = block_genome_size_without_gaps;
            block_snp_count = find_number_of_snps_in_block(current_start, next_end_position, snp_site_coords, branch_snp_sequence, number_of_branch_snps);
            block_genome_size_without_gaps = calculate_block_size_without_gaps(child_sequence, snp_locations, current_start, next_end_position, current_num_snps);
			  
            double next_block_likelihood = get_block_likelihood(branch_genome_size, number_of_branch_snps, block_genome_size_without_gaps, block_snp_count);
            if(next_block_likelihood*trimming_ratio <= current_block_likelihood && block_snp_count >= min_snps)
            {
                current_block_likelihood = next_block_likelihood;
                current_end = next_end_position;
                end_index--;
            }
            else
            {
                block_snp_count = previous_block_snp_count;
                block_genome_size_without_gaps = previous_block_genome_size_without_gaps;
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
	for(i=number_of_branch_snps_excluding_block; i< number_of_branch_snps; i++)
	{
		snp_site_coords[i] = 0;
	}
	
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

		print_block_details(block_file_pointer, candidate_blocks[0][smallest_index], candidate_blocks[1][smallest_index],  number_of_recombinations_in_window, current_node->taxon,  root->taxon, current_node->taxon_names, current_node->childNum, block_likelihooods[smallest_index]);
		print_gff_line(gff_file_pointer, candidate_blocks[0][smallest_index], candidate_blocks[1][smallest_index],  number_of_recombinations_in_window, current_node->taxon,  root->taxon, current_node->taxon_names, block_likelihooods[smallest_index]);
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

// calculate an approximate p value threshold corrected to multiple testing
double calculate_threshold(int branch_genome_size, int window_size, float uncorrected_p_value)
{
	return 1-(uncorrected_p_value/((branch_genome_size*1.0)/(window_size*1.0)));
}

int calculate_cutoff(int branch_genome_size, int window_size, int num_branch_snps, int min_snps, float uncorrected_p_value)
{
	double threshold = 0.0;
	int cutoff = 0;
	double pvalue = 0.0;
	double part1, part2, part3 = 0.0;
	
	threshold = calculate_threshold(branch_genome_size,
                                    window_size,
                                    uncorrected_p_value);
	
	while(pvalue <= threshold)
	{
		part1 = reduce_factorial(window_size,cutoff)-reduce_factorial(cutoff,cutoff);
		part2 = log10((num_branch_snps*1.0)/branch_genome_size)*cutoff;
		part3 = log10(1.0-((num_branch_snps*1.0)/branch_genome_size))*(window_size-cutoff);
		pvalue = pvalue + pow(10,(part1 + part2 + part3));
		cutoff++;
	}
	cutoff--;

    if (cutoff < min_snps)
    {
        cutoff = min_snps - 1;
    }
    
    // End if the SNP density of the branch is too high for the specified window size
    if (cutoff >= 2*(int)(window_size/2)) // Account for integer division/rounding in this condition
    {
        return 0; // In this case, it is impossible to call recombinations on the branch
    }
    
//    printf("Window size %i cutoff %i num_snps %i\n", window_size,cutoff,num_branch_snps);
    
	return cutoff;
}

int p_value_test(int branch_genome_size, int window_size, int num_branch_snps, int block_snp_count, int min_snps, float uncorrected_p_value)
{
	double threshold = 0.0;
	int cutoff = 0;
	double pvalue = 0.0;
	double part1, part2, part3 = 0.0;
	
	if( block_snp_count < min_snps)
	{
		return 0;	
	}
	
	threshold = uncorrected_p_value/branch_genome_size;
	
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

int calculate_genome_length_excluding_blocks_and_gaps(char * sequence, int length_of_sequence, int ** block_coordinates, int num_blocks)
{
    int genome_length = length_of_sequence;
    
    // Process list of recombinations
    int j = 0;
    int * filtered_start_coords;
    int * filtered_end_coords;
    int num_filtered_blocks = 0;
    filtered_start_coords = (int *) calloc(num_blocks,sizeof(int));
    filtered_end_coords = (int *) calloc(num_blocks,sizeof(int));
    for(j = 0; j<num_blocks; j++)
    {
        if(block_coordinates[0][j] != -1)
        {
          filtered_start_coords[num_filtered_blocks] = block_coordinates[0][j];
          filtered_end_coords[num_filtered_blocks] = block_coordinates[1][j];
          num_filtered_blocks++;
        }
    }
  
    int i = 0;
    for(i = 0; i<length_of_sequence; i++)
    {
        if(sequence[i] == 'N' || sequence[i] == '-')
        {
            genome_length--;
        }
        else
        {
            for(j = 0; j<num_filtered_blocks; j++)
            {
                if (i >= filtered_start_coords[j] && i <= filtered_end_coords[j])
                {
                    genome_length--;
                }
            }
        }
    }

    return genome_length;
}
