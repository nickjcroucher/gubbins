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

#ifndef _BRANCH_SEQUENCES_H_
#define _BRANCH_SEQUENCES_H_
#include "seqUtil.h"
#include "Newickform.h"
void generate_branch_sequences(newick_node *root, char ** node_sequences, char ** node_names, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns, int length_of_original_genome, int num_stored_nodes, FILE * block_file_pointer, FILE * gff_file_pointer,int min_snps, FILE * branch_snps_file_pointer, int window_min, int window_max, float uncorrected_p_value, float trimming_ratio, int extensive_search_flag);
void identify_recombinations(int number_of_branch_snps, int * branches_snp_sites,int length_of_original_genome);
double calculate_snp_density(int * branches_snp_sites, int number_of_branch_snps, int index);
void get_likelihood_for_windows(char * child_sequence, int current_num_snps, int * snp_site_coords, int branch_genome_size, int number_of_branch_snps, int * snp_locations, newick_node * current_node, FILE * block_file_pointer, newick_node *root, char * branch_snp_sequence, FILE * gff_file_pointer,int min_snps, int length_of_original_genome, char * original_sequence,int window_min, int window_max, float uncorrected_p_value, float trimming_ratio, int extensive_search_flag);
double get_block_likelihood(int branch_genome_size, int number_of_branch_snps, int block_genome_size_without_gaps, int number_of_block_snps);
int calculate_window_size(int branch_genome_size, int number_of_branch_snps,int window_min, int window_max, int min_snps, int window_factor);
double calculate_threshold(int branch_genome_size, int window_size, float uncorrected_p_value);
int p_value_test(int branch_genome_size, int window_size, int num_branch_snps, int block_snp_count, int min_snps, float uncorrected_p_value);
double reduce_factorial(int l, int i);
void fill_in_recombinations_with_gaps(newick_node *root, int * parent_recombinations, int parent_num_recombinations, int current_total_snps,int num_blocks, int ** current_block_coordinates,int length_of_original_genome,int * snp_locations, int number_of_snps);
int copy_and_concat_integer_arrays(int * array_1, int array_1_size, int * array_2, int array_2_size, int * output_array);
int copy_and_concat_2d_integer_arrays(int ** array_1, int array_1_size, int ** array_2, int array_2_size, int ** output_array);
double snp_density(int length_of_sequence, int number_of_snps);
int calculate_cutoff(int branch_genome_size, int window_size, int num_branch_snps, int min_snps, float uncorrected_p_value);
int get_smallest_log_likelihood(double * candidate_blocks, int number_of_candidate_blocks);
int exclude_snp_sites_in_block(int window_start_coordinate, int window_end_coordinate, int * snp_site_coords, int number_of_branch_snps);
int flag_smallest_log_likelihood_recombinations(int ** candidate_blocks, int number_of_candidate_blocks, int number_of_branch_snps, int * snp_site_coords, int * recombinations, int number_of_recombinations,newick_node * current_node, FILE * block_file_pointer, newick_node *root,int * snp_locations, int total_num_snps, FILE * gff_file_pointer, double * block_likelihooods);
int calculate_number_of_bases_in_recombations_excluding_gaps(int ** block_coordinates, int num_blocks,char * child_sequence, int * snp_locations,int current_total_snps);
void carry_unambiguous_gaps_up_tree(newick_node *root);
void move_blocks_inwards_while_likelihood_improves(int number_of_blocks,int ** block_coordinates, int min_snps, int * snp_site_coords,  int number_of_branch_snps,char * branch_snp_sequence, int * snp_locations, int branch_genome_size,char * child_sequence, int length_of_sequence, double * block_likelihoods, int cutoff_value, float trimming_ratio);
int get_blocks(int ** block_coordinates, int branch_genome_size,int * snp_site_coords,int number_of_branch_snps, int window_size, int cutoff, char * original_sequence, int * snp_locations, int number_of_snps);
int extend_lower_part_of_window(int starting_coord, int initial_min_coord, int genome_size, int8_t * gaps_in_original_genome_space);
int extend_upper_part_of_window(int starting_coord, int initial_max_coord, int genome_size, int8_t * gaps_in_original_genome_space);
int get_list_of_snp_indices_which_fall_in_downstream_recombinations(int ** current_block_coordinates,int num_blocks, int * snp_locations,int current_total_snps, int * snps_in_recombinations);

int calculate_genome_length_excluding_blocks_and_gaps(char * sequence, int length_of_sequence, int ** block_coordinates, int num_blocks);

#define MAX_SAMPLE_NAME_SIZE 1024

#endif
