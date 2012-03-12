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
char *generate_branch_sequences(newick_node *root, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns, char * reference_bases, char * leaf_sequence, int length_of_original_genome,  FILE * block_file_pointer);
void identify_recombinations(int number_of_branch_snps, int * branches_snp_sites,int length_of_original_genome);
double calculate_snp_density(int * branches_snp_sites, int number_of_branch_snps, int index);
void get_likelihood_for_windows(char * child_sequence, int length_of_sequence, int * snp_site_coords, int branch_genome_size, int number_of_branch_snps, int * snp_locations, newick_node *current_node, FILE * block_file_pointer, newick_node *root, char * branch_snp_sequence);
double get_block_likelihood(int branch_genome_size, int number_of_branch_snps, int block_genome_size_without_gaps, int number_of_block_snps);
int calculate_window_size(int branch_genome_size, int number_of_branch_snps);
int calculate_window_size(int branch_genome_size, int number_of_branch_snps);
double calculate_threshold(int branch_genome_size, int window_size);
int p_value_test(int branch_genome_size, int window_size, int num_branch_snps, int block_snp_count);
double reduce_factorial(int l, int i);
void fill_in_recombinations_with_reference_bases(newick_node *root, int * parent_recombinations, int parent_num_recombinations, char * reference_bases);
int copy_and_concat_integer_arrays(int * array_1, int array_1_size, int * array_2, int array_2_size, int * output_array);
double snp_density(int length_of_sequence, int number_of_snps);
int calculate_cutoff(int branch_genome_size, int window_size, int num_branch_snps);
int merge_adjacent_blocks(int ** block_coordinates, int number_of_blocks);
int get_smallest_log_likelihood(int ** candidate_blocks, int number_of_candidate_blocks);
int exclude_snp_sites_in_block(int window_start_coordinate, int window_end_coordinate, int * snp_site_coords, int number_of_branch_snps);
int flag_smallest_log_likelihood_recombinations(int ** candidate_blocks, int number_of_candidate_blocks, int number_of_branch_snps, int * snp_site_coords, int * recombinations, int number_of_recombinations,newick_node * current_node, FILE * block_file_pointer, newick_node *root );


#define MIN_SNPS_FOR_IDENTIFYING_RECOMBINATIONS 3
#define DEFAULT_SNP_DENSITY 0.000001
#define WINDOW_SNP_MODE_TARGET 10
#define MIN_WINDOW_SIZE 100
#define MAX_WINDOW_SIZE 10000
#define RANDOMNESS_DAMPNER 0.05
#define MAX_SAMPLE_NAME_SIZE 1024

#endif