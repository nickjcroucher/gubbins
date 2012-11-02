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

#ifndef _SNP_SEARCHING_H_
#define _SNP_SEARCHING_H_
int advance_window_start_to_next_snp(int window_start_coordinate, int * snp_locations, char * child_sequence, int number_of_branch_snps);
int find_starting_index(int window_start_coordinate, int * snp_locations, int start_index, int end_index);
int rewind_window_end_to_last_snp(int window_end_coordinate, int * snp_locations, char * child_sequence, int number_of_branch_snps);
int get_window_end_coordinates_excluding_gaps(int window_start_coordinate, int window_size, int * snp_locations, char * child_sequence, int number_of_snps);
int find_number_of_snps_in_block(int window_start_coordinate, int window_end_coordinate, int * snp_locations,  char * child_sequence, int number_of_snps);
int calculate_block_size_without_gaps(char * child_sequence, int * snp_locations, int starting_coordinate, int ending_coordinate,  int length_of_original_genome);
int calculate_size_of_genome_without_gaps(char * child_sequence, int start_index, int length_of_sequence,  int length_of_original_genome);
int calculate_number_of_snps_excluding_gaps(char * ancestor_sequence, char * child_sequence, int child_sequence_size, int * branch_snp_coords, int * snp_locations, char * branch_snp_sequence, char * branch_snp_ancestor_sequence);
int flag_recombinations_in_window(int window_start_coordinate, int window_end_coordinate, int number_of_snps, int * branch_snp_sites, int * recombinations, int number_of_recombinations,int * snp_locations, int total_num_snps);
int find_matching_coordinate_index(int window_start_coordinate, int * snp_sites, int number_of_snps, int starting_index);
#endif