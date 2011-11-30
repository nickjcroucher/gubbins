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
#include "snp_searching.h"

// Most of the methods in this file look the same, so should be DRYed out.

int advance_window_start_to_next_snp(int window_start_coordinate, int * snp_locations, char * child_sequence, int number_of_branch_snps)
{
	int i;
	int start_index = find_starting_index( window_start_coordinate, snp_locations,0, number_of_branch_snps);
	
	for(i = start_index; i < number_of_branch_snps; i++)
	{
		if(snp_locations[i]>= window_start_coordinate && child_sequence[i] != '-')
		{
			return snp_locations[i];
		}
	}
	return window_start_coordinate;
}

int find_starting_index(int window_start_coordinate, int * snp_locations, int start_index, int end_index)
{
	int current_index = 0;
	
	
	if(start_index == end_index  || start_index + 1 == end_index)
	{
		return start_index;
	}
	else if(end_index < start_index)
	{
		return end_index;
	}
	current_index = (int)((end_index-start_index)/2) + start_index;
	
	if( snp_locations[current_index] < window_start_coordinate)
	{
		start_index = current_index;
	}
	else if(  snp_locations[current_index] > window_start_coordinate)
	{
		end_index = current_index;
	}
	else
	{
		return current_index;
	}
	
	return find_starting_index(window_start_coordinate, snp_locations, start_index, end_index);
}

int rewind_window_end_to_last_snp(int window_end_coordinate, int * snp_locations, char * child_sequence, int number_of_branch_snps)
{
	int i;
	int end_index =  find_starting_index( window_end_coordinate, snp_locations, 0, number_of_branch_snps);
	if(end_index +1 <number_of_branch_snps )
	{
		end_index++;
	}
	if(end_index +1 <number_of_branch_snps )
	{
		end_index++;
	}
	
	for(i = end_index; i >= 0; i--)
	{
		if(snp_locations[i]< window_end_coordinate && child_sequence[i] != '-')
		{
			return (snp_locations[i] +1);
		}
		
	}
	return window_end_coordinate;
}

int get_window_end_coordinates_excluding_gaps(int window_start_coordinate, int window_size, int * snp_locations, char * child_sequence, int number_of_snps)
{
	int i;
	int window_end_coordinate = window_start_coordinate + window_size;
	int start_index = find_starting_index( window_start_coordinate, snp_locations, 0, number_of_snps);
	
	for(i = start_index; i < number_of_snps; i++)
	{
		if(snp_locations[i]>= window_start_coordinate && snp_locations[i] < window_end_coordinate)
		{
			if(child_sequence[i] == '-')
			{
				window_end_coordinate++;
			}
		}
		if(snp_locations[i] > window_end_coordinate)
		{
			break;	
		}
	}
	return window_end_coordinate;
}



int find_number_of_snps_in_block(int window_start_coordinate, int window_end_coordinate, int * snp_locations,  char * child_sequence, int number_of_snps)
{
	int i;
	int number_of_snps_in_block =0;
	int start_index = find_starting_index( window_start_coordinate, snp_locations,0, number_of_snps);
	
	for(i = start_index; i < number_of_snps; i++)
	{
		if(snp_locations[i]>= window_start_coordinate && snp_locations[i] < window_end_coordinate)
		{
			number_of_snps_in_block++;
		}
		if(snp_locations[i] > window_end_coordinate)
		{
			break;	
		}
	}
	return number_of_snps_in_block;
}

int calculate_block_size_without_gaps(char * child_sequence, int * snp_locations, int starting_coordinate, int ending_coordinate,  int length_of_original_genome)
{
	int i;
	int block_size = ending_coordinate - starting_coordinate;
	int start_index = find_starting_index( starting_coordinate, snp_locations,0, length_of_original_genome);
	
	for(i = start_index; i < length_of_original_genome ; i++)
	{
		if(snp_locations[i]< ending_coordinate && snp_locations[i]>= starting_coordinate)
		{
			if(child_sequence[i] == '-')
			{
				block_size--;
			}
		}
		
		if( snp_locations[i]> ending_coordinate)
		{
			break;	
		}
		
	}
	return block_size;
}


// take in a sequence, and calculate the size of the genome when gaps are excluded
int calculate_size_of_genome_without_gaps(char * child_sequence, int start_index, int length_of_sequence,  int length_of_original_genome)
{
	int i;
	int total_length_of_genome = length_of_original_genome;
	for(i = start_index; i< (start_index+length_of_sequence) && (i-start_index) < (total_length_of_genome); i++)
	{
		if(child_sequence[i] == '\0')
		{
			break;
		}
		
		if(child_sequence[i] == '-')
		{
			length_of_original_genome--;
		}
	}
	return length_of_original_genome;
}

int calculate_number_of_snps_excluding_gaps(char * ancestor_sequence, char * child_sequence, int child_sequence_size, int * branch_snp_coords, int * snp_locations)
{
	int i ;
	int number_of_branch_snp_sites = 0;
	
	for(i = 0; i< child_sequence_size; i++)
	{
		branch_snp_coords[i] = 0;
		if(ancestor_sequence[i] == '\0' || child_sequence[i] == '\0')
		{
			break;
		}
		
		// If there is a gap in the ancestor, and a base in the child, what happens?
		if(ancestor_sequence[i] != child_sequence[i]  && child_sequence[i] != '-')
		{
			branch_snp_coords[number_of_branch_snp_sites] = snp_locations[i];
			number_of_branch_snp_sites++;
		}
	}	
	branch_snp_coords = realloc(branch_snp_coords, number_of_branch_snp_sites*sizeof(int));
	return number_of_branch_snp_sites;
}


int flag_recombinations_in_window(int window_start_coordinate, int window_end_coordinate, int number_of_snps, int * snp_locations, int * recombinations, int number_of_recombinations)
{
	int i;
	int number_of_snps_in_block = 0;
	int start_index = find_starting_index( window_start_coordinate, snp_locations, 0, number_of_snps);
	
	for(i = start_index; i < number_of_snps; i++)
	{
		if(snp_locations[i]>= window_start_coordinate && snp_locations[i] < window_end_coordinate)
		{
			recombinations[number_of_recombinations+number_of_snps_in_block] = i;
			number_of_snps_in_block++;
		}
		if(snp_locations[i] > window_end_coordinate)
		{
			break;	
		}
	}
	return number_of_snps_in_block;
}


