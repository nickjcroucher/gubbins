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

char *generate_branch_sequences(newick_node *root, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns, char reference_bases, char * leaf_sequence, int length_of_original_genome)
{
	newick_child *child;
	int child_counter = 0;
	int current_branch =0;
	
	
	if (root->childNum == 0)
	{
		leaf_sequence = (char *) malloc(number_of_snps*sizeof(char));
		get_sequence_for_sample_name(leaf_sequence, root->taxon);
		
		return leaf_sequence;
	}
	else
	{
		child = root->child;
		char * child_sequences[root->childNum];
		newick_node * child_nodes[root->childNum];
		

		// generate pointers for each child seuqn
		
		while (child != NULL)
		{
			// recursion
			child_sequences[child_counter] = generate_branch_sequences(child->node, vcf_file_pointer, snp_locations, number_of_snps, column_names, number_of_columns, reference_bases, child_sequences[child_counter],length_of_original_genome);
			child_nodes[child_counter] = child->node;
			
			child = child->next;
			child_counter++;
		}
		
		leaf_sequence = (char *) malloc(number_of_snps*sizeof(char));
		if (root->taxon != NULL)
		{
			// this non leaf node has its own sequence
			get_sequence_for_sample_name(leaf_sequence, root->taxon);
		}
		else
		{
			// All child sequneces should be available use them to find the ancestor sequence
			leaf_sequence = calculate_ancestor_sequence(leaf_sequence, child_sequences, number_of_snps, root->childNum);
		}
		
		int * branches_snp_sites[root->childNum];
		
		for(current_branch = 0 ; current_branch< (root->childNum); current_branch++)
		{
			int number_of_branch_snps=0;
			branches_snp_sites[current_branch] = (int *) malloc(number_of_snps*sizeof(int));
			
			int branch_genome_size;
			branch_genome_size = calculate_size_of_genome_without_gaps(child_sequences[current_branch], 0,number_of_snps, length_of_original_genome);
			number_of_branch_snps = calculate_number_of_snps_excluding_gaps(leaf_sequence, child_sequences[current_branch], number_of_snps, branches_snp_sites[current_branch], snp_locations);
			
			get_likelihood_for_windows(child_sequences[current_branch], number_of_snps, branches_snp_sites[current_branch], branch_genome_size, number_of_branch_snps,snp_locations, child_nodes[current_branch]->recombinations);

		}
		
		return leaf_sequence;
	}
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
	realloc(branch_snp_coords, number_of_branch_snp_sites*sizeof(int));
	return number_of_branch_snp_sites;
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
	
	// Convert a double to an int????
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


// inefficient
int get_window_end_coordinates_excluding_gaps(int window_start_coordinate, int window_size, int * snp_locations, char * child_sequence, int number_of_snps)
{
	int i;
	int window_end_coordinate = window_start_coordinate + window_size;
	
	for(i = 0; i < number_of_snps; i++)
	{
		if(snp_locations[i]>= window_start_coordinate && snp_locations[i] < window_end_coordinate)
		{
			if(child_sequence[i] == '-')
			{
				window_end_coordinate++;
			}
		}
	}
	return window_end_coordinate;
}


// inefficient
int find_number_of_snps_in_block(int window_start_coordinate, int window_end_coordinate, int * snp_locations, char * child_sequence, int number_of_snps)
{
	int i;
	int number_of_snps_in_block =0;
	
	for(i = 0; i < number_of_snps; i++)
	{
		if(snp_locations[i]>= window_start_coordinate && snp_locations[i] < window_end_coordinate)
		{
			if(child_sequence[i] != '-')
			{
				number_of_snps_in_block++;
			}
		}
	}
	return number_of_snps_in_block;
}

void get_likelihood_for_windows(char * child_sequence, int length_of_sequence, int * snp_site_coords, int branch_genome_size, int number_of_branch_snps, int * snp_locations, int * recombinations)
{
	int i;
	int window_size;
	int window_start_coordinate = 0;
	int window_end_coordinate = 0;
	int number_of_snps_in_block;
	int block_genome_size_without_gaps;
	double block_likelihood;
	int cutoff;
	int number_of_recombinations = 0;
	
	// place to store coordinates of recombinations snps
	recombinations = (int *) malloc(length_of_sequence*sizeof(int));
	
	if(number_of_branch_snps == 0)
	{
		return;
	}
	
	window_size = calculate_window_size(branch_genome_size, number_of_branch_snps);
	// start at the coordinate of the first snp
	window_start_coordinate = snp_site_coords[0];
	
	cutoff =  calculate_cutoff(branch_genome_size, window_size, number_of_branch_snps);
	
	for(i = 0; i < ceil(branch_genome_size/window_size) && (window_start_coordinate < branch_genome_size); i++)
	{
		window_end_coordinate = get_window_end_coordinates_excluding_gaps(window_start_coordinate, window_size, snp_locations, child_sequence,length_of_sequence);		
		
		number_of_snps_in_block = find_number_of_snps_in_block(window_start_coordinate, window_end_coordinate, snp_locations, child_sequence, length_of_sequence);
		block_genome_size_without_gaps = window_end_coordinate - window_start_coordinate;
		if(number_of_snps_in_block >= MIN_SNPS_FOR_IDENTIFYING_RECOMBINATIONS)
		{
			block_likelihood =  get_block_likelihood(branch_genome_size, number_of_branch_snps, block_genome_size_without_gaps, number_of_snps_in_block);

			// If you've found a recombination, extend out and mark the whole window, then feed this back
			// Fix me
			if(block_likelihood > 9000)
			{
				number_of_recombinations += flag_recombinations_in_window(window_start_coordinate, window_end_coordinate,length_of_sequence, snp_locations, recombinations, number_of_recombinations);
			}
			
			//printf("cutoff: %d\tN: %d\tC: %d\twindowsize: %d\tstart %d\tn: %d\tc: %d\tLH: %f\n",cutoff,branch_genome_size,number_of_branch_snps,window_size,window_start_coordinate,block_genome_size_without_gaps,number_of_snps_in_block, block_likelihood);
		}
		window_start_coordinate = window_end_coordinate;
	}
	realloc(recombinations, number_of_recombinations*sizeof(int));
	
}

// Inefficient
int flag_recombinations_in_window(int window_start_coordinate, int window_end_coordinate, int length_of_sequence, int * snp_locations, int * recombinations, int number_of_recombinations)
{
	int i;
	int number_of_snps_in_block = 0;
	
	for(i = 0; i < length_of_sequence; i++)
	{
		if(snp_locations[i]>= window_start_coordinate && snp_locations[i] < window_end_coordinate)
		{
			recombinations[number_of_recombinations+number_of_snps_in_block] = snp_locations[i];
			number_of_snps_in_block++;
		}
	}
	return number_of_snps_in_block;
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
		part4=log(((branch_genome_size-block_genome_size_without_gaps)-(number_of_branch_snps-number_of_block_snps)*1.0 )/((branch_genome_size-block_genome_size_without_gaps)*1.0)) * ((branch_genome_size-block_genome_size_without_gaps)-(number_of_branch_snps-number_of_block_snps));
	}
	
	return (part1+part2+part3+part4)*-1;
}

// Get the total genome length
// for each branch sequence, go through the snp sites. if its a - reduce the genome length, otherwise increment the snps. The final genome length = N
// and the number of snps is = C

// create sliding windows containing 10 snps (which are not gaps)
// adjust the coordinates of the snps so that gaps are eliminated. keep a lookup table of coordinates to gapless coordinates
// then for each window do the same to get n and c



















