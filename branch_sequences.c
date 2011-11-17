#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

		// generate pointers for each child seuqn
		
		while (child != NULL)
		{
			// recursion
			child_sequences[child_counter] = generate_branch_sequences(child->node, vcf_file_pointer, snp_locations, number_of_snps, column_names, number_of_columns, reference_bases, child_sequences[child_counter],length_of_original_genome);
			
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
			
			number_of_branch_snps = find_branch_snp_sites(leaf_sequence, child_sequences[current_branch], snp_locations,number_of_snps, branches_snp_sites[current_branch]);

			identify_recombinations(number_of_branch_snps, branches_snp_sites[current_branch]);
		}
		
		return leaf_sequence;
	}
}

int find_branch_snp_sites(char * ancestor_sequence, char * child_sequence, int * snp_locations, int number_of_snps, int * branch_snp_sites)
{
	int i ;
	int number_of_branch_snp_sites = 0;
	
	for(i = 0; i< number_of_snps; i++)
	{
		branch_snp_sites[i] = 0;
		
		if(ancestor_sequence[i] == '\0' || child_sequence[i] == '\0')
		{
			break;
		}
		if(ancestor_sequence[i] != child_sequence[i])
		{
			
			branch_snp_sites[number_of_branch_snp_sites] = snp_locations[i];
			number_of_branch_snp_sites++;
		}
	}
	realloc(branch_snp_sites,number_of_branch_snp_sites*sizeof(int));
	
	return number_of_branch_snp_sites;
}


void identify_recombinations(int number_of_branch_snps, int * branches_snp_sites)
{
	int i;
	if(number_of_branch_snps < MIN_SNPS_FOR_IDENTIFYING_RECOMBINATIONS)
	{
		return;	
	}
	
	for(i=0; i< number_of_branch_snps; i++)
	{
		double density;
		density = calculate_snp_density(branches_snp_sites, number_of_branch_snps, i);
	}
}

// Assume that the branches_snps_sites array is sorted from smallest to largest
double calculate_snp_density(int * branches_snp_sites, int number_of_branch_snps, int index)
{
	int starting_index = index;
	int ending_index = index;
	if( index > 0)
	{
		starting_index = index-1;
	}
	if(index<(number_of_branch_snps - 1) )
	{
		ending_index = index+1;
	}
	
	if(ending_index == starting_index)
	{
		return DEFAULT_SNP_DENSITY;	
	}
	int index_interval = ending_index - starting_index;
	int coordinate_interval = branches_snp_sites[ending_index] - branches_snp_sites[starting_index];
	if(coordinate_interval > MAX_WINDOW)
	{
		return DEFAULT_SNP_DENSITY;
	}
		
	return ((index_interval)*1.0)/((coordinate_interval)*1.0);
}



// Get the total genome length
// for each branch sequence, go through the snp sites. if its a - reduce the genome length, otherwise increment the snps. The final genome length = N
// and the number of snps is = C

// create sliding windows containing 10 snps (which are not gaps)
// adjust the coordinates of the snps so that gaps are eliminated. keep a lookup table of coordinates to gapless coordinates
// then for each window do the same to get n and c



/*
def get_block_likelihood(start, end, binsnps, N, C):



n=0.0
c=0.0
for x in binsnps[start:end+1]:
if x==0:
n+=1
elif x==1:
c+=1
n+=1

#print start, end, c, n, C, N
part1=math.log((c/n),10)*c
if n-c==0:
part2=0
else:
part2=math.log((((n-c)/n)),10)*(n-c)
if C-c==0:
part3=0
else:
part3=math.log((((C-c)/(N-n))),10)*(C-c)
if ((N-n)-(C-c))==0:
part4=0
else:
part4=math.log(((((N-n)-(C-c))/(N-n))),10)*((N-n)-(C-c))

likelihood=(part1+part2+part3+part4)*-1

#print start, end, c, n, C, N, likelihood

return likelihood

*/


















