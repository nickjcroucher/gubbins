#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seqUtil.h"
#include "Newickform.h"
#include "branch_sequences.h"
#include "gubbins.h"
#include "parse_vcf.h"

char *generate_branch_sequences(newick_node *root, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns, char reference_bases, char * leaf_sequence)
{
	newick_child *child;
	int child_counter = 0;
	int number_of_children = root->childNum;
	int i;
	int j;
	
	if (number_of_children == 0)
	{
		leaf_sequence = (char *) malloc(number_of_snps*sizeof(char));
		get_sequence_from_column_in_vcf(vcf_file_pointer,  leaf_sequence, number_of_snps,  column_number_for_column_name(column_names, root->taxon, number_of_columns));
		
		return leaf_sequence;
	}
	else
	{
		child = root->child;
		char * child_sequences[number_of_children];

		// generate pointers for each child seuqn
		
		while (child != NULL)
		{
			// recursion
			child_sequences[child_counter] = generate_branch_sequences(child->node, vcf_file_pointer, snp_locations, number_of_snps, column_names, number_of_columns, reference_bases, child_sequences[child_counter]);
			
			child = child->next;
			child_counter++;
		}
		
		leaf_sequence = (char *) malloc(number_of_snps*sizeof(char));
		if (root->taxon != NULL)
		{
			// this non leaf node has its own sequence
			get_sequence_from_column_in_vcf(vcf_file_pointer, leaf_sequence, number_of_snps,  column_number_for_column_name(column_names, root->taxon, number_of_columns));
		}
		else
		{
			// All child sequneces should be available use them to find the ancestor sequence
			leaf_sequence = calculate_ancestor_sequence(leaf_sequence, child_sequences, number_of_snps, number_of_children);
		}
		
		int * branches_snp_sites[number_of_children];
		
		for(i = 0 ; i< number_of_children; i++)
		{
			int number_of_branch_snps=0;
			number_of_branch_snps = find_branch_snp_sites(leaf_sequence, child_sequences[i], snp_locations,number_of_snps, branches_snp_sites[i]);
			for(j = 0; j < number_of_branch_snps; j++)
			{
				printf("%d\t",branches_snp_sites[i][j]);
			}
			printf("\n");
		}
		
		return leaf_sequence;
		//free(child_sequences);
	}
}

int find_branch_snp_sites(char * ancestor_sequence, char * child_sequence, int * snp_locations, int number_of_snps, int * branch_snp_sites)
{
	int i ;
	int number_of_branch_snp_sites = 0;
	branch_snp_sites[number_of_snps];
	
	for(i = 0; i< number_of_snps; i++)
	{
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
	
	return number_of_branch_snp_sites;
}


























