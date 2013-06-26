/*
 *  Yu-Wei Wu  http://yuweibioinfo.blogspot.com/2008/10/newick-tree-parser-in-c-make-use-of.html
 *  Copyright (C) 2011  Yu-Wei Wu
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

#define __NEWICKFORM_C__

#include "seqUtil.h"
#include "Newickform.h"
#include "branch_sequences.h"
#include "gff_file.h"


#define STR_OUT	"out"

newick_node* build_newick_tree(char * filename, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns,int length_of_original_genome,int min_snps)
{
	int iLen, iMaxLen;
	char *pcTreeStr;
	char *pcInputFile;
	char *pcOutputFile;
	char acStrArray[256];
	newick_node *root;
	
	FILE *f;
	
	// Initialize memory management procedure
	seqMemInit();
	
	pcInputFile = NULL;
	pcOutputFile = NULL;
	
	pcInputFile = filename;

	// Open tree file
	f = fopen(pcInputFile, "r+");
	
	// Read in tree string
	pcTreeStr = NULL;
	iLen = 0;
	iMaxLen = 0;
	while (1)
	{
		memset(acStrArray, '\0', 256);
		fgets(acStrArray, 255, f);
		if (acStrArray[0] == '\0' && feof(f))
		{
			break;
		}
		inputString(acStrArray, &pcTreeStr, &iLen, &iMaxLen);
	}
	fclose(f);
	
	// Parse tree string
	root = parseTree(pcTreeStr);
	
	// output tab file
  FILE * block_file_pointer;
  char block_file_name[MAX_FILENAME_SIZE];
  strcpy(block_file_name, filename);
	block_file_pointer = fopen(strcat(block_file_name,".tab"), "w");
	
	// output tab file
  FILE * branch_snps_file_pointer;
  char branch_snps_file_name[MAX_FILENAME_SIZE];
  strcpy(branch_snps_file_name, filename);
	branch_snps_file_pointer = fopen(strcat(branch_snps_file_name,".branch_snps.tab"), "w");
	
	// output gff file
	FILE * gff_file_pointer;
  char gff_file_name[MAX_FILENAME_SIZE];
  strcpy(gff_file_name, filename);
	gff_file_pointer = fopen(strcat(gff_file_name,".gff"), "w");
	print_gff_header(gff_file_pointer,length_of_original_genome);
	
	char * root_sequence;

	carry_unambiguous_gaps_up_tree(root);
	root_sequence = generate_branch_sequences(root, vcf_file_pointer, snp_locations, number_of_snps, column_names, number_of_columns,root_sequence, length_of_original_genome, block_file_pointer,gff_file_pointer,min_snps,branch_snps_file_pointer);
	free(root_sequence);
	int * parent_recombinations;
	fill_in_recombinations_with_gaps(root, parent_recombinations, 0, 0,0,root->block_coordinates,length_of_original_genome,snp_locations);

	fclose(block_file_pointer);
	fclose(gff_file_pointer);
	fclose(branch_snps_file_pointer);
	return root;
}

char * strip_quotes(char *taxon)
{
	int i = 0;
	int target_i =0;
	char cleaned_taxon[MAX_FILENAME_SIZE];
	while(taxon[i] != '\0')
	{
		if(taxon[i] != '\'')
		{
			cleaned_taxon[target_i] = taxon[i];
			target_i++;
		}
		i++;
	}
	cleaned_taxon[target_i] = '\0';
	strcpy(taxon,cleaned_taxon);
	return taxon;
}

newick_node* parseTree(char *str)
{
	newick_node *node;
	newick_child *child;
	char *pcCurrent;
	char *pcStart;
	char *pcColon = NULL;
	char cTemp;
	int iCount;

	pcStart = str;
	
	if (*pcStart != '(')
	{
		// Leaf node. Separate taxon name from distance. If distance not exist then take care of taxon name only
		pcCurrent = str;
		while (*pcCurrent != '\0')
		{
			if (*pcCurrent == ':')
			{
				pcColon = pcCurrent;
			}
			pcCurrent++;
		}
		node = (newick_node*)seqMalloc(sizeof(newick_node));
		if (pcColon == NULL)
		{
			// Taxon only
			node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
		}
		else
		{
			// Taxon
			*pcColon = '\0';
			node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
			*pcColon = ':';
			// Distance
			pcColon++;
			node->dist = (float)atof(pcColon);
		}
		node->taxon = strip_quotes(node->taxon);
		node->number_of_blocks = 0;
		node->childNum = 0;
	}
	else
	{
		// Create node
		node = (newick_node*)seqMalloc(sizeof(newick_node));
		child = NULL;
		// Search for all child nodes
		// Find all ',' until corresponding ')' is encountered
		iCount = 0;
		pcStart++;
		pcCurrent = pcStart;
		while (iCount >= 0)
		{
			switch (*pcCurrent)
			{
				case '(':
					// Find corresponding ')' by counting
					pcStart = pcCurrent;
					pcCurrent++;
					iCount++;
					while (iCount > 0)
					{
						if (*pcCurrent == '(')
						{
							iCount++;
						}
						else if (*pcCurrent == ')')
						{
							iCount--;
						}
						pcCurrent++;
					}
					while (*pcCurrent != ',' && *pcCurrent != ')')
					{
						pcCurrent++;
					}
 					cTemp = *pcCurrent;
					*pcCurrent = '\0';
					// Create a child node
					if (child == NULL)
					{
						node->child = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum = 1;
						child = node->child;
					}
					else
					{
						child->next = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum++;
						child = child->next;
					}
					child->node = parseTree(pcStart);
					*pcCurrent = cTemp;
					if (*pcCurrent != ')')
					{
						pcCurrent++;
					}
				break;

				case ')':
					// End of tihs tree. Go to next part to retrieve distance
					iCount--;
				break;

				case ',':
					// Impossible separation since according to the algorithm, this symbol will never encountered.
					// Currently don't handle this and don't create any node
				break;

				default:
					// leaf node encountered
					pcStart = pcCurrent;
					while (*pcCurrent != ',' && *pcCurrent != ')')
					{
						pcCurrent++;
					}
					cTemp = *pcCurrent;
					*pcCurrent = '\0';
					// Create a child node
					if (child == NULL)
					{
						node->child = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum = 1;
						child = node->child;
					}
					else
					{
						child->next = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum++;
						child = child->next;
					}
					child->node = parseTree(pcStart);
					*pcCurrent = cTemp;
					if (*pcCurrent != ')')
					{
						pcCurrent++;
					}
				break;
			}
		}

		// If start at ':', then the internal node has no name.
		pcCurrent++;
		if (*pcCurrent == ':')
		{
			pcStart = pcCurrent + 1;
			while (*pcCurrent != '\0' && *pcCurrent != ';')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->dist = (float)atof(pcStart);
			*pcCurrent = cTemp;
		}
		else if (*pcCurrent != ';' && *pcCurrent != '\0')
		{
			// Find ':' to retrieve distance, if any.
			// At this time *pcCurrent should equal to ')'
			pcStart = pcCurrent;
			while (*pcCurrent != ':' && *pcCurrent != '\0' && *pcCurrent != ';')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->taxon = seqMalloc(strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
			node->taxon = strip_quotes(node->taxon);
			
			*pcCurrent = cTemp;
			pcCurrent++;
			pcStart = pcCurrent;
			while (*pcCurrent != '\0' && *pcCurrent != ';')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->dist = (float)atof(pcStart);
			*pcCurrent = cTemp;
		}
	}

	node->number_of_blocks = 0;
	node->total_bases_removed_excluding_gaps = 0;
	node->block_coordinates =  (int **) malloc((3)*sizeof(int *));	
	node->block_coordinates[0] = (int*) malloc((3)*sizeof(int ));
	node->block_coordinates[1] = (int*) malloc((3)*sizeof(int ));

	return node;
}



void print_tree(newick_node *root, FILE * outputfile)
{
	newick_child *child;
	if (root->childNum == 0)
	{
		fprintf(outputfile,"%s:%0.6f", root->taxon, root->dist);
	}
	else
	{
		child = root->child;
		fprintf(outputfile,"(");
		while (child != NULL)
		{
			print_tree(child->node,outputfile);
			if (child->next != NULL)
			{
				fprintf(outputfile,",");
			}
			child = child->next;
		}
		if (root->taxon != NULL)
		{
			fprintf(outputfile,")%s:%0.6f", root->taxon, root->dist);
		}
		else
		{
			fprintf(outputfile,"):%0.6f", root->dist);
		}
	}
	fflush(outputfile);
}

