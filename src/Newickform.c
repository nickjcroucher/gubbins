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

#include <pthread.h>
#include "stdio.h"
#include "seqUtil.h"
#include "Newickform.h"
#include "branch_sequences.h"
#include "gff_file.h"
#include "string_cat.h"


#define STR_OUT	"out"

// Define thread function
void* threadFunction(void* arg) {
  
    // Extract thread data
    struct ThreadData* data = (struct ThreadData*)arg;
    
    if (data->num_nodes_to_process > -1)
    {
      for (int node_index = data->start_node; node_index < data->start_node+data->num_nodes_to_process; ++node_index)
      {
        // Generate branch sequences and identify recombinations
        generate_branch_sequences(data->nodes[node_index],
                                   data->vcf_file_pointer,
                                   data->snp_locations,
                                   data->number_of_snps,
                                   data->column_names,
                                   data->number_of_columns,
                                   data->length_of_original_genome,
                                   data->block_file_pointer,
                                   data->gff_file_pointer,
                                   data->min_snps,
                                   data->branch_snps_file_pointer,
                                   data->window_min,
                                   data->window_max,
                                   data->uncorrected_p_value,
                                   data->trimming_ratio,
                                   data->extensive_search_flag,
                                   node_index);
      }
    }

    // Exit thread
    pthread_exit(NULL);
}

// Function to extract nodes relevant for each depth
void get_job_nodes(newick_node** jobNodeArray, newick_node** nodeArray, int* node_depths, int depth, int num_nodes)
{
  int j = 0;
  for (int i = 0; i < num_nodes; ++i)
  {
    if (node_depths[i] == depth)
    {
      jobNodeArray[j] = nodeArray[i]; // TO DO convert to pointer
      ++j;
    }
  }
}

// Function to count number of jobs to run at a particular depth
int get_job_counts(int *node_depths, int depth, int num_nodes)
{
  int count = 0;
  for (int i = 0; i < num_nodes; ++i)
  {
    if (node_depths[i] == depth)
    {
      count++;
    }
  }
  return count;
}

// Function to create an array of nodes
void fill_nodeArray(newick_node *root,newick_node** nodeArray, int num_nodes)
{
  if (root->childNum != 0)
  {
    newick_child *child;
    child = root->child;
    while (child != NULL)
    {
      fill_nodeArray(child->node,nodeArray, num_nodes);
      child = child->next;
    }
  }
  for (int i = 0; i < num_nodes; ++i)
  {
    if (nodeArray[i]->taxon == NULL)
    {
      nodeArray[i] = root;
      break;
    }
  }
}

// Function to count the total number of tree nodes
int count_tree_nodes(newick_node* root) {
    if (root == NULL) return 0;
    int count = 1; // Count the root node
    // Recursively count nodes in the child subtree
    if (root->childNum != 0)
    {
      newick_child *child;
      child = root->child;
      while (child != NULL)
      {
        count += count_tree_nodes(child->node);
        child = child->next;
      }
    }
    return count;
}

// Function to find the maximum distance from an internal node to the tips
int max_distance_to_tips(newick_node *root) {
    
    // Initialize the maximum distance to 0
    int max_distance = 0;
  
    // Distance is non-zero if not leaf node
    if (root->childNum != 0)
    {
      // Traverse each child of the internal node
      newick_child* child = root->child;
      while (child != NULL)
      {
          // Recursively find the maximum distance from the child to the tips
          int child_distance = max_distance_to_tips(child->node);
          // Add one for the distance from the child to the parent
          child_distance++;
          // Update the maximum distance if the distance from this child is greater
          if (child_distance > max_distance) {
              max_distance = child_distance;
          }
          child = child->next;
      }
    }

    return max_distance;
}

newick_node* build_newick_tree(char * filename, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns,int length_of_original_genome,int min_snps, int window_min, int window_max, float uncorrected_p_value, float trimming_ratio, int extensive_search_flag, int num_threads)
{
	int iLen, iMaxLen;
	char *pcTreeStr;
	char *pcInputFile;
	char *pcOutputFile;
	char acStrArray[256];
	newick_node *root;
  char *returnchar;
	
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
		returnchar = fgets(acStrArray, 255, f);
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
  char block_file_name[MAX_FILENAME_SIZE] = {""};
  char block_file_extension[5]= {".tab"};
	memcpy(block_file_name, filename, size_of_string(filename) +1);
	concat_strings_created_with_malloc(block_file_name,block_file_extension);
	block_file_pointer = fopen(block_file_name, "w");
	
	// output tab file
  FILE * branch_snps_file_pointer;
  char branch_snps_file_name[MAX_FILENAME_SIZE]= {""};
  char branchtab_extension[18]= {".branch_snps.tab"};
	memcpy(branch_snps_file_name, filename, size_of_string(filename) +1);
	concat_strings_created_with_malloc(branch_snps_file_name,branchtab_extension);
	branch_snps_file_pointer = fopen(branch_snps_file_name, "w");
	
	// output gff file
	FILE * gff_file_pointer;
  char gff_file_name[MAX_FILENAME_SIZE]= {""};
  memcpy(gff_file_name, filename, size_of_string(filename) +1);
  char gff_extension[5]= {".gff"};
	concat_strings_created_with_malloc(gff_file_name,gff_extension);
	gff_file_pointer = fopen(gff_file_name, "w");
	print_gff_header(gff_file_pointer,length_of_original_genome);
	
  // iterate through tree to get list of nodes
  int num_nodes = 0;
  newick_node * current_node = root;
  num_nodes = count_tree_nodes(current_node);
  newick_node** nodeArray = malloc(num_nodes * sizeof(newick_node*));
  for (int i = 0; i < num_nodes; ++i)
  {
    nodeArray[i] = (newick_node*)seqMalloc(sizeof(newick_node));
  }
  fill_nodeArray(current_node,nodeArray,num_nodes);
  
  // get depths of each node from tips
  int max_depth = max_distance_to_tips(root);
  int *node_depths = malloc(num_nodes * sizeof(int));
  for (int i = 0; i < num_nodes; ++i)
  {
    node_depths[i] = max_distance_to_tips(nodeArray[i]);
  }
  
  // Initiate multithreading
  pthread_t threads[num_threads];
  struct ThreadData threadData[num_threads];

  // iterate through depths and identify batches of analyses to be run
  for (int depth = 0; depth <= max_depth; ++depth)
  {
    
    // Identify number of nodes at the current depth
    int num_jobs = get_job_counts(node_depths,depth,num_nodes);
    newick_node** jobNodeArray = malloc(num_jobs * sizeof(newick_node*));
    // Allocate memory for each element of jobNodeArray
    for (int i = 0; i < num_jobs; i++) {
        jobNodeArray[i] = (newick_node*)seqMalloc(sizeof(newick_node));
    }
    get_job_nodes(jobNodeArray,nodeArray,node_depths,depth,num_nodes);

    // Divide jobNodeArray among threads
    int numJobsPerThread = num_jobs / num_threads;
    int remainder = num_jobs % num_threads;
    
    // Create and execute threads
    for (int i = 0; i < num_threads; ++i) {

        // Calculate start and end indices for current thread
        int startIndex = i * numJobsPerThread + (i < remainder ? i : remainder);
        int endIndex = startIndex + numJobsPerThread + (i < remainder ? 1 : 0) - 1;
      
        // Set thread data
        threadData[i].nodes = jobNodeArray;
        threadData[i].start_node = startIndex;
        threadData[i].num_nodes_to_process = endIndex - startIndex + 1; // Number of nodes for this thread
        threadData[i].vcf_file_pointer = vcf_file_pointer;
        threadData[i].snp_locations = snp_locations;
        threadData[i].number_of_snps = number_of_snps;
        threadData[i].column_names = column_names;
        threadData[i].number_of_columns = number_of_columns;
        threadData[i].length_of_original_genome = length_of_original_genome;
        threadData[i].block_file_pointer = block_file_pointer;
        threadData[i].gff_file_pointer = gff_file_pointer;
        threadData[i].min_snps = min_snps;
        threadData[i].branch_snps_file_pointer = branch_snps_file_pointer;
        threadData[i].window_min = window_min;
        threadData[i].window_max = window_max;
        threadData[i].uncorrected_p_value = uncorrected_p_value;
        threadData[i].trimming_ratio = trimming_ratio;
        threadData[i].extensive_search_flag = extensive_search_flag;
        threadData[i].thread_index = i;

        // Create thread
        if (pthread_create(&threads[i], NULL, threadFunction, (void*)&threadData[i]) != 0) {
            perror("pthread_create");
            exit(EXIT_FAILURE);
        }

    }

    // Join threads
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }
    
    // Free jobNodeArray
    free(jobNodeArray);
    
  }
  
  int * parent_recombinations = NULL;
	fill_in_recombinations_with_gaps(root, parent_recombinations, 0, 0,0,root->block_coordinates,length_of_original_genome,snp_locations,number_of_snps);

  // Free arrays
  free(nodeArray);
  free(node_depths);
  
	fclose(block_file_pointer);
	fclose(gff_file_pointer);
	fclose(branch_snps_file_pointer);
	return root;
}

char * strip_quotes(char *taxon)
{
	int i = 0;
	int target_i =0;
	char cleaned_taxon[MAX_FILENAME_SIZE] = {""};
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
	memcpy(taxon, cleaned_taxon, size_of_string(cleaned_taxon) +1);
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
			node->taxon = (char*)seqMalloc((int)strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
		}
		else
		{
			// Taxon
			*pcColon = '\0';
			node->taxon = (char*)seqMalloc((int)strlen(pcStart) + 1);
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
			node->taxon = seqMalloc((int)strlen(pcStart) + 1);
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
	node->block_coordinates =  (int **) calloc((3),sizeof(int *));	
	node->block_coordinates[0] = (int*) calloc((3),sizeof(int ));
	node->block_coordinates[1] = (int*) calloc((3),sizeof(int ));

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

