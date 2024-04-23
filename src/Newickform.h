#ifndef __NEWICKFORM_H__
#define __NEWICKFORM_H__

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

typedef struct newick_child
{
	struct newick_node *node;
	struct newick_child *next;
} newick_child;

typedef struct newick_node
{
	char *taxon;
	char *taxon_names;
	char *seq;
	
	float dist;
	int childNum;
	int * recombinations;
	int num_recombinations;
  int number_of_snps;
  int current_node_id;
  int number_of_blocks;
	int total_bases_removed_excluding_gaps;
  int ** block_coordinates;
  
	struct newick_child *child;
	struct newick_node *parent;
} newick_node;

// Define the structure to hold thread arguments
struct rec_ThreadData {
    newick_node** nodes; // Nodes to be processed by all threads
    int start_node;             // Index of starting node for this thread
    int num_nodes_to_process;    // Number of nodes to process by this thread
    char** node_sequences;      // Pointer to the array of node sequences
    char** node_names;          // Pointer to the array of node names
    FILE* vcf_file_pointer;     // Pointer to the VCF file
    int* snp_locations;         // Pointer to the array of SNP locations
    int number_of_snps;         // Number of SNPs
    char** column_names;        // Pointer to the array of column names
    int number_of_columns;      // Number of columns
    int length_of_original_genome;  // Length of the original genome
    int num_stored_nodes;       // Number of stored nodes
    FILE* block_file_pointer;   // Pointer to the block file
    FILE* gff_file_pointer;     // Pointer to the GFF file
    int min_snps;               // Minimum number of SNPs
    FILE* branch_snps_file_pointer; // Pointer to the branch SNPs file
    int window_min;             // Minimum window size
    int window_max;             // Maximum window size
    float uncorrected_p_value;  // Uncorrected p-value
    float trimming_ratio;       // Trimming ratio
    int extensive_search_flag;  // Extensive search flag
    int thread_index;           // Thread index
};


#define MAX_FILENAME_SIZE 1024

#ifdef __NEWICKFORM_C__
newick_node* parseTree(char *str);
newick_node* build_newick_tree(char * filename, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns,  int length_of_original_genome,int min_snps, int window_min, int window_max, float uncorrected_p_value, float trimming_ratio, int extensive_search_flag, int num_threads);
void print_tree(newick_node *root, FILE * outputfile);
char* strip_quotes(char *taxon);
#else
extern newick_node* parseTree(char *str);
extern newick_node* build_newick_tree(char * filename, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns, int length_of_original_genome,int min_snps, int window_min, int window_max, float uncorrected_p_value, float trimming_ratio, int extensive_search_flag, int num_threads);
void* rec_threadFunction(void* arg);
extern void print_tree(newick_node *root, FILE * outputfile);
void fill_nodeArray(newick_node *root, newick_node** nodeArray);
int count_tree_nodes(newick_node* root);
void get_job_nodes(newick_node** jobNodeArray,newick_node** nodeArray,int* node_depths,int depth,int num_nodes);
void get_job_node_indices(int* jobNodeIndexArray, newick_node** nodeArray, int* node_depths, int depth, int num_nodes);
void get_job_counts(int *node_depths, int depth, int num_nodes);
extern char* strip_quotes(char *taxon);
#endif

#endif

