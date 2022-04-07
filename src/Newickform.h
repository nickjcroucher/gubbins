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

#define MAX_FILENAME_SIZE 1024

#ifdef __NEWICKFORM_C__
newick_node* parseTree(char *str);
newick_node* build_newick_tree(char * filename, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns,  int length_of_original_genome,int min_snps, int window_min, int window_max, float uncorrected_p_value, float trimming_ratio);
void print_tree(newick_node *root, FILE * outputfile);
char* strip_quotes(char *taxon);
#else
extern newick_node* parseTree(char *str);
extern newick_node* build_newick_tree(char * filename, FILE *vcf_file_pointer,int * snp_locations, int number_of_snps, char** column_names, int number_of_columns, int length_of_original_genome,int min_snps, int window_min, int window_max, float uncorrected_p_value, float trimming_ratio);
extern void print_tree(newick_node *root, FILE * outputfile);
extern char* strip_quotes(char *taxon);
#endif

#endif

