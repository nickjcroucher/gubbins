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

#ifndef _BLOCK_TAB_FILE_H_
#define _BLOCK_TAB_FILE_H_

void print_block_details(FILE * block_file_pointer, int start_coordinate, int end_coordinate, int number_of_snps, char * current_node_id, char * parent_node_id, char * taxon_names);
void print_branch_snp_details(FILE * branch_snps_file_pointer, char * current_node_id, char * parent_node_id, int * branches_snp_sites, int number_of_branch_snps, char * branch_snp_sequence, char * branch_snp_ancestor_sequence,char * taxon_names);
#endif