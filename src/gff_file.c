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
#include "gff_file.h"

void print_gff_header(FILE * gff_file_pointer, int genome_length)
{
	fprintf(gff_file_pointer, "##gff-version 3\n");
	fprintf(gff_file_pointer, "##sequence-region SEQUENCE 1 %d\n", genome_length);
	fflush(gff_file_pointer);
}

void print_gff_line(FILE * gff_file_pointer, int start_coordinate, int end_coordinate, int number_of_snps, char * current_node_id, char * parent_node_id, char * taxon_names, double  neg_log_likelihood)
{
	fprintf(gff_file_pointer, "SEQUENCE\tGUBBINS\tCDS\t");
  fprintf(gff_file_pointer, "%d\t",start_coordinate);
  fprintf(gff_file_pointer, "%d\t",end_coordinate);
	fprintf(gff_file_pointer, "0.000\t.\t0\t");
	
	fprintf(gff_file_pointer, "node=\"%s->%s\";", parent_node_id, current_node_id );
	fprintf(gff_file_pointer, "neg_log_likelihood=\"%f\"", neg_log_likelihood);
	fprintf(gff_file_pointer, "taxa=\"%s\";", taxon_names);
	fprintf(gff_file_pointer, "snp_count=\"%d\"", number_of_snps);
	fprintf(gff_file_pointer, "\n");
	
  fflush(gff_file_pointer);
}
