/*
 *  Wellcome Trust Sanger Institute
 *  Copyright (C) 2012  Wellcome Trust Sanger Institute
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

#ifndef _GFF_FILE_H_
#define _GFF_FILE_H_
 void print_gff_header(FILE * gff_file_pointer, int genome_length);
 void print_gff_line(FILE * gff_file_pointer, int start_coordinate, int end_coordinate, int number_of_snps, char * current_node_id, char * parent_node_id, char * taxon_names, double  neg_log_likelihood);
#endif

