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

#ifndef _PARSE_PHYLIP_H_
#define _PARSE_PHYLIP_H_

void load_sequences_from_phylib(FILE * phylip_file_pointer);
int get_number_of_samples_from_phylip(char * phylip_string);
int get_number_of_snps_from_phylip(char * phylip_string);
void load_sequences_from_phylib_file(char phylip_filename[]);

#define MAX_READ_BUFFER 1048576
#define MAX_SAMPLE_NAME_SIZE 1024

#endif




