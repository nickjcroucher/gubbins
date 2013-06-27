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
#include <regex.h>
#include "fasta_of_snp_sites.h"
#include "parse_phylip.h"
#include "string_cat.h"

void create_fasta_of_snp_sites(char filename[], int number_of_snps, char ** bases_for_snps, char ** sequence_names, int number_of_samples,int internal_nodes[])
{
	FILE *fasta_file_pointer;
	int sample_counter;
	int snp_counter; 
	char * base_filename;
	
	base_filename = (char *) malloc(1024*sizeof(char));
	memcpy(base_filename, filename, 1024*sizeof(char));
	
	concat_strings_created_with_malloc(base_filename,".snp_sites.aln");
	fasta_file_pointer = fopen(base_filename, "w");
	
	for(sample_counter=0; sample_counter< number_of_samples; sample_counter++)
	{
		if(internal_nodes[sample_counter] == 1)
		{
			continue;
		}
		fprintf( fasta_file_pointer, ">%s\n", sequence_names[sample_counter]);
		for(snp_counter=0; snp_counter< number_of_snps; snp_counter++)
		{
			if(snp_counter > 0 && snp_counter % FASTA_LINE_LENGTH == 0)
			{
				fprintf( fasta_file_pointer, "\n");
			}
			fprintf( fasta_file_pointer, "%c", bases_for_snps[snp_counter][sample_counter]);
		}
		fprintf( fasta_file_pointer, "\n");
	}
  fclose(fasta_file_pointer);
	free(base_filename);
}