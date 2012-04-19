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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parse_phylip.h"
#include "tree_statistics.h"

void create_tree_statistics_file(char filename[], sample_statistics ** statistics_for_samples, int number_of_samples)
{
	FILE *file_pointer;
	int sample_counter;
	char * base_filename;
	
	base_filename = (char *) malloc(MAX_FILE_NAME_SIZE*sizeof(char));
	strcpy(base_filename, filename);
	
	file_pointer = fopen(strcat(base_filename,".stats"), "w");
	
	for(sample_counter=0; sample_counter< number_of_samples; sample_counter++)
	{
		sample_statistics * sample_details = ((sample_statistics *) statistics_for_samples[sample_counter]);
		
		fprintf( file_pointer, "%s\t", sample_details->sample_name);
    fprintf( file_pointer, "%i\t", sample_details->number_of_recombinations);
    fprintf( file_pointer, "%i\t", sample_details->number_of_snps);
    fprintf( file_pointer, "%i\t", sample_details->genome_length_without_gaps);
		fprintf( file_pointer, "%f", recombination_to_mutation_ratio(sample_details->number_of_recombinations, sample_details->number_of_snps));
		
		fprintf( file_pointer, "\n");
	}
  fclose(file_pointer);
}

float recombination_to_mutation_ratio(int number_of_recombinations, int number_of_snps)
{
	return (number_of_recombinations*1.0)/(number_of_snps*1.0);
}
