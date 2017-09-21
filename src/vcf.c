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
#include "vcf.h"
#include "alignment_file.h"
#include "snp_sites.h"
#include "parse_phylip.h"
#include "string_cat.h"


void create_vcf_file(char filename[], int snp_locations[],int number_of_snps, char ** bases_for_snps, char ** sequence_names, int number_of_samples,int internal_nodes[], int offset, int length_of_original_genome)
{
	FILE *vcf_file_pointer;
	char * base_filename;
	base_filename = (char *) calloc((1024 +1),sizeof(char));
	memcpy(base_filename, filename, (1024+1)*sizeof(char));
	char extension[5] = {".vcf"};
	concat_strings_created_with_malloc(base_filename,extension);
	vcf_file_pointer=fopen(base_filename, "w");
	output_vcf_header(vcf_file_pointer,sequence_names, number_of_samples,internal_nodes,length_of_original_genome);
	output_vcf_snps(vcf_file_pointer, bases_for_snps, snp_locations, number_of_snps, number_of_samples,internal_nodes,offset);
    fclose(vcf_file_pointer);
	free(base_filename);
}

void output_vcf_snps(FILE * vcf_file_pointer, char ** bases_for_snps, int * snp_locations, int number_of_snps, int number_of_samples,int internal_nodes[], int offset)
{
	int i;
	for(i=0; i < number_of_snps; i++)
	{
		output_vcf_row(vcf_file_pointer, bases_for_snps[i], snp_locations[i], number_of_samples,internal_nodes, offset);
	}
}

void output_vcf_header( FILE * vcf_file_pointer, char ** sequence_names, int number_of_samples,int internal_nodes[], int length_of_original_genome)
{
	int i;
	fprintf( vcf_file_pointer, "##fileformat=VCFv4.2\n" );
	fprintf( vcf_file_pointer, "##contig=<ID=1,length=%d>\n",length_of_original_genome );
	fprintf( vcf_file_pointer, "##FORMAT=<ID=AB,Number=1,Type=String,Description=\"Alt Base\">\n" );
	
	fprintf( vcf_file_pointer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" );
	
	for(i=0; i<number_of_samples; i++)
	{
		fprintf( vcf_file_pointer, "\t");
		if(internal_nodes[i] == 1)
		{
			continue;
		}
		fprintf( vcf_file_pointer, "%s",  sequence_names[i]);
	}
	fprintf( vcf_file_pointer, "\n");
}

void output_vcf_row(FILE * vcf_file_pointer, char * bases_for_snp, int snp_location, int number_of_samples,int internal_nodes[], int offset)
{
	char reference_base =  bases_for_snp[0];
	char alt_bases[30];
	if(reference_base == '\0')
	{
		return;	
	}
	
	// Chromosome
	fprintf( vcf_file_pointer, "1\t");
	
	// Position
	fprintf( vcf_file_pointer, "%d\t", (int) snp_location + offset);	
	
	//ID
	fprintf( vcf_file_pointer, ".\t");
	
	// REF
	fprintf( vcf_file_pointer, "%c\t", (char) reference_base );
	
	// ALT
	// Need to look through list and find unique characters
	
	alternative_bases(reference_base, bases_for_snp, alt_bases, number_of_samples);
	fprintf( vcf_file_pointer, "%s\t", alt_bases);
	
	// QUAL
	fprintf( vcf_file_pointer, ".\t");
	
	// FILTER
	fprintf( vcf_file_pointer, ".\t");
	
	// INFO
	fprintf( vcf_file_pointer, ".\t");
	
	// FORMAT
	fprintf( vcf_file_pointer, "GT\t");
	
	// Bases for each sample
	output_vcf_row_samples_bases(vcf_file_pointer, reference_base, alt_bases, bases_for_snp, number_of_samples,internal_nodes );
	
	fprintf( vcf_file_pointer, "\n");	
}


void alternative_bases(char reference_base, char * bases_for_snp, char alt_bases[], int number_of_samples)
{
	int i;
	int num_alt_bases = 0;
	for(i=0; i< number_of_samples; i++ )
	{
		if((bases_for_snp[i] != reference_base)  )
		{
			if(check_if_char_in_string(alt_bases, bases_for_snp[i], num_alt_bases) == 0)
			{
				alt_bases[num_alt_bases] = bases_for_snp[i];
				num_alt_bases++;
				alt_bases[num_alt_bases] = ',';
				num_alt_bases++;
			}
		}
	}
	if(num_alt_bases > 0 && alt_bases[num_alt_bases-1] == ',')
	{
		alt_bases[num_alt_bases-1] = '\0';
	}
	else
	{
		alt_bases[num_alt_bases] = '\0';
	}
}

int check_if_char_in_string(char search_string[], char target_char, int search_string_length)
{
	int i;
	for(i=0; i < search_string_length ; i++ )
	{
		if(search_string[i] == target_char)
		{
			return 1;
		}
	}
	return 0;
}

// One indexed. String must be null terminated
int check_where_char_in_string(char search_string[], char target_char)
{
	int i;
	while(search_string[i] != '\0')
	{
		if(search_string[i] == target_char)
		{
			return i+1;
		}
      i++;
	}
	return 0;
}


void output_vcf_row_samples_bases(FILE * vcf_file_pointer, char reference_base, char alt_bases[], char * bases_for_snp, int number_of_samples,int internal_nodes[])
{
	int i;
	
	for(i=0; i < number_of_samples ; i++ )
	{
		if(internal_nodes[i] == 1)
		{
			continue;
		}
		if(bases_for_snp[i] == reference_base)
		{
			fprintf( vcf_file_pointer, "%c", (char) '0' );
		}
		else
		{
         fprintf( vcf_file_pointer, "%c", (char) check_where_char_in_string(alt_bases, bases_for_snp[i]));
		}
		if(i+1 != number_of_samples)
		{
			fprintf( vcf_file_pointer, "\t");
		}
	}
}




