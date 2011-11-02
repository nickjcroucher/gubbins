/*

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "vcf.h"
#include "alignment_file.h"
#include "snp_sites.h"






void create_vcf_file(char filename[],  FILE * alignment_file_pointer, int snp_locations[], int length_of_genome, int number_of_snps)
{
	FILE *vcf_file_pointer;
	int number_of_samples;
	
	// TODO chunk up to reduce memory usage
	
	vcf_file_pointer=fopen(strcat(filename,".vcf"), "w");
	output_vcf_header(vcf_file_pointer);
	
	// store values for each snp location
	rewind(alignment_file_pointer);
	
	number_of_samples = count_lines_in_file(alignment_file_pointer)/2;
	char* bases_for_snps[number_of_snps];
	
	for(int i = 0; i < number_of_snps; ++i)
	{
		bases_for_snps[i] = malloc(number_of_samples*sizeof(char));
	}
	
	get_bases_for_each_snp(alignment_file_pointer, snp_locations, bases_for_snps, length_of_genome);
	output_vcf_snps(vcf_file_pointer, bases_for_snps, snp_locations);
}

void output_vcf_snps(FILE * vcf_file_pointer, char ** bases_for_snps, int * snp_locations)
{
	int i;
	for(i=0; bases_for_snps[i]; i++)
	{

		output_vcf_row(vcf_file_pointer, bases_for_snps[i], snp_locations[i]);
	}
}

void output_vcf_header( FILE * vcf_file_pointer)
{
	// populate the sample names
	fprintf( vcf_file_pointer, "##fileformat=VCFv4.1\n" );	
	fprintf( vcf_file_pointer, "##INFO=<ID=AB,Number=1,Type=String,Description=\"Alt Base\">\n" );
	fprintf( vcf_file_pointer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n" );
}

void output_vcf_row(FILE * vcf_file_pointer, char * bases_for_snp, int snp_location)
{
	char reference_base =  bases_for_snp[0];
	
	// Chromosome
	fprintf( vcf_file_pointer, "1\t");
	
	// Position
	fprintf( vcf_file_pointer, "%d\t", (int) snp_location );	
	
	//ID
	fprintf( vcf_file_pointer, ".\t");
	
	// REF
	fprintf( vcf_file_pointer, "%c\t", (char) reference_base );
	
	// ALT
	// Need to look through list and find unique characters
	fprintf( vcf_file_pointer, ".\t");
	
	// QUAL
	fprintf( vcf_file_pointer, ".\t");
	
	// FILTER
	fprintf( vcf_file_pointer, ".\t");
	
	// FORMAT
	fprintf( vcf_file_pointer, "AB\t");
	
	// INFO
	fprintf( vcf_file_pointer, ".\t");
	
	// Bases for each sample
	output_vcf_row_samples_bases(vcf_file_pointer,reference_base, bases_for_snp );
	
	fprintf( vcf_file_pointer, "\n");	
}


void output_vcf_row_samples_bases(FILE * vcf_file_pointer, char reference_base, char * bases_for_snp)
{
	int i;
	for(i=0; bases_for_snp[i]; i++ )
	{
		if(bases_for_snp[i] == reference_base)
		{
			fprintf( vcf_file_pointer, "." );	
		}
		else
		{
			fprintf( vcf_file_pointer, "%c", (char) bases_for_snp[i] );	
		}
		fprintf( vcf_file_pointer, "\t");
	}
}



