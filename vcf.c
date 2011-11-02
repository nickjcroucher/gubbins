/*

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "vcf.h"
#include "alignment_file.h"
#include "snp_sites.h"



void output_vcf_header( FILE * vcf_file_pointer)
{
	fprintf( vcf_file_pointer, "##fileformat=VCFv4.1\n" );	
	fprintf( vcf_file_pointer, "##INFO=<ID=AB,Number=1,Type=String,Description=\"Alt Base\">\n" );
	fprintf( vcf_file_pointer, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\n" );
}


void create_vcf_file(char filename[],  FILE * alignment_file_pointer, int snp_locations[], int length_of_genome)
{
	FILE *vcf_file_pointer;
	int number_of_samples;
	int number_of_snps;
	number_of_snps	= *snp_locations;
	char* bases_for_snps[number_of_snps];
	
	// TODO chunk up to reduce memory usage
	
	vcf_file_pointer=fopen(strcat(filename,".vcf"), "w");
	output_vcf_header(vcf_file_pointer);
	
	// store values for each snp location
	rewind(alignment_file_pointer);
	
	number_of_samples = count_lines_in_file(alignment_file_pointer)/2;
	
	for(int i = 0; i < number_of_snps; ++i)
	{
		bases_for_snps[i] = malloc(number_of_samples*sizeof(char));
	}
	
	get_bases_for_each_snp(alignment_file_pointer, snp_locations, bases_for_snps,length_of_genome);
	output_vcf_snps(vcf_file_pointer,bases_for_snps);
	
	
}

void output_vcf_snps(FILE * vcf_file_pointer, char ** bases_for_snps)
{
	int i;
	for(i=0; bases_for_snps[i]; i++)
	{
		fprintf( vcf_file_pointer, bases_for_snps[i] );	
		//output_vcf_row(vcf_file_pointer, bases_for_snps[i]);
	}
}

void output_vcf_row(FILE * vcf_file_pointer, char * bases_for_snp)
{
	// TODO fix the formatting to match the VCF spec
	fprintf( vcf_file_pointer, bases_for_snp );	
}



