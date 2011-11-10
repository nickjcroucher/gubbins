#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gubbins.h"
#include "parse_vcf.h"

// get reference sequence from VCF, and store snp locations
// given a sample name extract the sequences from the vcf
// compare two sequences to get pseudo sequnece and fill in with difference from reference sequence

void run_gubbins(char vcf_filename[], char tree_filename[])
{
	FILE *vcf_file_pointer;
	vcf_file_pointer=fopen(vcf_filename, "r");
	
	char * reference_bases;
	int * snp_locations;
	int number_of_snps;
	int number_of_columns;
	int i;
	int reference_column_number;
	
	number_of_columns = get_number_of_columns_from_file(vcf_file_pointer);
	char* column_names[number_of_columns];
	for(i = 0; i < number_of_columns; i++)
	{
		column_names[i] = malloc(100*sizeof(char));
	}
	get_column_names(vcf_file_pointer, column_names, number_of_columns);
	
	number_of_snps  = get_number_of_snps(vcf_file_pointer);
	snp_locations   = (int *)  malloc(number_of_snps*sizeof(int));
	reference_bases = (char *) malloc(number_of_snps*sizeof(char));
	// get reference sequence from VCF
	reference_column_number = column_number_for_column_name(column_names, "REF", number_of_columns);
	get_sequence_from_column_in_vcf(vcf_file_pointer, snp_locations, reference_bases, number_of_snps,reference_column_number);
	
	printf("%s\n", reference_bases);
	
	
	
}


















