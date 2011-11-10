#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vcf.h"
#include "parse_vcf.h"
#include "alignment_file.h"

void get_sequence_from_column_in_vcf(FILE * vcf_file_pointer, int snp_locations[], char * sequence_bases, int number_of_snps, int column_number)
{
	rewind(vcf_file_pointer);
	char szBuffer[100000] = {0};  
	int reference_index = 0;
	char result[10000] = {0};  
		
	do{
		strcpy(szBuffer,""); 
		// check the first character of the line to see if its in the header
		read_line(szBuffer, vcf_file_pointer);
		
		if(szBuffer[0] == '\0')
		{
			break;
		}
		
		if(szBuffer[0] != '#')
		{
			char * returned_result ;
			returned_result  = split_string_and_return_specific_index(result, szBuffer, column_number);
			sequence_bases[reference_index] = returned_result[0];
			reference_index++;
		}
		
	}while(szBuffer[0] != '\0');
}

char *split_string_and_return_specific_index(char * result, char * input_string, int token_index)
{
	char delims[] = "\t";
    int i = 0;	
	char tokenised_string[100000] = {0}; 
	strcpy(tokenised_string,input_string);
	
	// split the line and go through the tokens to find the reference char
	result = strtok( tokenised_string, delims );
	if(token_index ==0)
	{
		return result;	
	}
	
	while( result != NULL ) {
		result = strtok( NULL, delims );
		if(i == token_index -1)
		{
			return result;	
		}
		
		i++;
	}
	return "";
}


int get_number_of_snps(FILE * vcf_file_pointer)
{
	rewind(vcf_file_pointer);
	
	int i = 0;
	int length_of_line =0;
	char szBuffer[2] = {0};  
	
	do{
		// check the first character of the line to see if its in the header
		fgets(szBuffer, sizeof(szBuffer), vcf_file_pointer);
		if(szBuffer[0] != '#')
		{
			i++;
		}
		length_of_line = line_length(vcf_file_pointer);
		
	}while(length_of_line != 0);
		
	return (i> 0) ? i-1 : 0;	
}

// Assumes that all column headers have something in them
int get_number_of_columns(char * column_header)
{
	int number_of_columns = 0;
	char result[100] = {0}; 
	char * returned_result;
	
	do
	{
		returned_result = split_string_and_return_specific_index( result, column_header, number_of_columns);
		if(returned_result == NULL ||  returned_result[0] == '\n')
		{
			break;	
		}
		
		number_of_columns++;
	}
	while( ! (strcmp(returned_result,"") == 0) );
	
	return number_of_columns;
}

// Assumes that all column headers have something in them
int get_number_of_columns_from_file(FILE * vcf_file_pointer)
{
	rewind(vcf_file_pointer);
	char szBuffer[100000] = {0};
	char result[100] = {0};
	char * returned_result;
	
	do{
		strcpy(szBuffer,""); 
		// check the first character of the line to see if its in the header
		read_line(szBuffer, vcf_file_pointer);
		
		if(szBuffer[0] == '\0' || szBuffer[0] != '#')
		{
			break;
		}
		
		//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  _S_pneumoniae_Spanis    _3948_7_10
		returned_result = split_string_and_return_specific_index( result, szBuffer, 0);

		if(strcmp(returned_result, "#CHROM")==0)
		{
			return get_number_of_columns(szBuffer);
			
		}
		
	}while(szBuffer[0] != '\0');
	return 0;
}


void get_column_names(FILE * vcf_file_pointer, char ** column_names, int number_of_columns)
{
	rewind(vcf_file_pointer);
	char szBuffer[100000] = {0};
	char result[100] = {0};  
	int i;
	char * returned_result;
	
	do{
		strcpy(szBuffer,""); 
		// check the first character of the line to see if its in the header
		read_line(szBuffer, vcf_file_pointer);
		
		if(szBuffer[0] == '\0' || szBuffer[0] != '#')
		{
			break;
		}
		
		//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  _S_pneumoniae_Spanis    _3948_7_10
		returned_result = split_string_and_return_specific_index( result, szBuffer, 0);
		if(strcmp(returned_result, "#CHROM")==0)
		{
			for(i = 0; i< number_of_columns; i++)
			{
				returned_result = split_string_and_return_specific_index( result, szBuffer, i);
				strcpy(column_names[i], returned_result);
			}
		}
		
	}while(szBuffer[0] != '\0');
}

// Assume the sample names are unique
int column_number_for_column_name(char ** column_names, char * column_name, int number_of_columns)
{
	int i;
	for(i = 0; i< number_of_columns; i++)
	{
		if(strcmp(column_names[i], column_name) == 0)
		{
			return 	i;
		}
	}
	
	return -1;
}














