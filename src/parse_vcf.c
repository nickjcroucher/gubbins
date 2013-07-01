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
#include <sys/types.h>
#include <regex.h>
#include "vcf.h"
#include "parse_vcf.h"
#include "alignment_file.h"

int * column_data;

void get_integers_from_column_in_vcf(FILE * vcf_file_pointer, int * integer_values, int number_of_snps, int column_number)
{
	rewind(vcf_file_pointer);
	char * szBuffer;
	szBuffer = (char *) calloc(MAX_READ_BUFFER,sizeof(char));
	int reference_index = 0;
	char result[1000] = {0};  
	
	do{
		szBuffer[0] = '\0';
 		// check the first character of the line to see if its in the header
  	szBuffer = read_line(szBuffer, vcf_file_pointer);
		
		if(szBuffer[0] == '\0')
		{
			break;
		}
		
		if(szBuffer[0] != '#')
		{
			split_string_and_return_specific_index(result, szBuffer, column_number,100000);
			integer_values[reference_index] = atoi(result);
			reference_index++;
		}
		
	}while(szBuffer[0] != '\0');
	free(szBuffer);
}


void get_sequence_from_column_in_vcf(FILE * vcf_file_pointer, char * sequence_bases, int number_of_snps, int column_number)
{
	rewind(vcf_file_pointer);
	char * szBuffer;
	szBuffer = (char *) calloc(MAX_READ_BUFFER,sizeof(char));
	int reference_index = 0;
	char result[1000] = {0};  
		
	do{
		szBuffer[0] = '\0';
		// check the first character of the line to see if its in the header
		szBuffer = read_line(szBuffer, vcf_file_pointer);
		
		if(szBuffer[0] == '\0')
		{
			break;
		}
		
		if(szBuffer[0] != '#')
		{
			split_string_and_return_specific_index(result, szBuffer, column_number, 1000);
			sequence_bases[reference_index] = result[0];
			reference_index++;
		}
		
	}while(szBuffer[0] != '\0');
	
	sequence_bases[reference_index] = '\0';
}

void split_string_and_return_specific_index(char * result, char * input_string, int token_index, int input_string_length)
{
	int i;
	int tab_counter = 0;
	int result_counter = 0; 
	result[0] = '\0';
	
	for(i = 0 ; i< input_string_length; i++)
	{
		if(input_string[i] == '\0' || input_string[i] == '\n')
		{
			result[result_counter] = '\0';
			break;	
		}
		
		if(input_string[i] == '\t')
		{
			tab_counter++;
		}
		else if(tab_counter == token_index)
		{
			result[result_counter] = input_string[i];
			result_counter++;
		}
		else if(tab_counter > token_index)
		{
			result[result_counter] = '\0';
			break;	
		}
	}
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
	
	do
	{
		split_string_and_return_specific_index( result, column_header, number_of_columns,100000);
		if(result == NULL ||  result[0] == '\n')
		{
			break;	
		}
		
		number_of_columns++;
	}
	while( ! (strcmp(result,"") == 0) );
	
	return number_of_columns;
}

// Assumes that all column headers have something in them
int get_number_of_columns_from_file(FILE * vcf_file_pointer)
{
	rewind(vcf_file_pointer);
	char result[100] = {0};
	
	char * szBuffer;
	szBuffer = (char *) calloc(MAX_READ_BUFFER,sizeof(char));
	
	do{
		szBuffer[0] = '\0';
		// check the first character of the line to see if its in the header
		szBuffer = read_line(szBuffer, vcf_file_pointer);
		if(szBuffer[0] == '\0' || szBuffer[0] != '#')
		{
			break;
		}
		
		//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  _S_pneumoniae_Spanis    _3948_7_10
		split_string_and_return_specific_index( result, szBuffer, 0,100000);
		if(strcmp(result, "#CHROM")==0)
		{
			int number_of_columns =  get_number_of_columns(szBuffer);
			free(szBuffer);
			return number_of_columns;
		}
		
	}while(szBuffer[0] != '\0');
	free(szBuffer);
	return 0;
}


void get_column_names(FILE * vcf_file_pointer, char ** column_names, int number_of_columns)
{
	rewind(vcf_file_pointer);
	char * szBuffer;
	szBuffer = (char *) calloc(MAX_READ_BUFFER,sizeof(char));
	char result[100] = {0};  
	int i;
	
	do{
		szBuffer[0] = '\0';
		// check the first character of the line to see if its in the header
		szBuffer = read_line(szBuffer, vcf_file_pointer);
		
		if(szBuffer[0] == '\0' || szBuffer[0] != '#')
		{
			break;
		}
		
		//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  _S_pneumoniae_Spanis    _3948_7_10
		split_string_and_return_specific_index( result, szBuffer, 0,100000);
		if(strcmp(result, "#CHROM")==0)
		{
			for(i = 0; i< number_of_columns; i++)
			{
				split_string_and_return_specific_index( result, szBuffer, i,100000);
				memcpy(column_names[i], result, (strlen(result)+1)*sizeof(char));
			}
		}
		
	}while(szBuffer[0] != '\0');
	free(szBuffer);
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

