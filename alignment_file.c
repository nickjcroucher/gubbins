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

// Given a file handle, return the length of the current line
int line_length(FILE * alignment_file_pointer)
{
	char szBuffer[4194304] = {0};  
	char *pcRes         = NULL; 
	int  length_of_line    = 0;    
	
	while((pcRes = fgets(szBuffer, sizeof(szBuffer), alignment_file_pointer))  != NULL){
		length_of_line = strlen(szBuffer) - 1;
		if((szBuffer)[length_of_line] == '\n'){
			break;
		}
	}
	return length_of_line;
}

void advance_to_sequence(FILE * alignment_file_pointer)
{
	// Skip first line since its a comment, ToDo make this better by doing a regex on the line
	line_length(alignment_file_pointer);
}

void advance_to_sequence_name(FILE * alignment_file_pointer)
{
	// Skip sequence line, TODO make this work properly
	line_length(alignment_file_pointer);
}



int validate_alignment_file(FILE * alignment_file_pointer)
{
	return 1;
}




int genome_length(FILE * alignment_file_pointer)
{
	int length_of_genome;
	
	advance_to_sequence(alignment_file_pointer);
	
	length_of_genome = line_length(alignment_file_pointer);
	rewind(alignment_file_pointer);
	return length_of_genome;
}

int read_line(char sequence[], FILE * pFilePtr)
{
    char szBuffer[4194304] = {0};   
    char *pcRes         = NULL;  
    int   lineLength    = 0; 
	
    while((pcRes = fgets(szBuffer, sizeof(szBuffer), pFilePtr))  != NULL){
        //append string to line buffer
		
        strcat(sequence, szBuffer);
        strcpy(szBuffer, "");
        lineLength = strlen(sequence) - 1;
        //if end of line character is found then exit from loop
		
        if((sequence)[lineLength] == '\n'){
            break;
        }
    }
    return 1;
}

int count_lines_in_file(FILE * alignment_file_pointer)
{
	rewind(alignment_file_pointer);
	int i = 0;
	int length_of_line =0;
	
	do{
		length_of_line = line_length(alignment_file_pointer);
		i++;
	}while(length_of_line != 0);
	
	return i;	
}


void get_sample_names_for_header(FILE * alignment_file_pointer, char ** sequence_names, int number_of_samples)
{
	rewind(alignment_file_pointer);
	int i = 0;
	// remove this hardcoding and figure out number of lines in the file
	char * sequence_name;
	char filtered_sequence_name[20];
	int name_counter;
	
	do{
		sequence_name = (char *) malloc(1000*sizeof(char));
		read_line(sequence_name, alignment_file_pointer);
		advance_to_sequence_name(alignment_file_pointer);
		
		if(sequence_name[0] == '\0')
		{
			break;
		}
		
		for(name_counter=0; name_counter < number_of_samples; name_counter++)
		{
			if((sequence_name[name_counter+1] == '\0') || (sequence_name[name_counter+1] == '\n') || (sequence_name[name_counter+1] == '\r') || (name_counter >= 20))
			{
				filtered_sequence_name[name_counter]  = '\0';
				break;
			}
			
			if((sequence_name[name_counter+1] == '\t') || (sequence_name[name_counter+1] == ' ') || sequence_name[name_counter+1] == '>' )
			{
				filtered_sequence_name[name_counter] = '_';
			}
			else
			{
				filtered_sequence_name[name_counter] = sequence_name[name_counter+1];
			}
		}
		//TODO clean up the sample name before use
		strcpy(sequence_names[i], filtered_sequence_name);
		
		i++;
	}while(sequence_name[0] != '\0');
	free(sequence_name);
}



