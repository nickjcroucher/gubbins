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


#ifndef _PARSE_VCF_H_
#define _PARSE_VCF_H_

void get_sequence_from_column_in_vcf(FILE * vcf_file_pointer, char * sequence_bases, int number_of_snps, int column_number);
int get_number_of_snps(FILE * vcf_file_pointer);
int get_number_of_samples(FILE * vcf_file_pointer);
void split_string_and_return_specific_index(char * result, char * input_string, int token_index, int input_string_length);
int get_number_of_columns(char * column_header);
int get_number_of_columns_from_file(FILE * vcf_file_pointer);
void get_column_names(FILE * vcf_file_pointer, char ** column_names, int number_of_columns);
int column_number_for_column_name(char ** column_names, char * column_name, int number_of_columns);

void get_integers_from_column_in_vcf(FILE * vcf_file_pointer, int * integer_values, int number_of_snps, int column_number);


#endif