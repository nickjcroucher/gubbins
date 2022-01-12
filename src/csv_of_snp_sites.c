/*
 *  Imperial College London
 *  Copyright (C) 2022  Imperial College London
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
#include "csv_of_snp_sites.h"
#include "parse_phylip.h"
#include "string_cat.h"

void create_csv_of_snp_sites(char filename[], int number_of_snps, char ** bases_for_snps, int* snp_location, char ** sequence_names, int number_of_samples,int internal_nodes[]) {
    
    // Patterns CSV file
    FILE* patterns_file_pointer;
    char* patterns_base_filename;
    patterns_base_filename = (char*) calloc(1024,sizeof(char));
    memcpy(patterns_base_filename, filename, 1024*sizeof(char));
    char patterns_extension[19] = {".base_patterns.csv"};
    concat_strings_created_with_malloc(patterns_base_filename,patterns_extension);
    patterns_file_pointer = fopen(patterns_base_filename, "w");
    
    // Positions CSV file
    FILE* positions_file_pointer;
    char* positions_base_filename;
    positions_base_filename = (char*) calloc(1024,sizeof(char));
    memcpy(positions_base_filename, filename, 1024*sizeof(char));
    char positions_extension[20] = {".base_positions.csv"};
    concat_strings_created_with_malloc(positions_base_filename,positions_extension);
    positions_file_pointer = fopen(positions_base_filename, "w");

    // Avoid any large data structures to reduce peak memory
    int i = 0;
    for (i = 0; i < number_of_snps; i++)
    {
        int j = 0;
        for (j = 0; j < number_of_snps; j++)
        {
            if (strcmp(bases_for_snps[i],bases_for_snps[j]) == 0) {
                // If the pattern has not been observed earlier in the list,
                // print to file
                if (j == i)
                {
                    fprintf(patterns_file_pointer, "%s\n", bases_for_snps[i]);
                    fprintf(positions_file_pointer, "%i", snp_location[i]);
                }
                // Stop iterating if base pattern has already been seen
                // earlier in the list
                if (j < i)
                {
                    break;
                }
                // Print out other locations at which the same string
                // is observed on the first observation of the string
                else if (j > i)
                {
                    fprintf(positions_file_pointer, ",%i", snp_location[j]);
                }
            }
        }
        if (j >= i)
        {
            fprintf(positions_file_pointer, "\n");
        }
    }
    
    // Tidy up memory
    fclose(patterns_file_pointer);
    fclose(positions_file_pointer);
    free(patterns_base_filename);
    free(positions_base_filename);
    
}
