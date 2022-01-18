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

typedef struct {
    int index;
    char* pattern;
} indexed_pattern;

int qrcmp(const void *x, const void *y) {
    const indexed_pattern pattern_x = *(indexed_pattern *)x;
    const indexed_pattern pattern_y = *(indexed_pattern *)y;
    return (strcmp(pattern_x.pattern,pattern_y.pattern));
}

void create_csv_of_snp_sites(char filename[], int number_of_snps, char ** bases_for_snps) {
    
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
    
    // Indices run consecutively, rather than using SNP locations in whole genome alignment
    // This is because pyjar reconstructs only the polymorphic sites, not the whole sequences
    indexed_pattern* base_pattern_indices = malloc(number_of_snps * sizeof(indexed_pattern));
    int i = 0;
    for (i = 0; i < number_of_snps; i++)
    {
        base_pattern_indices[i].pattern = bases_for_snps[i];
        base_pattern_indices[i].index = i;
    }

    // Sort the base patterns
    qsort(base_pattern_indices, number_of_snps, sizeof(indexed_pattern), qrcmp);
    
    // Avoid any large data structures to reduce peak memory
    int j = 0;
    int8_t first = 1;
    for (j = 0; j < number_of_snps; j++)
    {
        // Identify if value is unique using sorted array
        if (j == 0 || strcmp(base_pattern_indices[j-1].pattern,base_pattern_indices[j].pattern) != 0)
        {
            // Set as first index
            first = 1;
            // End previous line if required
            if (j > 0)
            {
                fprintf(positions_file_pointer, "\n");
            }
            // Print base pattern to file
            fprintf(patterns_file_pointer, "%s\n", base_pattern_indices[j].pattern);
        }
        if (first == 1)
        {
            fprintf(positions_file_pointer, "%i", base_pattern_indices[j].index);
        }
        else
        {
            fprintf(positions_file_pointer, ",%i", base_pattern_indices[j].index);
        }
        first = 0;
    }
    // Final end of line
    fprintf(positions_file_pointer, "\n");
    
    // Tidy up memory
    fclose(patterns_file_pointer);
    fclose(positions_file_pointer);
    free(base_pattern_indices);
    free(patterns_base_filename);
    free(positions_base_filename);
    
}
