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

struct cmpargs {
    int *arr;
    char **words;
} args;

// String comparison code
int qcmp(const void *x, const void *y) {
    
    const char *pp1 = *(const char**)x;
    const char *pp2 = *(const char**)y;
    printf("PP1 is %s\n",pp1);
    int compval = strcmp(pp1,pp2);
//    return strncmp(pp1, pp2, 14);
    return compval;
    
//    const char *(x_str) = *(const char **)x;
//    const char *(y_str) = *(const char **)y;
//    return strcmp(x_str, y_str);
    
//    const int ix = *(const int *)x;
//    const int iy = *(const int *)y;
//    return ix - iy;
//    const char *(x_str) = *(const char **)bases[*ix];
//    const char *(y_str) = *(const char **)bases[*iy];
//    return strcmp(x_str, y_str);
//    const char *(x_str) = *(const char **)bases[*(const int*)x];
//    const char *(y_str) = *(const char **)bases[*(const int*)y];
//    return strcmp(x_str, y_str);
}

//int cmpfunc(const void * a, const void * b, void *_args) {
//    struct cmpargs args = _args;
//
//    int *a1 = a;
//    int *a2 = b;
//
//    int idx1 = a - args->arr;
//    int idx2 = b - args->arr;
//}

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

    // Copy pointers for sorting
    char** sorted_base_patterns = malloc(number_of_snps * sizeof(char*));
    int i = 0;
    for (i = 0; i < number_of_snps; i++)
    {
        sorted_base_patterns[i] = bases_for_snps[i];
    }
    
    // Create structure for sorting
    args.arr = sorted_base_patterns;
    args.words = bases_for_snps;

    // Sort the pointers to strings
    int stringLen = sizeof(bases_for_snps) / sizeof(char*);
    
    qsort(sorted_base_patterns, number_of_snps, sizeof(char*), qcmp);
//    qsort_r(sorted_base_patterns, number_of_snps, sizeof(int*), bases_for_snps, qcmp);
//    qsort_r(sorted_base_patterns, number_of_snps, sizeof(int), &args, qcmp);
    
    // Avoid any large data structures to reduce peak memory
    int j = 0;
    for (j = 0; j < number_of_snps; j++)
    {
        // Identify if value is unique using sorted array
        if (j == 0 || strcmp(sorted_base_patterns[j-1],sorted_base_patterns[j]) != 0)
        {
            // Print base pattern to file
            fprintf(patterns_file_pointer, "%s\n", sorted_base_patterns[j]);
            // Print each position at which the pattern is observed
            int k = 0;
            for (k = 0; k < number_of_snps; k++)
            {
                if (strcmp(bases_for_snps[j],bases_for_snps[k]) == 0)
                {
                    fprintf(positions_file_pointer, "%i\t", snp_location[k]);
                }
            }
            fprintf(positions_file_pointer, "\n");
        }
    }

    // Tidy up memory
    fclose(patterns_file_pointer);
    fclose(positions_file_pointer);
    free(sorted_base_patterns);
    free(patterns_base_filename);
    free(positions_base_filename);
    
}
