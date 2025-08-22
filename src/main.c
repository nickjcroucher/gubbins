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
#include <getopt.h>
#include <unistd.h>
#include "snp_sites.h"
#include "gubbins.h"
#include "../config.h"
#include "string_cat.h"

#define MAX_FILENAME_SIZE 1024
const char* program_name;

// Assumptions:
// The sequences in the multi fasta alignment file are the same length
// You are only interested in SNPs, INDELS are ignored
// The first sequence is chosen as the reference sequence
// If there is an indel in the reference sequence, the first normal base found in another strain is used.

void print_usage(FILE* stream, int exit_code)
{
  fprintf (stream, "This program is not supposed to be directly run. Use run_gubbins.py instead\n");
  fprintf (stream, "Usage:  %s [options] alignment_file\n", program_name);
  fprintf (stream, "Version: %s\n", PACKAGE_VERSION);
  fprintf (stream,
           "  -r    detect recombinations mode\n"
           "  -t    Newick tree file\n"
           "  -v    VCF file\n"
           "  -f    Original Multifasta file\n"
           "  -n    Number of threads to use for recombination detection\n"
           "  -m    Min SNPs for identifying a recombination block\n"
           "  -a    Min window size\n"
           "  -b    Max window size\n"
           "  -p    p value for detecting recombinations\n"
           "  -i    p value ratio for trimming recombinations\n"
           "  -s    scale branch lengths without correcting for invariant sites\n"
           "  -h    Display this usage information.\n\n"
);
  exit (exit_code);
}

int check_file_exists_or_exit(char * filename)
{
  if (access( filename, F_OK ) != -1 )
  {
    return 1;
  }
  else
  {
    printf("Error: File '%s' does not exist\n",filename);
    print_usage(stderr, EXIT_FAILURE);
    return 0;
  }
}

//int main (argc, argv) int argc; char **argv;
int main (int argc, char ** argv)
{
    int c;
    char multi_fasta_filename[MAX_FILENAME_SIZE] = {""};
    char vcf_filename[MAX_FILENAME_SIZE] = {""};
    char tree_filename[MAX_FILENAME_SIZE] = {""};
    char original_multi_fasta_filename[MAX_FILENAME_SIZE] = {""};

    int recombination_flag = 0 ;
    int min_snps = 3;
    int window_min = 100;
    int window_max = 10000;
    float uncorrected_p_value = 0.05;
    float trimming_ratio = 1.0;
    int extensive_search_flag = 0;
    int scaling_flag = 0;
    int num_threads = 1;
    program_name = argv[0];
  
    while (1)
    {
        static struct option long_options[] =
        {
            {"help",                no_argument,       0, 'h'},
            {"recombination",       no_argument,       0, 'r'},
            {"vcf",                 required_argument, 0, 'v'},
            {"tree",                required_argument, 0, 't'},
            {"original_multifasta", required_argument, 0, 'f'},
            {"min_snps",            required_argument, 0, 'm'},
            {"window_min",          required_argument, 0, 'a'},
            {"window_max",          required_argument, 0, 'b'},
            {"p_value",             required_argument, 0, 'p'},
            {"trimming_ratio",      required_argument, 0, 'i'},
            {"extended_search",     required_argument, 0, 'x'},
            {"ncpu",                required_argument, 0, 'n'},
            {"scaling",             required_argument, 0, 's'},

            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        c = getopt_long (argc, argv, "hrxsv:f:t:m:a:b:p:i:n:",
                           long_options, &option_index);
        /* Detect the end of the options. */
        if (c == -1)
        break;

        switch (c)
        {
            case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
                printf ("option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
                printf ("\n");
                break;
            case 'h':
                print_usage(stdout, EXIT_SUCCESS);
            case 'r':
                recombination_flag = 1;
                break;
            case 'x':
                extensive_search_flag = 1;
                break;
            case 's':
                scaling_flag = 1;
                break;
            case 'f':
                memcpy(original_multi_fasta_filename, optarg, size_of_string(optarg) +1);
                break;
            case 'v':
                memcpy(vcf_filename, optarg, size_of_string(optarg) +1);
                break;
            case 'm':
                min_snps = atoi(optarg);
                break;
            case 'a':
                window_min = atoi(optarg);
                break;
            case 'b':
                window_max = atoi(optarg);
                break;
            case 'p':
                uncorrected_p_value = atof(optarg);
                break;
            case 'i':
                trimming_ratio = atof(optarg);
                break;
            case 't':
                memcpy(tree_filename, optarg, size_of_string(optarg) +1);
                break;
            case 'n':
                num_threads = atoi(optarg);
                break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default:
                abort ();
        }
    }

    /* Print any remaining command line arguments (not options). */
    if (optind < argc)
    {
        memcpy(multi_fasta_filename, argv[optind], size_of_string(argv[optind]) +1);
        optind++;
    }

    check_file_exists_or_exit(multi_fasta_filename);

    if(recombination_flag == 1)
    {
        check_file_exists_or_exit(vcf_filename);
        check_file_exists_or_exit(tree_filename);
        check_file_exists_or_exit(original_multi_fasta_filename);
        run_gubbins(vcf_filename,
                    tree_filename,
                    multi_fasta_filename,
                    min_snps,
                    original_multi_fasta_filename,
                    window_min,
                    window_max,
                    uncorrected_p_value,
                    trimming_ratio,
                    extensive_search_flag,
                    scaling_flag,
                    num_threads);
    }
    else
    {
        generate_snp_sites(multi_fasta_filename, 0, ".gaps");
        generate_snp_sites(multi_fasta_filename, 1, "");
    }

    exit(EXIT_SUCCESS);
}
  


