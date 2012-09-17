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

#define MAX_FILENAME_SIZE 1024
const char* program_name;

// Assumptions:
// The sequences in the multi fasta alignment file are the same length
// Your only interested in SNPs, INDELS are ignored
// The first sequence is chosen as the reference sequence
// If there is an indel in the reference sequence, the first normal base found in another strain is used.

int check_file_exists_or_exit(char * filename)
{
  if( access( filename, F_OK ) != -1 ) {
		return 1;
  } else {
		printf("Error: File '%s' doesnt exist\n",filename);
		print_usage(stderr, EXIT_FAILURE);
  }
}

void print_usage(FILE* stream, int exit_code)
{
  fprintf (stream, "Usage:  %s [options] alignment_file\n", program_name);
  fprintf (stream, "Version: %s\n", PACKAGE_VERSION);
  fprintf (stream,
           "  -r    detect recombinations mode\n"
           "  -t    Newick tree file\n"
           "  -p    Phylip file\n"
           "  -v    VCF file\n"
           "  -m    Min SNPs for identifying a recombination block\n"
           "  -h    Display this usage information.\n\n"
);

  fprintf (stream, "Step 1: Detect SNP sites (generates inputs files for step 2)\n");
  fprintf (stream, "gubbins alignment_file\n\n", program_name);
  fprintf (stream, "Step 2: Detect recombinations\n");
  fprintf (stream, "gubbins -r -v vcf_file -t newick_tree -p phylip_file -m 10 alignment_file\n\n", program_name);
  exit (exit_code);
}

int main (argc, argv) int argc; char **argv;
{
  int c;
  char multi_fasta_filename[MAX_FILENAME_SIZE];
  char vcf_filename[MAX_FILENAME_SIZE];
  char tree_filename[MAX_FILENAME_SIZE];
  char phylip_filename[MAX_FILENAME_SIZE];
  int recombination_flag = 0 ;
	int min_snps = 3;
  program_name = argv[0];
  
  while (1)
    {
      static struct option long_options[] =
        {
					{"help",          no_argument,       0, 'h'},
          {"recombination", no_argument,       0, 'r'},
          {"vcf",           required_argument, 0, 'v'},
          {"tree",          required_argument, 0, 't'},
          {"min_snps",      required_argument, 0, 'm'},
          {"phylip",        required_argument, 0, 'p'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;
      c = getopt_long (argc, argv, "hrv:t:p:m:",
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
        case 'v':
          strcpy(vcf_filename,optarg);
          break;
	      case 'm':
	        min_snps = atoi(optarg);
	        break;
        case 't':
          strcpy(tree_filename,optarg);
          break;
        case 'p':
          strcpy(phylip_filename,optarg);
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
        strcpy(multi_fasta_filename,argv[optind++]);
    }

	
		check_file_exists_or_exit(multi_fasta_filename);
  
    if(recombination_flag == 1)
    {
			check_file_exists_or_exit(vcf_filename);
			check_file_exists_or_exit(tree_filename);
			check_file_exists_or_exit(phylip_filename);
      run_gubbins(vcf_filename,tree_filename,phylip_filename,multi_fasta_filename, min_snps);
    }
    else
    {
      generate_snp_sites(multi_fasta_filename, 0, ".gaps");
			generate_snp_sites(multi_fasta_filename, 1, "");
    }

    exit(EXIT_SUCCESS);
}
  


