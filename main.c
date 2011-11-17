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
#include "snp_sites.h"
#include "gubbins.h"

#define MAX_FILENAME_SIZE 200

static void print_usage()
{
	puts("Find SNP sites");
	puts("./multi_fasta_to_vcf -s file.aln ");
	puts("Run Gubbins");
	puts("./multi_fasta_to_vcf -r file.aln file.aln.vcf file.aln.vcf.tre file.aln.phylip");
	puts("");
	puts("Input:");
	puts("\tfile.aln\t:file containting a multi fasta alignment");
	puts("");
	puts("Output:");
	puts("\tfile.aln.vcf\t:VCF file containing all SNP sites");
}


// Assumptions:
// The sequences in the multi fasta alignment file are the same length
// Your only interested in SNPs, INDELS are ignored
// The first sequence is chosen as the reference sequence
// If there is an indel in the reference sequence, the first normal base found in another strain is used.


int main (int argc, const char * argv[]) {
	char multi_fasta_filename[MAX_FILENAME_SIZE];
	char vcf_filename[MAX_FILENAME_SIZE];
	char tree_filename[MAX_FILENAME_SIZE];
	char phylip_filename[MAX_FILENAME_SIZE];

	
	if(strcmp(argv[1], "--help") == 0)
    {
		print_usage();
		return 0;
    }
	else if(strcmp(argv[1], "-s") == 0)
    {
		strcpy(multi_fasta_filename,argv[2]);
		generate_snp_sites(multi_fasta_filename);
		return 0;
    }
	else if(strcmp(argv[1], "-r") == 0)
    {
		strcpy(multi_fasta_filename,argv[2]);
		strcpy(vcf_filename,argv[3]);
		strcpy(tree_filename,argv[4]);
		strcpy(phylip_filename,argv[5]);
		
        run_gubbins(vcf_filename,tree_filename,phylip_filename,multi_fasta_filename);
		return 0;
    }
	
	return 0;
}


