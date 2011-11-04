#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "snp_sites.h"


static void print_usage()
{
	puts("Usage:");
	puts("./multi_fasta_to_vcf file.aln ");
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
	if(strcmp(argv[1], "--help") == 0)
    {
		print_usage();
		return 0;
    }
	
	generate_snp_sites(argv[1]);
	
	
	return 0;
}


