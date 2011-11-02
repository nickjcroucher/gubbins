#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "vcf.h"
#include "alignment_file.h"
#include "snp_sites.h"


static void print_usage()
{
	puts("Usage:");
	puts("./multi_fasta_to_vcf file.aln ");
	puts("");
	puts("\tfile.aln\t\t\t: file containting a multi fasta alignment");
	puts("");
	puts("Output:");
	puts("\tfile.aln.vcf\t\t:VCF file containing all SNP sites");
}


// first pass read in first genome and store in array
// read in subsequent lines, if base difference mark the coord
// no need to check already marked base
// could speed up by comparing x strings at a time
// or by having very large input read buffer


int main (int argc, const char * argv[]) {
	if(strcmp(argv[1], "--help") == 0)
    {
		print_usage();
		return 0;
    }
	
	generate_snp_sites(argv[1]);
	
	
	return 0;
}


