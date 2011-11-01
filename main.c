#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>

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


int build_reference_sequence(char reference_sequence[], FILE * alignment_file_pointer)
{
	int i;
	
	rewind(alignment_file_pointer);
    advance_to_sequence(alignment_file_pointer);
	
    read_line(reference_sequence, alignment_file_pointer);
	
	for(i = 0; reference_sequence[ i ]; i++)
	{
		reference_sequence[i] = toupper(reference_sequence[i]);
	}
	
	return 1;
}


int detect_snps(char reference_sequence[], FILE * alignment_file_pointer, int length_of_genome)
{
	char * comparison_sequence;
	int i;
	int number_of_snps = 0;
	
	do{
		comparison_sequence = (char *) malloc(length_of_genome*sizeof(char));
		advance_to_sequence(alignment_file_pointer);
		read_line(comparison_sequence, alignment_file_pointer);
		
		if(comparison_sequence[0] == '\0')
		{
			break;
		}
	
		// Set the reference base to * if 
		for(i = 0; reference_sequence[ i ]; i++)
		{
			if(reference_sequence[i] != '*' && comparison_sequence[i] != '-' && reference_sequence[i] != toupper(comparison_sequence[i]))
			{
				reference_sequence[i] = '*';
				number_of_snps++;
			}
		}
	}while(comparison_sequence[0] != '\0');

	free(comparison_sequence);
	return number_of_snps;
}

void build_snp_locations(int snp_locations[], char reference_sequence[])
{
	int i;
	int snp_counter = 0;
	
	for(i = 0; reference_sequence[i]; i++)
    {
		if(reference_sequence[i] == '*')
		{
			snp_locations[snp_counter] = i;
			snp_counter++;
		}
	}
}

void output_vcf_header( FILE * vcf_file_pointer)
{
	fprintf( vcf_file_pointer, "##fileformat=VCFv4.1\n" );	
	fprintf( vcf_file_pointer, "##INFO=<ID=AB,Number=1,Type=String,Description=\"Alt Base\">\n" );
	fprintf( vcf_file_pointer, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\n" );
}


void get_bases_for_each_snp(FILE * alignment_file_pointer, int snp_locations[], char bases_for_snps[][], int length_of_genome)
{
	int i;
	int sequence_number = 0;
	
	do{
		advance_to_sequence(alignment_file_pointer);
		read_line(comparison_sequence, alignment_file_pointer);

		*****************
		if(comparison_sequence[0] == '\0')
		{
			break;
		}
		
		for(i = 0; snp_locations[i]; i++)
		{
		  	bases_for_snps[i][sequence_number] =  comparison_sequence[i];
		}
		
		sequence_number++;
	}while(comparison_sequence[0] != '\0');
	
	free(comparison_sequence);
}


void create_vcf_file(char filename[],  FILE * alignment_file_pointer, int snp_locations[], int length_of_genome)
{
	FILE *vcf_file_pointer;
	int number_of_samples;
	char * bases_for_snps;
	
	// TODO chunk up to reduce memory usage
	
	vcf_file_pointer=fopen(strcat(filename,".vcf"), "w");
	output_vcf_header(vcf_file_pointer);
	
	// store values for each snp location
	rewind(alignment_file_pointer);
	
	number_of_samples = count_lines_in_file(alignment_file_pointer)/2;
	bases_for_snps =(char *) malloc(number_of_samples*(sizeof(snp_locations)/sizeof(*snp_locations))*sizeof(char));
	
	get_bases_for_each_snp(alignment_file_pointer, snp_locations, bases_for_snps,length_of_genome);
	
	
}


int generate_snp_sites(char filename[])
{
	FILE *alignment_file_pointer;
	int length_of_genome;
	char * reference_sequence;
	int number_of_snps;
	int * snp_locations;

	alignment_file_pointer=fopen(filename, "r");
	
	if(validate_alignment_file(alignment_file_pointer) == 0)
	{
		return 0;
	}
	
	length_of_genome = genome_length(alignment_file_pointer);
	reference_sequence = (char *) malloc(length_of_genome*sizeof(char));
	
	build_reference_sequence(reference_sequence,alignment_file_pointer);
	number_of_snps = detect_snps(reference_sequence, alignment_file_pointer, length_of_genome);
	
	snp_locations = (int *) malloc(number_of_snps*sizeof(int));
	build_snp_locations(snp_locations, reference_sequence);
	free(reference_sequence);
	
	//printf("Number of SNPs: %d\n", number_of_snps);
	create_vcf_file(filename, alignment_file_pointer, snp_locations,length_of_genome);
	
	free(snp_locations);
	return 1;
}



void get_sample_names_for_header(FILE * alignment_file_pointer, char sequence_names[])
{
	rewind(alignment_file_pointer);
	int i = 0;
	// remove this hardcoding and figure out number of lines in the file
	char * sequence_name;
	
	do{
		sequence_name = (char *) malloc(500*sizeof(char));
		read_line(sequence_name, alignment_file_pointer);
		advance_to_sequence_name(alignment_file_pointer);
		
		if(sequence_name[0] == '\0')
		{
			break;
		}
		
		//TODO clean up the sample name before use
		strcpy(sequence_names[i],sequence_name);
		i++;
	}while(sequence_name[0] != '\0');
	free(sequence_name);
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


