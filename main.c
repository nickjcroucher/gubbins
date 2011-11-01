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


void detect_snps(char reference_sequence[], FILE * alignment_file_pointer, int length_of_genome)
{
	char * comparison_sequence;
	int i;
	int s;
	
	do{
		s = 0;
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
			}
			if(reference_sequence[i] == '*')
			{
				s++;	
			}
		}
		printf("Snps: %d", s);
	}while(comparison_sequence[0] != '\0');
	printf("Snps: %d", s);
	
	
	free(comparison_sequence);
}


int generate_snp_sites(char filename[])
{
	FILE *alignment_file_pointer;
	int length_of_genome;
	char * reference_sequence;

	alignment_file_pointer=fopen(filename, "r");
	
	if(validate_alignment_file(alignment_file_pointer) == 0)
	{
		return 0;
	}
	
	length_of_genome = genome_length(alignment_file_pointer);
	reference_sequence = (char *) malloc(length_of_genome*sizeof(char));
	
	build_reference_sequence(reference_sequence,alignment_file_pointer);
	detect_snps(reference_sequence, alignment_file_pointer, length_of_genome);
	
	//create_list_of_snps(reference_sequence);
	
	free(reference_sequence);
	
	return 1;
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


