#ifndef _ALIGNMENT_FILE_H_
#define _ALIGNMENT_FILE_H_

int line_length(FILE * alignment_file_pointer);
void advance_to_sequence(FILE * alignment_file_pointer);
void advance_to_sequence_name(FILE * alignment_file_pointer);
int validate_alignment_file(FILE * alignment_file_pointer);
int genome_length(FILE * alignment_file_pointer);
int read_line(char sequence[], FILE * pFilePtr);
int count_lines_in_file(FILE * alignment_file_pointer);
void get_sample_names_for_header(FILE * alignment_file_pointer, char ** sequence_names);

#endif