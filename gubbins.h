
#ifndef _GUBBINS_H_
#define _GUBBINS_H_

void run_gubbins(char vcf_filename[], char tree_filename[]);
void calculate_ancestor_sequence(char * ancestor_sequence, char * reference_sequence, char ** child_sequences, int sequence_length, int number_of_child_sequences);

#endif

