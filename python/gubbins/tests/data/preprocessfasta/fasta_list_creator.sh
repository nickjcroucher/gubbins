#!/bin/bash 
set -e
# Remove any older ska_fasta_list.txt file

DATA_DIR=$1

if [ -f ska_fasta_list.txt ]
then
	rm ska_fasta_list.txt
fi

ls -d $DATA_DIR/sequence_t* | while read line
do
FASTA_NAME=$(basename $line | sed 's/\..*$//g')
printf "%s\t%s\n" $FASTA_NAME $line >> ska_fasta_list.txt
done

