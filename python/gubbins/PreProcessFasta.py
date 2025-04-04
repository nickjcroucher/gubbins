#!/usr/bin/env python
# encoding: utf-8

import sys
import hashlib
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict


class PreProcessFasta(object):

    def __init__(self, input_filename, verbose=False, filter_percentage=25):
        self.input_filename = input_filename
        self.verbose = verbose
        self.filter_percentage = filter_percentage

    def hash_sequences(self):
        sequence_hash_to_taxa = defaultdict(list)
        with open(self.input_filename) as input_handle:
            alignments = AlignIO.parse(input_handle, "fasta")
            for alignment in alignments:
                for record in alignment:
                    sequence_hash = hashlib.md5()
                    sequence_hash.update(str(record.seq).encode('utf-8'))
                    hash_of_sequence = sequence_hash.digest()
                    sequence_hash_to_taxa[hash_of_sequence].append(record.id)

                    if self.verbose:
                        print("Sample " + str(record.id) + " has a hash of " + str(hash_of_sequence))
        input_handle.close()
        return sequence_hash_to_taxa

    def calculate_sequences_missing_data_percentage(self):
        sequences_to_missing_data = {}
        with open(self.input_filename) as input_handle:
            alignments = AlignIO.parse(input_handle, "fasta")
            for alignment in alignments:
                for record in alignment:
                    number_of_gaps = 0
                    number_of_gaps += record.seq.count('n')
                    number_of_gaps += record.seq.count('N')
                    number_of_gaps += record.seq.count('-')
                    sequence_length = len(record.seq)

                    if sequence_length == 0:
                        sequences_to_missing_data[record.id] = 100
                        if self.verbose:
                            print("Sample " + str(record.id) + " has no sequence ")
                    else:
                        per_missing_data = number_of_gaps*100/sequence_length
                        sequences_to_missing_data[record.id] = per_missing_data
                        if self.verbose:
                            print("Sample " + str(record.id) + " has missing data percentage of " +
                                  str(per_missing_data))

        input_handle.close()
        return sequences_to_missing_data

    def taxa_missing_too_much_data(self):
        taxa_to_remove = []
        for taxa, percentage_missing in self.calculate_sequences_missing_data_percentage().items():
            if percentage_missing > self.filter_percentage:
                taxa_to_remove.append(taxa)
                print("Excluded sequence " + taxa + " because it had " + str(percentage_missing) +
                      " percentage missing data while a maximum of " + str(self.filter_percentage) + " is allowed")

        return taxa_to_remove

    def taxa_of_duplicate_sequences(self):
        taxa_to_remove = []
        for sequence_hash, taxa in sorted(self.hash_sequences().items()):
            if len(taxa) > 1:
                taxon_to_keep = taxa.pop()
                for taxon in taxa:
                    print("Sequences in " + taxon + " and " + taxon_to_keep + " are identical, removing " + taxon +
                          " from analysis")
                    taxa_to_remove.append(taxon)

        return taxa_to_remove

    def remove_duplicate_sequences_and_sequences_missing_too_much_data(self, output_filename,
                                                                       remove_identical_sequences=None):

        if not remove_identical_sequences:
            taxa_to_remove = self.taxa_missing_too_much_data()
        else:
            taxa_to_remove = self.taxa_of_duplicate_sequences() + self.taxa_missing_too_much_data()

        with open(self.input_filename) as input_handle:
            alignments = AlignIO.parse(input_handle, "fasta")
            output_alignments = []
            number_of_included_alignments = 0
            for alignment in alignments:
                for record in alignment:
                    if record.id not in taxa_to_remove:
                        output_alignments.append(record)
                        number_of_included_alignments += 1
            if number_of_included_alignments <= 1:
                sys.exit("Not enough sequences are left after removing duplicates.Please check you input data.")

        with open(output_filename, "w") as output_handle:
            AlignIO.write(MultipleSeqAlignment(output_alignments), output_handle, "fasta")

        return taxa_to_remove

    def get_sequence_names(self):
        sequence_names = []
        with open(self.input_filename) as input_handle:
            alignments = AlignIO.parse(input_handle, "fasta")
            for alignment in alignments:
                for record in alignment:
                    sequence_names.append(record.id)
        return sequence_names
    
    def get_alignment_length(self, format = "fasta"):
        alignment = AlignIO.read(self.input_filename, format)
        return alignment.get_alignment_length()

    def get_alignment_information(self, format = "fasta"):
        sequence_names = []
        base_frequencies = [0,0,0,0]
        alignment = AlignIO.read(self.input_filename, format)
        alignment_length = alignment.get_alignment_length()
        
        for record in alignment:
            sequence_names.append(record.id)
            for index,base in enumerate(['A','C','G','T']):
                base_frequencies[index] += (record.seq.count(base) + record.seq.count(base.lower()))
        
        base_frequencies = [b/(alignment_length*len(sequence_names)) for b in base_frequencies]
        
        return alignment_length, sequence_names, base_frequencies
