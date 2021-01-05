#! /usr/bin/env python3
# encoding: utf-8

"""
Tests the validation of the input fasta file. It must be a multifasta alignment with sensible data
"""

import os
import unittest
from gubbins import common
from gubbins.ValidateFastaAlignment import ValidateFastaAlignment

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestValidateInputFastaFile(unittest.TestCase):

    def test_valid_fasta_file(self):
        validate_fasta = ValidateFastaAlignment('gubbins/tests/data/multiple_recombinations.aln')
        self.assertTrue(validate_fasta.is_input_fasta_file_valid())

    def test_does_each_sequence_have_a_name(self):
        validate_fasta = ValidateFastaAlignment('gubbins/tests/data/sequence_without_a_name.fa')
        self.assertFalse(validate_fasta.does_each_sequence_have_a_name_and_genomic_data())
        self.assertFalse(validate_fasta.is_input_fasta_file_valid())

    def test_does_the_sequence_have_sensible_characters(self):
        validate_fasta = ValidateFastaAlignment('gubbins/tests/data/sequence_with_odd_chars.fa')
        self.assertFalse(validate_fasta.does_each_sequence_have_a_name_and_genomic_data())
        self.assertFalse(validate_fasta.is_input_fasta_file_valid())

    def test_are_all_sequences_the_same_length(self):
        validate_fasta = ValidateFastaAlignment('gubbins/tests/data/valid_alignment.aln')
        self.assertTrue(validate_fasta.does_each_sequence_have_the_same_length())
        self.assertTrue(validate_fasta.is_input_fasta_file_valid())

        validate_fasta = ValidateFastaAlignment('gubbins/tests/data/sequences_of_different_lengths.fa')
        self.assertFalse(validate_fasta.does_each_sequence_have_the_same_length())
        self.assertFalse(validate_fasta.is_input_fasta_file_valid())

    def test_are_all_sequence_names_unique(self):
        validate_fasta = ValidateFastaAlignment('gubbins/tests/data/all_unique_sequence_names.fa')
        self.assertTrue(validate_fasta.are_sequence_names_unique())
        self.assertTrue(validate_fasta.is_input_fasta_file_valid())

        validate_fasta = ValidateFastaAlignment('gubbins/tests/data/non_unique_sequence_names.fa')
        self.assertFalse(validate_fasta.are_sequence_names_unique())
        self.assertFalse(validate_fasta.is_input_fasta_file_valid())

if __name__ == "__main__":
    unittest.main()
