#! /usr/bin/env python3
# encoding: utf-8

"""
Tests alignment manipulation and conversion with no external application dependencies.
"""

import unittest
import filecmp
import os
from gubbins import common

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAlignmentMethods(unittest.TestCase):

    def test_number_of_sequences_in_alignment(self):
        assert common.number_of_sequences_in_alignment(os.path.join(data_dir, 'small_alignment.aln')) == 5

    def test_get_sequence_names_from_alignment(self):
        assert common.get_sequence_names_from_alignment(os.path.join(data_dir, 'small_alignment.aln')) == \
               ['sequence1', 'sequence2', 'sequence3', 'sequence4', 'sequence5']

    def test_reinsert_gaps_into_fasta_file(self):
        common.reinsert_gaps_into_fasta_file(os.path.join(data_dir, 'gaps_to_be_reinserted.aln'),
                                             os.path.join(data_dir, 'gaps_to_be_reinserted.vcf'),
                                             os.path.join(data_dir, 'gaps_to_be_reinserted.aln.actual'))
        assert filecmp.cmp(os.path.join(data_dir, 'gaps_to_be_reinserted.aln.actual'),
                           os.path.join(data_dir, 'gaps_to_be_reinserted.aln.expected'))
        os.remove(os.path.join(data_dir, 'gaps_to_be_reinserted.aln.actual'))

    def test_reconvert_fasta_file(self):
        common.reconvert_fasta_file(os.path.join(data_dir, 'alignment_with_too_much_missing_data.aln'),
                                    os.path.join(data_dir, 'reconvert_fasta_file.aln.actual'))
        assert filecmp.cmp(os.path.join(data_dir, 'reconvert_fasta_file.aln.actual'),
                           os.path.join(data_dir, 'reconvert_fasta_file.aln.expected'))
        os.remove(os.path.join(data_dir, 'reconvert_fasta_file.aln.actual'))

    def test_concatenate_fasta_files(self):
        common.concatenate_fasta_files([os.path.join(data_dir, 'small_alignment.aln'),
                                        os.path.join(data_dir, 'further_alignment.aln')],
                                       os.path.join(data_dir, 'concatenate_fasta_files.aln.actual'))
        assert filecmp.cmp(os.path.join(data_dir, 'concatenate_fasta_files.aln.actual'),
                           os.path.join(data_dir, 'concatenate_fasta_files.aln.expected'))
        os.remove(os.path.join(data_dir, 'concatenate_fasta_files.aln.actual'))


if __name__ == "__main__":
    unittest.main()
