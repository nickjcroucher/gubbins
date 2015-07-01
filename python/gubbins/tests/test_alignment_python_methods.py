#! /usr/bin/env python3
# encoding: utf-8

"""
Tests alignement manipulation and conversion with no external application dependancies.
"""

import unittest
import re
import filecmp
import os
from gubbins import common

class TestAlignmentPythonMethods(unittest.TestCase):

  def test_number_of_sequences_in_alignment(self):
    assert common.GubbinsCommon.number_of_sequences_in_alignment('gubbins/tests/data/small_alignment.aln') == 5

  def test_get_sequence_names_from_alignment(self):
    assert  common.GubbinsCommon.get_sequence_names_from_alignment('gubbins/tests/data/small_alignment.aln') == ['sequence1','sequence2','sequence3','sequence4','sequence5']

  def test_filter_out_alignments_with_too_much_missing_data(self):
    common.GubbinsCommon.filter_out_alignments_with_too_much_missing_data('gubbins/tests/data/alignment_with_too_much_missing_data.aln', 'gubbins/tests/data/alignment_with_too_much_missing_data.aln.actual', 5,0)
    assert filecmp.cmp('gubbins/tests/data/alignment_with_too_much_missing_data.aln.actual','gubbins/tests/data/alignment_with_too_much_missing_data.aln.expected')
    os.remove('gubbins/tests/data/alignment_with_too_much_missing_data.aln.actual')

  def test_reinsert_gaps_into_fasta_file(self):
    common.GubbinsCommon.reinsert_gaps_into_fasta_file('gubbins/tests/data/gaps_to_be_reinserted.aln', 'gubbins/tests/data/gaps_to_be_reinserted.vcf', 'gubbins/tests/data/gaps_to_be_reinserted.aln.actual')
    assert filecmp.cmp('gubbins/tests/data/gaps_to_be_reinserted.aln.actual','gubbins/tests/data/gaps_to_be_reinserted.aln.expected')
    os.remove('gubbins/tests/data/gaps_to_be_reinserted.aln.actual')

  def test_reconvert_fasta_file(self):
    common.GubbinsCommon.reconvert_fasta_file('gubbins/tests/data/alignment_with_too_much_missing_data.aln', 'gubbins/tests/data/reconvert_fasta_file.aln.actual')
    assert filecmp.cmp('gubbins/tests/data/reconvert_fasta_file.aln.actual','gubbins/tests/data/reconvert_fasta_file.aln.expected')
    os.remove('gubbins/tests/data/reconvert_fasta_file.aln.actual')

if __name__ == "__main__":
  unittest.main()

