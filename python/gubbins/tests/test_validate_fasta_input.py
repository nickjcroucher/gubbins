#! /usr/bin/env python
# encoding: utf-8

"""
Tests the validation of the input fasta file. It must be a multifasta alignment with sensible data
"""

import unittest
import re
import os
from nose.tools import *
from gubbins import common

class TestValidateInputFastaFile(unittest.TestCase):

  def test_input_file_exists(self):
    assert common.GubbinsCommon.does_file_exist('non_existant_file', 'Alignment File') == 0
    
  def test_does_each_sequence_have_a_name(self):
    assert common.GubbinsCommon.does_each_sequence_have_a_name_and_genomic_data('gubbins/tests/data/sequence_without_a_name.fa') == 0

  def test_does_the_sequence_have_sensible_characters(self):
    assert common.GubbinsCommon.does_each_sequence_have_a_name_and_genomic_data('gubbins/tests/data/multiple_recombinations.aln') == 1
    assert common.GubbinsCommon.does_each_sequence_have_a_name_and_genomic_data('gubbins/tests/data/sequence_with_odd_chars.fa') == 0
     
  def test_are_all_sequences_the_same_length(self):
    assert common.GubbinsCommon.does_each_sequence_have_the_same_length('gubbins/tests/data/valid_alignment.aln') == 1
    assert common.GubbinsCommon.does_each_sequence_have_the_same_length('gubbins/tests/data/sequences_of_different_lengths.fa') == 0
  
  def test_are_all_sequence_names_unique(self):
    assert common.GubbinsCommon.are_sequence_names_unique('gubbins/tests/data/all_unique_sequence_names.fa') == 1
    assert common.GubbinsCommon.are_sequence_names_unique('gubbins/tests/data/non_unique_sequence_names.fa') == 0

if __name__ == "__main__":
  unittest.main()

