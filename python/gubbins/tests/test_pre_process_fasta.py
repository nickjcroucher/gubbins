#! /usr/bin/env python3
# encoding: utf-8

"""
Tests for preprocessing the input FASTA file
"""

import unittest
import re
import os
import subprocess
import pprint
from gubbins.PreProcessFasta import PreProcessFasta


class TestPreProcessFasta(unittest.TestCase):
      
  def test_input_file_with_no_duplicate_sequences(self):   
      preprocessfasta = PreProcessFasta('gubbins/tests/data/preprocessfasta/no_duplicates.aln')
      self.assertEqual(preprocessfasta._hash_sequences(), 
       {b"\x840\x89L\xfe\xb5J6%\xf1\x8f\xe2O\xce'.": ['sample1'],
        b'\x9c\xe6\x8b\xf7\xae\xe2\x1f\xf5j\xcfu\xf4\xfdO\x8b\xec': ['sample4'],
        b'\xc3n`\xf5t\x00\x1e\xf3\xde\nU\x1f\x95\x0b\xdb9': ['sample2'],
        b'\xdf\xedM\xdf1\xf2PPO\xc1\xf54T?\xdb\xa2': ['sample3']})

  def test_input_file_with_one_duplicate_sequences(self):   
      preprocessfasta = PreProcessFasta('gubbins/tests/data/preprocessfasta/one_duplicate.aln')
      self.assertEqual(preprocessfasta._hash_sequences(), 
       {b"\x840\x89L\xfe\xb5J6%\xf1\x8f\xe2O\xce'.": ['sample1', 'sample3'],
        b'\x9c\xe6\x8b\xf7\xae\xe2\x1f\xf5j\xcfu\xf4\xfdO\x8b\xec': ['sample4'],
        b'\xc3n`\xf5t\x00\x1e\xf3\xde\nU\x1f\x95\x0b\xdb9': ['sample2']})

  def test_input_file_with_multiple_duplicate_sequences(self):   
      preprocessfasta = PreProcessFasta('gubbins/tests/data/preprocessfasta/multiple_duplicates.aln')
      self.assertEqual(preprocessfasta._hash_sequences(), 
       {b"\x840\x89L\xfe\xb5J6%\xf1\x8f\xe2O\xce'.": ['sample1', 'sample3'],
        b'\x9c\xe6\x8b\xf7\xae\xe2\x1f\xf5j\xcfu\xf4\xfdO\x8b\xec': ['sample2', 'sample4']})

  def test_input_file_with_all_duplicate_sequences(self):   
      preprocessfasta = PreProcessFasta('gubbins/tests/data/preprocessfasta/all_same_sequence.aln')
      self.assertEqual(preprocessfasta._hash_sequences(), 
       {b"\x840\x89L\xfe\xb5J6%\xf1\x8f\xe2O\xce'.": ['sample1',
                                                      'sample2',
                                                      'sample3',
                                                      'sample4']})

if __name__ == "__main__":
  unittest.main()
