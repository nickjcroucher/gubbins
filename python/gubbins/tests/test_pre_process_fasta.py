# encoding: utf-8

"""
Tests for preprocessing the input FASTA file
"""

import unittest
import re
import os
import subprocess
import filecmp
import pprint
from gubbins.PreProcessFasta import PreProcessFasta

class TestPreProcessFasta(unittest.TestCase):
      
  def test_input_file_with_no_duplicate_sequences(self):   
      preprocessfasta = PreProcessFasta('gubbins/tests/data/preprocessfasta/no_duplicates.aln')
      self.assertEqual(preprocessfasta.hash_sequences(), 
       {b"\x840\x89L\xfe\xb5J6%\xf1\x8f\xe2O\xce'.": ['sample1'],
        b'\x9c\xe6\x8b\xf7\xae\xe2\x1f\xf5j\xcfu\xf4\xfdO\x8b\xec': ['sample4'],
        b'\xc3n`\xf5t\x00\x1e\xf3\xde\nU\x1f\x95\x0b\xdb9': ['sample2'],
        b'\xdf\xedM\xdf1\xf2PPO\xc1\xf54T?\xdb\xa2': ['sample3']})
      self.assertEqual(preprocessfasta.calculate_sequences_missing_data_percentage(), {'sample1': 0.0,
 'sample2': 0.0,
 'sample3': 0.0,
 'sample4': 0.0})
        
      self.assertEqual(preprocessfasta.taxa_of_duplicate_sequences(),[])

      preprocessfasta.remove_duplicate_sequences_and_sequences_missing_too_much_data('output.aln')
      self.assertTrue(filecmp.cmp('output.aln', 'gubbins/tests/data/preprocessfasta/no_duplicates.aln'))
      self.cleanup()

  def test_input_file_with_one_duplicate_sequences(self):   
      preprocessfasta = PreProcessFasta('gubbins/tests/data/preprocessfasta/one_duplicate.aln')
      self.assertEqual(preprocessfasta.hash_sequences(), 
       {b"\x840\x89L\xfe\xb5J6%\xf1\x8f\xe2O\xce'.": ['sample1', 'sample3'],
        b'\x9c\xe6\x8b\xf7\xae\xe2\x1f\xf5j\xcfu\xf4\xfdO\x8b\xec': ['sample4'],
        b'\xc3n`\xf5t\x00\x1e\xf3\xde\nU\x1f\x95\x0b\xdb9': ['sample2']})
        
      self.assertEqual(preprocessfasta.taxa_of_duplicate_sequences(),['sample1'])
      
      preprocessfasta.remove_duplicate_sequences_and_sequences_missing_too_much_data('output.aln')
      self.assertTrue(filecmp.cmp('output.aln', 'gubbins/tests/data/preprocessfasta/expected_one_duplicate.aln'))
      self.cleanup()
 
  def test_input_file_with_multiple_duplicate_sequences(self):   
      preprocessfasta = PreProcessFasta('gubbins/tests/data/preprocessfasta/multiple_duplicates.aln')
      self.assertEqual(preprocessfasta.hash_sequences(), 
       {b"\x840\x89L\xfe\xb5J6%\xf1\x8f\xe2O\xce'.": ['sample1', 'sample3'],
        b'\x9c\xe6\x8b\xf7\xae\xe2\x1f\xf5j\xcfu\xf4\xfdO\x8b\xec': ['sample2', 'sample4']})
        
      self.assertEqual(preprocessfasta.taxa_of_duplicate_sequences(),['sample1','sample2'])
      
      preprocessfasta.remove_duplicate_sequences_and_sequences_missing_too_much_data('output.aln')
      self.assertTrue(filecmp.cmp('output.aln', 'gubbins/tests/data/preprocessfasta/expected_multiple_duplicates.aln'))
      self.cleanup()
 
  def test_input_file_with_all_duplicate_sequences(self):   
      preprocessfasta = PreProcessFasta('gubbins/tests/data/preprocessfasta/all_same_sequence.aln')
      self.assertEqual(preprocessfasta.hash_sequences(), 
       {b"\x840\x89L\xfe\xb5J6%\xf1\x8f\xe2O\xce'.": ['sample1',
                                                      'sample2',
                                                      'sample3',
                                                      'sample4']})
                                                      
      self.assertEqual(preprocessfasta.taxa_of_duplicate_sequences(),['sample1',
                                                      'sample2',
                                                      'sample3'])
      self.cleanup()
                                                      
  def test_filter_out_alignments_with_too_much_missing_data(self):
    preprocessfasta = PreProcessFasta('gubbins/tests/data/preprocessfasta/missing_data.aln', False, 5)
    preprocessfasta.remove_duplicate_sequences_and_sequences_missing_too_much_data('output.aln')
    self.assertTrue(filecmp.cmp('output.aln','gubbins/tests/data/preprocessfasta/expected_missing_data.aln'))
    self.cleanup()           
      
  def cleanup(self):
      for file_to_delete in ['output.aln']:
          if os.path.exists(file_to_delete):
              os.remove(file_to_delete) 

if __name__ == "__main__":
  unittest.main()
