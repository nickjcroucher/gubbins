#! /usr/bin/env python
# encoding: utf-8

"""
Tests if the same recombinations are see in multiple iterations
"""

import unittest
import re
import os
from gubbins import common

class TestConvergingRecombinations(unittest.TestCase):

  def test_reading_embl_file(self):
    assert common.GubbinsCommon.extract_recombinations_from_embl('gubbins/tests/data/small_recombination.embl') == { 'sequence1': [[1,3],[5,7],[9,11]], 'sequence2': [[1,3],[9,11]], 'sequence3': [[1,3]] }
  

  def test_two_files_are_different(self):
    assert common.GubbinsCommon.have_recombinations_been_seen_before('gubbins/tests/data/small_recombination.embl', ['gubbins/tests/data/small_recombination_different.embl']) == 0
    
  def test_two_files_are_same(self):
      assert common.GubbinsCommon.have_recombinations_been_seen_before('gubbins/tests/data/small_recombination.embl', ['gubbins/tests/data/small_recombination.embl']) == 1
  
  def test_multiple_files_have_one_match(self):
          assert common.GubbinsCommon.have_recombinations_been_seen_before('gubbins/tests/data/small_recombination.embl', ['gubbins/tests/data/small_recombination_different.embl','gubbins/tests/data/small_recombination.embl','gubbins/tests/data/small_recombination_different.embl']) == 1

if __name__ == "__main__":
  unittest.main()

