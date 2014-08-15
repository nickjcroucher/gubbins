#! /usr/bin/env python
# encoding: utf-8

"""
Tests the validation of the starting newick tree.
"""

import unittest
import re
import os
from gubbins import common

class TestValidationOfStartingTree(unittest.TestCase):

  def test_are_all_names_unique(self):
    assert common.GubbinsCommon.is_starting_tree_valid('gubbins/tests/data/duplicate_names_in_tree.tre') == 0

  def test_is_it_a_valid_newick_tree(self):
    assert common.GubbinsCommon.is_starting_tree_valid('gubbins/tests/data/invalid_newick_tree.tre') == 0 
    assert common.GubbinsCommon.is_starting_tree_valid('gubbins/tests/data/valid_newick_tree.tre')  == 1

  def test_do_the_names_match_the_fasta_file(self):
    assert common.GubbinsCommon.do_the_names_match_the_fasta_file('gubbins/tests/data/valid_newick_tree.tre', 'gubbins/tests/data/valid_newick_tree.aln' ) == 1


if __name__ == "__main__":
  unittest.main()

