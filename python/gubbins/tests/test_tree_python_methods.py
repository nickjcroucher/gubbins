#! /usr/bin/env python
# encoding: utf-8

"""
Test manipulating and interogating trees using python modules with no external application dependancies.
"""

import unittest
import re
import shutil
import os
import difflib
from gubbins import common

class TestTreePythonMethods(unittest.TestCase):

  def test_robinson_foulds_distance(self):
    # two tree with different distances
    assert common.GubbinsCommon.robinson_foulds_distance('gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree2.tre') == 17.263494
    # two trees with same distance
    assert common.GubbinsCommon.robinson_foulds_distance('gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree1.tre') == 0
  
  def test_has_tree_been_seen_before(self):
    # same content so the tree has been seen before
    assert common.GubbinsCommon.has_tree_been_seen_before(['gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree1_dup.tre']) == 1
    # different trees
    assert common.GubbinsCommon.has_tree_been_seen_before(['gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree2.tre']) == 0

    #def test_reroot_tree(self):
    #  assert 1 == 0
    #
  def test_reroot_tree_with_outgroup(self):
    shutil.copyfile('gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree1.tre.copy')
    common.GubbinsCommon.reroot_tree_with_outgroup('gubbins/tests/data/robinson_foulds_distance_tree1.tre.copy', 'sequence_4')
    actual_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.copy', 'U').readlines()
    expected_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4', 'U').readlines()
    assert actual_file_content == expected_file_content
    os.remove('gubbins/tests/data/robinson_foulds_distance_tree1.tre.copy')
    
    
    #
    #def test_split_all_non_bi_nodes(self):
    #  assert 1 == 0
    #
    #def test_split_child_nodes(self):
    #  assert 1 == 0
    #
    #def test_reroot_tree_at_midpoint(self):
    #  assert 1 == 0
    #
    #def test_filter_out_removed_taxa_from_tree_and_return_new_file(self):
    #  assert 1 == 0
    #def test_create_pairwise_newick_tree(self):
    #  assert 1 == 0
    #
    
    
if __name__ == "__main__":
  unittest.main()