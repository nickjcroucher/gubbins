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
import tempfile
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

  def test_reroot_tree(self):
    shutil.copyfile('gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_actual')
    common.GubbinsCommon.reroot_tree('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_actual', 'sequence_4')
    actual_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_actual', 'U').readlines()
    expected_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4', 'U').readlines()
    assert actual_file_content == expected_file_content
    os.remove('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_actual')
    
    shutil.copyfile('gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual')
    common.GubbinsCommon.reroot_tree('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual','')
    actual_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual', 'U').readlines()
    expected_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_expected', 'U').readlines()
    assert actual_file_content == expected_file_content
    os.remove('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual')
    
  def test_reroot_tree_with_outgroup(self):
    shutil.copyfile('gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_actual')
    common.GubbinsCommon.reroot_tree_with_outgroup('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_actual', 'sequence_4')
    actual_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_actual', 'U').readlines()
    expected_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4', 'U').readlines()
    assert actual_file_content == expected_file_content
    os.remove('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_actual')
    
  def test_split_all_non_bi_nodes(self):
    # best way to access it is via reroot_tree_at_midpoint because it outputs to a file
    shutil.copyfile('gubbins/tests/data/non_bi_tree.tre','gubbins/tests/data/non_bi_tree.tre.actual')
    common.GubbinsCommon.reroot_tree_at_midpoint('gubbins/tests/data/non_bi_tree.tre.actual')
    actual_file_content = open('gubbins/tests/data/non_bi_tree.tre.actual', 'U').readlines()
    expected_file_content = open('gubbins/tests/data/non_bi_tree.tre.expected', 'U').readlines()
    assert actual_file_content == expected_file_content
    os.remove('gubbins/tests/data/non_bi_tree.tre.actual')
    
  def test_reroot_tree_at_midpoint(self):
    shutil.copyfile('gubbins/tests/data/robinson_foulds_distance_tree1.tre','gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual')
    common.GubbinsCommon.reroot_tree_at_midpoint('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual')
    actual_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual', 'U').readlines()
    expected_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_expected', 'U').readlines()
    assert actual_file_content == expected_file_content
    os.remove('gubbins/tests/data/robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual')

  def test_filter_out_removed_taxa_from_tree_and_return_new_file(self):
    temp_working_dir = tempfile.mkdtemp(dir=os.getcwd())
    common.GubbinsCommon.filter_out_removed_taxa_from_tree_and_return_new_file('gubbins/tests/data/robinson_foulds_distance_tree1.tre', temp_working_dir, ['sequence_1','sequence_2','sequence_3','sequence_4','sequence_5'])    
    actual_file_content = open(temp_working_dir + '/robinson_foulds_distance_tree1.tre', 'U').readlines()
    expected_file_content = open('gubbins/tests/data/robinson_foulds_distance_tree1.tre.filter_out_removed_taxa_from_tree_expected', 'U').readlines()
    assert actual_file_content == expected_file_content
    os.remove(temp_working_dir + '/robinson_foulds_distance_tree1.tre')
    os.removedirs(temp_working_dir)
    
  def test_create_pairwise_newick_tree(self):
    common.GubbinsCommon.create_pairwise_newick_tree(['sequence_2','sequence_3'], 'gubbins/tests/data/pairwise_newick_tree.actual')
    actual_file_content = open('gubbins/tests/data/pairwise_newick_tree.actual', 'U').readlines()
    expected_file_content = open('gubbins/tests/data/pairwise_newick_tree.expected', 'U').readlines()
    assert actual_file_content == expected_file_content
    os.remove('gubbins/tests/data/pairwise_newick_tree.actual')
    
if __name__ == "__main__":
  unittest.main()