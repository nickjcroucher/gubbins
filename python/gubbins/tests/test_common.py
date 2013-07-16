#! /usr/bin/env python
# encoding: utf-8

"""
Tests common.
"""

import unittest
import re
from gubbins import common

class TestCommon(unittest.TestCase):

  def test_which(self):
    # the location of ls varies depending on OS so just check end
    assert re.match('.*/ls$', common.GubbinsCommon.which('ls')) != None
    # Strip parameters
    assert re.match('.*/ls$', common.GubbinsCommon.which('ls -alrt')) != None
    assert common.GubbinsCommon.which('non_existant_program') == None
    
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
  
  #def test_parse_and_run(self):
  #  assert 1 == 0
  #def test_reroot_tree(self):
  #  assert 1 == 0
  #  
  #def test_reroot_tree_with_outgroup(self):
  #  assert 1 == 0
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
  #def test_raxml_base_name(self):
  #  assert 1 == 0
  #  
  #def test_raxml_current_tree_name(self):
  #  assert 1 == 0
  #  
  #def test_raxml_previous_tree_name(self):
  #  assert 1 == 0
  #  
  #def test_raxml_previous_tree(self):
  #  assert 1 == 0
  #  
  #def test_raxml_tree_building_command(self):
  #  assert 1 == 0
  #  
  #def test_raxml_gubbins_command(self):
  #  assert 1 == 0
  #  
  #def test_raxml_regex_for_file_deletions(self):
  #  assert 1 == 0
  #  
  #def test_fasttree_regex_for_file_deletions(self):
  #  assert 1 == 0
  #  
  #def test_starting_files_regex(self):
  #  assert 1 == 0
  #  
  #def test_fasttree_current_tree_name(self):
  #  assert 1 == 0
  #  
  #def test_fasttree_previous_tree_name(self):
  #  assert 1 == 0
  #  
  #def test_fasttree_tree_building_command(self):
  #  assert 1 == 0
  #  
  #def test_fasttree_gubbins_command(self):
  #  assert 1 == 0
  #  
  #def test_fasttree_fastml_command(self):
  #  assert 1 == 0
  #  
  #def test_raxml_fastml_command(self):
  #  assert 1 == 0
  #  
  #def test_generate_fastml_command(self):
  #  assert 1 == 0
  #  
  #def test_number_of_sequences_in_alignment(self):
  #  assert 1 == 0
  #  
  #def test_get_sequence_names_from_alignment(self):
  #  assert 1 == 0
  #  
  #def test_filter_out_alignments_with_too_much_missing_data(self):
  #  assert 1 == 0
  #  
  #def test_filter_out_removed_taxa_from_tree_and_return_new_file(self):
  #  assert 1 == 0
  #  
  #def test_reinsert_gaps_into_fasta_file(self):
  #  assert 1 == 0
  #  
  #def test_reconvert_fasta_file(self):
  #  assert 1 == 0
  #  
  #def test_pairwise_comparison(self):
  #  assert 1 == 0
  #  
  #def test_create_pairwise_newick_tree(self):
  #  assert 1 == 0
  #  
  #def test_delete_files_based_on_list_of_regexes(self):
  #  assert 1 == 0
  #  
  #def test_use_bundled_exec(self):
  #  assert 1 == 0
    
        
if __name__ == "__main__":
    unittest.main()

