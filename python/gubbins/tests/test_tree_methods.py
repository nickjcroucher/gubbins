#! /usr/bin/env python3
# encoding: utf-8

"""
Test manipulating and interogating trees using python modules with no external application dependancies.
"""

import unittest
import shutil
import os
import filecmp
from gubbins import common, treebuilders

modules_dir = os.path.dirname(os.path.abspath(treebuilders.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestTreeMethods(unittest.TestCase):

    def initialise_Star(self):
        algorithm = treebuilders.Star()
        assert algorithm.name == 'Star'

    def initialise_RapidNJ(self):
        algorithm = treebuilders.RapidNJ()
        assert algorithm.name == 'RapidNJ'

    def initialise_FastTree(self):
        algorithm = treebuilders.FastTree()
        assert algorithm.name == 'FastTree'

    def initialise_IQTree(self):
        algorithm = treebuilders.IQTree()
        assert algorithm.name == 'IQTree'

    def initialise_RAxML(self):
        algorithm = treebuilders.RAxML()
        assert algorithm.name == 'RAxML'

    def initialise_RAxMLNG(self):
        algorithm = treebuilders.RAxMLNG()
        assert algorithm.name == 'RAxMLNG'

    def initialise_VeryFastTree(self):
        algorithm = treebuilders.VeryFastTree()
        assert algorithm.name == 'VeryFastTree'

    def test_robinson_foulds_distance(self):
        # two tree with different distances
        assert common.robinson_foulds_distance(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                                               os.path.join(data_dir, 'robinson_foulds_distance_tree2.tre')) == 17.263494
        # two trees with same distance
        assert common.robinson_foulds_distance(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                                               os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre')) == 0

    def test_has_tree_been_seen_before(self):
        # same content so the tree has been seen before
        assert common.has_tree_been_seen_before([os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                                                 os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                                                 os.path.join(data_dir, 'robinson_foulds_distance_tree1_dup.tre')],
                                                'weighted_robinson_foulds') == 1
        # different trees
        assert common.has_tree_been_seen_before([os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                                                 os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                                                 os.path.join(data_dir, 'robinson_foulds_distance_tree2.tre')],
                                                'weighted_robinson_foulds') == 0

    def test_root_tree(self):
        common.root_tree(os.path.join(data_dir, 'unrooted_tree.newick'), os.path.join(data_dir, 'actual_rooted_tree.newick'))
        assert filecmp.cmp(os.path.join(data_dir, 'actual_rooted_tree.newick'),
                           os.path.join(data_dir, 'expected_rooted_tree.newick'), shallow=False)
        os.remove(os.path.join(data_dir, 'actual_rooted_tree.newick'))

    def test_reroot_tree(self):
        shutil.copyfile(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                        os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_test'))
        common.reroot_tree(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_test'),
                           'sequence_4')
        assert filecmp.cmp(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_test'),
                           os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_expected'))
        os.remove(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_test'))

        shutil.copyfile(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                        os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_test'))
        common.reroot_tree(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_test'), '')
        assert filecmp.cmp(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_test'),
                           os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_expected'))
        os.remove(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_test'))

    def test_reroot_tree_with_outgroup(self):
        shutil.copyfile(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                        os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_list_test'))
        common.reroot_tree_with_outgroup(
            os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_list_test'), ['sequence_4'])
        assert filecmp.cmp(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_list_test'),
                           os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_expected'))
        os.remove(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_list_test'))

    def test_reroot_tree_with_outgroups(self):
        shutil.copyfile(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                        os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_2_test'))
        common.reroot_tree_with_outgroup(
            os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_2_test'),
            ['sequence_4', 'sequence_2'])
        assert filecmp.cmp(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_2_test'),
                           os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_2_expected'))
        os.remove(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_at_sequence_4_2_test'))

    def test_reroot_tree_with_outgroups_all_in_one_clade(self):
        outgroups = ['A', 'B']
        expected_monophyletic_outgroup = ['A', 'B']
        expected_output_file = os.path.join(data_dir, 'expected_reroot_tree_with_outgroups_all_in_one_clade.tre')
        self.reroot_tree_check(outgroups, expected_output_file, expected_monophyletic_outgroup)

    def test_reroot_tree_with_outgroups_all_in_one_clade_large(self):
        outgroups = ['A', 'B', 'C']
        expected_monophyletic_outgroup = ['A', 'B', 'C']
        expected_output_file = os.path.join(data_dir, 'expected_test_reroot_tree_with_outgroups_all_in_one_clade_large.tre')
        self.reroot_tree_check(outgroups, expected_output_file, expected_monophyletic_outgroup)

    def test_reroot_tree_with_outgroups_all_in_different_clade(self):
        outgroups = ['A', 'D']
        expected_monophyletic_outgroup = ['A']
        expected_output_file = os.path.join(data_dir, 'expected_reroot_tree_with_outgroups_all_in_different_clade.tre')
        self.reroot_tree_check(outgroups, expected_output_file, expected_monophyletic_outgroup)

    def test_reroot_tree_with_outgroups_with_two_mixed_clades(self):
        outgroups = ['A', 'B', 'C', 'D']
        expected_monophyletic_outgroup = ['A']
        expected_output_file = os.path.join(data_dir, 'expected_reroot_tree_with_outgroups_with_two_mixed_clades.tre')
        self.reroot_tree_check(outgroups, expected_output_file, expected_monophyletic_outgroup)

    def reroot_tree_check(self, outgroups, expected_output_file, expected_monophyletic_outgroup):
        tmp_tree_file_name = expected_output_file.replace('expected_','generated_')
        shutil.copyfile(os.path.join(data_dir, 'outgroups_input.tre'), tmp_tree_file_name)
        assert expected_monophyletic_outgroup == common.get_monophyletic_outgroup(tmp_tree_file_name, outgroups)
        common.reroot_tree_with_outgroup(tmp_tree_file_name, outgroups)
        assert filecmp.cmp(tmp_tree_file_name, expected_output_file)
        os.remove(tmp_tree_file_name)

    def test_split_all_non_bi_nodes(self):
        # best way to access it is via reroot_tree_at_midpoint because it outputs to a file
        shutil.copyfile(os.path.join(data_dir, 'non_bi_tree.tre'), os.path.join(data_dir, 'non_bi_tree.tre.actual'))
        common.reroot_tree_at_midpoint(os.path.join(data_dir, 'non_bi_tree.tre.actual'))
        assert filecmp.cmp(os.path.join(data_dir, 'non_bi_tree.tre.actual'), os.path.join(data_dir, 'non_bi_tree.tre.expected'))
        os.remove(os.path.join(data_dir, 'non_bi_tree.tre.actual'))

    def test_reroot_tree_at_midpoint(self):
        shutil.copyfile(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
                        os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual'))
        common.reroot_tree_at_midpoint(
            os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual'))
        assert filecmp.cmp(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual'),
                           os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_expected'))
        os.remove(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.reroot_tree_at_midpoint_actual'))

    def test_filter_out_removed_taxa_from_tree(self):
        common.filter_out_removed_taxa_from_tree(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre'),
             os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.filter_out_removed_taxa_from_tree_actual'),
             ['sequence_1', 'sequence_2', 'sequence_3', 'sequence_4', 'sequence_5'])
        assert filecmp.cmp(
            os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.filter_out_removed_taxa_from_tree_actual'),
            os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.filter_out_removed_taxa_from_tree_expected'))
        os.remove(os.path.join(data_dir, 'robinson_foulds_distance_tree1.tre.filter_out_removed_taxa_from_tree_actual'))

    def test_internal_node_taxons_removed_when_used_as_starting_tree(self):
        common.filter_out_removed_taxa_from_tree(os.path.join(data_dir, 'tree_with_internal_nodes.tre'),
                                                 os.path.join(data_dir, 'tree_with_internal_nodes.tre_actual'), [])
        assert filecmp.cmp(os.path.join(data_dir, 'tree_with_internal_nodes.tre_actual'),
                           os.path.join(data_dir, 'tree_with_internal_nodes.tre_expected'))
        os.remove(os.path.join(data_dir, 'tree_with_internal_nodes.tre_actual'))

    def test_transfer_internal_node_labels_to_tree(self):
        reconstructor = treebuilders.IQTree(1, model = 'GTRGAMMA', internal_node_prefix = 'internal_')
        common.transfer_internal_node_labels_to_tree(os.path.join(data_dir, 'source_tree.tre'),
                                                     os.path.join(data_dir, 'destination_tree.tre'),
                                                     os.path.join(data_dir, 'renamed_output_tree'), reconstructor)
        assert os.path.exists(os.path.join(data_dir, 'renamed_output_tree'))
        assert filecmp.cmp(os.path.join(data_dir, 'renamed_output_tree'),
                           os.path.join(data_dir, 'expected_renamed_output_tree'), shallow=False)
        os.remove(os.path.join(data_dir, 'renamed_output_tree'))

    def test_remove_internal_node_labels(self):
        common.remove_internal_node_labels_from_tree(os.path.join(data_dir, 'final_tree_with_internal_labels.tre'),
                                                     'final_tree_with_internal_labels.tre')
        assert os.path.exists('final_tree_with_internal_labels.tre')
        assert filecmp.cmp('final_tree_with_internal_labels.tre',
                           os.path.join(data_dir, 'expected_final_tree_without_internal_labels.tre'))
        os.remove('final_tree_with_internal_labels.tre')


if __name__ == "__main__":
    unittest.main()
