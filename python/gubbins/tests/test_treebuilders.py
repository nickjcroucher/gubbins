#! /usr/bin/env python3
# encoding: utf-8

"""
Tests the construction of file names and commands with no dependencies on external applications.
"""

import unittest.mock
import os
import filecmp
from gubbins import treebuilders, utils

modules_dir = os.path.dirname(os.path.abspath(treebuilders.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestStringConstruction(unittest.TestCase):

    def test_fasttree_treebuilding_command(self):
        fasttree = treebuilders.FastTree()
        fasttree.executable = 'DDD'
        fasttree.tree_building_parameters = ['EEE']
        assert fasttree.tree_building_command('AAA', '', 'CCC') == 'DDD EEE -out CCC.tre AAA > /dev/null 2>&1'
        fasttree.verbose = True
        assert fasttree.tree_building_command('AAA', 'BBB', 'CCC') == 'DDD EEE -intree BBB -out CCC.tre AAA'

    def test_iqtree_treebuilding_command(self):
        iqtree = treebuilders.IQTree(8)
        iqtree.executable = 'DDD'
        iqtree.tree_building_parameters = ['EEE']
        assert iqtree.tree_building_command('AAA', '', 'CCC') == 'DDD EEE -s AAA -pre CCC -nt 8 -quiet'
        iqtree.verbose = True
        iqtree.threads = 0
        assert iqtree.tree_building_command('AAA', 'BBB', 'CCC') == 'DDD EEE -s AAA -pre CCC -nt AUTO -t BBB'

    def test_iqtree_sequence_reconstruction_command(self):
        iqtree = treebuilders.IQTree(8)
        iqtree.executable = 'DDD'
        iqtree.internal_sequence_reconstruction_parameters = ['EEE']
        assert iqtree.internal_sequence_reconstruction_command('AAA', 'BBB', 'CCC') \
               == 'DDD EEE -s AAA -pre CCC -nt 8 -te BBB -quiet'

    def test_iqtree_convert_raw_ancestral_states_to_fasta(self):
        iqtree = treebuilders.IQTree(8)
        iqtree.convert_raw_ancestral_states_to_fasta(os.path.join(data_dir, 'iqtree_ancestral.state'),
                                                     os.path.join(data_dir, 'iqtree_ancestral.fasta.actual'))
        assert filecmp.cmp(os.path.join(data_dir, 'iqtree_ancestral.fasta.actual'),
                           os.path.join(data_dir, 'iqtree_ancestral.fasta.expected'))
        os.remove(os.path.join(data_dir, 'iqtree_ancestral.fasta.actual'))

    def test_iqtree_replace_internal_node_label(self):
        iqtree = treebuilders.IQTree(8, 'AAA')
        assert iqtree.replace_internal_node_label('Node20') == 'AAA20'

    def test_raxml_treebuilding_command(self):
        raxml = treebuilders.RAxML(8)
        raxml.executable = 'DDD'
        raxml.tree_building_parameters = ['EEE']
        assert raxml.tree_building_command('AAA', '', 'CCC') == 'DDD EEE -s AAA -n CCC -T 8 > /dev/null 2>&1'
        raxml.verbose = True
        raxml.threads = 1
        assert raxml.tree_building_command('AAA', 'BBB', 'CCC') == 'DDD EEE -s AAA -n CCC -t BBB'

    def test_raxml_sequence_reconstruction_command(self):
        raxml = treebuilders.RAxML(8)
        raxml.executable = 'DDD'
        raxml.internal_sequence_reconstruction_parameters = ['EEE']
        assert raxml.internal_sequence_reconstruction_command('AAA', 'BBB', 'CCC') \
               == 'DDD EEE -s AAA -n CCC -T 8 -t BBB > /dev/null 2>&1'

    @unittest.mock.patch('gubbins.utils.choose_executable_based_on_processor')
    def test_raxml_select_executable_based_on_threads(self, mock_choice):
        self.assertIs(mock_choice, utils.choose_executable_based_on_processor)

        raxml = treebuilders.RAxML(8)
        raxml.single_threaded_executables = ['BBB', 'CCC']
        raxml.multi_threaded_executables = ['DDD', 'EEE']

        mock_choice.return_value = 'AAA'
        assert raxml.select_executable_based_on_threads() == 'AAA'
        mock_choice.assert_called_with(['DDD', 'EEE'])

        raxml.threads = 1
        mock_choice.reset_mock()
        mock_choice.return_value = None
        assert raxml.select_executable_based_on_threads() is None
        mock_choice.assert_any_call(['BBB', 'CCC'])
        mock_choice.assert_any_call(['DDD', 'EEE'])

    def test_raxml_convert_raw_ancestral_states_to_fasta(self):
        raxml = treebuilders.RAxML(8)
        raxml.convert_raw_ancestral_states_to_fasta(os.path.join(data_dir, 'raxml_ancestral.state'),
                                                    os.path.join(data_dir, 'raxml_ancestral.fasta.actual'))
        assert filecmp.cmp(os.path.join(data_dir, 'raxml_ancestral.fasta.actual'),
                           os.path.join(data_dir, 'raxml_ancestral.fasta.expected'))
        os.remove(os.path.join(data_dir, 'raxml_ancestral.fasta.actual'))

    def test_raxml_replace_internal_node_label(self):
        raxml = treebuilders.RAxML(8, internal_node_prefix='AAA')
        assert raxml.replace_internal_node_label('Node20') == 'AAANode20'


if __name__ == "__main__":
    unittest.main()
