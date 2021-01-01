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

    def test_iqtree_convert_raw_ancestral_states_to_fasta(self):
        iqtree = treebuilders.IQTree(1, model = 'GTRGAMMA')
        iqtree.convert_raw_ancestral_states_to_fasta(os.path.join(data_dir, 'iqtree_ancestral.state'),
                                                     os.path.join(data_dir, 'iqtree_ancestral.fasta.actual'))
        assert filecmp.cmp(os.path.join(data_dir, 'iqtree_ancestral.fasta.actual'),
                           os.path.join(data_dir, 'iqtree_ancestral.fasta.expected'))
        os.remove(os.path.join(data_dir, 'iqtree_ancestral.fasta.actual'))

    def test_iqtree_replace_internal_node_label(self):
        iqtree = treebuilders.IQTree(1, model = 'GTRGAMMA', internal_node_prefix='AAA')
        assert iqtree.replace_internal_node_label('Node20') == 'AAA20'

    def test_raxml_select_executable_based_on_threads(self):

        raxml_st = treebuilders.RAxML(1, model = 'GTRGAMMA')
        raxml_mt = treebuilders.RAxML(8, model = 'GTRGAMMA')
        assert raxml_st.select_executable_based_on_threads() != raxml_mt.select_executable_based_on_threads()

    def test_raxml_convert_raw_ancestral_states_to_fasta(self):
        raxml = treebuilders.RAxML(8)
        raxml.convert_raw_ancestral_states_to_fasta(os.path.join(data_dir, 'raxml_ancestral.state'),
                                                    os.path.join(data_dir, 'raxml_ancestral.fasta.actual'))
        assert filecmp.cmp(os.path.join(data_dir, 'raxml_ancestral.fasta.actual'),
                           os.path.join(data_dir, 'raxml_ancestral.fasta.expected'))
        os.remove(os.path.join(data_dir, 'raxml_ancestral.fasta.actual'))

    def test_raxml_replace_internal_node_label(self):
        raxml = treebuilders.RAxML(1, internal_node_prefix='AAA')
        assert raxml.replace_internal_node_label('Node20') == 'AAANode20'


if __name__ == "__main__":
    unittest.main()
