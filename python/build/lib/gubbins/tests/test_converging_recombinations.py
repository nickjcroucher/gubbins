#! /usr/bin/env python3
# encoding: utf-8

"""
Tests if the same recombinations are seen in multiple iterations
"""

import unittest
import os
from gubbins import common

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestConvergingRecombinations(unittest.TestCase):

    def test_reading_embl_file(self):
        assert common.extract_recombinations_from_embl(os.path.join(data_dir, 'small_recombination.embl')) \
               == {'sequence1': [[1, 3], [5, 7], [9, 11]], 'sequence2': [[1, 3], [9, 11]], 'sequence3': [[1, 3]]}

    def test_two_files_are_different(self):
        assert not common.have_recombinations_been_seen_before(os.path.join(data_dir, 'small_recombination.embl'),
                             [os.path.join(data_dir, 'small_recombination_different.embl')])

    def test_two_files_are_same(self):
        assert common.have_recombinations_been_seen_before(os.path.join(data_dir, 'small_recombination.embl'),
                                                           [os.path.join(data_dir, 'small_recombination.embl')])

    def test_multiple_files_have_one_match(self):
        assert common.have_recombinations_been_seen_before(os.path.join(data_dir, 'small_recombination.embl'),
            [os.path.join(data_dir, 'small_recombination_different.embl'),
             os.path.join(data_dir, 'small_recombination.embl'),
             os.path.join(data_dir, 'small_recombination_different.embl')])

    def test_get_recombination_files(self):
        assert common.get_recombination_files(['AAA', 'BBB', 'CCC']) == ('CCC.tab', ['AAA.tab', 'BBB.tab'])


if __name__ == "__main__":
    unittest.main()
