#! /usr/bin/env python3
# encoding: utf-8

"""
Integration testing of external dependencies. Likely to be the most brittle tests, but the most important.
"""

import unittest
import os
import glob
import argparse
import pkg_resources
from gubbins import common

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestExternalDependencies(unittest.TestCase):

    def test_fasttree(self):
        self.assertEqual(2, 1+1)
        parser = self.default_arg_parse()
        common.parse_and_run(parser.parse_args(["--tree_builder", "fasttree",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_iqtree(self):
        parser = self.default_arg_parse()
        common.parse_and_run(parser.parse_args(["--tree_builder", "iqtree",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_raxml(self):
        parser = self.default_arg_parse()
        common.parse_and_run(parser.parse_args(["--tree_builder", "raxml",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_hybrid(self):
        parser = self.default_arg_parse()
        common.parse_and_run(parser.parse_args(["--tree_builder", "hybrid",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_rename_final_output(self):
        parser = self.default_arg_parse()
        common.parse_and_run(parser.parse_args(["--prefix", "different_prefix",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('different_prefix')
        self.cleanup('different_prefix')

    def test_cleanup(self):
        parser = self.default_arg_parse()
        common.parse_and_run(parser.parse_args([os.path.join(data_dir, 'multiple_recombinations.aln')]))
        assert not glob.glob('multiple_recombinations.aln.*')
        assert not glob.glob('multiple_recombinations.iteration*')
        self.cleanup('multiple_recombinations')

    @staticmethod
    def check_for_output_files(prefix):
        assert os.path.exists(prefix + '.summary_of_snp_distribution.vcf')
        assert os.path.exists(prefix + '.recombination_predictions.embl')
        assert os.path.exists(prefix + '.per_branch_statistics.csv')
        assert os.path.exists(prefix + '.filtered_polymorphic_sites.fasta')
        assert os.path.exists(prefix + '.filtered_polymorphic_sites.phylip')
        assert os.path.exists(prefix + '.recombination_predictions.gff')
        assert os.path.exists(prefix + '.branch_base_reconstruction.embl')
        assert os.path.exists(prefix + '.final_tree.tre')
        assert os.path.exists(prefix + '.node_labelled.final_tree.tre')

    @staticmethod
    def cleanup(prefix):
        regex_to_remove = prefix + ".*"
        for file in glob.glob(regex_to_remove):
            os.remove(file)

    @staticmethod
    def default_arg_parse():
        # Keep this up to date with run_gubbins.py!
        parser = argparse.ArgumentParser()
        parser.add_argument('alignment_filename', help='Multifasta alignment file')
        parser.add_argument('--outgroup', '-o', help='Outgroup name for rerooting. A list of comma separated names '
                                                     'can be used if they form a clade')
        parser.add_argument('--starting_tree', '-s', help='Starting tree')
        parser.add_argument('--use_time_stamp', '-u', help='Use a time stamp in file names', action='store_true')
        parser.add_argument('--verbose', '-v', help='Turn on debugging', action='store_true')
        parser.add_argument('--no_cleanup', '-n', help='Dont cleanup intermediate files', action='store_true')
        parser.add_argument('--tree_builder', '-t', help='Application to use for tree building', default="raxml",
                            choices=['raxml', 'fasttree', 'hybrid', 'iqtree'])
        parser.add_argument('--iterations', '-i', help='Maximum No. of iterations', type=int, default=5)
        parser.add_argument('--min_snps', '-m', help='Min SNPs to identify a recombination block', type=int, default=3)
        parser.add_argument('--filter_percentage', '-f', help='Filter out taxa with more than this percentage of gaps',
                            type=int, default=25)
        parser.add_argument('--prefix', '-p', help='Add a prefix to the final output filenames')
        parser.add_argument('--threads', '-c', help='Number of threads to run with RAXML, but only if a PTHREADS '
                                                    'version is available', type=int, default=1)
        parser.add_argument('--converge_method', '-z', help='Criteria to use to know when to halt iterations',
                            default='weighted_robinson_foulds', choices=['weighted_robinson_foulds', 'robinson_foulds',
                                                                         'recombination'])
        parser.add_argument('--version', action='version',
                            version=str(pkg_resources.get_distribution("gubbins").version))
        parser.add_argument('--min_window_size', '-a', help='Minimum window size', type=int, default=100)
        parser.add_argument('--max_window_size', '-b', help='Maximum window size', type=int, default=10000)
        parser.add_argument('--raxml_model', '-r', help='RAxML model', default='GTRCAT', choices=['GTRGAMMA', 'GTRCAT'])
        parser.add_argument('--remove_identical_sequences', '-d', help='Remove identical sequences',
                            action='store_true')
        return parser


if __name__ == "__main__":
    unittest.main(buffer=True)
