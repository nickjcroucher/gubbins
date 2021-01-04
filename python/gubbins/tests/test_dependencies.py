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
from gubbins import common, run_gubbins

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestExternalDependencies(unittest.TestCase):

    # Test individual tree builders
    def test_fasttree(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "fasttree",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_iqtree(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_raxml(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_raxmlng(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxmlng",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_rapidnj(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "rapidnj",
                                                "--model","JC",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    # Test initial star tree
    def test_starting_star_fasttree(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                "--tree-builder", "fasttree",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_starting_star_iqtree(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                "--tree-builder", "iqtree",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_starting_star_raxml(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                "--tree-builder", "raxml",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_starting_star_raxmlng(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                "--tree-builder", "raxmlng",
                                                "--verbose",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_starting_star_rapidnj(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                "--tree-builder", "rapidnj",
                                                "--model","JC",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        
    # Test sequence reconstruction
    def test_raxml_seq_recon(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                "--seq-recon", "raxml",
                                                "--mar",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_iqtree_seq_recon(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                "--seq-recon", "iqtree",
                                                "--mar",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_raxmlng_seq_recon(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                "--seq-recon", "raxmlng",
                                                "--mar",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    # Test model fitting
    def test_fasttree_model_fit(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                "--model-fitter", "fasttree",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_iqtree_model_fit(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                "--model-fitter", "iqtree",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_raxml_model_fit(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                "--model-fitter", "raxml",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_raxmlng_model_fit(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                "--model-fitter", "raxmlng",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    # Test renaming of final output
    def test_rename_final_output(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--prefix", "different_prefix",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('different_prefix')
        self.cleanup('different_prefix')

    def test_cleanup(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args([os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.cleanup('multiple_recombinations')
        assert not glob.glob('multiple_recombinations.aln.*')
        assert not glob.glob('multiple_recombinations.iteration*')

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

if __name__ == "__main__":
    unittest.main(buffer=True)
