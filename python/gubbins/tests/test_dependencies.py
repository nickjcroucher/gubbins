#! /usr/bin/env python3
# encoding: utf-8

"""
Integration testing of external dependencies. Likely to be the most brittle tests, but the most important.
"""

import unittest
import os
import sys
import glob
import argparse
import pkg_resources
from gubbins import common, run_gubbins

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
working_dir = os.path.join(modules_dir, 'tests')

class TestExternalDependencies(unittest.TestCase):

    # Test individual tree builders
    def test_fasttree(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "fasttree",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_iqtree(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                "--verbose", "--iterations", "3",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxml(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxmlng(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "raxmlng",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_rapidnj(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "rapidnj",
                                                    "--model","JC",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test initial star tree
    def test_starting_star_fasttree(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "fasttree",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_starting_star_iqtree(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "iqtree",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_starting_star_raxml(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_starting_star_raxmlng(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "raxmlng",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_starting_star_rapidnj(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "rapidnj",
                                                    "--model","JC",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test sequence reconstruction
    def test_raxml_seq_recon(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--seq-recon", "raxml",
                                                    "--mar",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_iqtree_seq_recon(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--seq-recon", "iqtree",
                                                    "--mar",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxmlng_seq_recon(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--seq-recon", "raxmlng",
                                                    "--mar",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test model fitting
    def test_fasttree_model_fit(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--model-fitter", "fasttree",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_iqtree_model_fit(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--model-fitter", "iqtree",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxml_model_fit(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--model-fitter", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxmlng_model_fit(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--model-fitter", "raxmlng",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('multiple_recombinations')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test renaming of final output
    def test_rename_final_output(self):
        exit_code = 0
        parser = run_gubbins.parse_input_args()
        try:
            common.parse_and_run(parser.parse_args(["--prefix", "different_prefix",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
            self.check_for_output_files('different_prefix')
        except:
            sys.stderr.write("Test failed")
            exit_code = 1
        self.cleanup('different_prefix')
        assert exit_code == 0

    def test_cleanup(self):
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree", "--verbose", os.path.join(data_dir, 'multiple_recombinations.aln')]))
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
        os.chdir(working_dir)
        regex_to_remove = prefix + ".*"
        for file in glob.glob(regex_to_remove):
            os.remove(file)
        try:
            tmp_to_remove = "tmp*/*"
            for file in glob.glob(tmp_to_remove):
                os.remove(file)
            if os.isdir("tmp*"):
                os.rmdir("tmp*")
        except:
            sys.stderr.write("No tmp directory to remove\n")

if __name__ == "__main__":
    unittest.main(buffer=True)
