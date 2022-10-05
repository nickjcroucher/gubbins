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
import shutil
from gubbins import common, run_gubbins

unittest.TestLoader.sortTestMethodsUsing = None
modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
working_dir = os.path.join(modules_dir, 'tests')

class TestExternalDependencies(unittest.TestCase):

    # Test pairwise comparisons
    def test_pairwise(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--pairwise",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'pairwise.aln')]))
        exit_code = self.check_for_output_files('pairwise')
        self.cleanup('pairwise')
        assert exit_code == 0

    def test_pairwise_altered_rec_options(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--pairwise",
                                                    "--threads", "1",
                                                    "--p-value", "0.01",
                                                    "--trimming-ratio", "2.0",
                                                    os.path.join(data_dir, 'pairwise.aln')]))
        exit_code = self.check_for_output_files('pairwise')
        self.cleanup('pairwise')
        assert exit_code == 0

    def test_pairwise_altered_recon(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--pairwise",
                                                    "--threads", "1",
                                                    "--mar",
                                                    os.path.join(data_dir, 'pairwise.aln')]))
        exit_code = self.check_for_output_files('pairwise')
        self.cleanup('pairwise')
        assert exit_code == 0

    def test_for_outgroup_present(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "fasttree",
                                                    "--verbose", "--iterations", "3",
                                                    "--outgroup", "sequence_10",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test individual tree builders
    def test_fasttree(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "fasttree",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test resuming a default analysis
    def test_fasttree_resume(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "fasttree",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    "--no-cleanup",
                                                    "--prefix", "original",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        common.parse_and_run(parser.parse_args(["--tree-builder", "fasttree",
                                            "--verbose", "--iterations", "3",
                                            "--resume", "multiple_recombinations.iteration_1.tre",
                                            "--threads", "1",
                                            os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('original')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_iqtree(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                "--verbose", "--iterations", "3",
                                                "--threads", "1",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_iqtree_best_model(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                "--best-model",
                                                "--verbose", "--iterations", "3",
                                                "--threads", "1",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_iqtree_fast(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree-fast",
                                                "--verbose", "--iterations", "3",
                                                "--threads", "1",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_iqtree_with_dates(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                "--verbose", "--iterations", "3",
                                                "--date",os.path.join(data_dir, 'taxon.times'),
                                                "--threads", "1",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_iqtree_with_dates_and_outgroup(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                "--verbose", "--iterations", "3",
                                                "--date", os.path.join(data_dir, 'taxon.times'),
                                                "--outgroup", "sequence_10",
                                                "--threads", "1",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test resuming a default analysis
    def test_iqtree_resume(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    "--no-cleanup",
                                                    "--prefix", "original",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                            "--verbose", "--iterations", "3",
                                            "--resume", "multiple_recombinations.iteration_1.tre",
                                            "--threads", "1",
                                            os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('original')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test resuming a default analysis
    def test_iqtree_fast_resume(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree-fast",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    "--no-cleanup",
                                                    "--prefix", "original",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree-fast",
                                            "--verbose", "--iterations", "3",
                                            "--resume", "multiple_recombinations.iteration_1.tre",
                                            "--threads", "1",
                                            os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('original')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxml(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test resuming a default analysis
    def test_raxml_resume(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    "--no-cleanup",
                                                    "--prefix", "original",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                            "--verbose", "--iterations", "3",
                                            "--resume", "multiple_recombinations.iteration_1.tre",
                                            "--threads", "1",
                                            os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('original')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxml_quiet(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxmlng(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxmlng",
                                                    "--model","GTR",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test resuming a default analysis
    def test_raxmlng_resume(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxmlng",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    "--no-cleanup",
                                                    "--prefix", "original",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxmlng",
                                            "--verbose", "--iterations", "3",
                                            "--resume", "multiple_recombinations.iteration_1.tre",
                                            "--threads", "1",
                                            os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('original')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_rapidnj(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "rapidnj",
                                                    "--model","JC",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        # Copy file for subsequent tests
        shutil.copyfile(os.path.join(data_dir, 'multiple_recombinations.recombination_predictions.embl'),
                        os.path.join(data_dir, 'new_rapidnj_jc_output.recombination_predictions.embl'))
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test resuming a default analysis
    def test_rapidnj_resume(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "rapidnj",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    "--model","JC",
                                                    "--no-cleanup",
                                                    "--prefix", "original",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        common.parse_and_run(parser.parse_args(["--tree-builder", "rapidnj",
                                            "--verbose", "--iterations", "3",
                                            "--resume", "multiple_recombinations.iteration_1.tre",
                                            "--model","JC",
                                            "--threads", "1",
                                            os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('original')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def check_rapidnj_consistency(self):
        new_file = os.path.join(data_dir,'new_rapidnj_jc_output.recombination_predictions.embl')
        reference_file = os.path.join(data_dir,'ref_rapidnj_jc_output.recombination_predictions.embl')
        assert common.have_recombinations_been_seen_before(reference_file,[new_file])

    def test_defined_starting_tree(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--starting-tree",
                                                os.path.join(data_dir, 'starting_tree.tre'),
                                                "--tree-builder", "iqtree",
                                                "--verbose", "--iterations", "3",
                                                "--threads", "1",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_defined_mislabelled_starting_tree(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--starting-tree",
                                                os.path.join(data_dir, 'mislabelled.starting_tree.tre'),
                                                "--tree-builder", "iqtree",
                                                "--verbose", "--iterations", "3",
                                                "--threads", "1",
                                                os.path.join(data_dir, 'mislabelled.multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('mislabelled.multiple_recombinations')
        self.cleanup('mislabelled.multiple_recombinations')
        assert exit_code == 0

    # Test initial star tree
    def test_starting_star_fasttree(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "fasttree",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_starting_star_iqtree(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "iqtree",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_starting_star_iqtree_fast(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "iqtree-fast",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_starting_star_raxml(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_starting_star_raxmlng(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "raxmlng",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_starting_star_rapidnj(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--first-tree-builder","star",
                                                    "--tree-builder", "rapidnj",
                                                    "--model","JC",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test sequence reconstruction
    def test_raxml_seq_recon(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--seq-recon", "raxml",
                                                    "--mar",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_iqtree_seq_recon(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--seq-recon", "iqtree",
                                                    "--mar",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxmlng_seq_recon(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--seq-recon", "raxmlng",
                                                    "--mar",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test model fitting
    def test_fasttree_model_fit(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--model-fitter", "fasttree",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_iqtree_model_fit(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--model-fitter", "iqtree",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxml_model_fit(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--model-fitter", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_raxmlng_model_fit(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--model-fitter", "raxmlng",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test different model structures across iterations
    def test_different_model_selections(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--model-fitter", "iqtree",
                                                    "--first-model", "JC",
                                                    "--model", "GTRGAMMA",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test bootstrapping
    def test_raxml_bootstrapping(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    "--bootstrap","10",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    def test_iqtree_bootstrapping(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                    "--verbose", "--iterations", "3",
                                                    "--bootstrap","1000",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    def test_rapidnj_bootstrapping(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "rapidnj",
                                                    "--model","JC",
                                                    "--verbose", "--iterations", "3",
                                                    "--bootstrap","10",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    def test_raxmlng_bootstrapping(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxmlng",
                                                    "--verbose", "--iterations", "3",
                                                    "--bootstrap","100",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    def test_fasttree_bootstrapping(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "fasttree",
                                                    "--verbose", "--iterations", "3",
                                                    "--bootstrap","10",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    def test_fasttree_bootstrapping_with_raxml(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "fasttree",
                                                    "--model-fitter", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    "--bootstrap","10",
                                                    "--threads", "1",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    def test_raxmlng_transfer_bootstrapping(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxmlng",
                                                    "--verbose", "--iterations", "3",
                                                    "--bootstrap","100",
                                                    "--threads", "1",
                                                    "--transfer-bootstrap",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    def test_raxml_transfer_bootstrapping(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    "--bootstrap","10",
                                                    "--threads", "1",
                                                    "--transfer-bootstrap",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    def test_iqtree_transfer_bootstrapping(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                    "--verbose", "--iterations", "3",
                                                    "--bootstrap","1000",
                                                    "--threads", "1",
                                                    "--transfer-bootstrap",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    # Test SH tests
    def test_iqtree_sh_test(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    "--sh-test",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    def test_raxml_sh_test(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    "--sh-test",
                                                    os.path.join(data_dir, 'bootstrapping_test.aln')]))
        exit_code = self.check_for_output_files('bootstrapping_test')
        self.cleanup('bootstrapping_test')
        assert exit_code == 0

    # Test convergence
    def test_converge_on_rec(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    "--converge-method", "recombination",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    def test_converge_on_unweighted_rf(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                    "--verbose", "--iterations", "3",
                                                    "--threads", "1",
                                                    "--converge-method", "robinson_foulds",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')
        assert exit_code == 0

    # Test renaming of final output
    def test_rename_final_output(self):
        exit_code = 1
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--prefix", "different_prefix",
                                                    "--verbose", "--iterations", "3",
                                                    os.path.join(data_dir, 'multiple_recombinations.aln')]))
        exit_code = self.check_for_output_files('different_prefix')
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
        return 0

    @staticmethod
    def cleanup(prefix):
        #os.chdir(working_dir)
        regex_to_remove = prefix + ".*"
        for file in glob.glob(regex_to_remove):
            os.remove(file)
        tmp_to_remove = "./tmp*/*"
        for file in glob.glob(tmp_to_remove):
            os.remove(file)
        for dir in glob.glob("./tmp*"):
            if os.path.isdir(dir):
                os.rmdir(dir)

if __name__ == "__main__":
    unittest.main(buffer=True)
