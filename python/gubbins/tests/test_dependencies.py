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
        common.parse_and_run(parser.parse_args(["--tree-builder", "fasttree",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_iqtree(self):
        parser = self.default_arg_parse()
        common.parse_and_run(parser.parse_args(["--tree-builder", "iqtree",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_raxml(self):
        parser = self.default_arg_parse()
        common.parse_and_run(parser.parse_args(["--tree-builder", "raxml",
                                                os.path.join(data_dir, 'multiple_recombinations.aln')]))
        self.check_for_output_files('multiple_recombinations')
        self.cleanup('multiple_recombinations')

    def test_hybrid(self):
        parser = self.default_arg_parse()
        common.parse_and_run(parser.parse_args(["--tree-builder", "hybrid",
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
        parser = argparse.ArgumentParser(
            description='Croucher N. J., Page A. J., Connor T. R., Delaney A. J., Keane J. A., Bentley S. D., Parkhill J., '
                        'Harris S.R. "Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome '
                        'sequences using Gubbins". Nucleic Acids Res. 2015 Feb 18;43(3):e15. doi: 10.1093/nar/gku1196.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        ioGroup = parser.add_argument_group('Input and output options')
        ioGroup.add_argument('alignment_filename',        help='Multifasta alignment file')
        ioGroup.add_argument('--prefix',            '-p', help='Add a prefix to the final output filenames')
        ioGroup.add_argument('--starting-tree',     '-s', help='Starting tree')
        ioGroup.add_argument('--use-time-stamp',    '-u', help='Use a time stamp in file names', action='store_true')
        ioGroup.add_argument('--version',                 action='version',
                                                          version=str(pkg_resources.get_distribution("gubbins").version))
                            
        dataGroup = parser.add_argument_group('Data processing options')
        dataGroup.add_argument('--pairwise',              help='Compare two sequences (without using a tree)',
                                                          default = False, action = 'store_true') # fasttree model fit, star phylogeny, one iteration
        dataGroup.add_argument('--filter-percentage','-f',
                                                          help='Filter out taxa with more than this percentage of gaps',
                                                          type=float, default=25.0)
        dataGroup.add_argument('--remove-identical-sequences',
                                                    '-d', help='Remove identical sequences', action='store_true')
        dataGroup.add_argument('--threads',         '-c', help='Number of threads to run with RAXML, but only if a PTHREADS '
                                                          'version is available', type=int,  default=1)
        dataGroup.add_argument('--verbose',         '-v', help='Turn on debugging', action='store_true')
        dataGroup.add_argument('--no-cleanup',      '-n', help='Do not cleanup intermediate files', action='store_true')

        treeGroup = parser.add_argument_group('Tree building options')
        treeGroup.add_argument('--tree-builder',    '-t', help='Application to use for tree building',
                                                          default='raxml',
                                                          choices=['raxml', 'raxmlng', 'iqtree', 'fasttree', 'hybrid', 'rapidnj'])
        treeGroup.add_argument('--tree-args',             help='Quoted string of further arguments passed to tree building algorithm'
                                                          ' (start string with a space if there is a risk of being interpreted as a flag)',
                                                          default = None)
        treeGroup.add_argument('--first-tree-builder',    help='Application to use for building the first tree',
                                                          default=None,
                                                          choices=['raxml', 'raxmlng', 'iqtree', 'fasttree', 'rapidnj', 'star'])
        treeGroup.add_argument('--first-tree-args',       help='Further arguments passed to first tree building algorithm',
                                                          default = None)
        treeGroup.add_argument('--outgroup',        '-o', help='Outgroup name for rerooting. A list of comma separated '
                                                          'names can be used if they form a clade')
        treeGroup.add_argument('--bootstrap',       '-#', help='Number of bootstrap replicates to perform with final alignment '
                                                          '[default = 0]', type = int, default = 0)
        treeGroup.add_argument('--transfer-bootstrap',    help='Calculate bootstrap supporting transfer bootstrap expectation '
                                                          '[default = False]', default = False, action = 'store_true')
        treeGroup.add_argument('--sh-test',               help='Perform an SH test of node likelihoods', default = False,
                                                          action = 'store_true')
                                                              
        modelGroup = parser.add_argument_group('Nucleotide substitution model options')
        modelGroup.add_argument('--model-fitter',   '-F', help='Application to use for model fitting [default = same as'
                                                          ' tree builder if possible, else raxml]',
                                                          default = None,
                                                          choices=['raxml', 'raxmlng', 'iqtree', 'fasttree', None])
        modelGroup.add_argument('--model',          '-M', help='Nucleotide substitution model (not all available for all'
                                                          'tree building algorithms',
                                                          default='GTRGAMMA',
                                                          choices=['JC','K2P','HKY','GTR','GTRGAMMA','GTRCAT'])
        modelGroup.add_argument('--model-args',           help='Quoted string of further arguments passed to model fitting algorithm'
                                                          ' (start string with a space if there is a risk of being interpreted as a flag)'
                                                          '(default = same as --tree-builder-args)',
                                                          default=None)
        modelGroup.add_argument('--custom-model',         help='String corresponding to a substitution model for the selected tree'
                                                          ' building algorithm [default = None]', default = None)
        modelGroup.add_argument('--first-model-fitter',   help='Application to use for model fitting in first iteration'
                                                          ' [default = same as tree builder if possible, else raxml]',
                                                          default = None,
                                                          choices=['raxml', 'raxmlng', 'iqtree', 'fasttree', None])
        modelGroup.add_argument('--first-model',          help='Nucleotide substitution model used for first tree',
                                                          default=None,
                                                          choices=['JC','K2P','HKY','GTR','GTRGAMMA','GTRCAT'])
        modelGroup.add_argument('--first-model-args',     help='Further arguments passed to model fitting algorithm used in first'
                                                          'iteration (default = same as --first-tree-builder-args)',
                                                          default=None)
        modelGroup.add_argument('--custom-first-model',   help='String corresponding to a substitution model for the selected tree'
                                                          ' building algorithm for the first iteration [default = None]',
                                                          default = None)

        reconGroup = parser.add_argument_group('Ancestral sequence reconstruction options')
        reconGroup.add_argument('--mar',                  help='Use marginal, rather than joint, ancestral reconstruction',
                                                          action='store_true')
        reconGroup.add_argument('--seq-recon',            help='Algorithm to use for marginal reconstruction [default = '
                                                          'same as tree builder if possible, else raxml]',
                                                          default=None,
                                                          choices=['raxml', 'raxmlng', 'iqtree', None])
        reconGroup.add_argument('--seq-recon-args',       help='Further arguments passed to sequence reconstruction algorithm'
                                                          ' (start string with a space if there is a risk of being interpreted as a flag)',
                                                          default=None)
                                                                
        gubbinsGroup = parser.add_argument_group('Recombination detection options')
        gubbinsGroup.add_argument('--min-snps',     '-m', help='Min SNPs to identify a recombination block',
                                                          type=int,
                                                          default = 3)
        gubbinsGroup.add_argument('--min-window-size','-a',
                                                          help='Minimum window size', type=int, default=100)
        gubbinsGroup.add_argument('--max-window-size','-b',
                                                          help='Maximum window size', type=int, default=10000)

        stopGroup = parser.add_argument_group('Algorithm stop options')
        stopGroup.add_argument('--iterations',      '-i', help='Maximum No. of iterations', type=int, default=5)
        stopGroup.add_argument('--converge-method', '-z', help='Criteria to use to know when to halt iterations',
                            default='weighted_robinson_foulds', choices=['weighted_robinson_foulds', 'robinson_foulds',
                                                                         'recombination'])
        return parser

if __name__ == "__main__":
    unittest.main(buffer=True)
