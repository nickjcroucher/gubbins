#!/usr/bin/env python
# encoding: utf-8
#
# Wellcome Trust Sanger Institute
# Copyright (C) 2012  Wellcome Trust Sanger Institute
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

import sys
import argparse
from gubbins.__init__ import version
import gubbins.common

def parse_input_args():

    parser = argparse.ArgumentParser(
        description='Croucher N. J., Page A. J., Connor T. R., Delaney A. J., Keane J. A., Bentley S. D., Parkhill J., '
                    'Harris S.R. "Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome '
                    'sequences using Gubbins". Nucleic Acids Res. 2015 Feb 18;43(3):e15. doi: 10.1093/nar/gku1196.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ioGroup = parser.add_argument_group('Input and output options')
    ioGroup.add_argument('alignment_filename',        help='Multifasta alignment file')
    ioGroup.add_argument('--prefix',            '-p', help='Add a prefix to the final output filenames')
    ioGroup.add_argument('--starting-tree',     '-s', help='Starting tree')
    ioGroup.add_argument('--date',              '-D', help='Two-column text file in which the second column is the'
                                                      ' date of isolation in YYYY,YYYY-MM or YYYY-MM-DD format',
                                                      default = None)
    ioGroup.add_argument('--use-time-stamp',    '-u', help='Use a time stamp in file names', action='store_true')
    ioGroup.add_argument('--version',                 action='version',
                                                      version = version())
    ioGroup.add_argument('--threads',           '-c', help='Number of threads to use for parallelisation',
                                                      type=int,  default=1)
    ioGroup.add_argument('--verbose',           '-v', help='Turn on debugging', action='store_true')
    ioGroup.add_argument('--no-cleanup',        '-n', help='Do not cleanup intermediate files', action='store_true')

    dataGroup = parser.add_argument_group('Data processing options')
    dataGroup.add_argument('--pairwise',              help='Compare two sequences (without using a tree)',
                                                      default = False, action = 'store_true') # fasttree model fit, star phylogeny, one iteration
    dataGroup.add_argument('--filter-percentage','-f',
                                                      help='Filter out taxa with more than this percentage of gaps',
                                                      type=float, default=25.0)
    dataGroup.add_argument('--remove-identical-sequences',
                                                '-d', help='Remove identical sequences', action='store_true')

    treeGroup = parser.add_argument_group('Tree building options')
    treeGroup.add_argument('--tree-builder',    '-t', help='Application to use for tree building',
                                                      default='raxml',
                                                      choices=['raxml', 'raxmlng', 'iqtree', 'iqtree-fast', 'fasttree', 'hybrid', 'rapidnj'])
    treeGroup.add_argument('--tree-args',             help='Quoted string of further arguments passed to tree building algorithm'
                                                      ' (start string with a space if there is a risk of being interpreted as a flag)',
                                                      default = None)
    treeGroup.add_argument('--first-tree-builder',    help='Application to use for building the first tree',
                                                      default=None,
                                                      choices=['raxml', 'raxmlng', 'iqtree', 'iqtree-fast', 'fasttree', 'rapidnj', 'star'])
    treeGroup.add_argument('--first-tree-args',       help='Further arguments passed to first tree building algorithm',
                                                      default = None)
    treeGroup.add_argument('--outgroup',        '-o', help='Outgroup name for rerooting. A list of comma separated '
                                                      'names can be used if they form a clade',
                                                      default = None)
    treeGroup.add_argument('--bootstrap',       '-#', help='Number of bootstrap replicates to perform with final alignment',
                                                      type = int, default = 0)
    treeGroup.add_argument('--transfer-bootstrap',    help='Calculate bootstrap supporting transfer bootstrap expectation',
                                                      default = False, action = 'store_true')
    treeGroup.add_argument('--sh-test',               help='Perform an SH test of node likelihoods', default = False,
                                                      action = 'store_true')
    treeGroup.add_argument('--seed',                  help='Set seed for reproducibility of analysis',
                                                      default = None, type = int)
                                                          
    modelGroup = parser.add_argument_group('Nucleotide substitution model options')
    modelGroup.add_argument('--model',          '-M', help='Nucleotide substitution model (not all available for all '
                                                      'tree building algorithms)',
                                                      default=None,
                                                      choices=['JC','K2P','HKY','GTR','GTRGAMMA','GTRCAT'])
    modelGroup.add_argument('--first-model',          help='Nucleotide substitution model used for first tree',
                                                      default=None,
                                                      choices=['JC','K2P','HKY','GTR','GTRGAMMA','GTRCAT'])
    modelGroup.add_argument('--best-model',           help='Automatically select best substitution model using iqtree in later iterations',
                                                      default = False, action = 'store_true')
    modelGroup.add_argument('--custom-model',         help='String corresponding to a substitution model for the selected tree'
                                                      ' building algorithm', default = None)
    modelGroup.add_argument('--custom-first-model',   help='String corresponding to a substitution model for the selected tree'
                                                      ' building algorithm for the first iteration',
                                                      default = None)

    reconGroup = parser.add_argument_group('Ancestral sequence reconstruction options')
    reconGroup.add_argument('--model-fitter',   '-F', help='Application to use for model fitting for joint ancestral state'
                                                      ' reconstruction [if unspecified: same as tree builder if possible'
                                                      ', else iqtree]',
                                                      default = None,
                                                      choices=['raxml', 'raxmlng', 'iqtree', 'fasttree', None])
    reconGroup.add_argument('--recon-model',    '-R', help='Nucleotide substitution model used for ancestral state reconstruction'
                                                      ' (not all available for all tree building algorithms)',
                                                      default='GTRGAMMA',
                                                      choices=['JC','K2P','HKY','GTR','GTRGAMMA','GTRCAT'])
    reconGroup.add_argument('--custom-recon-model',   help='String corresponding to a substitution model for the selected '
                                                      ' model fitting algorithm',
                                                      default=None)
    reconGroup.add_argument('--recon-with-dates',     help='Use isolate date information in ancestral joint sequence'
                                                      ' reconstruction',
                                                      default=False, action='store_true')
    reconGroup.add_argument('--model-fitter-args',    help='Further arguments passed to model fitting algorithm',
                                                      default=None)
    reconGroup.add_argument('--mar',                  help='Use marginal, rather than joint, ancestral reconstruction',
                                                      action='store_true')
    reconGroup.add_argument('--seq-recon',            help='Algorithm to use for marginal reconstruction [if unspecified: '
                                                      'same as tree builder if possible, else iqtree; requires --mar flag]',
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
    gubbinsGroup.add_argument('--p-value',
                                                      help='Uncorrected p value used to identify recombinations',
                                                      type=float,
                                                      default=0.05)
    gubbinsGroup.add_argument('--trimming-ratio',
                                                      help='Ratio of log probabilities used to trim recombinations',
                                                      type=float,
                                                      default=1.0)
    gubbinsGroup.add_argument('--extensive-search',
                                                      help='Undertake slower, more thorough, search for recombination',
                                                      action='store_true',
                                                      default=False)

    stopGroup = parser.add_argument_group('Algorithm start/stop options')
    stopGroup.add_argument('--iterations',      '-i',
                                                        help='Maximum No. of iterations',
                                                        type=int,
                                                        default=5)
    stopGroup.add_argument('--converge-method', '-z',
                                                        help='Criteria to use to know when to halt iterations',
                                                        default='weighted_robinson_foulds',
                                                        choices=['weighted_robinson_foulds', 'robinson_foulds',
                                                                     'recombination'])
    stopGroup.add_argument('--resume',
                                                        help='Intermediate tree from previous run (must include'
                                                        ' "iteration_X" in file name)',
                                                        default=None)
    return parser


def main():
    parser = parse_input_args()
    gubbins.common.parse_and_run(parser.parse_args(), parser.description)

if __name__ == '__main__':
    main()
    sys.exit(0)
