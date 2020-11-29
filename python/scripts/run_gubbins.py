#!/usr/bin/env python3
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

import argparse
import pkg_resources
import gubbins.common


def main():

    parser = argparse.ArgumentParser(
        description='Croucher N. J., Page A. J., Connor T. R., Delaney A. J., Keane J. A., Bentley S. D., Parkhill J., '
                    'Harris S.R. "Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome '
                    'sequences using Gubbins". Nucleic Acids Res. 2015 Feb 18;43(3):e15. doi: 10.1093/nar/gku1196.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ioGroup = parser.add_argument_group('Input and output options')
    ioGroup.add_argument('alignment_filename',        help='Multifasta alignment file')
    ioGroup.add_argument('--prefix',            '-p', help='Add a prefix to the final output filenames')
    ioGroup.add_argument('--starting_tree',     '-s', help='Starting tree')
    ioGroup.add_argument('--use_time_stamp',    '-u', help='Use a time stamp in file names', action='store_true')
    ioGroup.add_argument('--version',                 action='version',
                        version=str(pkg_resources.get_distribution("gubbins").version))
                        
    dataGroup = parser.add_argument_group('Data processing options')
    dataGroup.add_argument('--pairwise',        '-2', help='Compare two sequences (without using a tree)',
                                                        default = False, action = 'store_true') # fasttree model fit, star phylogeny, one iteration
    dataGroup.add_argument('--filter_percentage',  '-f', help='Filter out taxa with more than this percentage of gaps',
                        type=int, default=25)
    dataGroup.add_argument('--remove_identical_sequences', '-d', help='Remove identical sequences', action='store_true')
    dataGroup.add_argument('--threads',            '-c', help='Number of threads to run with RAXML, but only if a PTHREADS '
                                                         'version is available', type=int,  default=1)
    dataGroup.add_argument('--verbose',            '-v', help='Turn on debugging', action='store_true')
    dataGroup.add_argument('--no_cleanup',         '-n', help="Don't cleanup intermediate files", action='store_true')

    treeGroup = parser.add_argument_group('Tree building options')
    treeGroup.add_argument('--first-tree-builder', '-1', help='Application to use for building the first tree',
                                                                default=None,
                                                                choices=['raxml', 'iqtree', 'fasttree', 'rapidnj', 'star'])
    treeGroup.add_argument('--tree-builder',       '-t', help='Application to use for tree building',
                                                            default="raxml",
                                                            choices=['raxml', 'iqtree', 'fasttree', 'hybrid', 'rapidnj'])
    treeGroup.add_argument('--outgroup',           '-o', help='Outgroup name for rerooting. A list of comma separated '
                                                          'names can be used if they form a clade')
                                                          
    modelGroup = parser.add_argument_group('Nucleotide substitution model options')
    modelGroup.add_argument('--model',             '-g', help='Nucleotide substitution model (GTRCAT not available for iqtree;'
                                                        ' GTR not available for raxml)',
                                                         default='GTRGAMMA',
                                                         choices=['GTR' ,'GTRGAMMA', 'GTRCAT'])
    modelGroup.add_argument('--model-fitter',      '-r', help='Application to use for model fitting [default = same as'
                                                         ' tree builder if possible, else raxml]',
                                                         default="raxml",
                                                         choices=['raxml', 'iqtree', 'fasttree'])
    modelGroup.add_argument('--mar',               '-M', help='Use marginal ancestral reconstruction', action='store_true')
    modelGroup.add_argument('--sequence-recon',    '-q', help='Application to use for marginal reconstruction [default = '
                                                            'same as tree builder if possible, else raxml]',
                                                            default="raxml",
                                                            choices=['raxml', 'iqtree'])
    
    gubbinsGroup = parser.add_argument_group('Recombination detection options')
    gubbinsGroup.add_argument('--min_snps',          '-m', help='Min SNPs to identify a recombination block', type=int,
                        default=3)
    gubbinsGroup.add_argument('--min_window_size',   '-a', help='Minimum window size', type=int, default=100)
    gubbinsGroup.add_argument('--max_window_size',   '-b', help='Maximum window size', type=int, default=10000)

    stopGroup = parser.add_argument_group('Algorithm stop options')
    stopGroup.add_argument('--iterations',        '-i', help='Maximum No. of iterations', type=int, default=5)
    stopGroup.add_argument('--converge_method',   '-z', help='Criteria to use to know when to halt iterations',
                        default='weighted_robinson_foulds', choices=['weighted_robinson_foulds', 'robinson_foulds',
                                                                     'recombination'])

    gubbins.common.parse_and_run(parser.parse_args(), parser.description)
