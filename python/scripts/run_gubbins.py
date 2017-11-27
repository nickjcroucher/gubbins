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

import sys
sys.path.append(".")
sys.path.append("..")
import argparse
import pkg_resources
from gubbins import common



parser = argparse.ArgumentParser(description='Croucher N. J., Page A. J., Connor T. R., Delaney A. J., Keane J. A., Bentley S. D., Parkhill J., Harris S.R. "Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins". Nucleic Acids Res. 2015 Feb 18;43(3):e15. doi: 10.1093/nar/gku1196 .', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('alignment_filename',       help='Multifasta alignment file')
parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting. A list of comma separated names can be used if they form a clade')
parser.add_argument('--starting_tree',    '-s', help='Starting tree')
parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names', default = 0)
parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging', default = 0)
parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files', default = 0)
parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building.', default = "raxml", choices=['raxml', 'fastree', 'hybrid'])
parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations', type=int,  default = 5)
parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block', type=int,  default = 3)
parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps', type=int,  default = 25)
parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')
parser.add_argument('--threads',          '-c', help='Number of threads to run with RAXML, but only if a PTHREADS version is available', type=int,  default = 1)
parser.add_argument('--converge_method',  '-z', help='Criteria to use to know when to halt iterations.',  default = 'weighted_robinson_foulds', choices=['weighted_robinson_foulds', 'robinson_foulds', 'recombination'])
parser.add_argument('--version',                action='version', version=str(pkg_resources.get_distribution("gubbins").version))
parser.add_argument('--min_window_size',  '-a', help='Minimum window size', type=int,  default = 100)
parser.add_argument('--max_window_size',  '-b', help='Maximum window size', type=int,  default = 10000)
parser.add_argument('--raxml_model',      '-r', help='RAxML model.',  default = 'GTRCAT', choices=['GTRGAMMA', 'GTRCAT'])
parser.add_argument('--remove_identical_sequences', '-d', action='count', help='Remove identical sequences', default = 0)

gubbins_runner  = common.GubbinsCommon(parser.parse_args())
gubbins_runner.parse_and_run()
