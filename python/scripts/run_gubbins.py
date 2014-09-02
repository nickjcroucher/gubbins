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
sys.path.append(".")
import argparse
from gubbins import common

parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
parser.add_argument('alignment_filename',       help='Multifasta alignment file')
parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
parser.add_argument('--starting_tree',    '-s', help='Starting tree')
parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging')
parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 25)
parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')

gubbins_runner  = common.GubbinsCommon(parser.parse_args())
gubbins_runner.parse_and_run()
