#! python

# encoding: utf-8
# Wellcome Trust Sanger Institute and Imperial College London
# Copyright (C) 2020  Wellcome Trust Sanger Institute and Imperial College London
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

# Generic imports
import os
import sys
import argparse
import re
import math

# Biopython imports
from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq

# command line parsing
def get_options():

    parser = argparse.ArgumentParser(description='Extract all the unique alleles at recombinant loci')

    # input options
    parser.add_argument('--aln',
                        help = 'Input alignment (FASTA format)',
                        required = True)
    parser.add_argument('--gff',
                        help = 'GFF of recombinant regions detected by Gubbins',
                        required = True)
    parser.add_argument('--out-dir',
                        help = 'Output directory',
                        required = True)
    parser.add_argument('--start',
                        help = 'Start of region of interest',
                        default = 1,
                        required = False)
    parser.add_argument('--end',
                        help = 'End of region of interest',
                        default = math.inf,
                        required = False)
    parser.add_argument('--terminal-only',
                        help = 'Only extract recombinations on terminal branches',
                        default = False,
                        action = 'store_true')

    return parser.parse_args()

# main code
if __name__ == "__main__":

    # Get command line options
    args = get_options()

    # Create output directory
    if not os.path.isdir(args.out_dir):
      os.mkdir(args.out_dir)

    # Read recombinant regions from GFF
    rec_start = []
    rec_end = []
    with open(args.gff,'r') as gff_file:
      for line in gff_file.readlines():
        if not line.startswith('##'):
          # Get coordinates
          info = line.rstrip().split('\t')
          taxon_pattern = re.compile('taxa="([^"]*)"')
          taxon_set = set(taxon_pattern.search(info[8]).group(1).split())
          if (len(taxon_set) == 1 or not args.terminal_only):
            rec_start.append(int(info[3]))
            rec_end.append(int(info[4]))
    
    # Read in alignment and identify recombinations
    alignment = AlignIO.read(args.aln,'fasta')
    for (start,end) in zip(rec_start,rec_end):
      if start >= args.start and end <= args.end:
        out_fn = 'locus_' + str(start) + '_' + str(end) + '.aln'
        with open(os.path.join(args.out_dir,out_fn),'w') as out_file:
          seen_seqs = []
          rec_locus_alignment = alignment[:,(start-1):end]
          for taxon in rec_locus_alignment:
            if taxon.seq not in seen_seqs:
              out_file.write('>' + taxon.id + '\n' + str(taxon.seq) + '\n')
              seen_seqs.append(taxon.seq)
