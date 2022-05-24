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
import sys
import argparse
import re
# Biopython imports
from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.Seq import MutableSeq

# command line parsing
def get_options():

    parser = argparse.ArgumentParser(description='Mask recombinant regions detected by'
                                                    ' Gubbins from the input alignment',
                                     prog='mask_gubbins_aln')

    # input options
    parser.add_argument('--aln',
                        help = 'Input alignment (FASTA format)',
                        required = True)
    parser.add_argument('--gff',
                        help = 'GFF of recombinant regions detected by Gubbins',
                        required = True)
    parser.add_argument('--out',
                        help = 'Output file name',
                        required = True)
    parser.add_argument('--out-fmt',
                        help = 'Format of output alignment',
                        default = 'fasta')
    parser.add_argument('--missing-char',
                        help = 'Character used to replace recombinant sequence',
                        default = '-')

    return parser.parse_args()

# main code
if __name__ == "__main__":

    # Get command line options
    args = get_options()
    
    # Read in alignment
    overall_taxon_list = []
    alignment = AlignIO.read(args.aln,'fasta')
    for taxon in alignment:
        overall_taxon_list.append(taxon.id)
        taxon.seq = MutableSeq(taxon.seq)
    overall_taxon_set = set(overall_taxon_list)
    
    # Read recombinant regions from GFF
    taxon_pattern = re.compile('taxa="([^"]*)"')
    with open(args.gff,'r') as gff_file:
        for line in gff_file.readlines():
            if not line.startswith('##'):
                info = line.rstrip().split('\t')
                start = int(info[3])-1
                end = int(info[4])-1
                taxon_set = set(taxon_pattern.search(info[8]).group(1).split())
                if taxon_set.issubset(overall_taxon_set):
                    for taxon in alignment:
                        if taxon.id in taxon_set:
                            taxon.seq[start:end+1] = args.missing_char*(end-start+1)

    # Write out masked sequence
    with open(args.out,'w') as out_aln:
        AlignIO.write(alignment, out_aln, args.out_fmt)
