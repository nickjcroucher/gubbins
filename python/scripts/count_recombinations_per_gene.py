#! python

# encoding: utf-8
# Wellcome Trust Sanger Institute and Imperial College London
# Copyright (C) 2023  Wellcome Trust Sanger Institute and Imperial College London
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

# command line parsing
def get_options():

    parser = argparse.ArgumentParser(description='Mask recombinant regions detected by'
                                                    ' Gubbins from the input alignment',
                                     prog='mask_gubbins_aln')

    # input options
    parser.add_argument('--rec-gff',
                        help = 'GFF of recombinant regions detected by Gubbins',
                        required = True)
    parser.add_argument('--anno-gff',
                        help = 'GFF of annotation corresponding to the input alignment',
                        required = True)
    parser.add_argument('--out',
                        help = 'Output file name',
                        required = True)
    return parser.parse_args()

# main code
if __name__ == "__main__":

    # Get command line options
    args = get_options()
    
    # Read recombinant regions from GFF
    rec_start = []
    rec_end = []
    rec_affected = []
    taxon_pattern = re.compile('taxa="([^"]*)"')
    with open(args.rec_gff,'r') as gff_file:
        for line in gff_file.readlines():
            if not line.startswith('##'):
                # Calculate stats
                info = line.rstrip().split('\t')
                start = int(info[3])
                end = int(info[4])
                taxon_set = set(taxon_pattern.search(info[8]).group(1).split())
                # Record stats
                rec_start.append(start)
                rec_end.append(end)
                rec_affected.append(len(taxon_set))

    # Read annotation from GFF
    cds_start = []
    cds_end = []
    cds_name = []
    cds_index = []
    contig_starts = {}
    cumulative_length = 0
    with open(args.anno_gff,'r') as gff_file:
        for line in gff_file.readlines():
            if line.startswith('##sequence-region'):
                info = line.rstrip().split(' ')
                contig_starts[info[1]] = int(cumulative_length)
                cumulative_length = contig_starts[info[1]] + int(info[3])
            elif not line.startswith('##'):
                info = line.rstrip().split('\t')
                if len(info) >= 7:
                  if info[2] == 'CDS':
                      name = '-'
                      index = None
                      cds_data = info[8].replace('"','').split(';')
                      for datum in cds_data:
                          qualifier, value = datum.split('=')
                          if qualifier == 'locus_tag':
                              index = value
                          elif qualifier == 'gene':
                              name = value
                      if index is not None:
                          contig_start = contig_starts[info[0]]
                          cds_start.append(int(info[3]) + contig_start)
                          cds_end.append(int(info[4]) + contig_start)
                          cds_index.append(index)
                          cds_name.append(name)

    # Run checks on whether genes have been detected
    if len(cds_index) == 0:
        sys.stderr.write('No genes detected in annotation\n')
        sys.exit()
    elif not (len(cds_start) == len(cds_end) and len(cds_start) == len(cds_index) and len(cds_start) == len(cds_name)):
        sys.stderr.write('Error with extraction of information on annotation\n')
        sys.exit()
    elif not (len(rec_start) == len(rec_end) and len(rec_start) == len(rec_affected)):
        sys.stderr.write('Error with extraction of information on recombination\n')
        sys.exit()

    # Write out summary statistics
    with open(args.out,'w') as out_file:
        out_file.write('CDS\tGeneName\tStart\tEnd\tNumRec\tNumAffectedTaxa\n')
        for cnum,cindex in enumerate(cds_index):
            num_rec = 0
            num_affected_taxa = 0
            cstart = cds_start[cnum]
            cend = cds_end[cnum]
            for rnum,rstart in enumerate(rec_start):
                rend = rec_end[rnum]
                if (rstart < cend and rend > cend) or (rstart < cstart and rend > cend):
                    num_rec = num_rec + 1
                    num_affected_taxa = num_affected_taxa + rec_affected[rnum]
            out_file.write(cds_index[cnum] + '\t' + cds_name[cnum] + '\t' + str(cds_start[cnum]) + '\t' + str(cds_end[cnum]) + '\t' + str(num_rec) + '\t' + str(num_affected_taxa) + '\n')
