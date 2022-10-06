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
# Phylogenetic imports
import dendropy
# Biopython imports
from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq

# command line parsing
def get_options():

    parser = argparse.ArgumentParser(description='Extract a clade from a Gubbins output',
                                     prog='extract_clade')

    # input options
    parser.add_argument('--list',
                        help = 'List of sequences to extract',
                        required = True)
    parser.add_argument('--aln',
                        help = 'Input alignment (FASTA format)',
                        required = True)
    parser.add_argument('--gff',
                        help = 'GFF of recombinant regions detected by Gubbins',
                        required = True)
    parser.add_argument('--tree',
                        help = 'Final tree generated by Gubbins',
                        required = True)
    parser.add_argument('--out',
                        help = 'Output file prefix',
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
    
    # Parse list of input sequences
    subset = set()
    # Read in FASTA assemblies
    with open(args.list,'r') as seq_list:
        for line in seq_list.readlines():
            subset.add(line.strip().split()[0])
    
    # Extract from alignment
    output_aln_name = args.out + '.aln'
    names_in_alignment = set()
    sequences_to_print = []
    alignment = AlignIO.read(args.aln,'fasta')
    for taxon in alignment:
        names_in_alignment.add(taxon.id)
        if taxon.id in subset:
            sequences_to_print.append(taxon)

    # Check subset sequences are found in alignment
    not_found_in_dataset = subset - names_in_alignment
    if len(not_found_in_dataset) > 0:
        sys.stderr.write('Sequences in subset missing from alignment: ' + \
                            str(not_found_in_dataset) + '\n')
        sys.exit(1)
    
    # Prune from the tree
    output_tree_name = args.out + '.tree'
    tree = dendropy.Tree.get(path = args.tree,
                             schema = 'newick',
                             preserve_underscores = True)
    tree.retain_taxa_with_labels(subset)
    tree.write_to_path(output_tree_name,
                        'newick')

    # Identify relevant recombination blocks
    recombination_starts = []
    recombination_ends = []
    recombination_taxa = []
    output_gff_name = args.out + '.gff'
    taxon_pattern = re.compile('taxa="([^"]*)"')
    with open(args.gff,'r') as in_gff, open(output_gff_name,'w') as out_gff:
        for line in in_gff.readlines():
            if line.startswith('##'):
                out_gff.write(line)
            else:
                info = line.rstrip().split('\t')
                taxon_set = set(taxon_pattern.search(info[8]).group(1).split())
                if not taxon_set.isdisjoint(subset):
                    out_gff.write(line)
                    if len(taxon_set) < len(subset):
                        recombination_starts.append(int(info[3]))
                        recombination_ends.append(int(info[4]))
                        recombination_taxa.append(frozenset(taxon_set))

    # Write out masked alignment for clade
    with open(output_aln_name,'w') as out_aln:
        for taxon in sequences_to_print:
            for index,taxa in enumerate(recombination_taxa):
                if taxon.id in taxa:
                    seq_string = str(taxon.seq)
                    replacement_string = args.missing_char*(recombination_ends[index] - recombination_starts[index] + 1)
                    seq_string = seq_string[:recombination_starts[index]-1] + replacement_string + seq_string[recombination_ends[index]:]
                    taxon.seq = Seq(seq_string)
            SeqIO.write(taxon, out_aln, args.out_fmt)
