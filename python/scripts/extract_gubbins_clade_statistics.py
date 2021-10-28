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
                                     prog='extract_gubbins_clade')

    # input options
    parser.add_argument('--clades',
                        help = 'Two column file assigning isolates (first column) to clades (second column)',
                        required = True)
    parser.add_argument('--gff',
                        help = 'recombination prediction GFF file output by Gubbins',
                        required = True)
    parser.add_argument('--snps',
                        help = 'branch base reconstruction EMBL file output by Gubbins',
                        required = True)
    parser.add_argument('--exclude-regions',
                        help = 'Two column file specifying start and end of regions to be excluded',
                        required = False)
    parser.add_argument('--tree',
                        help = 'Labelled tree output by Gubbins',
                        required = True)
    parser.add_argument('--out',
                        help = 'Output file prefix; suffix is "_clades.csv"',
                        required = True)

    return parser.parse_args()

# main code
if __name__ == "__main__":

    # Get command line options
    args = get_options()
    
    # Parse clades
    clades = {}
    clade_names = set()
    with open(args.clades,'r') as clade_list:
        for line in clade_list.readlines():
            info = line.strip().split()
            if len(info) == 2:
                clades[info[0]] = info[1]
                clade_names.add(info[1])
            else:
                sys.stderr.write('Line needs two columns: ' + line + '\n')
    
    # Store SNP information
    node_snps = {}
    snp_total = 0
    with open(args.snps,'r') as snp_file:
        pos = 0
        for line in snp_file.readlines():
            info = line.strip().split()
            if info[1] == 'variation':
                pos = int(info[2])
                snp_total += 1
            if info[1].startswith('/node='):
                node = info[1].replace('"','').split('->')
                if node[1] in node_snps:
                    node_snps[node[1]].append(pos)
                else:
                    node_snps[node[1]] = [pos]

    # Store recombination information
    node_rec_starts = {}
    node_rec_ends = {}
    with open(args.gff,'r') as gff_file:
        for line in gff_file.readlines():
            if not line.startswith('##'):
                info = line.rstrip().split('\t')
                start = int(info[3])
                end = int(info[4])
                node = info[8].split(';')[0].replace('"','').split('->')[1]
                if node not in node_rec_starts:
                    node_rec_starts[node] = [start]
                    node_rec_ends[node] = [end]
                else:
                    node_rec_starts[node].append(start)
                    node_rec_ends[node].append(end)
    
    # Divide SNPs into recombinant and non-recombinant
    rec_snps = {node:0 for node in node_snps}
    pm_snps = {node:0 for node in node_snps}
    for node in node_snps:
        for p in node_snps[node]:
            rec_snp = False
            if node in node_rec_starts:
                for s,e in zip(node_rec_starts[node],node_rec_ends[node]):
                    if p >= s and p <= e:
                        rec_snp = True
                        break
            if rec_snp:
                rec_snps[node] += 1
            else:
                pm_snps[node] += 1
    
    # Parse tree
    info_labels = ['total_snps','rec_snps','mutation_snps','recombinations']
    tree_info_labels = ['n_taxa','n_branches','branch_length']
    tree = dendropy.Tree.get(path = args.tree,
                             schema = 'newick',
                             preserve_underscores = True,
                             rooting='force-rooted')
    
    # Calculate statistics per clade
    with open(args.out + '_clades.csv','w') as out_file:
        out_file.write('Clade,')
        out_file.write(','.join(info_labels + tree_info_labels))
        out_file.write('\n')
        for clade_name in clade_names:
            out_file.write(clade_name + ',')
            clade_members = [sequence for sequence in clades if clades[sequence] == clade_name]
            clade_tree = tree.clone(depth = 1)
            clade_tree.retain_taxa_with_labels(clade_members)
            clade_info = {label:0 for label in info_labels + tree_info_labels}
            for node in clade_tree.preorder_node_iter():
                if node != clade_tree.seed_node:
                    clade_info['n_branches'] += 1
                    clade_info['branch_length'] += node.edge_length
                    if node.is_leaf():
                        clade_info['n_taxa'] += 1
                        node_label_string = node.taxon.label
                    else:
                        node_label_string = node.label
                    if node_label_string in node_snps:
                        clade_info['total_snps'] += len(node_snps[node_label_string])
                        clade_info['rec_snps'] += rec_snps[node_label_string]
                        clade_info['mutation_snps'] += pm_snps[node_label_string]
                    if node_label_string in node_rec_starts:
                        clade_info['recombinations'] += len(node_rec_starts[node_label_string])
                    
            out_file.write(','.join([str(clade_info[label]) for label in info_labels + tree_info_labels]))
            out_file.write('\n')
