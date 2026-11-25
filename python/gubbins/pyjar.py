#!/usr/bin/env python
# encoding: utf-8
#
# pyjar written by Simon Harris
# code modified from https://github.com/simonrharris/pyjar
# pyjar is free software, licensed under GPLv3.

import shutil
from scipy import linalg
import numpy
import dendropy
import sys
import os
import time
from Bio import AlignIO
from math import log, exp
from functools import partial
import numba
from numba import jit, njit, prange, types, from_dtype
from math import ceil
import collections
import datetime
import multiprocessing
try:
    from multiprocessing import Pool, shared_memory
    from multiprocessing.managers import SharedMemoryManager
    NumpyShared = collections.namedtuple('NumpyShared', ('name', 'shape', 'dtype'))
except ImportError as e:
    sys.stderr.write("This version of Gubbins requires the multiprocessing library and python v3.8 or higher for memory management\n")
    sys.exit(201)

from gubbins.utils import generate_shared_mem_array

###########################
# Python-native functions #
###########################

class Pyjar:
    def __init__(self, model: str):
        """Initialises the object"""
        self.citation = "https://doi.org/10.1093/oxfordjournals.molbev.a026369"
        self.process = "Sequence reconstructor"
        self.version = "1.0"
        self.algorithm = "pyjar"
        self.executable = "pyjar"
        self.model = model


class FlattenedBasePositions:
    """Holds CSR-flattened base pattern positions while supporting sequence indexing/iteration."""
    __slots__ = ('positions_flat', 'position_offsets')

    def __init__(self, positions_flat, position_offsets):
        self.positions_flat = positions_flat
        self.position_offsets = position_offsets

    def __len__(self):
        return len(self.position_offsets) - 1

    def __getitem__(self, idx):
        return self.positions_flat[self.position_offsets[idx]:self.position_offsets[idx + 1]]

    def __iter__(self):
        offsets = self.position_offsets
        flat = self.positions_flat
        for i in range(len(offsets) - 1):
            yield flat[offsets[i]:offsets[i + 1]]


# Split a list into chunks for multiprocessing
# from https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length/37414115#37414115
def chunks(l, k):
    n = len(l)
    return [l[i * (n // k) + min(i, n % k):(i+1) * (n // k) + min(i+1, n % k)] for i in range(k)]

#Read the tree file and root
def read_tree(treefile):
    if not os.path.isfile(treefile):
        print("Error: tree file does not exist")
        sys.exit(204)
    t=dendropy.Tree.get(path=treefile,
                        schema="newick",
                        preserve_underscores=True,
                        rooting="force-rooted")
    return t

#Calculate Pij from Q matrix and branch length
def calculate_pij(branch_length,rate_matrix):
    if branch_length==0:
        pij = numpy.full((4,4), -numpy.inf, dtype = numpy.float32)
        numpy.fill_diagonal(pij, 0.0)
    else:
        pij = numpy.array(numpy.log(linalg.expm(numpy.multiply(branch_length,rate_matrix))), dtype = numpy.float32) # modified
    return pij.flatten()

# Create an instanteous rate matrix
def create_rate_matrix(f, r):
    #convert f and r to Q matrix
    rm=numpy.array([[0, f[0]*r[1], f[0]*r[2], f[0]*r[3]],[f[1]*r[0], 0, f[1]*r[3],f[1]*r[4]],[f[2]*r[1], f[2]*r[3], 0, f[2]*r[5]],[f[3]*r[2], f[3]*r[4], f[3]*r[5], 0]])

    rm[0][0]=numpy.sum(rm[0])*-1
    rm[1][1]=numpy.sum(rm[1])*-1
    rm[2][2]=numpy.sum(rm[2])*-1
    rm[3][3]=numpy.sum(rm[3])*-1

    return rm

# Read the info file from the selected phylogenetic software
# to get rates and frequencies
def read_info(infofile, type = 'raxml'):

    if not os.path.isfile(infofile):
        print("Error: model information file " + infofile + " does not exist")
        sys.exit(205)
    
    if type not in ['raxml', 'raxmlng', 'iqtree','fasttree']:
        sys.stderr.write('Only able to parse GTR-type models from raxml, iqtree or fasttree')
        sys.exit(206)
    
    r=[-1.0] * 6 # initialise rates
    f=[-1.0] * 4 # initialise frequencies
    
    with open(infofile, "r") as info_file:
        for line in info_file:
            line=line.strip()
            if type == 'raxml':
                if "freq pi" in line:
                    words=line.split()
                    if "pi(A)" in line:
                        f[0] = float(words[2])
                    elif "pi(C)" in line:
                        f[1] = float(words[2])
                    elif "pi(G)" in line:
                        f[2] = float(words[2])
                    elif "pi(T)" in line:
                        f[3] = float(words[2])
                elif "Base frequencies:" in line:
                    words=line.split()
                    f=[float(words[2]), float(words[3]), float(words[4]), float(words[5])]
                elif "<->" in line:
                    # order is ac ag at cg ct gt
                    words=line.split()
                    if "A <-> C" in line:
                        r[0] = float(words[4])
                    elif "A <-> G" in line:
                        r[1] = float(words[4])
                    elif "A <-> T" in line:
                        r[2] = float(words[4])
                    elif "C <-> G" in line:
                        r[3] = float(words[4])
                    elif "C <-> T" in line:
                        r[4] = float(words[4])
                    elif "G <-> T" in line:
                        r[5] = float(words[4])
                elif "alpha[0]:" in line:
                    # order is ac ag at cg ct gt
                    words=line.split()
                    r=[float(words[9]), float(words[10]), float(words[11]), float(words[12]), float(words[13]), float(words[14])]
            elif type == 'raxmlng':
                sep_by_braces = line.replace('{','}').split('}')
                if sep_by_braces[0] == "GTR":
                    r = [float(rate) for rate in sep_by_braces[1].split('/')]
                    f = [float(rate) for rate in sep_by_braces[3].split('/')]
                elif sep_by_braces[0] == "K80":
                    sep_rates = [float(rate) for rate in sep_by_braces[1].split('/')]
                    r = [sep_rates[0], sep_rates[1], sep_rates[0], sep_rates[0], sep_rates[1], sep_rates[0]]
                    f = [0.25,0.25,0.25,0.25]
                elif sep_by_braces[0] == "HKY":
                    sep_rates = [float(rate) for rate in sep_by_braces[1].split('/')]
                    r = [sep_rates[0], sep_rates[1], sep_rates[0], sep_rates[0], sep_rates[1], sep_rates[0]]
                    f = [float(rate) for rate in sep_by_braces[3].split('/')]
                elif line.startswith("JC"):
                    f = [0.25,0.25,0.25,0.25]
                    r = [1.0,1.0,1.0,1.0,1.0,1.0]
            elif type == 'iqtree':
                # Base frequencies
                words=line.split()
                if "pi(A)" in line:
                    f[0] = float(words[2])
                elif "pi(C)" in line:
                    f[1] = float(words[2])
                elif "pi(G)" in line:
                    f[2] = float(words[2])
                elif "pi(T)" in line:
                    f[3] = float(words[2])
                # Rate estimates
                if "A-C:" in line:
                    r[0] = float(words[1])
                elif "A-G:" in line:
                    r[1] = float(words[1])
                elif "A-T:" in line:
                    r[2] = float(words[1])
                elif "C-G:" in line:
                    r[3] = float(words[1])
                elif "C-T:" in line:
                    r[4] = float(words[1])
                elif "G-T:" in line:
                    r[5] = float(words[1])
            
                if line.startswith('Base frequencies:'):
                    words=line.split()
                    f=[float(words[3]), float(words[5]), float(words[7]), float(words[9])]
                elif line.startswith('Rate parameters:'):
                    words=line.split()
                    # order is ac ag at cg ct gt
                    r=[float(words[3]), float(words[5]), float(words[7]), float(words[9]), float(words[11]), float(words[13])]
            elif type == 'fasttree':
                if line.startswith('GTRFreq'):
                    words=line.split()
                    f=[float(words[1]), float(words[2]), float(words[3]), float(words[4])]
                elif line.startswith('GTRRates'):
                    words=line.split()
                    r=[float(words[1]), float(words[2]), float(words[3]), float(words[4]), float(words[5]), float(words[6])]
                elif 'ML Model: Jukes-Cantor' in line: # Jukes-Cantor model for Fasttree
                    f = [0.25,0.25,0.25,0.25]
                    r = [1.0,1.0,1.0,1.0,1.0,1.0]

    # Check frequencies and rates have been extracted correctly
    if -1.0 in f or -1.0 in r:
        sys.stderr.write('Problem with extracting model parameters - frequencies are ' + str(f) + ' and rates are ' + str(r))
        sys.exit(207)

    return numpy.array(f, dtype = numpy.float32), numpy.array(r, dtype = numpy.float32)

# Convert arrays of variable-length arrays to a square matrix
# for compatibility with numba
# from https://stackoverflow.com/questions/32037893/numpy-fix-array-with-rows-of-different-lengths-by-filling-the-empty-elements-wi
def convert_to_square_numpy_array(data):
    # Get lengths of each row of data
    lens = numpy.array([len(i) for i in data])

    # Mask of valid places in each row
    mask = numpy.arange(lens.max()) < lens[:,None]

    # Setup output array and put elements from data into masked positions
    out = numpy.full(mask.shape, -1, dtype = numpy.int32)
    out[mask] = numpy.concatenate(data)
    return out

# Get the unique base patterns within the numpy array
# Based on https://stackoverflow.com/questions/21888406/getting-the-indexes-to-the-duplicate-columns-of-a-numpy-array
def get_unique_columns(data):
    dt = numpy.dtype((numpy.void, data.dtype.itemsize * data.shape[0]))
    dataf = numpy.asfortranarray(data).view(dt)
    u,uind = numpy.unique(dataf, return_inverse=True)
    u = u.view(data.dtype).reshape(-1,data.shape[0]).T
    return (u,uind)

##########################
# JIT-compiled functions #
##########################

# Fast byte-buffer parsers for base_patterns.csv and base_positions.csv
#######################################################################
@njit(numba.uint8[:, :](numba.uint8[:], numba.int64, numba.int64), cache=True)
def parse_patterns_bytes(buf, n_rows, n_cols):
    out = numpy.empty((n_rows, n_cols), dtype=numpy.uint8)
    lut = numpy.full(256, 5, dtype=numpy.uint8)
    lut[65] = 0; lut[97] = 0   # A, a
    lut[67] = 1; lut[99] = 1   # C, c
    lut[71] = 2; lut[103] = 2  # G, g
    lut[84] = 3; lut[116] = 3  # T, t
    lut[45] = 4                # -
    lut[78] = 5; lut[110] = 5  # N, n
    row = 0
    col = 0
    for i in range(buf.size):
        b = buf[i]
        if b == 10:
            if col > 0:
                row += 1
                col = 0
        elif b != 13:
            out[row, col] = lut[b]
            col += 1
    return out


@njit(numba.types.Tuple((numba.int64, numba.int64))(numba.uint8[:]), cache=True)
def count_positions_and_lines(buf):
    n_lines = 0
    n_nums = 0
    in_num = False
    for i in range(buf.size):
        b = buf[i]
        if 48 <= b <= 57:
            if not in_num:
                in_num = True
                n_nums += 1
        elif b == 10:
            if in_num:
                n_lines += 1
                in_num = False
        else:
            in_num = False
    if in_num:
        n_lines += 1
    return n_lines, n_nums


@njit(numba.types.Tuple((numba.int32[:], numba.int32[:], numba.int64))(numba.uint8[:], numba.int64, numba.int64), cache=True)
def parse_positions_bytes(buf, n_lines, n_nums):
    positions_flat = numpy.empty(n_nums, dtype=numpy.int32)
    offsets = numpy.empty(n_lines + 1, dtype=numpy.int32)
    offsets[0] = 0
    line_idx = 0
    num_idx = 0
    val = numpy.int32(0)
    in_num = False
    max_val = numpy.int32(0)
    for i in range(buf.size):
        b = buf[i]
        if 48 <= b <= 57:
            val = val * numpy.int32(10) + numpy.int32(b - 48)
            in_num = True
        elif b == 44:
            if in_num:
                positions_flat[num_idx] = val
                if val > max_val:
                    max_val = val
                num_idx += 1
                val = numpy.int32(0)
                in_num = False
        elif b == 10:
            if in_num:
                positions_flat[num_idx] = val
                if val > max_val:
                    max_val = val
                num_idx += 1
                val = numpy.int32(0)
                in_num = False
                line_idx += 1
                offsets[line_idx] = num_idx
    if in_num:
        positions_flat[num_idx] = val
        if val > max_val:
            max_val = val
        num_idx += 1
        line_idx += 1
        offsets[line_idx] = num_idx
    return positions_flat, offsets, int(max_val) + 1


# Convert bases to integers
###########################
@njit(numba.void(numba.typeof(numpy.dtype('U1'))[:],
                numba.uint8[:]),
                cache = True)
def seq_to_int(seq,out_seq):
    for i,b in enumerate(seq):
        if b == 'A':
            out_seq[i] = 0
        elif b == 'C':
            out_seq[i] = 1
        elif b == 'G':
            out_seq[i] = 2
        elif b == 'T':
            out_seq[i] = 3
        elif b == '-':
            out_seq[i] = 4
        elif b == 'N':
            out_seq[i] = 5
        else:
            print('Unable to process character',b)

# Calculate most likely base given bases in descendents
#######################################################
@njit(numba.void(numba.float32[:,:],
                numba.uint8[:,:],
                numba.float32[:,::1],
                numba.int32,
                numba.int32[:,:],
                numba.uint8[::1]),
                cache=True)
def find_most_likely_base_given_descendents(Lmat, Cmat, pij, node_index, child_nodes, column_base_indices):
    #2a. Lz(i) = maxj Pij(tz) x Lx(j) x Ly(j)
    #2b. Cz(i) = the value of j attaining the above maximum.
    max_children = child_nodes.shape[1]
    n_bases = column_base_indices.size
    for b_end in range(n_bases):
        end_index = column_base_indices[b_end]
        c = numpy.float32(0.0)
        for k in range(max_children):
            ch = child_nodes[node_index, k]
            if ch < 0:
                break
            c += Lmat[ch, end_index]
        for b_start in range(n_bases):
            start_index = column_base_indices[b_start]
            j = pij[start_index, end_index] + c
            if j > Lmat[node_index, start_index]:
                Lmat[node_index, start_index] = j
                Cmat[node_index, start_index] = end_index

# Calculate the most likely base at the root node
#################################################
@njit(numba.void(numba.float32[:,:],
                numba.uint8[:,:],
                numba.float32[:],
                numba.int32,
                numba.int32[:],
                numba.uint8[:]),
                cache=True)
def calculate_root_likelihood(Lmat, Cmat, base_frequencies, node_index, child_node_indices, column_base_indices):
    n_children = child_node_indices.size
    n_bases = column_base_indices.size
    for b_end in range(n_bases):
        end_index = column_base_indices[b_end]
        c = numpy.float32(0.0)
        for k in range(n_children):
            ch = child_node_indices[k]
            if ch < 0:
                break
            c += Lmat[ch, end_index]
        j = log(base_frequencies[end_index]) + c
        for b_start in range(n_bases):
            start_index = column_base_indices[b_start]
            if j > Lmat[node_index, start_index]:
                Lmat[node_index, start_index] = j
                Cmat[node_index, start_index] = end_index

# Fill in matrices given known or unknown base in sequence
##########################################################
@njit(numba.void(numba.float32[:,:],
                numba.uint8[:,:],
                numba.float32[:,:],
                numba.int32,
                numba.uint8),
                cache=True)
def process_leaf(Lmat, Cmat, pij, node_index, taxon_base_index):
    if taxon_base_index < 4:
        #1a. Let j be the amino acid at y. Set, for each amino acid i: Cy(i)= j.
        #1b. Set for each amino acid i: Ly(i) = Pij(ty), where ty is the branch length between y and its father.
        for i in range(4):
            Cmat[node_index, i] = taxon_base_index
            Lmat[node_index, i] = pij[i, taxon_base_index]
    else:
        # Cmat stays as default when base is unknown
        for i in range(4):
            Lmat[node_index, i] = pij[i, i]

# Count the number of substitutions occurring on a branch
#########################################################
@njit(numba.void(numba.int32[:],
                numba.int32[:],
                numba.int32[:],
                numba.int32,
                numba.uint8[:],
                numba.int32[:]),
                cache=True)
def count_node_snps(node_snps,preordered_nodes,parent_nodes,seed_node,reconstructed_alleles,base_pattern_columns):
    # Note that preordered node list does not include the root
    n_cols = base_pattern_columns.size
    for i in range(preordered_nodes.size):
        node_index = preordered_nodes[i]
        parent_node_index = parent_nodes[node_index]
        a_node = reconstructed_alleles[node_index]
        a_parent = reconstructed_alleles[parent_node_index]
        if a_node < 4 and a_parent < 4 and a_node != a_parent:
            node_snps[node_index] += n_cols

# Reconstruct missing data at internal nodes
############################################
@njit(numba.void(numba.uint8[:],
                numba.int32[:],
                numba.int32[:],
                numba.int32[:],
                numba.uint8[:],
                numba.int32[:,:],
                numba.uint8[:]),
                cache=True)
def reconstruct_alleles(reconstructed_alleles,
                        postordered_nodes,
                        leaf_nodes,
                        node_index_to_aln_row,
                        column,
                        child_nodes,
                        reconstructed_base_indices):
    max_children = child_nodes.shape[1]
    for i in range(postordered_nodes.size):
        node_index = postordered_nodes[i]
        alignment_index = node_index_to_aln_row[node_index]
        if alignment_index >= 0:
            reconstructed_alleles[node_index] = column[alignment_index]
        else:
            has_child_base = False
            for k in range(max_children):
                child_taxon_index = child_nodes[node_index, k]
                if child_taxon_index < 0:
                    break
                if reconstructed_alleles[child_taxon_index] < 4:
                    has_child_base = True
                    break
            if has_child_base:
                reconstructed_alleles[node_index] = reconstructed_base_indices[node_index]
            else:
                reconstructed_alleles[node_index] = numpy.uint8(4)

# Transfer reconstructed alleles into alignment
###############################################
@njit(numba.void(numba.typeof(numpy.dtype('i1'))[:,:],
                numba.uint8[:],
                numba.typeof(numpy.dtype('i1'))[:],
                numba.int32[:],
                numba.int32[:]),
                cache=True)
def fill_out_aln(out_aln,reconstructed_alleles,ordered_bases,ancestral_node_order,base_pattern_columns):
    n_anc = ancestral_node_order.size
    n_cols = base_pattern_columns.size
    for index in range(n_anc):
        node_index = ancestral_node_order[index]
        base_index = reconstructed_alleles[node_index]
        base = ordered_bases[base_index]
        for c_i in range(n_cols):
            out_aln[base_pattern_columns[c_i], index] = base

# Reconstruct each base pattern
###############################
@njit(numba.void(numba.uint8[:],
                numba.int32[:],
                numba.float32[:,:],
                numba.uint8[:,:],
                numba.typeof(numpy.dtype('i1'))[:,:],
                numba.typeof(numpy.dtype('i1'))[:],
                numba.int32[:],
                numba.int32[:],
                numba.int32[:],
                numba.int32[:,:],
                numba.int32,
                numba.int32[:],
                numba.int32[:],
                numba.float32[:,:],
                numba.float32[:],
                numba.int32[:],
                numba.uint8[:],
                numba.int32[:]),
                cache=True)
def iterate_over_base_patterns(column,
                                base_pattern_columns,
                                Lmat,
                                Cmat,
                                tmp_out_aln,
                                ordered_bases,
                                postordered_nodes,
                                preordered_nodes,
                                parent_nodes,
                                child_nodes,
                                seed_node,
                                leaf_nodes,
                                ancestral_node_order,
                                node_pij,
                                base_frequencies,
                                node_index_to_aln_row,
                                reconstructed_base_indices,
                                node_snps):

    num_nodes = postordered_nodes.size
    for n in range(num_nodes):
        for b in range(4):
            Lmat[n, b] = -numpy.inf
            Cmat[n, b] = b

    # Count unknown bases and find unique observed bases in {0, 1, 2, 3}
    unknown_base_count = 0
    seen0 = False; seen1 = False; seen2 = False; seen3 = False
    for i in range(column.size):
        val = column[i]
        if val > 3:
            unknown_base_count += 1
        elif val == 0:
            seen0 = True
        elif val == 1:
            seen1 = True
        elif val == 2:
            seen2 = True
        elif val == 3:
            seen3 = True
    n_bases = (1 if seen0 else 0) + (1 if seen1 else 0) + (1 if seen2 else 0) + (1 if seen3 else 0)
    column_base_indices = numpy.empty(n_bases, dtype=numpy.uint8)
    idx = 0
    if seen0: column_base_indices[idx] = 0; idx += 1
    if seen1: column_base_indices[idx] = 1; idx += 1
    if seen2: column_base_indices[idx] = 2; idx += 1
    if seen3: column_base_indices[idx] = 3; idx += 1
      
    # Heuristic for speed: if all taxa are monomorphic, with a gap in only one sequence, then the ancestral states
    # will all be the observed base, as no ancestral node will have two child nodes with unknown bases at this site
    if unknown_base_count == 1 and n_bases == 1:
        # If site is monomorphic - replace entire column
        base_val = ordered_bases[column_base_indices[0]]
        n_anc = tmp_out_aln.shape[1]
        for c_i in range(base_pattern_columns.size):
            col_pos = base_pattern_columns[c_i]
            for a_i in range(n_anc):
                tmp_out_aln[col_pos, a_i] = base_val
    else:
        for i in range(num_nodes):
            node_index = postordered_nodes[i]
            if node_index == seed_node:
                continue
            pij = numpy.reshape(numpy.ascontiguousarray(node_pij[node_index, :]), (4, 4))
            alignment_index = node_index_to_aln_row[node_index]
            if alignment_index >= 0:
                taxon_base_index = column[alignment_index]
                process_leaf(Lmat,
                                Cmat,
                                pij,
                                node_index,
                                taxon_base_index)
            else:
                find_most_likely_base_given_descendents(Lmat,
                                                        Cmat,
                                                        pij,
                                                        node_index,
                                                        child_nodes,
                                                        column_base_indices)

        # Calculate likelihood of base at root node
        root_index = postordered_nodes[num_nodes - 1]
        child_node_indices = child_nodes[root_index, :]
        calculate_root_likelihood(Lmat, Cmat, base_frequencies, root_index, child_node_indices, column_base_indices)
        max_root_base_index = Cmat[root_index, numpy.argmax(Lmat[root_index, :])]
        reconstructed_base_indices[root_index] = max_root_base_index
        
        for i in range(preordered_nodes.size):
            node_index = preordered_nodes[i]
            parent_node_index = parent_nodes[node_index]
            reconstructed_base_indices[node_index] = Cmat[node_index, reconstructed_base_indices[parent_node_index]]

        reconstructed_alleles = numpy.full(num_nodes, 8, dtype = numpy.uint8)
        reconstruct_alleles(reconstructed_alleles,
                            postordered_nodes,
                            leaf_nodes,
                            node_index_to_aln_row,
                            column,
                            child_nodes,
                            reconstructed_base_indices
                            )

        fill_out_aln(tmp_out_aln,
                    reconstructed_alleles,
                    ordered_bases,
                    ancestral_node_order,
                    base_pattern_columns
                    )

        count_node_snps(node_snps,
                        preordered_nodes,
                        parent_nodes,
                        seed_node,
                        reconstructed_alleles,
                        base_pattern_columns,
                        )


# Allocation-free batch kernel operating directly on transposed ASCII buffer (n_ancestral, max_pos)
###################################################################################################
@njit(cache=True)
def _reconstruct_single_pattern_T(column,
                                  base_pattern_columns,
                                  Lmat,
                                  Cmat,
                                  out_aln_T,
                                  ordered_ascii_bases,
                                  postordered_nodes,
                                  preordered_nodes,
                                  parent_nodes,
                                  child_nodes,
                                  seed_node,
                                  ancestral_node_order,
                                  node_pij_3d,
                                  log_base_frequencies,
                                  node_index_to_aln_row,
                                  reconstructed_base_indices,
                                  reconstructed_alleles,
                                  column_base_indices_buf,
                                  node_snps):
    num_nodes = postordered_nodes.size
    max_children = child_nodes.shape[1]
    for n in range(num_nodes):
        for b in range(4):
            Lmat[n, b] = -numpy.inf
            Cmat[n, b] = b

    unknown_base_count = 0
    seen0 = False; seen1 = False; seen2 = False; seen3 = False
    for i in range(column.size):
        val = column[i]
        if val > 3:
            unknown_base_count += 1
        elif val == 0:
            seen0 = True
        elif val == 1:
            seen1 = True
        elif val == 2:
            seen2 = True
        elif val == 3:
            seen3 = True
    n_bases = 0
    if seen0: column_base_indices_buf[n_bases] = 0; n_bases += 1
    if seen1: column_base_indices_buf[n_bases] = 1; n_bases += 1
    if seen2: column_base_indices_buf[n_bases] = 2; n_bases += 1
    if seen3: column_base_indices_buf[n_bases] = 3; n_bases += 1

    n_cols = base_pattern_columns.size
    n_anc = ancestral_node_order.size

    if unknown_base_count == 1 and n_bases == 1:
        base_val = ordered_ascii_bases[column_base_indices_buf[0]]
        for a_i in range(n_anc):
            row = out_aln_T[a_i]
            for c_i in range(n_cols):
                row[base_pattern_columns[c_i]] = base_val
    else:
        for i in range(num_nodes):
            node_index = postordered_nodes[i]
            if node_index == seed_node:
                continue
            pij_node = node_pij_3d[node_index]
            aln_row = node_index_to_aln_row[node_index]
            if aln_row >= 0:
                process_leaf(Lmat, Cmat, pij_node, node_index, column[aln_row])
            else:
                for b_end in range(n_bases):
                    end_index = column_base_indices_buf[b_end]
                    c = numpy.float32(0.0)
                    for k in range(max_children):
                        ch = child_nodes[node_index, k]
                        if ch < 0:
                            break
                        c += Lmat[ch, end_index]
                    for b_start in range(n_bases):
                        start_index = column_base_indices_buf[b_start]
                        j = pij_node[start_index, end_index] + c
                        if j > Lmat[node_index, start_index]:
                            Lmat[node_index, start_index] = j
                            Cmat[node_index, start_index] = end_index

        root_idx = postordered_nodes[num_nodes - 1]
        for b_end in range(n_bases):
            end_index = column_base_indices_buf[b_end]
            c = numpy.float32(0.0)
            for k in range(max_children):
                ch = child_nodes[root_idx, k]
                if ch < 0:
                    break
                c += Lmat[ch, end_index]
            j = log_base_frequencies[end_index] + c
            for b_start in range(n_bases):
                start_index = column_base_indices_buf[b_start]
                if j > Lmat[root_idx, start_index]:
                    Lmat[root_idx, start_index] = j
                    Cmat[root_idx, start_index] = end_index

        best_b = 0
        best_l = Lmat[root_idx, 0]
        for b in range(1, 4):
            if Lmat[root_idx, b] > best_l:
                best_l = Lmat[root_idx, b]
                best_b = b
        reconstructed_base_indices[root_idx] = Cmat[root_idx, best_b]

        for i in range(preordered_nodes.size):
            node_index = preordered_nodes[i]
            parent_node_index = parent_nodes[node_index]
            reconstructed_base_indices[node_index] = Cmat[node_index, reconstructed_base_indices[parent_node_index]]

        for i in range(num_nodes):
            node_index = postordered_nodes[i]
            aln_row = node_index_to_aln_row[node_index]
            if aln_row >= 0:
                reconstructed_alleles[node_index] = column[aln_row]
            else:
                has_child_base = False
                for k in range(max_children):
                    ch = child_nodes[node_index, k]
                    if ch < 0:
                        break
                    if reconstructed_alleles[ch] < 4:
                        has_child_base = True
                        break
                if has_child_base:
                    reconstructed_alleles[node_index] = reconstructed_base_indices[node_index]
                else:
                    reconstructed_alleles[node_index] = numpy.uint8(4)

        for index in range(n_anc):
            node_index = ancestral_node_order[index]
            base = ordered_ascii_bases[reconstructed_alleles[node_index]]
            row = out_aln_T[index]
            for c_i in range(n_cols):
                row[base_pattern_columns[c_i]] = base

        for i in range(preordered_nodes.size):
            node_index = preordered_nodes[i]
            parent_node_index = parent_nodes[node_index]
            a_node = reconstructed_alleles[node_index]
            a_parent = reconstructed_alleles[parent_node_index]
            if a_node < 4 and a_parent < 4 and a_node != a_parent:
                node_snps[node_index] += n_cols


@njit(cache=True)
def reconstruct_all_patterns_T(base_patterns,
                               positions_flat,
                               position_offsets,
                               out_aln_T,
                               ordered_ascii_bases,
                               postordered_nodes,
                               preordered_nodes,
                               parent_nodes,
                               child_nodes,
                               seed_node,
                               ancestral_node_order,
                               node_pij_3d,
                               log_base_frequencies,
                               node_index_to_aln_row):
    n_patterns = base_patterns.shape[0]
    num_nodes = postordered_nodes.size
    Lmat = numpy.full((num_nodes, 4), -numpy.inf, dtype=numpy.float32)
    Cmat = numpy.zeros((num_nodes, 4), dtype=numpy.uint8)
    reconstructed_base_indices = numpy.full(num_nodes, 8, dtype=numpy.uint8)
    reconstructed_alleles = numpy.full(num_nodes, 8, dtype=numpy.uint8)
    column_base_indices_buf = numpy.zeros(4, dtype=numpy.uint8)
    node_snps = numpy.zeros(num_nodes, dtype=numpy.int32)
    for b in range(n_patterns):
        start = position_offsets[b]
        end = position_offsets[b + 1]
        _reconstruct_single_pattern_T(
            base_patterns[b],
            positions_flat[start:end],
            Lmat,
            Cmat,
            out_aln_T,
            ordered_ascii_bases,
            postordered_nodes,
            preordered_nodes,
            parent_nodes,
            child_nodes,
            seed_node,
            ancestral_node_order,
            node_pij_3d,
            log_base_frequencies,
            node_index_to_aln_row,
            reconstructed_base_indices,
            reconstructed_alleles,
            column_base_indices_buf,
            node_snps
        )
    return node_snps


@njit(parallel=True, cache=True)
def reconstruct_all_patterns_T_parallel(base_patterns,
                                        positions_flat,
                                        position_offsets,
                                        out_aln_T,
                                        ordered_ascii_bases,
                                        postordered_nodes,
                                        preordered_nodes,
                                        parent_nodes,
                                        child_nodes,
                                        seed_node,
                                        ancestral_node_order,
                                        node_pij_3d,
                                        log_base_frequencies,
                                        node_index_to_aln_row,
                                        n_threads):
    n_patterns = base_patterns.shape[0]
    num_nodes = postordered_nodes.size
    n_chunks = min(n_threads, n_patterns)
    if n_chunks < 1:
        n_chunks = 1
    node_snps_chunks = numpy.zeros((n_chunks, num_nodes), dtype=numpy.int32)
    for t in prange(n_chunks):
        chunk_start = (t * n_patterns) // n_chunks
        chunk_end = ((t + 1) * n_patterns) // n_chunks
        Lmat = numpy.full((num_nodes, 4), -numpy.inf, dtype=numpy.float32)
        Cmat = numpy.zeros((num_nodes, 4), dtype=numpy.uint8)
        reconstructed_base_indices = numpy.full(num_nodes, 8, dtype=numpy.uint8)
        reconstructed_alleles = numpy.full(num_nodes, 8, dtype=numpy.uint8)
        column_base_indices_buf = numpy.zeros(4, dtype=numpy.uint8)
        local_snps = node_snps_chunks[t]
        for b in range(chunk_start, chunk_end):
            start = position_offsets[b]
            end = position_offsets[b + 1]
            _reconstruct_single_pattern_T(
                base_patterns[b],
                positions_flat[start:end],
                Lmat,
                Cmat,
                out_aln_T,
                ordered_ascii_bases,
                postordered_nodes,
                preordered_nodes,
                parent_nodes,
                child_nodes,
                seed_node,
                ancestral_node_order,
                node_pij_3d,
                log_base_frequencies,
                node_index_to_aln_row,
                reconstructed_base_indices,
                reconstructed_alleles,
                column_base_indices_buf,
                local_snps
            )
    node_snps = numpy.zeros(num_nodes, dtype=numpy.int32)
    for t in range(n_chunks):
        for n in range(num_nodes):
            node_snps[n] += node_snps_chunks[t, n]
    return node_snps


# Convert integers to bases
###########################

@njit(numba.void(numba.int8[:],
                 numba.typeof(numpy.dtype('U1'))[:]),
                cache = True)
def int_to_seq(seq,out_seq):

    for i,b in enumerate(seq):
        if b == 0:
            out_seq[i] = 'A'
        elif b == 1:
            out_seq[i] = 'C'
        elif b == 2:
            out_seq[i] = 'G'
        elif b == 3:
            out_seq[i] = 'T'
        elif b == 4:
            out_seq[i] = '-'
        elif b == 5:
            out_seq[i] = 'N'
        else:
            print('Unable to process integer')

##################
# Main functions #
##################

####################################################
# Function for converting alignment to numpy array #
####################################################

def get_base_patterns(prefix, verbose, threads = 1):
    
    # Identify unique base patterns
    if verbose:
        print("Finding unique base patterns")
    t1=time.process_time()

    # Read in ordered sequence names
    sequence_names_fn = prefix + '.gaps.sequence_names.csv'
    sequence_names = []
    if not os.path.isfile(sequence_names_fn):
        sys.exit("Unable to open sequence names file " + sequence_names_fn + "\n")
    with open(sequence_names_fn,'r') as names_file:
        for line in names_file:
            sequence_names.append(line.rstrip())

    # Read in base patterns and base positions using compiled byte-buffer parsers
    base_patterns_fn = prefix + '.gaps.base_patterns.csv'
    if not os.path.isfile(base_patterns_fn):
        sys.exit("Unable to open base patterns file " + base_patterns_fn + "\n")
    base_positions_fn = prefix + '.gaps.base_positions.csv'
    if not os.path.isfile(base_positions_fn):
        sys.exit("Unable to open base positions file " + base_positions_fn + "\n")

    with open(base_positions_fn, 'rb') as positions_file:
        pos_bytes = numpy.fromfile(positions_file, dtype=numpy.uint8)
    n_lines, n_nums = count_positions_and_lines(pos_bytes)
    positions_flat, position_offsets, array_max = parse_positions_bytes(pos_bytes, n_lines, n_nums)

    with open(base_patterns_fn, 'rb') as pattern_file:
        pat_bytes = numpy.fromfile(pattern_file, dtype=numpy.uint8)
    vstacked_patterns = parse_patterns_bytes(pat_bytes, n_lines, len(sequence_names))
    array_of_position_arrays = FlattenedBasePositions(positions_flat, position_offsets)

    # Record timing
    t2=time.process_time()
    if verbose:
        print("Time taken to load unique base patterns:", t2-t1, "seconds")
        print("Unique base patterns: ", array_max)
    
    # Return output
    return sequence_names,vstacked_patterns,array_of_position_arrays,array_max

########################################################
# Function for reconstructing individual base patterns #
########################################################

def reconstruct_alignment_column(column_index,
                                base_pattern_positions,
                                tree = None,
                                preordered_nodes = None,
                                postordered_nodes = None,
                                leaf_nodes = None,
                                parent_nodes = None,
                                child_nodes = None,
                                seed_node = None,
                                node_pij = None,
                                node_index_to_aln_row = None,
                                ancestral_node_order = None,
                                base_patterns = None,
                                base_frequencies = None,
                                new_aln = None,
                                threads = 1,
                                verbose = False):
    
    # Load shared memory output alignment
    out_aln_shm = shared_memory.SharedMemory(name = new_aln.name)
    out_aln = numpy.ndarray(new_aln.shape, dtype = numpy.int8, buffer = out_aln_shm.buf)
    
    # Generate data structures for reconstructions
    num_nodes = postordered_nodes.size if postordered_nodes is not None else len(tree.nodes())
    Lmat = numpy.full((num_nodes,4), -numpy.inf, dtype = numpy.float32)
    Cmat = numpy.full((num_nodes,4), [0,1,2,3], dtype = numpy.uint8)
    reconstructed_base_indices = numpy.full(num_nodes, 8, dtype = numpy.uint8)

    ordered_bases = numpy.array([0,1,2,3,4], dtype=numpy.int8)
    node_snps = numpy.zeros(num_nodes, dtype = numpy.int32)

    # Load base pattern information
    base_patterns_shm = shared_memory.SharedMemory(name = base_patterns.name)
    base_patterns = numpy.ndarray(base_patterns.shape, dtype = base_patterns.dtype, buffer = base_patterns_shm.buf)
    
    column = base_patterns[column_index]

    iterate_over_base_patterns(column,
                                base_pattern_positions,
                                Lmat,
                                Cmat,
                                out_aln,
                                ordered_bases,
                                postordered_nodes,
                                preordered_nodes,
                                parent_nodes,
                                child_nodes,
                                seed_node,
                                leaf_nodes,
                                ancestral_node_order,
                                node_pij,
                                base_frequencies,
                                node_index_to_aln_row,
                                reconstructed_base_indices,
                                node_snps)

    # Close shared memory
    out_aln_shm.close()
    base_patterns_shm.close()

    return node_snps

##################################################
# Function for reconstructing complete alignment #
##################################################

def jar(sequence_names = None,
        base_patterns = None,
        base_pattern_positions = None,
        alignment_filename = None,
        tree_filename = None,
        info_filename = None,
        info_filetype = None,
        output_prefix = None,
        outgroup_name = None,
        threads = 1,
        verbose = False,
        mp_method = "spawn",
        max_pos = None):

    if verbose:
        prep_time = 0.0
        calc_time = 0.0
        prep_time_start = time.process_time()
    # Create a new alignment for the output containing all taxa in the input alignment
    alignment_sequence_names = {}
    for i, name in enumerate(sequence_names):
        alignment_sequence_names[name] = i
    
    # Read the tree
    if verbose:
        print("Reading tree file:", tree_filename)
    tree=read_tree(tree_filename)
    
    # Read the info file and get frequencies and rates
    if info_filename != "":
        if verbose:
            print("Reading info file:", info_filename)
        f,r = read_info(info_filename, type = info_filetype)
    else:
        if verbose:
            print("Using default JC rates and frequencies")
        f=numpy.array([0.25,0.25,0.25,0.25], dtype = numpy.float32)
        r=numpy.array([1.0,1.0,1.0,1.0,1.0,1.0], dtype = numpy.float32)
    
    if verbose:
        print("Frequencies:", ", ".join(map(str,f)))
        print("Rates:", ", ".join(map(str,r)))
    
    # Create rate matrix from f and r
    rm = create_rate_matrix(f,r)

    # Label internal nodes in tree and add these to the new alignment and calculate pij per non-root branch
    nodecounter=0
    num_nodes = len(tree.nodes())
    node_indices = {}
    child_nodes_array = numpy.empty(num_nodes, dtype=object)
    leaf_node_list = []
    node_labels = numpy.empty(num_nodes, dtype=object)
    node_pij_3d = numpy.full((num_nodes, 4, 4), -numpy.inf, dtype=numpy.float32)
    postordered_nodes = numpy.arange(num_nodes, dtype=numpy.int32)
    seed_node = -1
    seed_node_edge_truncation = True
    node_index_to_aln_row = numpy.full(num_nodes, -1, dtype=numpy.int32)
    ancestral_node_indices = {}

    for node_index,node in zip(postordered_nodes,tree.postorder_node_iter()):
        if node.taxon == None:
            nodecounter+=1
            nodename="Node_"+str(nodecounter)
            tree.taxon_namespace.add_taxon(dendropy.Taxon(nodename))
            node.taxon=tree.taxon_namespace.get_taxon(nodename)
            if nodename in alignment_sequence_names:
                print(nodename, "already in alignment. Quitting")
                sys.exit(209)
            ancestral_node_indices[node_index] = nodename
        else:
            node.taxon.label = node.taxon.label.strip("'")
            if node.taxon.label in alignment_sequence_names:
                node_index_to_aln_row[node_index] = alignment_sequence_names[node.taxon.label]
            else:
                sys.stderr.write('Unable to find ' + node.taxon.label + ' in alignment')
                sys.exit(1)
        if node.parent_node == None:
            seed_node = node_index
        elif node.parent_node == tree.seed_node and seed_node_edge_truncation:
            if outgroup_name is None or \
                            (outgroup_name is not None and node.taxon.label != outgroup_name):
                seed_node_edge_truncation = False
                node_pij_3d[node_index] = calculate_pij(node.edge_length/1e6, rm).reshape((4, 4))
            else:
                node_pij_3d[node_index] = calculate_pij(node.edge_length, rm).reshape((4, 4))
        else:
            node_pij_3d[node_index] = calculate_pij(node.edge_length, rm).reshape((4, 4))
        node_label = node.taxon.label
        node_indices[node_label] = node_index
        node_labels[node_index] = node_label
        if node.is_leaf():
            leaf_node_list.append(node_index)
            child_nodes_array[node_index] = numpy.full(1, -1, dtype=numpy.int32)
        else:
            child_nodes_array[node_index] = numpy.array([node_indices[child.taxon.label] for child in node.child_node_iter()],
                                                    dtype=numpy.int32)
    child_nodes = convert_to_square_numpy_array(child_nodes_array)

    # Store the preordered nodes and record parent node information
    parent_nodes = numpy.full(num_nodes, -1, dtype = numpy.int32)
    preordered_nodes = numpy.full(num_nodes-1, -1, dtype=numpy.int32)
    for node_count,node in enumerate(tree.preorder_node_iter()):
        if node.parent_node != None:
            node_index = node_indices[node.taxon.label]
            preordered_nodes[node_count-1] = node_index
            parent_nodes[node_index] = node_indices[node.parent_node.taxon.label]

    # Transposed contiguous ASCII alignment buffer (n_ancestral, max_pos), initialized to 78 (ASCII 'N')
    out_aln_T = numpy.full((len(ancestral_node_indices), max_pos), 78, dtype=numpy.uint8)
    ancestral_node_order = numpy.fromiter(ancestral_node_indices.keys(), dtype=numpy.int32)
    ordered_ascii_bases = numpy.array([65, 67, 71, 84, 45], dtype=numpy.uint8)
    log_f = numpy.log(f).astype(numpy.float32)

    if isinstance(base_pattern_positions, FlattenedBasePositions):
        positions_flat = base_pattern_positions.positions_flat
        position_offsets = base_pattern_positions.position_offsets
    else:
        lengths = numpy.array([len(p) for p in base_pattern_positions], dtype=numpy.int32)
        position_offsets = numpy.empty(len(lengths) + 1, dtype=numpy.int32)
        position_offsets[0] = 0
        numpy.cumsum(lengths, out=position_offsets[1:])
        positions_flat = numpy.concatenate(base_pattern_positions).astype(numpy.int32, copy=False)

    # Reconstruct each base position
    if verbose:
        print("Reconstructing sites on tree")
        prep_time_end = time.process_time()
        prep_time = prep_time_end - prep_time_start
        calc_time_start = time.process_time()

    if threads > 1:
        try:
            numba.set_num_threads(min(int(threads), numba.config.NUMBA_NUM_THREADS))
        except Exception:
            pass
        total_node_snps = reconstruct_all_patterns_T_parallel(
            base_patterns,
            positions_flat,
            position_offsets,
            out_aln_T,
            ordered_ascii_bases,
            postordered_nodes,
            preordered_nodes,
            parent_nodes,
            child_nodes,
            numpy.int32(seed_node),
            ancestral_node_order,
            node_pij_3d,
            log_f,
            node_index_to_aln_row,
            int(threads)
        )
    else:
        total_node_snps = reconstruct_all_patterns_T(
            base_patterns,
            positions_flat,
            position_offsets,
            out_aln_T,
            ordered_ascii_bases,
            postordered_nodes,
            preordered_nodes,
            parent_nodes,
            child_nodes,
            numpy.int32(seed_node),
            ancestral_node_order,
            node_pij_3d,
            log_f,
            node_index_to_aln_row
        )

    # Process outputs: write contiguous ASCII rows directly
    if verbose:
        print("Printing alignment with internal node sequences: ", output_prefix+".joint.aln")
    source = alignment_filename
    destination = output_prefix+".joint.aln"
    dest = shutil.copy(source, destination)
    with open(dest, "ab") as asr_output:
        for i, node_index in enumerate(ancestral_node_order):
            taxon = ancestral_node_indices[node_index]
            asr_output.write(('>' + taxon + '\n').encode('ascii'))
            asr_output.write(out_aln_T[i])
            asr_output.write(b"\n")

    # Set branch lengths to total reconstructed SNPs per branch
    for node in tree.preorder_node_iter():
        node_index = node_indices[node.taxon.label]
        node.edge_length = float(total_node_snps[node_index])

    ### TIMING
    if verbose:
        calc_time_end = time.process_time()
        calc_time = (calc_time_end - calc_time_start)

    # Print tree
    from gubbins.common import tree_as_string
    
    if verbose:
        print("Printing tree with internal nodes labelled: ", output_prefix+".joint.tre")
    with open(output_prefix+".joint.tre", "w") as tree_output:
    
        recon_tree = tree_as_string(tree,
                                    suppress_rooting=True,
                                    suppress_internal=False)
        print(recon_tree.replace('\'', ''),
              file = tree_output)

    if verbose:
        print("Done")
        print('Time for JAR preparation:\t' + str(prep_time))
        print('Time for JAR calculation:\t' + str(calc_time))
