#!/usr/bin/env python3

# pyjar written by Simon Harris
# code modified from https://github.com/simonrharris/pyjar
# pyjar is free software, licensed under GPLv3.

from asyncio import subprocess
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
from numba import jit, njit, types, from_dtype
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
        pij = numpy.full((4,4), numpy.NINF, dtype = numpy.float32)
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
    child_node_indices = child_nodes[node_index,:]
    child_node_indices = child_node_indices[child_node_indices > -1]
    for end_index in column_base_indices:
        c = Lmat[child_node_indices,end_index].sum()
        for start_index in column_base_indices:
            j = pij[start_index,end_index]+c
            if j > Lmat[node_index,start_index]:
                Lmat[node_index,start_index] = j
                Cmat[node_index,start_index] = end_index

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
    for end_index in column_base_indices:
        c = Lmat[child_node_indices,end_index].sum()
        for start_index in column_base_indices:
            j = log(base_frequencies[end_index]) + c
            if j > Lmat[node_index,start_index]:
                Lmat[node_index,start_index] = j
                Cmat[node_index,start_index] = end_index

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
        #1a. Let j be the amino acid at y. Set, for each amino acid i: Cy(i)= j. This implies that no matter what is the amino acid in the father of y, j is assigned to node y.
        Cmat[node_index,:] = taxon_base_index
        #1b. Set for each amino acid i: Ly(i) = Pij(ty), where ty is the branch length between y and its father.
        Lmat[node_index,:] = [pij[0][taxon_base_index],pij[1][taxon_base_index],pij[2][taxon_base_index],pij[3][taxon_base_index]]
    else:
        # Cmat stays as default when base is unknown
        Lmat[node_index,:] = [pij[0][0],pij[1][1],pij[2][2],pij[3][3]]

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
    for node_index in preordered_nodes:
        parent_node_index = parent_nodes[node_index]
        if reconstructed_alleles[node_index] < 4 and reconstructed_alleles[parent_node_index] < 4 \
          and reconstructed_alleles[node_index] != reconstructed_alleles[parent_node_index]:
            node_snps[node_index] += len(base_pattern_columns)

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
    for node_index in postordered_nodes:
        if node_index in leaf_nodes:
            alignment_index = node_index_to_aln_row[node_index]
            reconstructed_alleles[node_index] = column[alignment_index]
        else:
            has_child_base = False
            for child_taxon_index in child_nodes[node_index,:]:
                if child_taxon_index > -1:
                    if reconstructed_alleles[child_taxon_index] < 4:
                        has_child_base = True
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
    for index in numpy.arange(ancestral_node_order.size, dtype=numpy.int32):
        node_index = ancestral_node_order[index]
        base_index = reconstructed_alleles[node_index]
        base = ordered_bases[base_index]
        for column in base_pattern_columns:
            out_aln[column,index] = base

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

    Cmat_null = numpy.array([0,1,2,3], dtype = numpy.uint8)
        
    # Reset matrices
    Lmat.fill(numpy.NINF)
    Cmat[:] = Cmat_null

    # Count unknown bases
    unknown_base_count = numpy.count_nonzero(column > 3)
    column_base_indices = numpy.unique(column[numpy.where(column <= 3)])
      
    # Heuristic for speed: if all taxa are monomorphic, with a gap in only one sequence, then the ancestral states
    # will all be the observed base, as no ancestral node will have two child nodes with unknown bases at this site
    if unknown_base_count == 1 and column_base_indices.size == 1:
        # If site is monomorphic - replace entire column
        tmp_out_aln[base_pattern_columns,:] = ordered_bases[column_base_indices[0]]
    else:
        # Otherwise perform a full ML inference
        #1 For each OTU y perform the following:
        #Visit a nonroot internal node, z, which has not been visited yet, but both of whose sons, nodes x and y, have already been visited, i.e., Lx(j), Cx(j), Ly(j), and Cy(j) have already been defined for each j. Let tz be the length of the branch connecting node z and its father. For each amino acid i, compute Lz(i) and Cz(i) according to the following formulae:
        #Denote the three sons of the root by x, y, and z. For each amino acid k, compute the expression Pk x Lx(k) x Ly(k) x Lz(k). Reconstruct r by choosing the amino acid k maximizing this expression. The maximum value found is the likelihood of the best reconstruction.
        for node_index in postordered_nodes:
            if node_index == seed_node:
                continue
            #calculate the transistion matrix for the branch
            pij=numpy.reshape(node_pij[node_index,:].copy(),(4,4))
            if node_index in leaf_nodes:
                alignment_index = node_index_to_aln_row[node_index]
                taxon_base_index = column[alignment_index]
                process_leaf(Lmat,
                                Cmat,
                                pij,
                                node_index,
                                taxon_base_index)
            else:
                #2a. Lz(i) = maxj Pij(tz) x Lx(j) x Ly(j)
                #2b. Cz(i) = the value of j attaining the above maximum.
                find_most_likely_base_given_descendents(Lmat,
                                                        Cmat,
                                                        pij,
                                                        node_index,
                                                        child_nodes,
                                                        column_base_indices)

        # Calculate likelihood of base at root node
        child_node_indices = child_nodes[node_index,:]
        child_node_indices = child_node_indices[child_node_indices > -1]
        calculate_root_likelihood(Lmat, Cmat, base_frequencies, node_index, child_node_indices, column_base_indices)
        max_root_base_index = Cmat[node_index,numpy.argmax(Lmat[node_index,:])]
        reconstructed_base_indices[node_index] = max_root_base_index
        
        #Traverse the tree from the root in the direction of the OTUs, assigning to each node its most likely ancestral character as follows:
        # Note that preordered node list does not include the root
        for node_index in preordered_nodes:
            #5a. Visit an unreconstructed internal node x whose father y has already been reconstructed. Denote by i the reconstructed amino acid at node y.
            parent_node_index = parent_nodes[node_index]
            i = reconstructed_base_indices[parent_node_index]
            #5b. Reconstruct node x by choosing Cx(i).
            reconstructed_base_indices[node_index] = Cmat[node_index,i]

        # Put gaps back in and check that any ancestor with only gaps downstream is made a gap
        # store reconstructed alleles
        reconstructed_alleles = numpy.full(postordered_nodes.size, 8, dtype = numpy.uint8)
        reconstruct_alleles(reconstructed_alleles,
                            postordered_nodes,
                            leaf_nodes,
                            node_index_to_aln_row,
                            column,
                            child_nodes,
                            reconstructed_base_indices
                            )

        # If site is not monomorphic - replace specific entries
        fill_out_aln(tmp_out_aln,
                    reconstructed_alleles,
                    ordered_bases,
                    ancestral_node_order,
                    base_pattern_columns
                    )

        # enumerate the number of base subtitutions reconstructed occurring on each branch
        count_node_snps(node_snps,
                        preordered_nodes,
                        parent_nodes,
                        seed_node,
                        reconstructed_alleles,
                        base_pattern_columns,
                        )

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
    
    # Check njit function is compiled before multiprocessing
    try:
        seq_to_int()
    except:
        pass

    # Read in ordered sequence names
    sequence_names_fn = prefix + '.gaps.sequence_names.csv'
    sequence_names = []
    if not os.path.isfile(sequence_names_fn):
        sys.exit("Unable to open sequence names file " + sequence_names_fn + "\n")
    with open(sequence_names_fn,'r') as names_file:
        for line in names_file:
            sequence_names.append(line.rstrip())

    # Convert alignment to Numpy array
    codec = 'utf-32-le' if sys.byteorder == 'little' else 'utf-32-be'
    # Read in the unique base patterns
    base_patterns_fn = prefix + '.gaps.base_patterns.csv'
    if not os.path.isfile(base_patterns_fn):
        sys.exit("Unable to open base patterns file " + base_patterns_fn + "\n")
    array_of_pattern_arrays = []
    with open(base_patterns_fn,'r') as pattern_file:
        for line in pattern_file:
            unicode_seq = numpy.frombuffer(bytearray(line.rstrip(), codec), dtype = 'U1')
            out_array = numpy.ndarray(unicode_seq.shape, dtype = numpy.uint8)
            seq_to_int(unicode_seq,out_array)
            array_of_pattern_arrays.append(out_array)
    vstacked_patterns = numpy.vstack(array_of_pattern_arrays)
    # Read in base positions
    array_of_position_arrays = []
    base_positions_fn = prefix + '.gaps.base_positions.csv'
    if not os.path.isfile(base_positions_fn):
        sys.exit("Unable to open base positions file " + base_positions_fn + "\n")
    array_max = 0
    with open(base_positions_fn, 'r') as positions_file:
        for line in positions_file:
            subarray = numpy.asarray(list(map(int,line.rstrip().split(','))),
                                    dtype = numpy.int32)
            array_of_position_arrays.append(subarray)
            subarray_max = numpy.amax(subarray)
            if subarray_max > array_max:
                    array_max = subarray_max

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
    
    ### TIMING
    # if verbose:
    #     prep_time = 0.0
    #     calc_time = 0.0
    #     prep_time_start = time.process_time()

    # Load shared memory output alignment
    out_aln_shm = shared_memory.SharedMemory(name = new_aln.name)
    out_aln = numpy.ndarray(new_aln.shape, dtype = numpy.int8, buffer = out_aln_shm.buf)
    
    # Generate data structures for reconstructions
    num_nodes = len(tree.nodes())
    Lmat = numpy.full((num_nodes,4), numpy.NINF, dtype = numpy.float32)
    Cmat = numpy.full((num_nodes,4), [0,1,2,3], dtype = numpy.uint8)
    reconstructed_base_indices = numpy.full(num_nodes, 8, dtype = numpy.uint8)

    # Record SNPs reconstructed as occurring on each branch
    
    ordered_bases = numpy.array([0,1,2,3,4], dtype=numpy.int8)
    node_snps = numpy.zeros(num_nodes, dtype = numpy.int32)

    # Load base pattern information
    base_patterns_shm = shared_memory.SharedMemory(name = base_patterns.name)
    base_patterns = numpy.ndarray(base_patterns.shape, dtype = base_patterns.dtype, buffer = base_patterns_shm.buf)
    # Load base pattern position information
    
    # Extract information for iterations
    column = base_patterns[column_index]

    ### TIMING
    # if verbose:
    #     prep_time_end = time.process_time()
    #     prep_time = prep_time_end - prep_time_start
    #     calc_time_start = time.process_time()

    # Iterate over columns
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


    ### TIMING
    # if verbose:
    #     calc_time_end = time.process_time()
    #     calc_time = (calc_time_end - calc_time_start)

    # Close shared memory
    out_aln_shm.close()
    base_patterns_shm.close()
    

    ### TIMING
    # if verbose:
    #     print('Time for JAR preparation:\t' + str(prep_time))
    #     print('Time for JAR calculation:\t' + str(calc_time))

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
    node_pij = numpy.full((num_nodes,16), numpy.NINF, dtype=numpy.float32)
    postordered_nodes = numpy.arange(num_nodes, dtype=numpy.int32)
    seed_node = None
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
            # Set the length of one root-to-child branch to ~zero
            # as reconstruction should occur with rooting at a node
            # midpoint rooting causes problems at the root, especially w/JC69
            # With outgroup:
            # Prefer to truncate the ingroup branch to force the recombinations
            # onto the outgroup branch
            if outgroup_name is None or \
                            (outgroup_name is not None and node.taxon.label != outgroup_name):
                seed_node_edge_truncation = False
                node_pij[node_index,:]=calculate_pij(node.edge_length/1e6, rm)
            else:
                node_pij[node_index,:]=calculate_pij(node.edge_length, rm)
        else:
            node_pij[node_index,:]=calculate_pij(node.edge_length, rm)
        # Store information to avoid subsequent recalculation as
        # look up of taxon labels with dendropy is slower than native data structures
        node_label = node.taxon.label
        node_indices[node_label] = node_index
        node_labels[node_index] = node_label
        if node.is_leaf():
            leaf_node_list.append(node_index)
            child_nodes_array[node_index] = numpy.full(1, -1, dtype=numpy.int32) # Cannot leave array empty
        else:
            child_nodes_array[node_index] = numpy.array([node_indices[child.taxon.label] for child in node.child_node_iter()],
                                                    dtype=numpy.int32)
    leaf_nodes = numpy.array(leaf_node_list, dtype = numpy.int32)
    child_nodes = convert_to_square_numpy_array(child_nodes_array)

    # Store the preordered nodes and record parent node information
    parent_nodes = numpy.full(num_nodes, -1, dtype = numpy.int32)
    preordered_nodes = numpy.full(num_nodes-1, -1, dtype=numpy.int32)
    for node_count,node in enumerate(tree.preorder_node_iter()):
        if node.parent_node != None:
            node_index = node_indices[node.taxon.label]
            preordered_nodes[node_count-1] = node_index # Do not add root node to preordered nodes
            parent_nodes[node_index] = node_indices[node.parent_node.taxon.label]

        
    ## Switch this to a numpy.int8 data type use the default as 5 the corresponding value to N
    ## form the seq_to_int transformation o the sequence 
    new_aln_array = numpy.full((max_pos,len(ancestral_node_indices)), 5, dtype = numpy.int8)
    
    # Index names for reconstruction
    ancestral_node_order = numpy.fromiter(ancestral_node_indices.keys(), dtype=numpy.int32)

    # Compile functions prior to multiprocessing
    for func in [find_most_likely_base_given_descendents,
                process_leaf,
                calculate_root_likelihood,
                count_node_snps,
                reconstruct_alleles,
                fill_out_aln,
                iterate_over_base_patterns,
                int_to_seq]:
        try:
            func()
        except:
            pass

    # Reconstruct each base position
    if verbose:
        print("Reconstructing sites on tree")

    with SharedMemoryManager() as smm:
    
        # Convert alignment to shared memory numpy array
        new_aln_shared_array = generate_shared_mem_array(new_aln_array, smm)
        
        # Convert base patterns to shared memory numpy array
        base_patterns_shared_array = generate_shared_mem_array(base_patterns, smm)
        
        # Convert base pattern positions to shared memory numpy array
        bp_list = list(range(len(base_patterns)))
        npatterns = len(base_patterns)

        if verbose:
            prep_time_end = time.process_time()
            prep_time = prep_time_end - prep_time_start
            calc_time_start = time.process_time()


        if threads > 1:

            # Parallelise reconstructions across alignment columns using multiprocessing
            with multiprocessing.get_context(method=mp_method).Pool(processes = threads) as pool:
                reconstruction_results = pool.starmap(partial(
                                            reconstruct_alignment_column,
                                                tree = tree,
                                                preordered_nodes = preordered_nodes,
                                                postordered_nodes = postordered_nodes,
                                                leaf_nodes = leaf_nodes,
                                                parent_nodes = parent_nodes,
                                                child_nodes = child_nodes,
                                                seed_node = seed_node,
                                                node_pij = node_pij,
                                                node_index_to_aln_row = node_index_to_aln_row,
                                                ancestral_node_order = ancestral_node_order,
                                                base_patterns = base_patterns_shared_array,
                                                base_frequencies = f,
                                                new_aln = new_aln_shared_array,
                                                threads = threads,
                                                verbose = verbose),
                                            zip(bp_list, base_pattern_positions)
                                        )
            
            # Write out alignment while shared memory manager still active
            out_aln_shm = shared_memory.SharedMemory(name = new_aln_shared_array.name)
            out_aln = numpy.ndarray(new_aln_array.shape, dtype = 'i1', buffer = out_aln_shm.buf)
                     
            # Release pool nodes
            pool.join()
        
        else:
        
            # Run as a loop on a single core
            reconstruction_results = \
                [
                            reconstruct_alignment_column(
                                b,
                                p,
                                tree = tree,
                                preordered_nodes = preordered_nodes,
                                postordered_nodes = postordered_nodes,
                                leaf_nodes = leaf_nodes,
                                parent_nodes = parent_nodes,
                                child_nodes = child_nodes,
                                seed_node = seed_node,
                                node_pij = node_pij,
                                node_index_to_aln_row = node_index_to_aln_row,
                                ancestral_node_order = ancestral_node_order,
                                base_patterns = base_patterns_shared_array,
                                base_frequencies = f,
                                new_aln = new_aln_shared_array,
                                threads = threads,
                                verbose = verbose) for b,p in zip(bp_list, base_pattern_positions)
                ]
            
            # Extract final result
            out_aln_shm = shared_memory.SharedMemory(name = new_aln_shared_array.name)
            out_aln = numpy.ndarray(new_aln_array.shape, dtype = 'i1', buffer = out_aln_shm.buf)

        # Process outputs
        aln_line = numpy.full(len(out_aln[:,0]),"?",dtype="U1")
        if verbose:
            print("Printing alignment with internal node sequences: ", output_prefix+".joint.aln")
        source = alignment_filename
        destination = output_prefix+".joint.aln"
        dest = shutil.copy(source, destination)
        with open(dest, "a") as asr_output:
            for i,node_index in enumerate(ancestral_node_order):
                taxon = ancestral_node_indices[node_index]
                asr_output.write('>' + taxon + '\n')
                int_to_seq(out_aln[:,i], aln_line)
                asr_output.write(''.join(aln_line) + "\n")

        # Combine results for each base across the alignment
        for node in tree.preorder_node_iter():
            node.edge_length = 0.0 # reset lengths to convert to SNPs
            node_index = node_indices[node.taxon.label]
            for x in range(len(reconstruction_results)):
                try:
                    node.edge_length += reconstruction_results[x][node_index];
                except AttributeError:
                    continue

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

