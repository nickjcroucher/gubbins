#!/usr/bin/env python3

# pyjar written by Simon Harris
# code modified from https://github.com/simonrharris/pyjar
# pyjar is free software, licensed under GPLv3.

from scipy import linalg
import numpy
import dendropy
import sys
import os
import time
from Bio import AlignIO
from math import log, exp, ceil
from functools import partial
import numba
from numba import jit, njit, types, from_dtype
import collections
from memory_profiler import profile
import psutil
import datetime
import argparse
import pickle

fp = open("memory_log", "w+")

try:
    from multiprocessing import shared_memory
    from multiprocessing.managers import SharedMemoryManager
    import multiprocessing
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

# Read in sequence to enable conversion to integers with JIT function

def process_sequence(index_list,alignment ,codec = None,align_array = None):
    # Load shared memory output alignment
    out_aln_shm = shared_memory.SharedMemory(name = align_array.name)
    out_aln = numpy.ndarray(align_array.shape, dtype = numpy.uint8, buffer = out_aln_shm.buf)
    for seq,i in enumerate(index_list):
        # Add sequence
        unicode_seq = numpy.frombuffer(bytearray(str(alignment[seq]), codec), dtype = 'U1')
        seq_to_int(unicode_seq,out_aln[i])

# Function to read an alignment in various formats
def read_alignment(filename, file_type, verbose=False, list_out=False):

    if not os.path.isfile(filename):
        print("Error: alignment file " + filename + " does not exist")
        sys.exit(202)
    if verbose:
        print("Trying to open file " + filename + " as " + file_type)
    try:
        with open(filename,'r') as aln_in:
            alignmentObject = AlignIO.read(aln_in, file_type)
        if verbose:
            print("Alignment read successfully")
    except:
        print("Cannot open alignment file " + filename + " as " + file_type)
        sys.exit(203)

    ## Convert to list of lists of seq types
    if list_out:
        aln_list = []
        for aln in alignmentObject:
            aln_list.append(str(aln.seq))
        return aln_list
    else:
        return alignmentObject


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
@njit(numba.void(numba.typeof(numpy.dtype('U1'))[:,:],
                numba.uint8[:],
                numba.typeof(numpy.dtype('U1'))[:],
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

# Return positions of columns in alignment
##########################################
@njit(numba.int32[:](numba.int32[:],
                numba.int32[:,:],
                numba.int32),
                cache=True)
def get_columns(base_pattern_columns_padded,column_positions,column_index):
    base_pattern_columns_indices = numpy.argmax(base_pattern_columns_padded == -1)
    if base_pattern_columns_indices == 0:
        base_pattern_columns_indices = base_pattern_columns_padded.size
    base_pattern_columns = column_positions[column_index,0:base_pattern_columns_indices]
    return base_pattern_columns

# Reconstruct each base pattern
###############################
@njit(numba.void(numba.uint8[:,:],
                numba.int32[:,:],
                numba.float32[:,:],
                numba.uint8[:,:],
                numba.typeof(numpy.dtype('U1'))[:,:],
                numba.typeof(numpy.dtype('U1'))[:],
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
def iterate_over_base_patterns(columns,
                                column_positions,
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

    column_indices = numpy.arange(columns.shape[0], dtype = numpy.int32)
    Cmat_null = numpy.array([0,1,2,3], dtype = numpy.uint8)
    
    for column_index in column_indices:
    
        # Get column bases
        column = columns[column_index]
        
        # Get column positions
        # base_pattern_columns_padded = column_positions[column_index]
        # base_pattern_columns = get_columns(base_pattern_columns_padded,column_positions,column_index)
        base_pattern_columns = column_positions[column_index]
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

##################
# Main functions #
##################

####################################################
# Function for converting alignment to numpy array #
####################################################
@profile(stream = fp)
def get_base_patterns(alignment, verbose,
                      printero = "printer_output", fit_method = "spawn",
                      threads = 1, pickle_aln = False):

    if verbose:
        print("Finding unique base patterns")
    # Identify unique base patterns
    t1=time.process_time()
    # Convert alignment to Numpy array
    ntaxa = len(alignment)

    seq_length = len(alignment[0])
    ## Now to create the list of alignments
    ## Manipulate aln list into list of lists
    print_file = open(printero, "a")
    print_file.write("Creating list of lists for aln " + str(datetime.datetime.now()) + "\n")
    print_file.write("Starting mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.close()
    ntaxa_jumps = ceil(ntaxa  / threads)
    aln_list = [alignment[i: i+ntaxa_jumps] for i in range(0, len(alignment), ntaxa_jumps)]
    print_file = open(printero, "a")
    print_file.write("Finished creating list of lists " + str(datetime.datetime.now()) + "\n")
    print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()

    print_file = open(printero, "a")
    print_file.write("Deleting input aln " + str(datetime.datetime.now()) + "\n")
    print_file.write("Starting mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.close()
    del(alignment)
    print_file = open(printero, "a")
    print_file.write("Deleting input aln " + str(datetime.datetime.now()) + "\n")
    print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()

    print_file = open(printero, "a")
    print_file.write("Creating initial align array " + str(datetime.datetime.now()) + "\n")
    print_file.write("Starting mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3)+ "\n")
    print_file.close()
    align_array = numpy.full((ntaxa,seq_length), 8, dtype = numpy.uint8, order='F')
    print_file = open(printero, "a")
    print_file.write("Finished initial align array " + str(datetime.datetime.now()) + "\n")
    print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3)+ "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()

    # Check njit function is compiled before multiprocessing
    try:
        seq_to_int()
    except:
        pass
    # Convert alignment to Numpy array
    codec = 'utf-32-le' if sys.byteorder == 'little' else 'utf-32-be'
    ntaxa_range_list = list(range(ntaxa))

    ntaxa_range_indices = [ntaxa_range_list[i: i+ntaxa_jumps] for i in range(0, len(ntaxa_range_list), ntaxa_jumps)]
    #list(chunks(ntaxa_range_list,threads))
    if pickle_aln:
        with open("pickled_aln.txt","wb") as fh:
            pickle.dump(aln_list, fh)



    with SharedMemoryManager() as smm:
        print_file = open(printero, "a")
        print_file.write("Starting shared memory array " + str(datetime.datetime.now()) + "\n")
        print_file.write("Starting mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
        #print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
        print_file.close()
        align_array_shared = generate_shared_mem_array(align_array, smm)
        print_file = open(printero, "a")
        print_file.write("Created shared memory array " + str(datetime.datetime.now()) + "\n")
        print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3)+ "\n")
        print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
        print_file.close()
        with multiprocessing.get_context(fit_method).Pool() as pool:
            print_file = open(printero, "a")
            print_file.write("Starting process sequence job " + str(datetime.datetime.now()) + "\n")
            print_file.write("Starting mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
            #print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
            print_file.close()
            pool.starmap(partial(
                process_sequence,
                    codec = codec,
                    align_array = align_array_shared
                ),
                zip(ntaxa_range_indices, aln_list)
            )
            print_file = open(printero, "a")
            print_file.write("Finished process sequence job " + str(datetime.datetime.now()) + "\n")
            print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
            print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
            print_file.close()

        # Write out alignment while shared memory manager still active
        align_array_shm = shared_memory.SharedMemory(name = align_array_shared.name)
        align_array = numpy.ndarray(align_array.shape, dtype = numpy.uint8, buffer = align_array_shm.buf)

    # Get unique base patterns and their indices in the alignment

    print_file = open(printero, "a")
    print_file.write("Staring unique column names " + str(datetime.datetime.now()) + "\n")
    print_file.write("Start mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    #print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()
    if pickle_aln:
        with open("align_array.npy","wb") as f:
            numpy.save(f, align_array)

    base_pattern_bases_array, base_pattern_positions_array = get_unique_columns(align_array)
    print_file = open(printero, "a")
    print_file.write("End unique column names " + str(datetime.datetime.now()) + "\n")
    print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3)+ "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()

    base_pattern_positions_array_of_arrays = \
        [numpy.where(base_pattern_positions_array==x)[0] for x in range(base_pattern_bases_array.shape[1])]

    # Convert the array of arrays into an ndarray that can be saved to shared memory

    print_file = open(printero, "a")
    print_file.write("Staring conversion to square numpy array " + str(datetime.datetime.now()) + "\n")
    print_file.write("Start mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3)+ "\n")
    print_file.close()
    #square_base_pattern_positions_array = convert_to_square_numpy_array(base_pattern_positions_array_of_arrays)
    print_file = open(printero, "a")
    print_file.write("End conversion to square numpy array " + str(datetime.datetime.now()) + "\n")
    print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()
    # Finish
    t2=time.process_time()
    if verbose:
        print("Time taken to find unique base patterns:", t2-t1, "seconds")

        print("Unique base patterns:", str(base_pattern_bases_array.shape[1]))
    return base_pattern_bases_array.transpose(), base_pattern_positions_array_of_arrays#square_base_pattern_positions_array


########################################################
# Function for reconstructing individual base patterns #
########################################################

def reconstruct_alignment_column(column_indices,
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
                                #base_pattern_positions = None,
                                base_frequencies = None,
                                new_aln = None,
                                threads = 1,
                                verbose = False,
                                printero = "./printer_output"):

    
    ### TIMING
    if verbose:
        prep_time = 0.0
        calc_time = 0.0
        prep_time_start = time.process_time()

    # Load shared memory output alignment
    out_aln_shm = shared_memory.SharedMemory(name = new_aln.name)
    out_aln = numpy.ndarray(new_aln.shape, dtype = 'U1', buffer = out_aln_shm.buf)
    
    # Generate data structures for reconstructions
    num_nodes = len(tree.nodes())
    Lmat = numpy.full((num_nodes,4), numpy.NINF, dtype = numpy.float32)
    Cmat = numpy.full((num_nodes,4), [0,1,2,3], dtype = numpy.uint8)
    reconstructed_base_indices = numpy.full(num_nodes, 8, dtype = numpy.uint8)

    # Record SNPs reconstructed as occurring on each branch
    bases = frozenset(['A','C','G','T'])
    ordered_bases = numpy.array(['A','C','G','T','-'], dtype = 'U1')
    node_snps = numpy.zeros(num_nodes, dtype = numpy.int32)

    # Load base pattern information
    base_patterns_shm = shared_memory.SharedMemory(name = base_patterns.name)
    base_patterns = numpy.ndarray(base_patterns.shape, dtype = base_patterns.dtype, buffer = base_patterns_shm.buf)
    # Load base pattern position information
    # No Need to with the chunked up positions
    #base_pattern_positions_shm = shared_memory.SharedMemory(name = base_pattern_positions.name)
    #base_pattern_positions = numpy.ndarray(base_pattern_positions.shape, dtype = base_pattern_positions.dtype, buffer = base_pattern_positions_shm.buf)

    # Extract information for iterations
    if threads == 1:
        columns = base_patterns
        column_positions = base_pattern_positions
    else:
        columns = base_patterns[column_indices]
        print_file = open(printero, "a")
        print_file.write("Converting back into square array " + str(datetime.datetime.now()) + "\n")
        print_file.write("Start mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
        print_file.write("Process id: " + str(multiprocessing.current_process()) + "\n")
        print_file.close()
        column_positions = convert_to_square_numpy_array(base_pattern_positions)
        print_file = open(printero, "a")
        print_file.write("End conversion to square numpy array " + str(datetime.datetime.now()) + "\n")
        print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) +  "\n")
        print_file.write("Process id: " + str(multiprocessing.current_process()) + "\n")
        print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" +  "\n")
        print_file.close()


    ### TIMING
    if verbose:
        prep_time_end = time.process_time()
        prep_time = prep_time_end - prep_time_start
        calc_time_start = time.process_time()

    # Iterate over columns

    print_file = open(printero, "a")
    print_file.write("Starting iteration over base patterns " + str(datetime.datetime.now()) + "\n")
    print_file.write("Start mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.write("Process id: " + str(multiprocessing.current_process()) + "\n")
    print_file.close()
    mem_start = psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3
    tim_start = datetime.datetime.now()
    iterate_over_base_patterns(columns,
                               column_positions,
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
    mem_use = (psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) - mem_start
    tim_end = datetime.datetime.now() - tim_start

    print_file = open(printero, "a")
    print_file.write("End iteration over base patterns " + str(datetime.datetime.now()) + "\n")
    print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.write("Process id: " + str(multiprocessing.current_process()) + "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.write("End iteration time, mem " + str(multiprocessing.current_process()) + "\n")
    print_file.write("End iteration time === mem " + str(tim_end) + "(time) === " + str(mem_use) + " (GB)" + "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()
    ### TIMING
    if verbose:
        calc_time_end = time.process_time()
        calc_time = (calc_time_end - calc_time_start)

    # Close shared memory
    out_aln_shm.close()
    base_patterns_shm.close()
    #base_pattern_positions_shm.close()

    ### TIMING
    if verbose:
        print('Time for JAR preparation:\t' + str(prep_time))
        print('Time for JAR calculation:\t' + str(calc_time))

    return node_snps

##################################################
# Function for reconstructing complete alignment #
##################################################

def jar(alignment = None,
        base_patterns = None,
        base_pattern_positions = None,
        tree_filename = None,
        info_filename = None,
        info_filetype = None,
        output_prefix = None,
        threads = 1,
        verbose = False,
        printero = "./printer_output",
        mp_metho = "spawn"):

    # Lookup for each base
    mb={"A": 0, "C": 1, "G": 2, "T": 3}
    print_file = open(printero, "a")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.write("<><><><><><><><><><><><><><><><><><><>" + "\n")
    print_file.write("JAR function start " + str(datetime.datetime.now()) + "\n")
    print_file.write("Start mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3)+ "\n")
    print_file.write("<><><><><><><><><><><><><><><><><><><>" + "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()
    #square_base_pattern_positions_array = convert_to_square_numpy_array(base_pattern_positions_array_of_arrays)
    print_file = open(printero, "a")
    print_file.write("End conversion to square numpy array " + str(datetime.datetime.now()) + "\n")
    print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()

    # Create a new alignment for the output containing all taxa in the input alignment
    alignment_sequence_names = {}
    for i, x in enumerate(alignment):
        alignment_sequence_names[x.id] = i
    
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
    print_file = open(printero, "a")
    print_file.write("Staring tree labelling jar " + str(datetime.datetime.now()) + "\n")
    print_file.write("Start mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.close()

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
            seed_node_edge_truncation = False
            node_pij[node_index,:]=calculate_pij(node.edge_length/1e6, rm)
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
    print_file = open(printero, "a")
    print_file.write("End tree labelling jar " + str(datetime.datetime.now()) + "\n")
    print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()

    # Store the preordered nodes and record parent node information
    parent_nodes = numpy.full(num_nodes, -1, dtype = numpy.int32)
    preordered_nodes = numpy.full(num_nodes-1, -1, dtype=numpy.int32)
    for node_count,node in enumerate(tree.preorder_node_iter()):
        if node.parent_node != None:
            node_index = node_indices[node.taxon.label]
            preordered_nodes[node_count-1] = node_index # Do not add root node to preordered nodes
            parent_nodes[node_index] = node_indices[node.parent_node.taxon.label]

    # Create new empty array

    print_file = open(printero, "a")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.write("|/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\|" + "\n")
    print_file.write("Alignment array creation " + str(datetime.datetime.now()) + "\n")
    print_file.write("Start mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3)+ "\n")
    
    print_file.close()
    
    new_aln_array = numpy.full((len(alignment[0]),len(ancestral_node_indices)), '?', dtype = 'U1')
    print_file = open(printero, "a")
    print_file.write("End Alignment array creation " + str(datetime.datetime.now()) + "\n")
    print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
    print_file.write("Size of new_aln_array " + str(new_aln_array.shape) + "\n")
    print_file.write("Memory size of new aln_array (bytes): " + str(new_aln_array.nbytes) + " (bytes)" + "\n")
    print_file.write("|/-\|/-\||/-\|/-\||/-\|/-\||/-\|/-\|" + "\n")
    print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
    print_file.close()
    # Index names for reconstruction
    ancestral_node_order = numpy.fromiter(ancestral_node_indices.keys(), dtype=numpy.int32)

    # Compile functions prior to multiprocessing
    for func in [find_most_likely_base_given_descendents,
                process_leaf,
                calculate_root_likelihood,
                count_node_snps,
                reconstruct_alleles,
                fill_out_aln,
                iterate_over_base_patterns]:
        try:
            func()
        except:
            pass

    # Reconstruct each base position
    if verbose:
        print("Reconstructing sites on tree")

    with SharedMemoryManager() as smm:
    
        # Convert alignment to shared memory numpy array
        print_file = open(printero, "a")
        print_file.write("Start Alignment shared jar " + str(datetime.datetime.now()) + "\n")
        print_file.write("Start mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
        print_file.close()
        new_aln_shared_array = generate_shared_mem_array(new_aln_array, smm)
        print_file = open(printero, "a")
        print_file.write("End Alignment shared jar " + str(datetime.datetime.now()) + "\n")
        print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
        print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
        print_file.close()

        # Convert base patterns to shared memory numpy array
        print_file = open(printero, "a")
        print_file.write("Generating base patterns shared jar " + str(datetime.datetime.now()) + "\n")
        print_file.write("Start mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
        print_file.close()
        base_patterns_shared_array = generate_shared_mem_array(base_patterns, smm)
        print_file = open(printero, "a")
        print_file.write("End base patterns shared jar " + str(datetime.datetime.now()) + "\n")
        print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
        print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
        print_file.close()
        # Convert base pattern positions to shared memory numpy array
        #base_pattern_positions_shared_array = generate_shared_mem_array(base_pattern_positions, smm)

        # split list of sites into chunks per core
        bp_list = list(range(len(base_patterns)))
        #base_pattern_indices = list(chunks(bp_list,threads))
        npatterns = len(base_patterns)
        ntaxa_jumps = ceil(npatterns / threads)
        ntaxa_range_list = list(range(npatterns))
        base_pattern_indices = [bp_list[i: i + ntaxa_jumps] for i in range(0, len(bp_list), ntaxa_jumps)]

        base_positions = [base_pattern_positions[i:i + ntaxa_jumps] for i in range(0, len(base_pattern_positions), ntaxa_jumps)]

        # Parallelise reconstructions across alignment columns using multiprocessing
        print_file = open(printero, "a")
        print_file.write("Starting Alignment reconstruction " + str(datetime.datetime.now()) + "\n")
        print_file.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" + "\n")
        print_file.write("This is the dimension of the base_pattern_positions " + str(base_patterns.shape) + "\n")
        print_file.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" + "\n")
        print_file.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" + "\n")
        print_file.write("This is the dimension of the shared base_pattern_positions " + str(base_patterns_shared_array.shape) + "\n")
        print_file.write("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" + "\n")
        print_file.write("Start mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
        print_file.close()
        with multiprocessing.get_context(method=mp_metho).Pool(processes=threads) as pool:
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
                                            #base_pattern_positions = base_pattern_positions_shared_array,
                                            base_frequencies = f,
                                            new_aln = new_aln_shared_array,
                                            threads = threads,
                                            verbose = verbose,
                                            printero = printero),
                                        zip(base_pattern_indices, base_positions)
                                    )

        # Write out alignment while shared memory manager still active
        out_aln_shm = shared_memory.SharedMemory(name = new_aln_shared_array.name)
        out_aln = numpy.ndarray(new_aln_array.shape, dtype = 'U1', buffer = out_aln_shm.buf)
        print_file = open(printero, "a")
        print_file.write("End alignment reconstruction jar " + str(datetime.datetime.now()) + "\n")
        print_file.write("End mem usage (GB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3) + "\n")
        print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print_file.close()

        
        if verbose:
            print("Printing alignment with internal node sequences: ", output_prefix+".joint.aln")
        with open(output_prefix+".joint.aln", "w") as asr_output:
            for taxon in alignment:
                print(">" + taxon.id, file = asr_output)
                print(taxon.seq, file = asr_output)
            for i,node_index in enumerate(ancestral_node_order):
                taxon = ancestral_node_indices[node_index]
                asr_output.write('>' + taxon + '\n')
                asr_output.write(''.join(out_aln[:,i]) + '\n')

        # Release pool nodes
        pool.join()

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

def main_func(alignment, input_args):
    if input_args.base_patterns:
        base_pattern_bases_array, base_pattern_positions_array = get_base_patterns(alignment, verbose=True,
                                                                                   printero=input_args.print_file,
                                                                                   threads=input_args.threads,
                                                                                   fit_method=input_args.mp_method,
                                                                                   pickle_aln=input_args.pickle_aln)
        with open("base_patterns.npy", "wb") as f:
            numpy.save(f, base_pattern_bases_array)
        with open("base_positions.npy", "wb") as f:
            numpy.save(f, base_pattern_positions_array, allow_pickle=True)
    if input_args.jar:
        poly_aln = read_alignment(input_args.aln, "fasta", True)
        with open("base_patterns.npy", "rb") as f:
            base_pattern_bases_array = numpy.load(f, allow_pickle=True)
        with open("base_positions.npy", "rb") as f:
            base_pattern_positions_array = numpy.load(f, allow_pickle=True)


        jar(alignment=poly_aln,  # complete polymorphism alignment
            base_patterns=base_pattern_bases_array,  # array of unique base patterns in alignment
            base_pattern_positions=base_pattern_positions_array,
            # nparray of positions of unique base patterns in alignment
            tree_filename=input_args.tree,  # tree generated by model fit
            info_filename=input_args.info,  # file containing evolutionary model parameters
            info_filetype=input_args.model_fitter,
            # model fitter - format of file containing evolutionary model parameters
            output_prefix="temp_working_dir" + "/" + "ancestral_sequence_basename",  # output prefix
            threads=input_args.threads,  # number of cores to use
            verbose=True,
            printero=input_args.print_file,
            mp_metho=input_args.mp_method)


def get_args():
    parser = argparse.ArgumentParser(
        description='Debug the memory usage of the pyjar reconstruction get_base_patterns function ')
        

    parser.add_argument('--aln','-a',dest="aln",
                        help='Multifasta alignment file', required=True)
    parser.add_argument('--threads', '-t', dest="threads",
                        help='Number of threads to use (must be > 1)', type=int, required=True)
    parser.add_argument('--print-file', '-p', dest="print_file",
                        help='File to print debug statements to', type=str, required=True)
    parser.add_argument('--mp-method', '-m', dest="mp_method",
                        help='method to run the pool jobs with. Either spawn, fork or forkserver',
                         choices=['spawn','fork','forkserver'], required=True, type=str)
    parser.add_argument('--jar','-j', dest="jar",
                        default=False, action='store_true',
                        help="Whether or not to run the full jar function as well")
    parser.add_argument('--tree','-tr',dest="tree",
                        help="tree file for jar", type=str)
    parser.add_argument('--info','-i', dest="info",
                        help="log file for jar recon", type=str)
    parser.add_argument('--model-fitter','-mf',dest="model_fitter",
                        help="Name of tree model for jar",type=str)
    parser.add_argument('--base-patterns','-b',dest="base_patterns",
                        help="Run the base pattern reconstructions or not",
                        default=False, action="store_true")
    parser.add_argument('--pickle-aln','-pa', dest="pickle_aln",
                        help="Whether to write out the aln list of lists for further inspection",
                        default=False, action="store_true")


    return parser.parse_args()


if __name__ == '__main__':
    input_args = get_args()
    print("reading in the alignment")
    aln_read = read_alignment(input_args.aln, "fasta", True, True)
    print("running the gap inserter")
    main_func(aln_read, input_args)

