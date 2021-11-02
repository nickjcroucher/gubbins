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
from math import log, exp
from functools import partial
from scipy.sparse import csr_matrix
from numba import jit, njit, types, from_dtype
from numba.types import bool_
from numba.typed import Dict
import collections
try:
    from multiprocessing import Pool, shared_memory
    from multiprocessing.managers import SharedMemoryManager
    NumpyShared = collections.namedtuple('NumpyShared', ('name', 'shape', 'dtype'))
except ImportError as e:
    sys.stderr.write("This version of Gubbins requires the multiprocessing library and python v3.8 or higher for memory management\n")
    sys.exit(201)

from gubbins.utils import generate_shared_mem_array

####################################################
# Function to read an alignment in various formats #
####################################################

def read_alignment(filename, file_type, verbose=False):
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
    return alignmentObject

#Calculate Pij from Q matrix and branch length
def calculate_pij(branch_length,rate_matrix):
    if branch_length==0:
        pij = numpy.full((4,4), numpy.NINF, dtype = numpy.float32)
        numpy.fill_diagonal(pij, 0.0)
    else:
        pij = numpy.log(linalg.expm(numpy.multiply(branch_length,rate_matrix))) # modified
    return pij

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

# Read the RAxML info file to get rates and frequencies
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

def create_rate_matrix(f, r):
    #convert f and r to Q matrix
    rm=numpy.array([[0, f[0]*r[1], f[0]*r[2], f[0]*r[3]],[f[1]*r[0], 0, f[1]*r[3],f[1]*r[4]],[f[2]*r[1], f[2]*r[3], 0, f[2]*r[5]],[f[3]*r[2], f[3]*r[4], f[3]*r[5], 0]])
    
    rm[0][0]=numpy.sum(rm[0])*-1
    rm[1][1]=numpy.sum(rm[1])*-1
    rm[2][2]=numpy.sum(rm[2])*-1
    rm[3][3]=numpy.sum(rm[3])*-1
    
    return rm

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

@njit
def process_sequence(seq,seq_to_int):
    return numpy.array([seq_to_int[x] for x in seq], dtype = numpy.uint8)

# Based on https://stackoverflow.com/questions/21888406/getting-the-indexes-to-the-duplicate-columns-of-a-numpy-array
def unique_columns2(data):
    dt = numpy.dtype((numpy.void, data.dtype.itemsize * data.shape[0]))
    dataf = numpy.asfortranarray(data).view(dt)
    u,uind = numpy.unique(dataf, return_inverse=True)
    u = u.view(data.dtype).reshape(-1,data.shape[0]).T
    return (u,uind)

def get_base_patterns(alignment, verbose):
    if verbose:
        print("Finding unique base patterns")
    # Identify unique base patterns
    t1=time.process_time()
    # Set up data structures to convert characters to integers with numba
    np_unichar = numpy.dtype('<U5')
    UnicharType = from_dtype(np_unichar)
    seq_to_int = Dict.empty(
        key_type=UnicharType,
        value_type=types.uint8,
    )
    seq_to_int['A'] = numpy.uint8(0)
    seq_to_int['C'] = numpy.uint8(1)
    seq_to_int['G'] = numpy.uint8(2)
    seq_to_int['T'] = numpy.uint8(3)
    seq_to_int['-'] = numpy.uint8(4)
    seq_to_int['N'] = numpy.uint8(5)
    # Convert alignment to Numpy array
    align_array = numpy.array([process_sequence(numpy.array(record.seq, dtype = np_unichar),
                                                seq_to_int) for record in alignment],
                                dtype = numpy.uint8,
                                order='F')
    # Get unique base patterns and their indices in the alignment
    base_pattern_bases_array, base_pattern_positions_array = unique_columns2(align_array)
    base_pattern_positions_array_of_arrays = [numpy.where(base_pattern_positions_array==x)[0] for x in range(base_pattern_bases_array.shape[1])]
    # Convert the array of arrays into an ndarray that can be saved to shared memory
    square_base_pattern_positions_array = convert_to_square_numpy_array(base_pattern_positions_array_of_arrays)
    # Finish
    t2=time.process_time()
    if verbose:
        print("Time taken to find unique base patterns:", t2-t1, "seconds")
        print("Unique base patterns:", str(square_base_pattern_positions_array.shape[0]))
    return base_pattern_bases_array.transpose(), square_base_pattern_positions_array

@njit
def find_most_likely_base_given_descendents(Lmat, Cmat, pij, node_index, child_node_indices, column_base_indices):
    #2a. Lz(i) = maxj Pij(tz) x Lx(j) x Ly(j)
    #2b. Cz(i) = the value of j attaining the above maximum.
    for end_index in column_base_indices:
        c = Lmat[child_node_indices,end_index].sum()
        for start_index in column_base_indices:
            j = pij[start_index,end_index]+c
            if j > Lmat[node_index,start_index]:
                Lmat[node_index,start_index] = j
                Cmat[node_index,start_index] = end_index

@njit
def process_leaf(Lmat, Cmat, pij, node_index, taxon_base_index):
    if taxon_base_index < 4:
        #1a. Let j be the amino acid at y. Set, for each amino acid i: Cy(i)= j. This implies that no matter what is the amino acid in the father of y, j is assigned to node y.
        Cmat[node_index,:] = taxon_base_index
        #1b. Set for each amino acid i: Ly(i) = Pij(ty), where ty is the branch length between y and its father.
        Lmat[node_index,:] = [pij[0][taxon_base_index],pij[1][taxon_base_index],pij[2][taxon_base_index],pij[3][taxon_base_index]]
    else:
        # Cmat stays as default when base is unknown
        Lmat[node_index,:] = [pij[0][0],pij[1][1],pij[2][2],pij[3][3]]

@njit
def calculate_root_likelihood(Lmat, Cmat, base_frequencies, node_index, child_node_indices, column_base_indices):
    for end_index in column_base_indices:
        c = Lmat[child_node_indices,end_index].sum()
        for start_index in column_base_indices:
            j = log(base_frequencies[end_index]) + c
            if j > Lmat[node_index,start_index]:
                Lmat[node_index,start_index] = j
                Cmat[node_index,start_index] = end_index

def reconstruct_alignment_column(column_indices, tree = None, alignment_sequence_names = None, ancestral_node_indices = None, base_patterns = None, base_pattern_positions = None, base_matrix = None, base_frequencies = None, new_aln = None, threads = 1, verbose = False):
    
    ### TIMING
    if verbose:
        calc_time = 0.0
        conversion_time = 0.0
        writing_time = 0.0
        prep_time = 0.0
        prep_time_start = time.process_time()
        
    # Record SNPs reconstructed as occurring on each branch
    bases = frozenset(['A','C','G','T'])
    ordered_bases = numpy.array(['A','C','G','T','-'], dtype = numpy.unicode_)
    node_snps = {node.taxon.label:0 for node in tree.postorder_node_iter()}
    ancestrally_conserved = {b:list() for b in ['A','C','G','T']}
    ancestrally_variable = {b:{ancestral_node_indices[x]:list() for x in ancestral_node_indices} for b in range(5)}

    # Generate data structures for reconstructions
    num_nodes = len(tree.nodes())
    Lmat = numpy.full((num_nodes,4), numpy.NINF, dtype = numpy.float32)
    Cmat = numpy.full((num_nodes,4), [0,1,2,3], dtype = numpy.uint8)
    node_indices = {}
    reconstructed_base_indices = [-1]*num_nodes

    # Load base pattern information
    base_patterns_shm = shared_memory.SharedMemory(name = base_patterns.name)
    base_patterns = numpy.ndarray(base_patterns.shape, dtype = base_patterns.dtype, buffer = base_patterns_shm.buf)
    # Load base pattern position information
    base_pattern_positions_shm = shared_memory.SharedMemory(name = base_pattern_positions.name)
    base_pattern_positions = numpy.ndarray(base_pattern_positions.shape, dtype = base_pattern_positions.dtype, buffer = base_pattern_positions_shm.buf)
    
    # Extract information for iterations
    if threads == 1:
        columns = base_patterns
        column_positions = base_pattern_positions
    else:
        columns = base_patterns[column_indices]
        column_positions = base_pattern_positions[column_indices,:]

    ### TIMING
    if verbose:
        prep_time_end = time.process_time()
        prep_time = prep_time_end - prep_time_start
    
    # Iterate over columns
    for column,base_pattern_columns_padded in zip(columns,column_positions):
        
        # Start calculation time
        if verbose:
            calc_time_start = time.process_time()
        
        # Get column information
        base_pattern_columns = base_pattern_columns_padded[base_pattern_columns_padded > -1].tolist()
        unknown_base_count = numpy.count_nonzero(column > 4)
        column_base_indices = numpy.unique(column[numpy.where(column <= 3)])

        # Reset matrices
        Lmat.fill(numpy.NINF)
        Cmat[:] = [0,1,2,3]

        # Heuristic for speed: if all taxa are monomorphic, with a gap in only one sequence, then the ancestral states
        # will all be the observed base, as no ancestral node will have two child nodes with unknown bases at this site
        if unknown_base_count == 1 and column_base_indices.size == 1:
            
            ancestrally_conserved[ordered_bases[column_base_indices[0]]].extend(base_pattern_columns)
            
        else:
            # Otherwise perform a full ML inference
            #1 For each OTU y perform the following:

            #Visit a nonroot internal node, z, which has not been visited yet, but both of whose sons, nodes x and y, have already been visited, i.e., Lx(j), Cx(j), Ly(j), and Cy(j) have already been defined for each j. Let tz be the length of the branch connecting node z and its father. For each amino acid i, compute Lz(i) and Cz(i) according to the following formulae:
            
            #Denote the three sons of the root by x, y, and z. For each amino acid k, compute the expression Pk x Lx(k) x Ly(k) x Lz(k). Reconstruct r by choosing the amino acid k maximizing this expression. The maximum value found is the likelihood of the best reconstruction.
            
            # Create data structures for storing likelihoods
            
            for node_index,node in enumerate(tree.postorder_node_iter()):
                taxon=str(node.taxon.label)
                node_indices[taxon] = node_index
                if node.parent_node==None:
                    continue
                #calculate the transistion matrix for the branch
                pij=node.pij
                if node.is_leaf():
                    try:
                        alignment_index = alignment_sequence_names[taxon]
                        taxon_base_index = column[alignment_index]
                        process_leaf(Lmat, Cmat, pij, node_index, taxon_base_index)
                    except KeyError:
                        print("Cannot find", taxon, "in alignment")
                        sys.exit(208)
                else:
                    child_node_indices = [node_indices[child.taxon.label] for child in node.child_node_iter()]
                    #2a. Lz(i) = maxj Pij(tz) x Lx(j) x Ly(j)
                    #2b. Cz(i) = the value of j attaining the above maximum.
                    find_most_likely_base_given_descendents(Lmat, Cmat, pij, node_index, numpy.array(child_node_indices), column_base_indices)
            

            # Calculate likelihood of base at root node
            child_node_indices = [node_indices[child.taxon.label] for child in node.child_node_iter()]
            calculate_root_likelihood(Lmat, Cmat, base_frequencies, node_index, numpy.array(child_node_indices), column_base_indices)
            max_root_base_index = Cmat[node_index,numpy.argmax(Lmat[node_index,:])]
            reconstructed_base_indices[node_index] = max_root_base_index

            #Traverse the tree from the root in the direction of the OTUs, assigning to each node its most likely ancestral character as follows:
            for node in tree.preorder_node_iter():
                node_label = node.taxon.label
                node_index = node_indices[node_label]
                try:
                    #5a. Visit an unreconstructed internal node x whose father y has already been reconstructed. Denote by i the reconstructed amino acid at node y.
                    parent_node_index = node_indices[node.parent_node.taxon.label]
                    i = reconstructed_base_indices[parent_node_index]
                except AttributeError:
                    continue
                #5b. Reconstruct node x by choosing Cx(i).
                reconstructed_base_indices[node_index] = Cmat[node_index,i]

            ### TIMING
            if verbose:
                calc_time_end = time.process_time()
                calc_time += (calc_time_end - calc_time_start)
                conversion_time_start = time.process_time()
            
            # Put gaps back in and check that any ancestor with only gaps downstream is made a gap
            # store reconstructed alleles
            reconstructed_alleles = {}
            for node in tree.postorder_node_iter():
                node_label = node.taxon.label
                node_index = node_indices[node_label]
                if node.is_leaf():
                    alignment_index = alignment_sequence_names[node_label]
                    reconstructed_alleles[node_label] = column[alignment_index]
                else:
                    has_child_base = False
                    child_taxon_labels = [child.taxon.label for child in node.child_node_iter()for child in node.child_node_iter()]
                    for child_taxon_label in child_taxon_labels:
                        if reconstructed_alleles[child_taxon_label] < 4:
                            has_child_base = True
                    if has_child_base:
                        reconstructed_alleles[node_label] = reconstructed_base_indices[node_index]
                    else:
                        reconstructed_alleles[node_label] = 4
        
            # If site is monomorphic - replace whole column; else replace specific entries
            reconstructed_allele_set = set(reconstructed_alleles.values())
            if len(reconstructed_allele_set) == 1:
                ancestrally_conserved[reconstructed_allele_set.pop()].extend(base_pattern_columns)
            else:
                for taxon in ancestral_node_indices:
                    ancestrally_variable[reconstructed_alleles[taxon]][ancestral_node_indices[taxon]].extend(base_pattern_columns)
            
            # iterate through tree
            for node in tree.preorder_node_iter():
                try:
                    if node.r in bases and node.parent_node.r in bases and node.r!=node.parent_node.r:
                        node_snps[node.taxon.label] += len(base_pattern_columns)
                except AttributeError:
                    continue
            
            ### TIMING
            if verbose:
                conversion_time_end = time.process_time()
                conversion_time += (conversion_time_end - conversion_time_start)

    ### TIMING
    if verbose:
        writing_time_start = time.process_time()

    # combine results across columns to access shared memory object as few times as possible
    # load output alignment
    out_aln_shm = shared_memory.SharedMemory(name = new_aln.name)
    out_aln = numpy.ndarray(new_aln.shape, dtype = new_aln.dtype, buffer = out_aln_shm.buf)

    for b in ancestrally_conserved:
        if len(ancestrally_conserved[b]) > 0:
            out_aln[ancestrally_conserved[b],:] = b
    for b_index in ancestrally_variable:
        for index in ancestrally_variable[b_index]:
            if len(ancestrally_variable[b_index][index]) > 0:
                out_aln[ancestrally_variable[b_index][index],index] = ordered_bases[b_index]

    # Close shared memory
    out_aln_shm.close()
    base_patterns_shm.close()
    base_pattern_positions_shm.close()

    ### TIMING
    if verbose:
        writing_time_end = time.process_time()
        writing_time += (writing_time_end - writing_time_start)
        print('Time for JAR preparation:\t' + str(prep_time))
        print('Time for JAR calculation:\t' + str(calc_time))
        print('Time for JAR conversion:\t' + str(conversion_time))
        print('Time for JAR writing:\t' + str(writing_time))

    return node_snps

# from https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length/37414115#37414115
def chunks(l, k):
    n = len(l)
    return [l[i * (n // k) + min(i, n % k):(i+1) * (n // k) + min(i+1, n % k)] for i in range(k)]

def jar(alignment = None, base_patterns = None, base_pattern_positions = None, tree_filename = None, info_filename = None, info_filetype = None, output_prefix = None, threads = 1, verbose = False):

    # Lookup for each base
    mb={"A": 0, "C": 1, "G": 2, "T": 3}

    # Create a new alignment for the output containing all taxa in the input alignment
    alignment_sequence_names = {}
    ancestral_node_names = []
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
    nodecounter=0
    
    for node in tree.preorder_node_iter():
        if node.taxon == None:
            nodecounter+=1
            nodename="Node_"+str(nodecounter)
            tree.taxon_namespace.add_taxon(dendropy.Taxon(nodename))
            node.taxon=tree.taxon_namespace.get_taxon(nodename)
            if nodename in alignment_sequence_names:
                print(nodename, "already in alignment. Quitting")
                sys.exit(209)
            ancestral_node_names.append(nodename) # index for reconstruction
        else:
            node.taxon.label = node.taxon.label.strip("'")
        if node.parent_node != None:
            node.pij=calculate_pij(node.edge_length, rm)

    # Create new empty array
    new_aln_array = numpy.full((len(alignment[0]),len(ancestral_node_names)), '?', dtype = numpy.unicode_)

    # Index names for reconstruction
    ancestral_node_indices = {name:i for i,name in enumerate(ancestral_node_names)}

    # Reconstruct each base position
    if verbose:
        print("Reconstructing sites on tree")
    
    with SharedMemoryManager() as smm:
    
        # Convert alignment to shared memory numpy array
        new_aln_shared_array = generate_shared_mem_array(new_aln_array, smm)

        # Convert base patterns to shared memory numpy array
        base_patterns_shared_array = generate_shared_mem_array(base_patterns, smm)

        # Convert base pattern positions to shared memory numpy array
        base_pattern_positions_shared_array = generate_shared_mem_array(base_pattern_positions, smm)

        # split list of sites into chunks per core
        bp_list = list(range(len(base_patterns)))
        base_pattern_indices = list(chunks(bp_list,threads))

        # Parallelise reconstructions across alignment columns using multiprocessing
        with Pool(processes = threads) as pool:
            reconstruction_results = pool.map(partial(
                                        reconstruct_alignment_column,
                                            tree = tree,
                                            alignment_sequence_names = alignment_sequence_names,
                                            ancestral_node_indices = ancestral_node_indices,
                                            base_patterns = base_patterns_shared_array,
                                            base_pattern_positions = base_pattern_positions_shared_array,
                                            base_matrix = mb,
                                            base_frequencies = f,
                                            new_aln = new_aln_shared_array,
                                            threads = threads,
                                            verbose = verbose),
                                        base_pattern_indices
                                    )

        # Write out alignment while shared memory manager still active
        out_aln_shm = shared_memory.SharedMemory(name = new_aln_shared_array.name)
        out_aln = numpy.ndarray(new_aln_array.shape, dtype = new_aln_array.dtype, buffer = out_aln_shm.buf)
        
        if verbose:
            print("Printing alignment with internal node sequences: ", output_prefix+".joint.aln")
        with open(output_prefix+".joint.aln", "w") as asr_output:
            for taxon in alignment:
                print(">" + taxon.id, file = asr_output)
                print(taxon.seq, file=asr_output)
            for taxon in ancestral_node_indices:
                print(">" + taxon, file = asr_output)
                print(''.join(out_aln[:,ancestral_node_indices[taxon]]), file=asr_output)

        # Release pool nodes
        pool.join()

        # Combine results for each base across the alignment
        for node in tree.preorder_node_iter():
            node.edge_length = 0.0 # reset lengths to convert to SNPs
            for x in range(len(reconstruction_results)):
                try:
                    node.edge_length += reconstruction_results[x][node.taxon.label];
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
