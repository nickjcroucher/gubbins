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
import collections
try:
    from multiprocessing import Pool, shared_memory
    from multiprocessing.managers import SharedMemoryManager
    NumpyShared = collections.namedtuple('NumpyShared', ('name', 'shape', 'dtype'))
except ImportError as e:
    sys.stderr.write("This version of Gubbins requires python v3.8 or higher\n")
    sys.exit(0)

####################################################
# Function to read an alignment in various formats #
####################################################

def read_alignment(filename, file_type, verbose=False):
    if not os.path.isfile(filename):
        print("Error: alignment file " + filename + " does not exist")
        sys.exit()
    if verbose:
        print("Trying to open file " + filename + " as " + file_type)
    try:
        alignmentObject = AlignIO.read(open(filename), file_type)
        if verbose:
            print("Alignment read successfully")
    except:
        print("Cannot open alignment file " + filename + " as " + file_type)
        sys.exit()
    return alignmentObject

#Calculate Pij from Q matrix and branch length
def calculate_pij(branch_length,rate_matrix):
    if branch_length==0:
        return numpy.array([[1, 0, 0, 0,], [0, 1, 0, 0,], [0, 0, 1, 0,], [0, 0, 0, 1,]])
    else:
        return numpy.log(linalg.expm(numpy.multiply(branch_length,rate_matrix)))

#Read the tree file and root
def read_tree(treefile):
    if not os.path.isfile(treefile):
        print("Error: tree file does not exist")
        sys.exit()
    t=dendropy.Tree.get(path=treefile, schema="newick", preserve_underscores=True, rooting="force-rooted")
    return t

# Read the RAxML info file to get rates and frequencies
def read_info(infofile, type = 'raxml'):

    if not os.path.isfile(infofile):
        print("Error: model information file " + infofile + " does not exist")
        sys.exit()
    
    if type not in ['raxml','iqtree','fasttree']:
        sys.stderr.write('Only able to parse GTR-type models from raxml, iqtree ot fasttree')
        sys.exit()
    
    r=[-1.0] * 6 # initialise rates
    f=[-1.0] * 4 # initialise frequencies
    
    for line in open(infofile, "r"):
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

    # Check frequencies and rates have been extracted correctly
    if -1.0 in f or -1.0 in r:
        sys.stderr.write('Problem with extracting model parameters - frequencies are ' + str(f) + ' and rates are ' + str(r))
        sys.exit()

    return f, r

def create_rate_matrix(f, r):
    #convert f and r to Q matrix
    rm=numpy.array([[0, f[0]*r[1], f[0]*r[2], f[0]*r[3]],[f[1]*r[0], 0, f[1]*r[3],f[1]*r[4]],[f[2]*r[1], f[2]*r[3], 0, f[2]*r[5]],[f[3]*r[2], f[3]*r[4], f[3]*r[5], 0]])
    
    rm[0][0]=numpy.sum(rm[0])*-1
    rm[1][1]=numpy.sum(rm[1])*-1
    rm[2][2]=numpy.sum(rm[2])*-1
    rm[3][3]=numpy.sum(rm[3])*-1
    
    return rm

def get_base_patterns(alignment, verbose):
    if verbose:
        print("Finding unique base patterns")
    base_patterns={}
    t1=time.process_time()
    for x in range(len(alignment[0])):
        try:
            base_patterns[alignment[:,x]].append(x)
        except KeyError:
            base_patterns[alignment[:,x]]=[x]
    t2=time.process_time()
    if verbose:
        print("Time taken to find unique base patterns:", t2-t1, "seconds")
        print("Unique base patterns:", len(base_patterns))
    return base_patterns

def reconstruct_alignment_column(columns, tree = None, alignment_sequence_names = None, ancestral_node_indices = None, base_patterns = None, base_matrix = None, base_frequencies = None, new_aln = None, verbose = False):
    
    ### TIMING
    if verbose:
        calc_time = 0.0
        storage_time = 0.0
        writing_time = 0.0
        prep_time_start = time.process_time()
    
    # process bases for alignment column
    bases = frozenset(["A", "C", "G", "T"])

    # Record SNPs reconstructed as occurring on each branch
    node_snps = {node.taxon.label:0 for node in tree.postorder_node_iter()}
    ancestrally_conserved = {b:list() for b in ["A", "C", "G", "T", "-"]}
    ancestrally_variable = {b:{ancestral_node_indices[x]:list() for x in ancestral_node_indices} for b in ["A", "C", "G", "T", "-"]}

    ### TIMING
    if verbose:
        prep_time_end = time.process_time()
        prep_time = prep_time_end - prep_time_start

    # Iterate over columns
    for column in columns:
    
        ### TIMING
        if verbose:
            calc_time_start = time.process_time()
        
        base_pattern_columns = base_patterns[column]
        
        columnbases=set([])
        base={}
        unknown_base_count = 0
        for i, y in enumerate(column):
            base[alignment_sequence_names[i]]=y
            if y in bases:
                columnbases.add(y)
            else:
                unknown_base_count = unknown_base_count + 1
        
        # Heuristic for speed: if all taxa are monomorphic, with a gap in only one sequence, then the ancestral states
        # will all be the observed base, as no ancestral node will have two child nodes with unknown bases at this site
        if unknown_base_count == 1 and len(columnbases) == 1:
            
            ancestrally_conserved[columnbases.pop()].extend(base_patterns[column])
            
        else:
            # Otherwise perform a full ML inference
            #1 For each OTU y perform the following:

            #Visit a nonroot internal node, z, which has not been visited yet, but both of whose sons, nodes x and y, have already been visited, i.e., Lx(j), Cx(j), Ly(j), and Cy(j) have already been defined for each j. Let tz be the length of the branch connecting node z and its father. For each amino acid i, compute Lz(i) and Cz(i) according to the following formulae:
            
            #Denote the three sons of the root by x, y, and z. For each amino acid k, compute the expression Pk x Lx(k) x Ly(k) x Lz(k). Reconstruct r by choosing the amino acid k maximizing this expression. The maximum value found is the likelihood of the best reconstruction.
            for node in tree.postorder_node_iter():
                if node.parent_node==None:
                    continue
                #calculate the transistion matrix for the branch
                pij=node.pij
                
                if node.is_leaf():
                    taxon=str(node.taxon.label).strip("'")
                    try:
                        if base[taxon] in ["A", "C", "G", "T"]:
                            #1a. Let j be the amino acid at y. Set, for each amino acid i: Cy(i)= j. This implies that no matter what is the amino acid in the father of y, j is assigned to node y.
                            node.C={"A": base[taxon], "C": base[taxon], "G": base[taxon], "T": base[taxon]}
                        
                            #1b. Set for each amino acid i: Ly(i) = Pij(ty), where ty is the branch length between y and its father.
                            node.L={"A": pij[base_matrix["A"]][base_matrix[base[taxon]]], "C": pij[base_matrix["C"]][base_matrix[base[taxon]]], "G": pij[base_matrix["G"]][base_matrix[base[taxon]]], "T": pij[base_matrix["T"]][base_matrix[base[taxon]]]}
                        else:
                            
                            node.C={"A": "A", "C": "C", "G": "G", "T": "T"}
                            node.L={"A": pij[base_matrix["A"]][base_matrix["A"]], "C": pij[base_matrix["C"]][base_matrix["C"]], "G": pij[base_matrix["G"]][base_matrix["G"]], "T": pij[base_matrix["T"]][base_matrix["T"]]}
                        
                    except KeyError:
                        print("Cannot find", taxon, "in base")
                        sys.exit()
                
                else:
                    node.L={}
                    node.C={}
                    
                    #2a. Lz(i) = maxj Pij(tz) x Lx(j) x Ly(j)
                    #2b. Cz(i) = the value of j attaining the above maximum.
                    
                    for basenum in columnbases:
                        node.L[basenum]=float("-inf")
                        node.C[basenum]=None
                    
                    for end in columnbases:
                        c=0.0
                        for child in node.child_node_iter():
                            c+=child.L[end]
                        for start in columnbases:
                            j=pij[base_matrix[start],base_matrix[end]]+c
                            
                            if j>node.L[start]:
                                node.L[start]=j
                                node.C[start]=end

            node.L={}
            node.C={}
            for basenum in columnbases:
                node.L[basenum]=float("-inf")
                node.C[basenum]=None
            for end in columnbases:
                c=0
                for child in node.child_node_iter():
                    c+=child.L[end]
                for start in columnbases:
                    j=log(base_frequencies[base_matrix[end]])+c

                    if j>node.L[start]:
                        node.L[start]=j
                        node.C[start]=end
                
            max_root_base=None
            max_root_base_likelihood=float("-inf")
            for root_base in columnbases:
                if node.L[root_base]>max_root_base_likelihood:
                    max_root_base_likelihood=node.L[root_base]
                    max_root_base=node.C[root_base]
            node.r=max_root_base
            
            #Traverse the tree from the root in the direction of the OTUs, assigning to each node its most likely ancestral character as follows:
            for node in tree.preorder_node_iter():
            
                try:
                    #5a. Visit an unreconstructed internal node x whose father y has already been reconstructed. Denote by i the reconstructed amino acid at node y.
                    i=node.parent_node.r
                except AttributeError:
                    continue
                #5b. Reconstruct node x by choosing Cx(i).
                node.r=node.C[i]

            rootlens=[]
            for child in tree.seed_node.child_node_iter():
                rootlens.append([child.edge_length,child,child.r])
            rootlens.sort(key = lambda x: x[0])
            tree.seed_node.r=rootlens[-1][1].r
            
            ### TIMING
            if verbose:
                calc_time_end = time.process_time()
                calc_time += (calc_time_end - calc_time_start)
                storage_time_start = time.process_time()

            # Put gaps back in and check that any ancestor with only gaps downstream is made a gap
            # store reconstructed alleles
            reconstructed_alleles = {}
            for node in tree.postorder_node_iter():
                if node.is_leaf():
                    node.r=base[node.taxon.label]
                else:
                    has_child_base=False
                    for child in node.child_node_iter():
                        if child.r in bases:
                            has_child_base=True
                            break
                    if not has_child_base:
                        node.r = "-"
                    # Store reconstructed allele to determine how it should be inserted into the new alignment
                    reconstructed_alleles[node.taxon.label] = node.r
        
            # If site is monomorphic - replace whole column; else replace specific entries
            reconstructed_allele_set = set(reconstructed_alleles.values())
            if len(reconstructed_allele_set) == 1:
                ancestrally_conserved[reconstructed_allele_set.pop()].extend(base_patterns[column])
            else:
                for taxon in reconstructed_alleles:
                    ancestrally_variable[reconstructed_alleles[taxon]][ancestral_node_indices[taxon]].extend(base_patterns[column])
            
            # iterate through tree
            for node in tree.preorder_node_iter():
                try:
                    if node.r in bases and node.parent_node.r in bases and node.r!=node.parent_node.r:
                        node_snps[node.taxon.label] += len(base_pattern_columns)
                except AttributeError:
                    continue

    ### TIMING
    if verbose:
        storage_time_end = time.process_time()
        storage_time += (storage_time_end - storage_time_start)
        writing_time_start = time.process_time()

    # combine results across columns to access shared memory object as few times as possible
    # load output alignment
    out_aln_shm = shared_memory.SharedMemory(name = new_aln.name)
    out_aln = numpy.ndarray(new_aln.shape, dtype = new_aln.dtype, buffer = out_aln_shm.buf)

    for b in ["A", "C", "G", "T", "-"]:
        if len(ancestrally_conserved[b]) > 0:
            out_aln[ancestrally_conserved[b],:] = b
        for index in ancestrally_variable[b]:
            if len(ancestrally_variable[b][index]) > 0:
                out_aln[ancestrally_variable[b][index],index] = b

    out_aln_shm.close()

    ### TIMING
    if verbose:
        writing_time_end = time.process_time()
        writing_time += (writing_time_end - writing_time_start)
        print('Time for JAR preparation:\t' + str(prep_time))
        print('Time for JAR calculation:\t' + str(calc_time))
        print('Time for JAR storage:\t' + str(storage_time))
        print('Time for JAR writing:\t' + str(writing_time))

    return node_snps

# from https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length/37414115#37414115
def chunks(l, k):
    n = len(l)
    return [l[i * (n // k) + min(i, n % k):(i+1) * (n // k) + min(i+1, n % k)] for i in range(k)]

def jar(alignment = None, base_patterns = None, tree_filename = None, info_filename = None, info_filetype = None, output_prefix = None, threads = 1, verbose = False):

    # Lookup for each base
    mb={"A": 0, "C": 1, "G": 2, "T":3 }

    # Create a new alignment for the output containing all taxa in the input alignment
    alignment_sequence_names = []
    ancestral_node_names = []
    for i, x in enumerate(alignment):
        alignment_sequence_names.append(x.id)
    
    # Read the tree
    if verbose:
        print("Reading tree file:", tree_filename)
    tree=read_tree(tree_filename)
    
    # Read the info file and get frequencies and rates
    if info_filename!="":
        if verbose:
            print("Reading info file:", info_filename)
        f, r=read_info(info_filename, type = info_filetype)
    else:
        if verbose:
            print("Using default JC rates and frequencies")
        f=[0.25,0.25,0.25,0.25]
        r=[1.0,1.0,1.0,1.0,1.0,1.0]
    
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
                sys.exit()
            ancestral_node_names.append(nodename) # index for reconstruction
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
    
        # Declare global
        global all_base_patterns
        all_base_patterns = base_patterns
    
        # Convert alignment to shared memory numpy array
        new_aln_array_raw = smm.SharedMemory(size = new_aln_array.nbytes)
        new_aln_shared_array = numpy.ndarray(new_aln_array.shape, dtype = new_aln_array.dtype, buffer = new_aln_array_raw.buf)
        new_aln_shared_array[:] = new_aln_array[:]
        new_aln_shared_array = NumpyShared(name = new_aln_array_raw.name, shape = new_aln_array.shape, dtype = new_aln_array.dtype)

        # split list of sites into chunks per core
        bp_list = list(base_patterns.keys())
        base_pattern_lists = list(chunks(bp_list,threads))

        # Parallelise reconstructions across alignment columns using multiprocessing
        with Pool(processes = threads) as pool:
            reconstruction_results = pool.map(partial(
                                        reconstruct_alignment_column,
                                            tree = tree,
                                            alignment_sequence_names = alignment_sequence_names,
                                            ancestral_node_indices = ancestral_node_indices,
                                            base_patterns = base_patterns,
                                            base_matrix = mb,
                                            base_frequencies = f,
                                            new_aln = new_aln_shared_array,
                                            verbose = verbose),
                                        base_pattern_lists
                                    )

        # Write out alignment while shared memory manager still active
        out_aln_shm = shared_memory.SharedMemory(name = new_aln_shared_array.name)
        out_aln = numpy.ndarray(new_aln_array.shape, dtype = new_aln_array.dtype, buffer = out_aln_shm.buf)
        if verbose:
            print("Printing alignment with internal node sequences: ", output_prefix+".joint.aln")
        asr_output = open(output_prefix+".joint.aln", "w")
        for taxon in alignment:
            print(">"+taxon.id, file=asr_output)
            print(taxon.seq, file=asr_output)
        for taxon in ancestral_node_indices:
            print(">"+taxon, file=asr_output)
            print(''.join(out_aln[:,ancestral_node_indices[taxon]]), file=asr_output)
        asr_output.close()

        # Combine results for each base across the alignment
        for node in tree.preorder_node_iter():
            node.edge_length = 0.0 # reset lengths to convert to SNPs
            for x in range(len(reconstruction_results)):
                try:
                    node.edge_length += reconstruction_results[x][node.taxon.label];
                except AttributeError:
                    continue

        # Print tree
        if verbose:
            print("Printing tree with internal nodes labelled: ", output_prefix+".joint.tre")
        tree_output=open(output_prefix+".joint.tre", "w")
        print(tree.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True, suppress_internal_node_labels=True).replace("'",""), file=tree_output)
        tree_output.close()
        
    if verbose:
        print("Done")
