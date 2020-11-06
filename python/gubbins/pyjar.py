#!/usr/bin/env python3

from scipy import linalg
import numpy
import dendropy
import sys
import os
import time
from Bio import AlignIO
from math import log, exp
from argparse import ArgumentParser

##########################################
# Function to Get command line arguments #
##########################################


def getargv():

        usage = "usage: %(prog)s [options]"
        description = "A python implementation of the joint ancestral state reconstruction algorithm of Pupko et al."
        epilog = "For the method, please cite \"Tal Pupko, Itsik Pe, Ron Shamir, Dan Graur; A Fast Algorithm for Joint Reconstruction of Ancestral Amino Acid Sequences, Molecular Biology and Evolution, Volume 17, Issue 6, 1 June 2000, Pages 890-896\""
        parser = ArgumentParser(usage=usage, description=description, epilog=epilog)
        
        parser.add_argument("-a", "--alignment", action="store", dest="alignment", help="Input alignment file. Required.", default="", metavar="FILE", required=True)
        parser.add_argument("-i", "--info", action="store", dest="info", help="Input RAxML info file. Optional. By default a JC model will be applied.", default="", metavar="FILE")
        parser.add_argument("-t", "--tree", action="store", dest="tree", help="Input tree file. Required.", default="", metavar="FILE", required=True)
        parser.add_argument("-o", "--output_prefix", action="store", dest="prefix", help="Output prefix. Required.", default="", required=True)
        parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="More verbose output", default=False)
        

        return parser.parse_args()


####################################################
# Function to read an alignment in various formats #
####################################################

def read_alignment(filename, verbose=False):

    if not os.path.isfile(filename):
        print("Error: alignment file does not exist")
        sys.exit()

    filetype=["phylip", "phylip-relaxed", "fasta", "clustal", "nexus", "emboss", "stockholm", "fasta-m10", "ig"]

    if filename.split(".")[-1].lower() in filetype:
        guesstype=filename.split(".")[-1].lower()
    elif filename.split(".")[-1].lower() in ["phy"]:
        guesstype="phylip"
    elif filename.split(".")[-1].lower() in ["fna", "dna", "aa", "aln", "fas"]:
        guesstype="fasta"
    elif filename.split(".")[-1].lower() in ["nxs", "nex", "nexus"]:
        guesstype="nexus"
    else:
        guesstype=""
    
    readok=False
    
    
    if guesstype!="":
        if verbose:
            print("Guessing file is in "+guesstype+" format")
            print("Trying to open file "+filename+" as "+guesstype)
        try:
            alignmentObject = AlignIO.read(open(filename, "r"), guesstype)
        except:
            print("Cannot open alignment file as "+guesstype)
        else:
            readok=True
            
        filetype.remove(guesstype)

    x=0

    while readok==False and x<len(filetype):
        if verbose:
            print("Trying to open file "+filename+" as "+filetype[x])

        try:
            alignmentObject = AlignIO.read(open(filename), filetype[x])
        except:
            print("Cannot open alignment file "+filename+" as "+filetype[x])
        else:
            readok=True

        x=x+1

    if readok==False:
        print("Failed to open alignment")
        sys.exit()
    else:
        if verbose:
            print("Alignment read successfully")
        return alignmentObject

#Calculate Pij from Q matrix and branch length
def calculate_pij(branch_length,rate_matrix):
    #print(linalg.expm(numpy.multiply(branch_length,rate_matrix)))
    if branch_length==0:
        return numpy.array([[1, 0, 0, 0,], [0, 1, 0, 0,], [0, 0, 1, 0,], [0, 0, 0, 1,]])
    else:
        return numpy.log(linalg.expm(numpy.multiply(branch_length,rate_matrix)))
        

#Read the tree file and root
def read_tree(treefile):
    if not os.path.isfile(treefile):
        print("Error: alignment file does not exist")
        sys.exit()
    t=dendropy.Tree.get(path=treefile, schema="newick", preserve_underscores=True, rooting="force-rooted")
# not for gubbins
#    t.reroot_at_midpoint()
#    t.update_bipartitions()
    return t


#Read the RAxML info file to get rates and frequencies
def read_info(infofile):
    if not os.path.isfile(infofile):
        print("Error: alignment file does not exist")
        sys.exit()
    r=[]
    f=[]
    for line in open(infofile, "r"):
        line=line.strip()
        if "freq pi" in line:
            words=line.split()
            f.append(float(words[2]))
        elif "Base frequencies:" in line:
            words=line.split()
            f=[float(words[2]), float(words[3]), float(words[4]), float(words[5])]
        elif "<->" in line:
            # order is ac ag at cg ct gt
            words=line.split()
            r.append(float(words[4]))
        elif "alpha[0]:" in line:
            # order is ac ag at cg ct gt
            words=line.split()
            r=[float(words[9]), float(words[10]), float(words[11]), float(words[12]), float(words[13]), float(words[14])]
    return f, r

def create_rate_matrix(f, r):
    #convert f and r to Q matrix
    rm=numpy.array([[0, f[0]*r[1], f[0]*r[2], f[0]*r[3]],[f[1]*r[0], 0, f[1]*r[3],f[1]*r[4]],[f[2]*r[1], f[2]*r[3], 0, f[2]*r[5]],[f[3]*r[2], f[3]*r[4], f[3]*r[5], 0]])
    
    rm[0][0]=numpy.sum(rm[0])*-1
    rm[1][1]=numpy.sum(rm[1])*-1
    rm[2][2]=numpy.sum(rm[2])*-1
    rm[3][3]=numpy.sum(rm[3])*-1
    
    return rm

def jar(alignment_filename, tree_filename, info_filename, output_prefix, verbose=False):
    
    #Lookup for each base
    mb={"A": 0, "C": 1, "G": 2, "T":3 }
    
    # read the alignment
    if verbose:
        print("Reading alignment file:", alignment_filename)
    alignment=read_alignment(alignment_filename, verbose)
    if verbose:
        print("Alignment size:", len(alignment), "taxa and", len(alignment[0]), "sites")
    
    #Create a new alignment for the output containing all taxa in the input alignment
    new_alignment={}
    for i, x in enumerate(alignment):
        new_alignment[x.id]=list(str(x.seq))
    
    # read the tree
    if verbose:
        print("Reading tree file:", tree_filename)
    tree=read_tree(tree_filename)
    
    
    #read the info file and get frequencids and rates
    if info_filename!="":
        if verbose:
            print("Reading info file:", info_filename)
        f, r=read_info(info_filename)
    else:
        if verbose:
            print("Using default JC rates and frequencies")
        f=[0.25,0.25,0.25,0.25]
        r=[1.0,1.0,1.0,1.0,1.0,1.0]
    
    if verbose:
        print("Frequencies:", ", ".join(map(str,f)))
        print("Rates:", ", ".join(map(str,r)))
    
    #create rate matrix from f and r
    rm=create_rate_matrix(f,r)
    
    #label internal nodes in tree and add these to the new alignment and calculate pij per non-root branch
    nodecounter=0
    for node in tree.preorder_node_iter():
    
        if node.taxon==None:
            nodecounter+=1
            nodename="Node_"+str(nodecounter)
            tree.taxon_namespace.add_taxon(dendropy.Taxon(nodename))
            node.taxon=tree.taxon_namespace.get_taxon(nodename)
            if nodename in new_alignment:
                print(nodename, "already in alignment. Quitting")
                sys.exit()
            new_alignment[nodename]=["?"]*len(alignment[0])
        if node.parent_node!=None:
            node.pij=calculate_pij(node.edge_length, rm)
            
        node.snps=0;
    
    
    #Find unique base patterns to speed up calculations
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
    
    
    allbases=set(["A", "C", "G", "T"])
    
    if verbose:
        print("Reconstructing sites on tree")
    
    # For each base
    t=0
    timestart=time.process_time()
    onetime=0.0
    twotime=0.0
    threetime=0.0
    for x, column in enumerate(base_patterns):
        t=t+1
        if t==1000:
            timeend=time.process_time()
            if verbose:
                print("Reconstructed", x+1, "of", len(base_patterns), "base patterns at", (timeend-timestart)/(x+1)*1000, "seconds per 1000 patterns")
            t=0
            #break
            
        
        columnbases=set([])
        base={}
        for i, y in enumerate(column):
            base[alignment[i].id]=y
            if y in allbases:
                columnbases.add(y)
        
        #1 For each OTU y perform the following:
    
        #Visit a nonroot internal node, z, which has not been visited yet, but both of whose sons, nodes x and y, have already been visited, i.e., Lx(j), Cx(j), Ly(j), and Cy(j) have already been defined for each j. Let tz be the length of the branch connecting node z and its father. For each amino acid i, compute Lz(i) and Cz(i) according to the following formulae:
        
        #Denote the three sons of the root by x, y, and z. For each amino acid k, compute the expression Pk x Lx(k) x Ly(k) x Lz(k). Reconstruct r by choosing the amino acid k maximizing this expression. The maximum value found is the likelihood of the best reconstruction.
        t=time.process_time()
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
                        node.L={"A": pij[mb["A"]][mb[base[taxon]]], "C": pij[mb["C"]][mb[base[taxon]]], "G": pij[mb["G"]][mb[base[taxon]]], "T": pij[mb["T"]][mb[base[taxon]]]}
                    else:
                        
                        node.C={"A": "A", "C": "C", "G": "G", "T": "T"}
                        node.L={"A": pij[mb["A"]][mb["A"]], "C": pij[mb["C"]][mb["C"]], "G": pij[mb["G"]][mb["G"]], "T": pij[mb["T"]][mb["T"]]}
                    
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
                        j=pij[mb[start],mb[end]]+c
                        
                        
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
                j=log(f[mb[end]])+c
    
                if j>node.L[start]:
                    node.L[start]=j
                    node.C[start]=end
            
        
        max_root_base=None
        max_root_base_likelihood=float("-inf")
        for root_base in columnbases:
            #print max_root_base, max_root_base_likelihood, root_base, node.L[root_base]
            if node.L[root_base]>max_root_base_likelihood:
                max_root_base_likelihood=node.L[root_base]
                max_root_base=node.C[root_base]
        node.r=max_root_base
        
        
        twoend=time.process_time()
        twotime+=twoend-t
        
        #Traverse the tree from the root in the direction of the OTUs, assigning to each node its most likely ancestral character as follows:
        threestart=time.process_time()
        for node in tree.preorder_node_iter():
        
            try:
                #5a. Visit an unreconstructed internal node x whose father y has already been reconstructed. Denote by i the reconstructed amino acid at node y.
                i=node.parent_node.r
            except AttributeError:
                continue
            #5b. Reconstruct node x by choosing Cx(i).
            node.r=node.C[i]
            #new_alignment[node.taxon.label].append(node.r)
    
        rootlens=[]
        for child in tree.seed_node.child_node_iter():
            rootlens.append([child.edge_length,child,child.r])
        rootlens.sort()
        tree.seed_node.r=rootlens[-1][1].r

    
        threeend=time.process_time()
        threetime+=threeend-threestart
        
        #Put gaps back in and check that any ancestor with only gaps downstream is made a gap
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                node.r=base[node.taxon.label]
            else:
                has_child_base=False
                for child in node.child_node_iter():
                    if child.r in allbases:
                        has_child_base=True
                        break
                if not has_child_base:
                    node.r="-"
                for bp in base_patterns[column]:
                    new_alignment[node.taxon.label][bp]=node.r
        
        
        for node in tree.preorder_node_iter():
            try:
                if node.r in ["A", "C", "G", "T"] and node.parent_node.r in ["A", "C", "G", "T"] and node.r!=node.parent_node.r:
                    node.snps+=len(base_patterns[column])
            except AttributeError:
                continue

    #if verbose:
    #    print(onetime, twotime, threetime)
    
    for node in tree.preorder_node_iter():
        try:
            node.edge_length=node.snps;
        except AttributeError:
            continue
    
    if verbose:
        print("Printing alignment with internal node sequences: ", output_prefix+".joint.aln")
    asr_output = open(output_prefix+".joint.aln", "w")
    for taxon in new_alignment:
        print(">"+taxon, file=asr_output)
        print(''.join(new_alignment[taxon]), file=asr_output)
    asr_output.close()
    
    if verbose:
        print("Printing tree with internal nodes labelled: ", output_prefix+".joint.tre")
    tree_output=open(output_prefix+".joint.tre", "w")
    print(tree.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True, suppress_internal_node_labels=True).replace("'",""), file=tree_output)
    tree_output.close()
    
    if verbose:
        print("Done")



################
# Main program #
################   

def main():
    args=getargv()
    jar(args.alignment, args.tree, args.info, args.prefix, args.verbose)
            

if __name__ == "__main__":
    main()
