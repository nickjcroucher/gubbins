#!/usr/bin/env python-2.7
# encoding: utf-8
import sys
import argparse
import subprocess
import os
import time
from Bio import Phylo
import dendropy
from array import *
from Bio import SeqIO
from cStringIO import StringIO
import shutil


# config variables
RAXML_EXEC = 'raxmlHPC -f d  -m GTRGAMMA'
FASTTREE_EXEC = 'FastTree'
FASTTREE_PARAMS = '-gtr -gamma -nt'
GUBBINS_EXEC = 'gubbins'

def robinson_foulds_distance(input_tree_name,output_tree_name):
  input_tree  = dendropy.Tree.get_from_path(input_tree_name, 'newick')
  output_tree = dendropy.Tree.get_from_path(output_tree_name, 'newick')
  return input_tree.robinson_foulds_distance(output_tree)

def reroot_tree(tree_name, outgroup):
  if outgroup:
    reroot_tree_with_outgroup(tree_name, outgroup)
  else:
    reroot_tree_at_midpoint(tree_name)

def reroot_tree_with_outgroup(tree_name, outgroup):
  tree = Phylo.read(tree_name, 'newick')
  tree.root_with_outgroup({'name': outgroup})
  Phylo.write(tree, tree_name, 'newick')
    
def reroot_tree_at_midpoint(tree_name):
  tree  = dendropy.Tree.get_from_path(tree_name, 'newick',
            preserve_underscores=True)
  tree.reroot_at_midpoint(update_splits=False, delete_outdegree_one=False)
  tree.deroot()
  tree.write_to_path(
    tree_name, 
    'newick',
    taxon_set=None,
    suppress_leaf_taxon_labels=False,
    suppress_leaf_node_labels=True,
    suppress_internal_taxon_labels=False,
    suppress_internal_node_labels=False,
    suppress_rooting=True,
    suppress_edge_lengths=False,
    unquoted_underscores=False,
    preserve_spaces=False,
    store_tree_weights=False,
    suppress_annotations=True,
    annotations_as_nhx=False,
    suppress_item_comments=True,
    node_label_element_separator=' ',
    node_label_compose_func=None)

def raxml_current_tree_name(base_filename_without_ext,current_time, i):
  return "RAxML_result."+base_filename_without_ext+"."+str(current_time) +".iteration_"+str(i)
    
    
def raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i):
  previous_tree_name = base_filename
  
  if i> 1:
    previous_tree_name = "RAxML_result."+base_filename_without_ext+"."+str(current_time)+".iteration_"+ str(i-1)
  return previous_tree_name
  
  
def raxml_previous_tree(base_filename_without_ext, base_filename, current_time,i,previous_tree_name):
  previous_tree = ""

  if i> 1:
    previous_tree = "-t "+ previous_tree_name
  return previous_tree
  
  
def raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time, raxml_exec,previous_tree_name):
  previous_tree = raxml_previous_tree(base_filename_without_ext, base_filename, current_time,i,previous_tree_name)
  return raxml_exec+ " -s "+previous_tree_name+".phylip -n "+base_filename_without_ext+"."+str(current_time)+".iteration_"+str(i)+" "+previous_tree


def raxml_gubbins_command(base_filename_without_ext,starting_base_filename,current_time, i,alignment_filename,gubbins_exec,min_snps):
  current_tree_name = raxml_current_tree_name(base_filename_without_ext,current_time, i)
  return gubbins_exec+" -r -v "+starting_base_filename+".vcf -t "+str(current_tree_name)+" -p "+starting_base_filename+".phylip -m "+ str(min_snps)+" "+ alignment_filename
  
def fasttree_current_tree_name(base_filename, i):
  return base_filename+".iteration_"+str(i)

def fasttree_previous_tree_name(base_filename, i):
  return base_filename+".iteration_"+str(i-1)

def fasttree_tree_building_command(i, starting_tree, current_tree_name,starting_base_filename, previous_tree_name,fasttree_exec, fasttree_params ):
  current_tree_name = fasttree_current_tree_name(base_filename, i)
  
  input_tree = ""
  if starting_tree is not None:
    input_tree = " -intree " + starting_tree
  elif i > 1 :
    input_tree = " -intree "+ previous_tree_name

  return fasttree_exec+" "+ input_tree+" "+ fasttree_params+" "+ starting_base_filename+".snp_sites.aln   > "+current_tree_name

def  fasttree_gubbins_command(base_filename,starting_base_filename, i,alignment_filename,gubbins_exec,min_snps):
  current_tree_name = fasttree_current_tree_name(base_filename, i)
  return gubbins_exec+" -r -v "+starting_base_filename+".vcf -t "+str(current_tree_name)+" -p "+starting_base_filename+".phylip -m "+ str(min_snps)+" "+ alignment_filename

def  starting_tree_gubbins_command(base_filename,starting_base_filename, i,alignment_filename,gubbins_exec,min_snps,starting_tree):
  return gubbins_exec+" -r -v "+starting_base_filename+".vcf -t "+str(starting_tree)+" -p "+starting_base_filename+".phylip -m "+ str(min_snps)+" "+ alignment_filename


def number_of_sequences_in_alignment(filename):
  return len(get_sequence_names_from_alignment(filename))
  
def get_sequence_names_from_alignment(filename):
  sequence_names = []
  handle = open(filename, "rU")
  for record in SeqIO.parse(handle, "fasta") :
    sequence_names.append(record.id)
  handle.close()
  return sequence_names

def pairwise_comparison(filename,base_filename,gubbins_exec,alignment_filename):
  sequence_names = get_sequence_names_from_alignment(filename)
  create_pairwise_newick_tree(sequence_names, base_filename+".tre")
  subprocess.call(gubbins_exec+" -r -v "+base_filename+".vcf -t "+base_filename+".tre "+" -p "+base_filename+".phylip "+ alignment_filename, shell=True)
  
def create_pairwise_newick_tree(sequence_names, output_filename):
  tree = Phylo.read(StringIO('('+sequence_names[0]+','+sequence_names[1]+')'), "newick")
  Phylo.write(tree, output_filename, 'newick')

def which(program):
  def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
  fpath, fname = os.path.split(program)
  if fpath:
    if is_exe(program):
      return program
  else:
    for path in os.environ["PATH"].split(os.pathsep):
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return exe_file
 
  return None



parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
parser.add_argument('alignment_filename',       help='Multifasta alignment file')
parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
parser.add_argument('--starting_tree',    '-s', help='Starting tree')
parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging')
parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
parser.add_argument('--min_snps',       '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
args = parser.parse_args()

# check that all the external executable dependancies are available
if which(GUBBINS_EXEC) is None:
  print "gubbins is not in your path"
  sys.exit()
if (args.tree_builder == "raxml" or args.tree_builder == "hybrid") and which('raxmlHPC') is None:
  print "RAxML is not in your path"
  sys.exit()
if (args.tree_builder == "fasttree" or args.tree_builder == "hybrid") and which(FASTTREE_EXEC) is None:
  print "FastTree is not in your path"
  sys.exit()

current_time = int(time.time())
snps_time =0;
recombinations_time =0;
if args.verbose > 0:
  print current_time

# find all snp sites
if args.verbose > 0:
  print GUBBINS_EXEC + args.alignment_filename
subprocess.call([GUBBINS_EXEC, args.alignment_filename])
if args.verbose > 0:
  snps_time = int(time.time() - current_time);
  print "SNPs Time:"+ int(snps_time)

# get the base filename
(base_directory,base_filename) = os.path.split(args.alignment_filename)
(base_filename_without_ext,extension) = os.path.splitext(base_filename)
starting_base_filename = base_filename

# Perform pairwise comparison if there are only 2 sequences
number_of_sequences = number_of_sequences_in_alignment(args.alignment_filename)
if(number_of_sequences == 2):
  pairwise_comparison(args.alignment_filename,starting_base_filename,GUBBINS_EXEC,args.alignment_filename)
  if args.verbose > 0:
    print "Recombinations Time:"+ int(time.time() - current_time - snps_time)
    print "Total Time:"+ int(time.time() - current_time)
  sys.exit()



previous_robinson_foulds_distances = array('d',[])

tree_building_command = ""
gubbins_command       = ""
previous_tree_name    = ""
current_tree_name     = ""

for i in range(1, args.iterations+1):
    
  if args.tree_builder == "hybrid" :
    if i == 1:
      previous_tree_name    = fasttree_previous_tree_name(base_filename, i)
      current_tree_name     = fasttree_current_tree_name(base_filename, i)
      tree_building_command = fasttree_tree_building_command(i, args.starting_tree,current_tree_name,base_filename,previous_tree_name,FASTTREE_EXEC,'-fastest' )
      gubbins_command       = fasttree_gubbins_command(base_filename,starting_base_filename+".gaps", i,args.alignment_filename,GUBBINS_EXEC,args.min_snps)
      
    elif i == 2:
      previous_tree_name    = current_tree_name
      current_tree_name     = raxml_current_tree_name(base_filename_without_ext,current_time, i)
      tree_building_command = raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name)
      gubbins_command       = raxml_gubbins_command(base_filename_without_ext,starting_base_filename+".gaps",current_time, i,args.alignment_filename,GUBBINS_EXEC,args.min_snps)
    else:
      previous_tree_name    = raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i)
      current_tree_name     = raxml_current_tree_name(base_filename_without_ext,current_time, i)
      tree_building_command = raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name)
      gubbins_command       = raxml_gubbins_command(base_filename_without_ext,starting_base_filename+".gaps",current_time, i,args.alignment_filename,GUBBINS_EXEC,args.min_snps)
  
  elif args.tree_builder == "raxml":
    previous_tree_name    = raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i)
    current_tree_name     = raxml_current_tree_name(base_filename_without_ext,current_time, i)
    tree_building_command = raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name)
    gubbins_command       = raxml_gubbins_command(base_filename_without_ext,starting_base_filename+".gaps",current_time, i,args.alignment_filename,GUBBINS_EXEC,args.min_snps)
    
  elif args.tree_builder == "fasttree":
    previous_tree_name    = fasttree_previous_tree_name(base_filename, i)
    if i == 1:
      previous_tree_name  = base_filename
    current_tree_name     = fasttree_current_tree_name(base_filename, i)

    tree_building_command = fasttree_tree_building_command(i, args.starting_tree,current_tree_name,previous_tree_name,previous_tree_name,FASTTREE_EXEC,FASTTREE_PARAMS )
    gubbins_command       = fasttree_gubbins_command(base_filename,starting_base_filename+".gaps", i,args.alignment_filename,GUBBINS_EXEC,args.min_snps)
  
  if args.verbose > 0:
    print tree_building_command
    
    
  if args.starting_tree is not None and i == 1:
    shutil.copyfile(args.starting_tree, current_tree_name)
  else:
    subprocess.call(tree_building_command, shell=True)
    
  if args.verbose > 0:
    print int(time.time() - current_time)
  
  reroot_tree(str(current_tree_name), args.outgroup)

  if(os.path.exists("latest.tre")):
    os.remove("latest.tre")
  os.symlink(str(current_tree_name), "latest.tre")
 
  start_gubbins = int(time.time())
  if args.verbose > 0:
    print gubbins_command
  subprocess.call(gubbins_command, shell=True)
  if args.verbose > 0:
    recombinations_time = recombinations_time + int(time.time() - start_gubbins)
    print int(time.time() - current_time)
  
  # first iteration creates tree 1
  # 2nd iteration creates tree 2, and you can calculate first RF distance
  # 3rd iteration creates tree 3, and you can now compare RF distances with the previous iteration
  if i == 2:
    previous_robinson_foulds_distances.append(robinson_foulds_distance(previous_tree_name,current_tree_name))
  elif i > 2:
    current_robinson_foulds_distance  = robinson_foulds_distance(previous_tree_name,current_tree_name)
    if args.verbose > 0:
      print "RF Distance (previous, current): "+ str(previous_robinson_foulds_distances) +", "+ str(current_robinson_foulds_distance)
      
    try:
      previous_robinson_foulds_distances.index(current_robinson_foulds_distance)
      break
    except ValueError:
      previous_robinson_foulds_distances.append(current_robinson_foulds_distance)

if args.verbose > 0:
  print "Tree Building Time:"+ int(time.time() - current_time - snps_time - recombinations_time )
  print "Recombinations Time:"+ int(recombinations_time)
  print "Total Time:"+ int(time.time() - current_time)
  