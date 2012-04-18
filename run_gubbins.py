#!/usr/bin/env python
# encoding: utf-8
import sys
import argparse
import subprocess
import os
import time
from Bio import Phylo
import dendropy

def robinson_foulds_distance(input_tree_name,output_tree_name):
  input_tree  = dendropy.Tree.get_from_path(input_tree_name, 'newick')
  output_tree = dendropy.Tree.get_from_path(output_tree_name, 'newick')
  input_tree.robinson_foulds_distance(output_tree)

def reroot_tree_with_outgroup(tree_name, outgroup):
  if outgroup:
    tree = Phylo.read(tree_name, 'newick')
    print tree
    tree.root_with_outgroup({'name': outgroup})
    Phylo.write(tree, tree_name, 'newick')

# config variables
TREE_BUILDER_EXEC = 'raxmlHPC -f d  -m GTRGAMMA'
GUBBINS_EXEC = './gubbins'

parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
parser.add_argument('alignment_filename',       help='Multifasta alignment file')
parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
parser.add_argument('--starting_tree',    '-s', help='Starting tree')
parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', default = 5)
args = parser.parse_args()

# find all snp sites
subprocess.call([GUBBINS_EXEC, "-s", args.alignment_filename])

# get the base filename
(base_directory,base_filename) = os.path.split(args.alignment_filename)
(base_filename_without_ext,extension) = os.path.splitext(base_filename)
starting_base_filename = base_filename

current_time = int(time.time())

previous_robinson_foulds_distance = 0

for i in range(1, args.iterations+1):
  current_tree = "RAxML_result."+base_filename_without_ext+"."+str(current_time) +".iteration_"+str(i)
  previous_tree_name = base_filename
  previous_tree = ""

  if i> 1:
    previous_tree_name = "RAxML_result."+base_filename_without_ext+"."+str(current_time)+".iteration_"+ str(i-1)
    previous_tree = "-t "+ previous_tree_name
    base_filename = previous_tree_name

  tree_building_command = TREE_BUILDER_EXEC+ " -s "+previous_tree_name+".phylip -n "+base_filename_without_ext+"."+str(current_time)+".iteration_"+str(i)+" "+previous_tree
  gubbins_command = GUBBINS_EXEC+" -r "+ args.alignment_filename+" "+starting_base_filename+".vcf "+str(current_tree)+" "+starting_base_filename+".phylip"
  print tree_building_command
  print gubbins_command
  
  #subprocess.call(tree_building_command)
  reroot_tree_with_outgroup(str(current_tree), args.outgroup)
  #subprocess.call(gubbins_command)

  current_robinson_foulds_distance  = robinson_foulds_distance(previous_tree_name,current_tree)
  if current_robinson_foulds_distance == previous_robinson_foulds_distance:
    break
  else
    previous_robinson_foulds_distance = current_robinson_foulds_distance
  