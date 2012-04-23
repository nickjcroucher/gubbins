#!/usr/bin/env python
# encoding: utf-8
import sys
import argparse
import subprocess
import os
import time
from Bio import Phylo
import dendropy
from array import *


# config variables
RAXML_EXEC = 'raxmlHPC -f d  -m GTRGAMMA'
FASTTREE_EXEC = 'FastTree'
FASTTREE_PARAMS = '-gtr -gamma -nt'
GUBBINS_EXEC = './gubbins'

def robinson_foulds_distance(input_tree_name,output_tree_name):
  input_tree  = dendropy.Tree.get_from_path(input_tree_name, 'newick')
  output_tree = dendropy.Tree.get_from_path(output_tree_name, 'newick')
  return input_tree.robinson_foulds_distance(output_tree)


def reroot_tree_with_outgroup(tree_name, outgroup):
  if outgroup:
    tree = Phylo.read(tree_name, 'newick')
    tree.root_with_outgroup({'name': outgroup})
    Phylo.write(tree, tree_name, 'newick')
    
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


def raxml_gubbins_command(base_filename_without_ext,starting_base_filename,current_time, i,alignment_filename,gubbins_exec):
  current_tree_name = raxml_current_tree_name(base_filename_without_ext,current_time, i)
  return gubbins_exec+" -r "+ alignment_filename+" "+starting_base_filename+".vcf "+str(current_tree_name)+" "+starting_base_filename+".phylip"
  
def fasttree_current_tree_name(base_filename, i):
  return base_filename+".iteration_"+str(i)

def fasttree_previous_tree_name(base_filename, i):
  return base_filename+".iteration_"+str(i-1)

def fasttree_tree_building_command(i, starting_tree, current_tree_name,starting_base_filename, previous_tree_name,fasttree_exec, fasttree_params ):
  current_tree_name = fasttree_current_tree_name(base_filename, i)
  
  input_tree = ""
  if starting_tree is not None:
    input_tree = " -intree" + starting_tree
  elif i > 1 :
    input_tree = " -intree "+ previous_tree_name

  return fasttree_exec+" "+ input_tree+" "+ fasttree_params+" "+ starting_base_filename+".snp_sites.aln   > "+current_tree_name

def  fasttree_gubbins_command(base_filename,starting_base_filename, i,alignment_filename,gubbins_exec):
  current_tree_name = fasttree_current_tree_name(base_filename, i)
  return gubbins_exec+" -r "+ alignment_filename+" "+starting_base_filename+".vcf "+str(current_tree_name)+" "+starting_base_filename+".phylip"
 


parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
parser.add_argument('alignment_filename',       help='Multifasta alignment file')
parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
parser.add_argument('--starting_tree',    '-s', help='Starting tree')
parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging')
parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree), default is to use FastTree for 1st iteration and RAxML for rest', default = "hybrid")
parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
args = parser.parse_args()

# find all snp sites
subprocess.call([GUBBINS_EXEC, "-s", args.alignment_filename])

# get the base filename
(base_directory,base_filename) = os.path.split(args.alignment_filename)
(base_filename_without_ext,extension) = os.path.splitext(base_filename)
starting_base_filename = base_filename

current_time = int(time.time())

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
      tree_building_command = fasttree_tree_building_command(i, args.starting_tree,current_tree_name,args.alignment_filename,previous_tree_name,FASTTREE_EXEC,FASTTREE_PARAMS )
      gubbins_command       = fasttree_gubbins_command(base_filename,starting_base_filename, i,args.alignment_filename,GUBBINS_EXEC)
      
    elif i == 2:
      previous_tree_name    = current_tree_name
      current_tree_name     = raxml_current_tree_name(base_filename_without_ext,current_time, i)
      tree_building_command = raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name)
      gubbins_command       = raxml_gubbins_command(base_filename_without_ext,starting_base_filename,current_time, i,args.alignment_filename,GUBBINS_EXEC)
    else:
      previous_tree_name    = raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i)
      current_tree_name     = raxml_current_tree_name(base_filename_without_ext,current_time, i)
      tree_building_command = raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name)
      gubbins_command       = raxml_gubbins_command(base_filename_without_ext,starting_base_filename,current_time, i,args.alignment_filename,GUBBINS_EXEC)
  
  elif args.tree_builder == "raxml":
    previous_tree_name    = raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i)
    current_tree_name     = raxml_current_tree_name(base_filename_without_ext,current_time, i)
    tree_building_command = raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name)
    gubbins_command       = raxml_gubbins_command(base_filename_without_ext,starting_base_filename,current_time, i,args.alignment_filename,GUBBINS_EXEC)
    
  elif args.tree_builder == "fasttree":
    previous_tree_name    = fasttree_previous_tree_name(base_filename, i)
    current_tree_name     = fasttree_current_tree_name(base_filename, i)
    
    if i == 1:
      os.rename(base_filename+".vcf",           previous_tree_name+".vcf")
      os.rename(base_filename+".phylip",        previous_tree_name+".phylip")
      os.rename(base_filename+".snp_sites.aln", previous_tree_name+".snp_sites.aln")
    tree_building_command = fasttree_tree_building_command(i, args.starting_tree,current_tree_name,args.alignment_filename,previous_tree_name,FASTTREE_EXEC,FASTTREE_PARAMS )
    gubbins_command       = fasttree_gubbins_command(base_filename,starting_base_filename, i,args.alignment_filename,GUBBINS_EXEC)
  
  if args.verbose > 0:
    print tree_building_command
  subprocess.call(tree_building_command, shell=True)
  
  reroot_tree_with_outgroup(str(current_tree_name), args.outgroup)
  
  if args.verbose > 0:
    print gubbins_command
  subprocess.call(gubbins_command, shell=True)
  
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

  