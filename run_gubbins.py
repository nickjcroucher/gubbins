#!/usr/bin/env python-2.7
# encoding: utf-8
#
# Wellcome Trust Sanger Institute
# Copyright (C) 2012  Wellcome Trust Sanger Institute
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

# fastml
# disable timestamp by default to allow pipelineing
# check a file exists before trying to delete it


import sys
import argparse
import subprocess
import os
import time
import re
from Bio import Phylo
import dendropy
from Bio import SeqIO
from Bio import AlignIO
from cStringIO import StringIO
import shutil

# config variables
RAXML_EXEC = 'raxmlHPC -f d  -m GTRGAMMA'
FASTTREE_EXEC = 'FastTree'
FASTTREE_PARAMS = '-gtr -gamma -nt'
GUBBINS_EXEC = 'gubbins'
FASTML_EXEC = 'fastml -mg -qf -b '

# Todo
# extract code into modules
# catch exception errors when shelling out


def robinson_foulds_distance(input_tree_name,output_tree_name):
  input_tree  = dendropy.Tree.get_from_path(input_tree_name, 'newick')
  output_tree = dendropy.Tree.get_from_path(output_tree_name, 'newick')
  return input_tree.robinson_foulds_distance(output_tree)
  
def has_tree_been_seen_before(tree_file_names):
  if len(tree_file_names) <= 2:
    return 0

  for tree_file_name in tree_file_names:
    if tree_file_name is not tree_file_names[-1]:
      current_rf_distance = robinson_foulds_distance(tree_file_name,tree_file_names[-1])
      if(current_rf_distance == 0.0):
        return 1

  return 0

def reroot_tree(tree_name, outgroup):
  if outgroup:
    reroot_tree_with_outgroup(tree_name, outgroup)
  else:
    reroot_tree_at_midpoint(tree_name)

def reroot_tree_with_outgroup(tree_name, outgroup):
  tree = Phylo.read(tree_name, 'newick')
  tree.root_with_outgroup({'name': outgroup})
  Phylo.write(tree, tree_name, 'newick')
    
def split_all_non_bi_nodes(node):
  if node.is_leaf():
    return None
  elif len(node.child_nodes()) > 2:
    split_child_nodes(node)

  for child_node in node.child_nodes():
    split_all_non_bi_nodes(child_node)

  return None

def split_child_nodes(node):
  all_child_nodes = node.child_nodes()
  #skip over the first node
  first_child = all_child_nodes.pop()
  # create a placeholder node to hang everything else off
  new_child_node = node.new_child(edge_length=0)
  # move the subtree into the placeholder
  new_child_node.set_child_nodes(all_child_nodes)
  # probably not really nessisary
  node.set_child_nodes((first_child,new_child_node))
    
def reroot_tree_at_midpoint(tree_name):
  tree  = dendropy.Tree.get_from_path(tree_name, 'newick',
            preserve_underscores=True)
  split_all_non_bi_nodes(tree.seed_node)
    
  tree.reroot_at_midpoint(update_splits=True, delete_outdegree_one=False)
  tree.deroot()
  tree.update_splits()
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
    unquoted_underscores=True,
    preserve_spaces=False,
    store_tree_weights=False,
    suppress_annotations=True,
    annotations_as_nhx=False,
    suppress_item_comments=True,
    node_label_element_separator=' ',
    node_label_compose_func=None)

def raxml_base_name(base_filename_without_ext,current_time):
  return base_filename_without_ext+"."+str(current_time) +".iteration_"

def raxml_current_tree_name(base_filename_without_ext,current_time, i):
  return "RAxML_result."+raxml_base_name(base_filename_without_ext,current_time)+str(i)
    
    
def raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i):
  previous_tree_name = base_filename
  
  if i> 1:
    previous_tree_name = "RAxML_result."+raxml_base_name(base_filename_without_ext,current_time)+ str(i-1)
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
  return gubbins_exec+" -r -v "+starting_base_filename+".vcf -t "+str(current_tree_name)+" -m "+ str(min_snps)+" "+ starting_base_filename+".snp_sites.aln"
  
def raxml_regex_for_file_deletions(base_filename_without_ext,current_time,starting_base_filename, max_intermediate_iteration):
  regex_for_file_deletions = []
  # Can delete all of these files
  regex_for_file_deletions.append("^RAxML_(bestTree|info|log|parsimonyTree)."+raxml_base_name(base_filename_without_ext,current_time))
  
  regex_for_file_deletions.append(starting_files_regex)
  
  # loop over previous iterations and delete
  for file_iteration in range(1,max_intermediate_iteration):
    regex_for_file_deletions.append("^RAxML_result."+raxml_base_name(base_filename_without_ext,current_time)+str(file_iteration))

  return regex_for_file_deletions
  
def fasttree_regex_for_file_deletions(starting_base_filename, max_intermediate_iteration):
  regex_for_file_deletions = []
  regex_for_file_deletions.append(starting_files_regex(starting_base_filename))

  # loop over previous iterations and delete
  for file_iteration in range(1,max_intermediate_iteration):
    regex_for_file_deletions.append("^"+starting_base_filename+".iteration_"+str(file_iteration))

  return regex_for_file_deletions

def starting_files_regex(starting_base_filename):
  # starting file with gapped and ungapped snps
  return starting_base_filename+".(gaps|vcf|snp_sites|phylip)"

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
  return gubbins_exec+" -r -v "+starting_base_filename+".vcf -t "+str(current_tree_name)+" -m "+ str(min_snps)+" "+ starting_base_filename+".snp_sites.aln"

def fasttree_fastml_command(fastml_exec, alignment_filename, base_filename,i):
  current_tree_name = fasttree_current_tree_name(base_filename, i)
  return generate_fastml_command(fastml_exec, alignment_filename, current_tree_name)


def raxml_fastml_command(fastml_exec, alignment_filename, base_filename_without_ext,current_time, i):
  current_tree_name = raxml_current_tree_name(base_filename_without_ext,current_time, i)
  return generate_fastml_command(fastml_exec, alignment_filename, current_tree_name)

def generate_fastml_command(fastml_exec, alignment_filename, tree_filename):
  return (fastml_exec 
    + " -s " + alignment_filename 
    + " -t " + tree_filename 
    + " -x " + tree_filename + ".output_tree"
    + " -y " + tree_filename + ".ancestor.tre"
    + " -j " + tree_filename + ".seq.joint.txt"
    + " -k " + tree_filename + ".seq.marginal.txt"
    + " -d " + tree_filename + ".prob.joint.txt"
    + " -e " + tree_filename + ".prob.marginal.txt")


def number_of_sequences_in_alignment(filename):
  return len(get_sequence_names_from_alignment(filename))
  
def get_sequence_names_from_alignment(filename):
  sequence_names = []
  handle = open(filename, "rU")
  for record in SeqIO.parse(handle, "fasta") :
    sequence_names.append(record.id)
  handle.close()
  return sequence_names
  
  # reparsing a fasta file splits the lines which makes fastml work
def reconvert_fasta_file(input_filename, output_filename):
  input_handle = open(input_filename, "rU")
  output_handle = open(output_filename, "w+")
  alignments = AlignIO.parse(input_handle, "fasta")
  AlignIO.write(alignments, output_handle, "fasta")
  output_handle.close()
  input_handle.close()
  return

def pairwise_comparison(filename,base_filename,gubbins_exec,alignment_filename):
  sequence_names = get_sequence_names_from_alignment(filename)
  create_pairwise_newick_tree(sequence_names, base_filename+".tre")
  subprocess.check_call(gubbins_exec+" -r -v "+base_filename+".vcf -t "+base_filename+".tre "+ base_filename+".snp_sites.aln", shell=True)
  
def create_pairwise_newick_tree(sequence_names, output_filename):
  tree = Phylo.read(StringIO('('+sequence_names[0]+','+sequence_names[1]+')'), "newick")
  Phylo.write(tree, output_filename, 'newick')
  
def delete_files_based_on_list_of_regexes(directory_to_search, regex_for_file_deletions, verbose):
  for dirname, dirnames, filenames in os.walk(directory_to_search):
    for filename in filenames:
      for deletion_regex in regex_for_file_deletions:
        full_path_of_file_for_deletion = os.path.join(directory_to_search, filename)
        if(re.match(str(deletion_regex), filename) != None and os.path.exists(full_path_of_file_for_deletion)):
          if verbose > 0:
            print "Deleting file: "+ os.path.join(directory_to_search, filename)
          os.remove(full_path_of_file_for_deletion)

def which(program):
  executable = program.split(" ")
  program = executable[0]
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
parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging')
parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
args = parser.parse_args()

# check that all the external executable dependancies are available
if which(GUBBINS_EXEC) is None:
  print "gubbins is not in your path"
  sys.exit()
if which(FASTML_EXEC) is None:
  print "fastml is not in your path"
  sys.exit()
if (args.tree_builder == "raxml" or args.tree_builder == "hybrid") and which('raxmlHPC') is None:
  print "RAxML is not in your path"
  sys.exit()
if (args.tree_builder == "fasttree" or args.tree_builder == "hybrid") and which(FASTTREE_EXEC) is None:
  print "FastTree is not in your path"
  sys.exit()

current_time = ''
if args.use_time_stamp > 0:
  current_time = int(time.time())
  if args.verbose > 0:
    print current_time

# find all snp sites
if args.verbose > 0:
  print GUBBINS_EXEC + args.alignment_filename
subprocess.check_call([GUBBINS_EXEC, args.alignment_filename])
if args.verbose > 0:
  print int(time.time())

# get the base filename
(base_directory,base_filename) = os.path.split(args.alignment_filename)
(base_filename_without_ext,extension) = os.path.splitext(base_filename)
starting_base_filename = base_filename

reconvert_fasta_file(starting_base_filename+".gaps.snp_sites.aln",starting_base_filename+".start")

# Perform pairwise comparison if there are only 2 sequences
number_of_sequences = number_of_sequences_in_alignment(args.alignment_filename)
if(number_of_sequences == 2):
  pairwise_comparison(args.alignment_filename,starting_base_filename,GUBBINS_EXEC,args.alignment_filename)
  sys.exit()


latest_file_name = "latest_tree."+base_filename_without_ext+"."+str(current_time)+".tre"
tree_file_names = []

tree_building_command = ""
gubbins_command       = ""
previous_tree_name    = ""
current_tree_name     = ""
max_iteration = 1

for i in range(1, args.iterations+1):
  max_iteration += 1
  
  if args.tree_builder == "hybrid" :
    if i == 1:
      previous_tree_name    = fasttree_previous_tree_name(base_filename, i)
      current_tree_name     = fasttree_current_tree_name(base_filename, i)
      tree_building_command = fasttree_tree_building_command(i, args.starting_tree,current_tree_name,base_filename,previous_tree_name,FASTTREE_EXEC, FASTTREE_PARAMS )
      fastml_command        = fasttree_fastml_command(FASTML_EXEC, starting_base_filename+".gaps.snp_sites.aln", base_filename, i)
      gubbins_command       = fasttree_gubbins_command(base_filename,starting_base_filename+".gaps", i,args.alignment_filename,GUBBINS_EXEC,args.min_snps)

    elif i == 2:
      previous_tree_name    = current_tree_name
      current_tree_name     = raxml_current_tree_name(base_filename_without_ext,current_time, i)
      tree_building_command = raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name)
      fastml_command        = raxml_fastml_command(FASTML_EXEC, starting_base_filename+".gaps.snp_sites.aln", base_filename_without_ext,current_time, i)
      gubbins_command       = raxml_gubbins_command(base_filename_without_ext,starting_base_filename+".gaps",current_time, i,args.alignment_filename,GUBBINS_EXEC,args.min_snps)
    else:
      previous_tree_name    = raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i)
      current_tree_name     = raxml_current_tree_name(base_filename_without_ext,current_time, i)
      tree_building_command = raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name)
      fastml_command        = raxml_fastml_command(FASTML_EXEC, starting_base_filename+".gaps.snp_sites.aln", base_filename_without_ext,current_time, i)
      gubbins_command       = raxml_gubbins_command(base_filename_without_ext,starting_base_filename+".gaps",current_time, i,args.alignment_filename,GUBBINS_EXEC,args.min_snps)
  
  elif args.tree_builder == "raxml":
    previous_tree_name    = raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i)
    current_tree_name     = raxml_current_tree_name(base_filename_without_ext,current_time, i)
    tree_building_command = raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name)
    fastml_command        = raxml_fastml_command(FASTML_EXEC, starting_base_filename+".start", base_filename_without_ext,current_time, i)
    gubbins_command       = raxml_gubbins_command(base_filename_without_ext,starting_base_filename+".gaps",current_time, i,args.alignment_filename,GUBBINS_EXEC,args.min_snps)
    
  elif args.tree_builder == "fasttree":
    previous_tree_name    = fasttree_previous_tree_name(base_filename, i)
    if i == 1:
      previous_tree_name  = base_filename
    current_tree_name     = fasttree_current_tree_name(base_filename, i)

    tree_building_command = fasttree_tree_building_command(i, args.starting_tree,current_tree_name,previous_tree_name,previous_tree_name,FASTTREE_EXEC,FASTTREE_PARAMS )
    fastml_command        = fasttree_fastml_command(FASTML_EXEC, starting_base_filename+".gaps.snp_sites.aln", base_filename, i)
    gubbins_command       = fasttree_gubbins_command(base_filename,starting_base_filename+".gaps", i,args.alignment_filename,GUBBINS_EXEC,args.min_snps)
  
  if args.verbose > 0:
    print tree_building_command
    
    
  if args.starting_tree is not None and i == 1:
    shutil.copyfile(args.starting_tree, current_tree_name)
  else:
    subprocess.check_call(tree_building_command, shell=True)
    
  if args.verbose > 0:
    print int(time.time())
  
  reroot_tree(str(current_tree_name), args.outgroup)

  if(os.path.exists(latest_file_name)):
    os.remove(latest_file_name)
  os.symlink(str(current_tree_name), latest_file_name)
 
  if args.verbose > 0:
    print fastml_command
  subprocess.check_call(fastml_command, shell=True)
  shutil.copyfile(current_tree_name+'.output_tree',current_tree_name)
  reconvert_fasta_file(current_tree_name+'.seq.joint.txt',starting_base_filename+".gaps.snp_sites.aln")

  if args.verbose > 0:
    print int(time.time())
 
  if args.verbose > 0:
    print gubbins_command
  subprocess.check_call(gubbins_command, shell=True)
  if args.verbose > 0:
    print int(time.time())
  
  tree_file_names.append(current_tree_name)
  if i > 2:
    if has_tree_been_seen_before(tree_file_names):
      if args.verbose > 0:
        print "Tree observed before so stopping: "+ str(current_tree_name)
      break

# cleanup intermediate files
if args.no_cleanup == 0 or args.no_cleanup is None:
  max_intermediate_iteration  = max_iteration - 1
  
  raxml_regex_for_file_deletions = raxml_regex_for_file_deletions(base_filename_without_ext,current_time,starting_base_filename, max_intermediate_iteration)
  delete_files_based_on_list_of_regexes('.', raxml_regex_for_file_deletions, args.verbose)
  
  fasttree_regex_for_file_deletions = fasttree_regex_for_file_deletions(starting_base_filename, max_intermediate_iteration)
  delete_files_based_on_list_of_regexes('.', fasttree_regex_for_file_deletions, args.verbose)
  
