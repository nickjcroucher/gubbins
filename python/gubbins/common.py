# encoding: utf-8
# Wellcome Trust Sanger Institute
# Copyright (C) 2013  Wellcome Trust Sanger Institute
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

import sys
import argparse
import subprocess
import os
import time
import re
import tempfile
from collections import Counter
from Bio import Phylo
import dendropy
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from cStringIO import StringIO
import shutil

class GubbinsError(Exception):
  def __init__(self, value,message):
      self.value = value
      self.message = message
  def __str__(self):
      return self.message

class GubbinsCommon():
  "Methods used by Gubbins"
  def __init__(self, input_args):
    self.args = input_args
  
  @staticmethod
  def is_starting_tree_valid(starting_tree):
    try:
      Phylo.read(starting_tree, 'newick')
      tree  = dendropy.Tree.get_from_path(starting_tree, 'newick', preserve_underscores=True)
    except:
      print "Error with the input starting tree: Is it a valid Newick file?"
      return 0
    return 1
    
  @staticmethod
  def do_the_names_match_the_fasta_file(starting_tree, alignment_filename):
    input_handle  = open(alignment_filename, "rU")
    alignments = AlignIO.parse(input_handle, "fasta")
    sequence_names = {}
    for alignment in alignments:
        for record in alignment:
            sequence_names[record.name] = 1
    input_handle.close()
    
    tree = dendropy.Tree.get_from_path(starting_tree, 'newick', preserve_underscores=True)
    
    leaf_nodes = tree.leaf_nodes()
    for i,lf in enumerate(leaf_nodes):
      if not leaf_nodes[i].taxon.label in sequence_names:
        print "Error: A taxon referenced in the starting tree isnt found in the input fasta file"
        return 0

    return 1

  @staticmethod
  def does_file_exist(alignment_filename, file_type_msg):
    if(not os.path.exists(alignment_filename)):
      print GubbinsError('','Cannot access the input '+file_type_msg+'. Check its been entered correctly')
      return 0
    return 1

  def parse_and_run(self):
    # Default parameters
    RAXML_EXEC = 'raxmlHPC -f d -p 1 -m GTRGAMMA'
    FASTTREE_EXEC = 'FastTree'
    FASTTREE_PARAMS = '-nosupport -gtr -gamma -nt'
    GUBBINS_EXEC = 'gubbins'
    FASTML_EXEC = 'fastml -mg -qf -b '

    GUBBINS_BUNDLED_EXEC = '../src/gubbins'

    # check that all the external executable dependancies are available
    if GubbinsCommon.which(GUBBINS_EXEC) is None:
      GUBBINS_EXEC = GubbinsCommon.use_bundled_exec(GUBBINS_EXEC, GUBBINS_BUNDLED_EXEC)
      if GubbinsCommon.which(GUBBINS_EXEC) is None:
        sys.exit(GUBBINS_EXEC+" is not in your path")
    if GubbinsCommon.which(FASTML_EXEC) is None:
      sys.exit("fastml is not in your path")
    if (self.args.tree_builder == "raxml" or self.args.tree_builder == "hybrid") and GubbinsCommon.which(RAXML_EXEC) is None:
      sys.exit("RAxML is not in your path")
    if (self.args.tree_builder == "fasttree" or self.args.tree_builder == "hybrid") and GubbinsCommon.which(FASTTREE_EXEC) is None:
      sys.exit("FastTree is not in your path")

    if(GubbinsCommon.does_file_exist(self.args.alignment_filename, 'Alignment File') == 0 or GubbinsCommon.is_input_fasta_file_valid(self.args.alignment_filename) == 0 ):
       sys.exit("There is a problem with your input fasta file so nothing can be done until you fix it")
       
    if(self.args.starting_tree is not None and self.args.starting_tree != "" and (GubbinsCommon.does_file_exist(self.args.starting_tree, 'Starting Tree') == 0 or GubbinsCommon.is_input_starting_tree_valid(self.args.starting_tree) )):
       sys.exit( "The starting tree is invalid")
       
    if(self.args.starting_tree is not None and self.args.starting_tree != "" and GubbinsCommon.do_the_names_match_the_fasta_file(self.args.starting_tree,self.args.alignment_filename) == 0):
      sys.exit("The names in the starting tree dont match the names in the fasta file")

    current_time = ''
    if self.args.use_time_stamp > 0:
      current_time = str(int(time.time()))+'.'
      if self.args.verbose > 0:
        print current_time

    # get the base filename
    (base_directory,base_filename) = os.path.split(self.args.alignment_filename)
    (base_filename_without_ext,extension) = os.path.splitext(base_filename)
    starting_base_filename = base_filename

    # put a filtered copy into a temp directory and work from that
    temp_working_dir = tempfile.mkdtemp(dir=os.getcwd())
    taxa_removed = GubbinsCommon.filter_out_alignments_with_too_much_missing_data(self.args.alignment_filename, temp_working_dir+"/"+starting_base_filename, self.args.filter_percentage, self.args.verbose)
    self.args.alignment_filename = temp_working_dir+"/"+starting_base_filename

    # If a starting tree has been provided make sure that taxa filtered out in the previous step are removed from the tree
    self.args.starting_tree = GubbinsCommon.filter_out_removed_taxa_from_tree_and_return_new_file(self.args.starting_tree, temp_working_dir, taxa_removed)

    # get the base filename
    (base_directory,base_filename) = os.path.split(self.args.alignment_filename)
    (base_filename_without_ext,extension) = os.path.splitext(base_filename)
    starting_base_filename = base_filename

    # find all snp sites
    if self.args.verbose > 0:
      print GUBBINS_EXEC +" "+ self.args.alignment_filename
    subprocess.check_call([GUBBINS_EXEC, self.args.alignment_filename])
    if self.args.verbose > 0:
      print int(time.time())

    GubbinsCommon.reconvert_fasta_file(starting_base_filename+".gaps.snp_sites.aln",starting_base_filename+".start")

    # Perform pairwise comparison if there are only 2 sequences
    number_of_sequences = GubbinsCommon.number_of_sequences_in_alignment(self.args.alignment_filename)
    if(number_of_sequences == 2):
      GubbinsCommon.pairwise_comparison(self.args.alignment_filename,starting_base_filename,GUBBINS_EXEC,self.args.alignment_filename,FASTML_EXEC)
      sys.exit()

    latest_file_name = "latest_tree."+base_filename_without_ext+"."+str(current_time)+"tre"
    tree_file_names = []

    tree_building_command = ""
    gubbins_command       = ""
    previous_tree_name    = ""
    current_tree_name     = ""
    max_iteration = 1


    # cleanup RAxML intermediate files
    if self.args.no_cleanup == 0 or self.args.no_cleanup is None:
      raxml_files_to_delete = GubbinsCommon.raxml_regex_for_file_deletions(base_filename_without_ext,current_time,starting_base_filename, self.args.iterations)
      GubbinsCommon.delete_files_based_on_list_of_regexes('.', raxml_files_to_delete, self.args.verbose)

    for i in range(1, self.args.iterations+1):
      max_iteration += 1

      if self.args.tree_builder == "hybrid" :
        if i == 1:
          previous_tree_name    = GubbinsCommon.fasttree_previous_tree_name(base_filename, i)
          current_tree_name     = GubbinsCommon.fasttree_current_tree_name(base_filename, i)
          tree_building_command = GubbinsCommon.fasttree_tree_building_command(i, self.args.starting_tree,current_tree_name,base_filename,previous_tree_name,FASTTREE_EXEC, FASTTREE_PARAMS,base_filename )
          fastml_command        = GubbinsCommon.fasttree_fastml_command(FASTML_EXEC, starting_base_filename+".snp_sites.aln", base_filename, i)
          gubbins_command       = GubbinsCommon.fasttree_gubbins_command(base_filename,starting_base_filename+".gaps", i,self.args.alignment_filename,GUBBINS_EXEC,self.args.min_snps,self.args.alignment_filename)

        elif i == 2:
          previous_tree_name    = current_tree_name
          current_tree_name     = GubbinsCommon.raxml_current_tree_name(base_filename_without_ext,current_time, i)
          tree_building_command = GubbinsCommon.raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name, self.args.verbose)
          fastml_command        = GubbinsCommon.raxml_fastml_command(FASTML_EXEC, starting_base_filename+".snp_sites.aln", base_filename_without_ext,current_time, i)
          gubbins_command       = GubbinsCommon.raxml_gubbins_command(base_filename_without_ext,starting_base_filename+".gaps",current_time, i,self.args.alignment_filename,GUBBINS_EXEC,self.args.min_snps,self.args.alignment_filename)
        else:
          previous_tree_name    = GubbinsCommon.raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i)
          current_tree_name     = GubbinsCommon.raxml_current_tree_name(base_filename_without_ext,current_time, i)
          tree_building_command = GubbinsCommon.raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name, self.args.verbose)
          fastml_command        = GubbinsCommon.raxml_fastml_command(FASTML_EXEC, starting_base_filename+".snp_sites.aln", base_filename_without_ext,current_time, i)
          gubbins_command       = GubbinsCommon.raxml_gubbins_command(base_filename_without_ext,starting_base_filename+".gaps",current_time, i,self.args.alignment_filename,GUBBINS_EXEC,self.args.min_snps,self.args.alignment_filename)

      elif self.args.tree_builder == "raxml":
        previous_tree_name    = GubbinsCommon.raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i)
        current_tree_name     = GubbinsCommon.raxml_current_tree_name(base_filename_without_ext,current_time, i)
        tree_building_command = GubbinsCommon.raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time,RAXML_EXEC,previous_tree_name, self.args.verbose)
        fastml_command        = GubbinsCommon.raxml_fastml_command(FASTML_EXEC, starting_base_filename+".snp_sites.aln", base_filename_without_ext,current_time, i)
        gubbins_command       = GubbinsCommon.raxml_gubbins_command(base_filename_without_ext,starting_base_filename+".gaps",current_time, i,self.args.alignment_filename,GUBBINS_EXEC,self.args.min_snps,self.args.alignment_filename)

      elif self.args.tree_builder == "fasttree":
        previous_tree_name    = GubbinsCommon.fasttree_previous_tree_name(base_filename, i)
        if i == 1:
          previous_tree_name  = base_filename
        current_tree_name     = GubbinsCommon.fasttree_current_tree_name(base_filename, i)

        tree_building_command = GubbinsCommon.fasttree_tree_building_command(i, self.args.starting_tree,current_tree_name,previous_tree_name,previous_tree_name,FASTTREE_EXEC,FASTTREE_PARAMS,base_filename )
        fastml_command        = GubbinsCommon.fasttree_fastml_command(FASTML_EXEC, starting_base_filename+".snp_sites.aln", base_filename, i)
        gubbins_command       = GubbinsCommon.fasttree_gubbins_command(base_filename,starting_base_filename+".gaps", i,self.args.alignment_filename,GUBBINS_EXEC,self.args.min_snps,self.args.alignment_filename)

      if self.args.verbose > 0:
        print tree_building_command


      if self.args.starting_tree is not None and i == 1:
        shutil.copyfile(self.args.starting_tree, current_tree_name)
      else:
        subprocess.check_call(tree_building_command, shell=True)

      if self.args.verbose > 0:
        print int(time.time())

      GubbinsCommon.reroot_tree(str(current_tree_name), self.args.outgroup)

      fastml_command_suffix = ' > /dev/null 2>&1'
      if self.args.verbose > 0:
        print fastml_command
        fastml_command_suffix = ''

      subprocess.check_call(fastml_command+fastml_command_suffix, shell=True)
      shutil.copyfile(current_tree_name+'.output_tree',current_tree_name)
      shutil.copyfile(starting_base_filename+".start", starting_base_filename+".gaps.snp_sites.aln")
      GubbinsCommon.reinsert_gaps_into_fasta_file(current_tree_name+'.seq.joint.txt', starting_base_filename +".gaps.vcf", starting_base_filename+".gaps.snp_sites.aln")

      if self.args.verbose > 0:
        print int(time.time())

      if self.args.verbose > 0:
        print gubbins_command
      subprocess.check_call(gubbins_command, shell=True)
      if self.args.verbose > 0:
        print int(time.time())

      tree_file_names.append(current_tree_name)
      if i > 2:
        if GubbinsCommon.has_tree_been_seen_before(tree_file_names):
          if self.args.verbose > 0:
            print "Tree observed before so stopping: "+ str(current_tree_name)
          break

    # cleanup intermediate files
    if self.args.no_cleanup == 0 or self.args.no_cleanup is None:
      max_intermediate_iteration  = max_iteration - 1

      raxml_files_to_delete = GubbinsCommon.raxml_regex_for_file_deletions(base_filename_without_ext,current_time,starting_base_filename, max_intermediate_iteration)
      GubbinsCommon.delete_files_based_on_list_of_regexes('.', raxml_files_to_delete, self.args.verbose)

      fasttree_files_to_delete = GubbinsCommon.fasttree_regex_for_file_deletions(starting_base_filename, max_intermediate_iteration)
      GubbinsCommon.delete_files_based_on_list_of_regexes('.', fasttree_files_to_delete, self.args.verbose)
      shutil.rmtree(temp_working_dir)
      
      GubbinsCommon.delete_files_based_on_list_of_regexes('.', [GubbinsCommon.starting_files_regex("^" + starting_base_filename),"^log.txt"], self.args.verbose)
      
    output_filenames_to_final_filenames = {}
    if self.args.prefix is None:
      self.args.prefix = base_filename_without_ext
    if self.args.tree_builder == "hybrid" or self.args.tree_builder == "raxml":
      output_filenames_to_final_filenames = GubbinsCommon.translation_of_raxml_filenames_to_final_filenames(base_filename_without_ext,current_time, max_iteration - 1, self.args.prefix)
    else:
      output_filenames_to_final_filenames = GubbinsCommon.translation_of_fasttree_filenames_to_final_filenames(starting_base_filename, max_iteration - 1,  self.args.prefix)
    GubbinsCommon.rename_files(output_filenames_to_final_filenames)
  
  @staticmethod
  def robinson_foulds_distance(input_tree_name,output_tree_name):
    input_tree  = dendropy.Tree.get_from_path(input_tree_name, 'newick')
    output_tree = dendropy.Tree.get_from_path(output_tree_name, 'newick')
    return input_tree.robinson_foulds_distance(output_tree)
    
  @staticmethod
  def has_tree_been_seen_before(tree_file_names):
    if len(tree_file_names) <= 2:
      return 0

    for tree_file_name in tree_file_names:
      if tree_file_name is not tree_file_names[-1]:
        current_rf_distance = GubbinsCommon.robinson_foulds_distance(tree_file_name,tree_file_names[-1])
        if(current_rf_distance == 0.0):
          return 1

    return 0

  @staticmethod
  def reroot_tree(tree_name, outgroup):
    if outgroup:
      GubbinsCommon.reroot_tree_with_outgroup(tree_name, outgroup)
    else:
      GubbinsCommon.reroot_tree_at_midpoint(tree_name)
      
  @staticmethod
  def reroot_tree_with_outgroup(tree_name, outgroup):
    tree = Phylo.read(tree_name, 'newick')
    tree.root_with_outgroup({'name': outgroup})
    Phylo.write(tree, tree_name, 'newick')

    tree  = dendropy.Tree.get_from_path(tree_name, 'newick',
              preserve_underscores=True)
    tree.deroot()
    tree.update_splits()
    output_tree_string = tree.as_string(
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
      node_label_compose_func=None
      )
    output_file = open(tree_name, 'w+')
    output_file.write(output_tree_string.replace('\'', ''))
    output_file.closed

  @staticmethod
  def split_all_non_bi_nodes(node):
    if node.is_leaf():
      return None
    elif len(node.child_nodes()) > 2:
      GubbinsCommon.split_child_nodes(node)

    for child_node in node.child_nodes():
      GubbinsCommon.split_all_non_bi_nodes(child_node)

    return None

  @staticmethod
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

  @staticmethod
  def reroot_tree_at_midpoint(tree_name):
    tree  = dendropy.Tree.get_from_path(tree_name, 'newick',
              preserve_underscores=True)
    GubbinsCommon.split_all_non_bi_nodes(tree.seed_node)

    tree.reroot_at_midpoint(update_splits=True, delete_outdegree_one=False)
    tree.deroot()
    tree.update_splits()
    output_tree_string = tree.as_string(
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
      node_label_compose_func=None
      )
    output_file = open(tree_name, 'w+')
    output_file.write(output_tree_string.replace('\'', ''))
    output_file.closed

  @staticmethod
  def raxml_base_name(base_filename_without_ext,current_time):
    return base_filename_without_ext+"."+str(current_time) +"iteration_"

  @staticmethod
  def raxml_current_tree_name(base_filename_without_ext,current_time, i):
    return "RAxML_result."+GubbinsCommon.raxml_base_name(base_filename_without_ext,current_time)+str(i)

  @staticmethod
  def raxml_previous_tree_name(base_filename_without_ext,base_filename, current_time,i):
    previous_tree_name = base_filename

    if i> 1:
      previous_tree_name = "RAxML_result."+GubbinsCommon.raxml_base_name(base_filename_without_ext,current_time)+ str(i-1)
    return previous_tree_name

  @staticmethod
  def raxml_previous_tree(base_filename_without_ext, base_filename, current_time,i,previous_tree_name):
    previous_tree = ""

    if i> 1:
      previous_tree = "-t "+ previous_tree_name
    return previous_tree

  @staticmethod
  def raxml_tree_building_command(i,base_filename_without_ext,base_filename,current_time, raxml_exec,previous_tree_name, verbose):
    previous_tree = GubbinsCommon.raxml_previous_tree(base_filename_without_ext, base_filename, current_time,i,previous_tree_name)

    command_suffix = '> /dev/null 2>&1'
    if verbose > 0:
      command_suffix = ''

    return raxml_exec+ " -s "+previous_tree_name+".phylip -n "+base_filename_without_ext+"."+str(current_time)+"iteration_"+str(i)+" "+previous_tree+ command_suffix

  @staticmethod
  def raxml_gubbins_command(base_filename_without_ext,starting_base_filename,current_time, i,alignment_filename,gubbins_exec,min_snps, original_aln):
    current_tree_name = GubbinsCommon.raxml_current_tree_name(base_filename_without_ext,current_time, i)
    return gubbins_exec+" -r -v "+starting_base_filename+".vcf -f "+original_aln+" -t "+str(current_tree_name)+" -m "+ str(min_snps)+" "+ starting_base_filename+".snp_sites.aln"

  @staticmethod
  def raxml_regex_for_file_deletions(base_filename_without_ext,current_time,starting_base_filename, max_intermediate_iteration):
    regex_for_file_deletions = []
    # Can delete all of these files
    regex_for_file_deletions.append("^RAxML_(bestTree|info|log|parsimonyTree)."+GubbinsCommon.raxml_base_name(base_filename_without_ext,current_time))

    # loop over previous iterations and delete
    for file_iteration in range(1,max_intermediate_iteration):
      regex_for_file_deletions.append("^RAxML_result."+GubbinsCommon.raxml_base_name(base_filename_without_ext,current_time)+str(file_iteration)+'\.')
      regex_for_file_deletions.append("^RAxML_result."+GubbinsCommon.raxml_base_name(base_filename_without_ext,current_time)+str(file_iteration)+'$')

    regex_for_file_deletions.append("^RAxML_result."+GubbinsCommon.raxml_base_name(base_filename_without_ext,current_time)+str(max_intermediate_iteration)+'.ancestor.tre')
    regex_for_file_deletions.append("^RAxML_result."+GubbinsCommon.raxml_base_name(base_filename_without_ext,current_time)+str(max_intermediate_iteration)+'.seq.joint.txt')
    regex_for_file_deletions.append("^RAxML_result."+GubbinsCommon.raxml_base_name(base_filename_without_ext,current_time)+str(max_intermediate_iteration)+'.prob.joint.txt') 
    return regex_for_file_deletions

  @staticmethod
  def fasttree_regex_for_file_deletions(starting_base_filename, max_intermediate_iteration):
    regex_for_file_deletions = []

    # loop over previous iterations and delete
    for file_iteration in range(1,max_intermediate_iteration):
      regex_for_file_deletions.append("^"+starting_base_filename+".iteration_"+str(file_iteration)+'[$|\.]')

    return regex_for_file_deletions
  
  @staticmethod
  def rename_files(input_to_output_filenames):
    for input_file in input_to_output_filenames:
      shutil.move(input_file, input_to_output_filenames[input_file])

  @staticmethod
  def starting_files_regex(starting_base_filename):
    # starting file with gapped and ungapped snps
    return starting_base_filename+".(gaps|vcf|snp_sites|phylip|start)"

  @staticmethod
  def translation_of_fasttree_filenames_to_final_filenames(starting_base_filename, max_intermediate_iteration, output_prefix):
    return GubbinsCommon.translation_of_filenames_to_final_filenames(starting_base_filename+".iteration_"+str(max_intermediate_iteration), output_prefix)
    
  @staticmethod
  def translation_of_raxml_filenames_to_final_filenames(base_filename_without_ext,current_time, max_intermediate_iteration, output_prefix):
    return GubbinsCommon.translation_of_filenames_to_final_filenames("RAxML_result."+GubbinsCommon.raxml_base_name(base_filename_without_ext,current_time)+str(max_intermediate_iteration), output_prefix)
  
  @staticmethod
  def translation_of_filenames_to_final_filenames(input_prefix, output_prefix):
    input_names_to_output_names = {  
      str(input_prefix)+".vcf":             str(output_prefix)+".summary_of_snp_distribution.vcf"  ,
      str(input_prefix)+".branch_snps.tab": str(output_prefix)+".branch_base_reconstruction.embl"  ,
      str(input_prefix)+".tab":             str(output_prefix)+".recombination_predictions.embl"   ,
      str(input_prefix)+".gff":             str(output_prefix)+".recombination_predictions.gff"    ,
      str(input_prefix)+".stats":           str(output_prefix)+".per_branch_statistics.csv"        ,
      str(input_prefix)+".snp_sites.aln":   str(output_prefix)+".filtered_polymorphic_sites.fasta" ,
      str(input_prefix)+".phylip":          str(output_prefix)+".filtered_polymorphic_sites.phylip",
      str(input_prefix)+".output_tree":     str(output_prefix)+".node_labelled.tre",
      str(input_prefix):                    str(output_prefix)+".final_tree.tre"
    }
    return input_names_to_output_names

  @staticmethod
  def fasttree_current_tree_name(base_filename, i):
    return base_filename+".iteration_"+str(i)

  @staticmethod
  def fasttree_previous_tree_name(base_filename, i):
    return base_filename+".iteration_"+str(i-1)

  @staticmethod
  def fasttree_tree_building_command(i, starting_tree, current_tree_name,starting_base_filename, previous_tree_name,fasttree_exec, fasttree_params,base_filename ):
    current_tree_name = GubbinsCommon.fasttree_current_tree_name(base_filename, i)

    input_tree = ""
    if starting_tree is not None:
      input_tree = " -intree " + starting_tree
    elif i > 1 :
      input_tree = " -intree "+ previous_tree_name

    return fasttree_exec+" "+ input_tree+" "+ fasttree_params+" "+ starting_base_filename+".snp_sites.aln   > "+current_tree_name

  @staticmethod
  def  fasttree_gubbins_command(base_filename,starting_base_filename, i,alignment_filename,gubbins_exec,min_snps,original_aln):
    current_tree_name = GubbinsCommon.fasttree_current_tree_name(base_filename, i)
    return gubbins_exec+" -r -v "+starting_base_filename+".vcf -f "+original_aln+" -t "+str(current_tree_name)+" -m "+ str(min_snps)+" "+ starting_base_filename+".snp_sites.aln"

  @staticmethod
  def fasttree_fastml_command(fastml_exec, alignment_filename, base_filename,i):
    current_tree_name = GubbinsCommon.fasttree_current_tree_name(base_filename, i)
    return GubbinsCommon.generate_fastml_command(fastml_exec, alignment_filename, current_tree_name)

  @staticmethod
  def raxml_fastml_command(fastml_exec, alignment_filename, base_filename_without_ext,current_time, i):
    current_tree_name = GubbinsCommon.raxml_current_tree_name(base_filename_without_ext,current_time, i)
    return GubbinsCommon.generate_fastml_command(fastml_exec, alignment_filename, current_tree_name)

  @staticmethod
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

  @staticmethod
  def number_of_sequences_in_alignment(filename):
    return len(GubbinsCommon.get_sequence_names_from_alignment(filename))

  @staticmethod
  def get_sequence_names_from_alignment(filename):
    sequence_names = []
    handle = open(filename, "rU")
    for record in SeqIO.parse(handle, "fasta") :
      sequence_names.append(record.id)
    handle.close()
    return sequence_names

  @staticmethod
  def is_input_fasta_file_valid(input_filename):
    try:
      if GubbinsCommon.does_each_sequence_have_the_same_length(input_filename) == 0:
        print "Each sequence must be the same length"
        return 0
      if GubbinsCommon.are_sequence_names_unique(input_filename) == 0:
        print "All sequence names in the fasta file must be unique"
        return 0
      if GubbinsCommon.does_each_sequence_have_a_name_and_genomic_data(input_filename) == 0:
        print "Each sequence must have a name and some genomic data"
        return 0
    except:
      return 0

    return 1
    
  @staticmethod
  def is_input_starting_tree_valid(starting_tree):
    try:
      if GubbinsCommon.is_starting_tree_valid(starting_tree) ==0:
        return 0
      if GubbinsCommon.do_all_leaf_nodes_have_unique_names(starting_tree) == 0:
        return 0
    except:
      return 0

    return 1
  

  @staticmethod
  def does_each_sequence_have_a_name_and_genomic_data(input_filename):
    input_handle  = open(input_filename, "rU")
    alignments = AlignIO.parse(input_handle, "fasta")
    number_of_sequences = 0
    for alignment in alignments:
        for record in alignment:
            number_of_sequences +=1
            if record.name is None or record.name == "":
              print "Error with the input FASTA file: One of the sequence names is blank"
              return 0
            if record.seq is None or record.seq == "":
              print "Error with the input FASTA file: One of the sequences is empty"
              return 0
            if re.search('[^ACGTNacgtn-]', str(record.seq))  != None:
              print "Error with the input FASTA file: One of the sequences contains odd characters, only ACGTNacgtn- are permitted"
              return 0
    if number_of_sequences <= 3:
      print "Error with input FASTA file: you need more than 3 sequences to build a meaningful tree"
      return 0
    input_handle.close()
    return 1
    
    
  @staticmethod
  def does_each_sequence_have_the_same_length(input_filename):
    try:
      input_handle  = open(input_filename, "rU")
      alignments = AlignIO.parse(input_handle, "fasta")
      sequence_length = -1
      for alignment in alignments:
          for record in alignment:
             if sequence_length == -1:
               sequence_length = len(record.seq)
             elif sequence_length != len(record.seq):
               print "Error with the input FASTA file: The sequences dont have the same lengths this isnt an alignment: "+record.name
               return 0
      input_handle.close()
    except:
      print "Error with the input FASTA file: It is in the wrong format so check its an alignment"
      return 0
    return 1

  @staticmethod
  def are_sequence_names_unique(input_filename):
    input_handle  = open(input_filename, "rU")
    alignments = AlignIO.parse(input_handle, "fasta")
    sequence_names = []
    for alignment in alignments:
        for record in alignment:
            sequence_names.append(record.name)
            
    if [k for k,v in Counter(sequence_names).items() if v>1] != []:
      return 0
    input_handle.close()
    return 1

  @staticmethod
  def filter_out_alignments_with_too_much_missing_data(input_filename, output_filename, filter_percentage,verbose):
    input_handle  = open(input_filename, "rU")
    output_handle = open(output_filename, "w+")
    alignments = AlignIO.parse(input_handle, "fasta")
    output_alignments = []
    taxa_removed = []
    number_of_included_alignments = 0
    for alignment in alignments:
        for record in alignment:
          number_of_gaps = 0
          number_of_gaps += record.seq.count('n')
          number_of_gaps += record.seq.count('N')
          number_of_gaps += record.seq.count('-')
          sequence_length = len(record.seq)

          if sequence_length == 0:
            taxa_removed.append(record.id)
            print "Excluded sequence " + record.id + " because there werent enough bases in it"
          elif((number_of_gaps*100/sequence_length) <= filter_percentage):
            output_alignments.append(record)
            number_of_included_alignments += 1
          else:
            taxa_removed.append(record.id)
            print "Excluded sequence " + record.id + " because it had " + str(number_of_gaps*100/sequence_length) +" percentage gaps while a maximum of "+ str(filter_percentage) +" is allowed"

    if number_of_included_alignments <= 1:
      sys.exit("Too many sequences have been excluded so theres no data left to work with. Please increase the -f parameter")

    AlignIO.write(MultipleSeqAlignment(output_alignments), output_handle, "fasta")
    output_handle.close()
    input_handle.close()
    return taxa_removed

  @staticmethod
  def filter_out_removed_taxa_from_tree_and_return_new_file(starting_tree, temp_working_dir, taxa_removed):
    if starting_tree is None:
       return None

    (base_directory,base_filename) = os.path.split(starting_tree)
    temp_starting_tree = temp_working_dir + '/'+ base_filename

    tree  = dendropy.Tree.get_from_path(starting_tree, 'newick',
              preserve_underscores=True)

    tree.prune_taxa_with_labels(taxa_removed, update_splits=True, delete_outdegree_one=False)          
    tree.prune_leaves_without_taxa(update_splits=True, delete_outdegree_one=False)
    tree.deroot()
    output_tree_string = tree.as_string(
      'newick',
      taxon_set=None,
      suppress_leaf_taxon_labels=False,
      suppress_leaf_node_labels=True,
      suppress_internal_taxon_labels=True,
      suppress_internal_node_labels=True,
      suppress_rooting=True,
      suppress_edge_lengths=False,
      unquoted_underscores=True,
      preserve_spaces=False,
      store_tree_weights=False,
      suppress_annotations=True,
      annotations_as_nhx=False,
      suppress_item_comments=True,
      node_label_element_separator=' ',
      node_label_compose_func=None
      )
    output_file = open(temp_starting_tree, 'w+')
    output_file.write(output_tree_string.replace('\'', ''))
    output_file.closed

    return temp_starting_tree

  @staticmethod
  def reinsert_gaps_into_fasta_file(input_fasta_filename, input_vcf_file, output_fasta_filename):
    # find out where the gaps are located
    # PyVCF removed for performance reasons
    vcf_file = open(input_vcf_file, 'r')

    sample_names  = []
    gap_position = []
    gap_alt_base = []

    for vcf_line in vcf_file:
      if re.match('^#CHROM', vcf_line)  != None :
         sample_names = vcf_line.split('\t' )[9:-1]
      elif re.match('^\d', vcf_line)  != None :
        # If the alternate is only a gap it wont have a base in this column
        if  re.match('^([^\t]+\t){3}([ACGTacgt])\t([^ACGTacgt])\t', vcf_line)  != None:
          m = re.match('^([^\t]+\t){3}([ACGTacgt])\t([^ACGTacgt])\t', vcf_line) 
          gap_position.append(1)
          gap_alt_base.append(m.group(2))
        elif re.match('^([^\t]+\t){3}([^ACGTacgt])\t([ACGTacgt])\t', vcf_line)  != None:
          # sometimes the ref can be a gap only 
          m = re.match('^([^\t]+\t){3}([^ACGTacgt])\t([ACGTacgt])\t', vcf_line) 
          gap_position.append(1)
          gap_alt_base.append(m.group(3))
        else:
          gap_position.append(0)
          gap_alt_base.append('-')

    gapped_alignments = []
    # interleave gap only and snp bases
    input_handle = open(input_fasta_filename, "rU")
    alignments = AlignIO.parse(input_handle, "fasta")
    for alignment in alignments:
      for record in alignment:
        inserted_gaps = []
        if record.id in sample_names:
          # only apply to internal nodes
          continue
        gap_index = 0
        for input_base in record.seq:
          while gap_index < len(gap_position) and gap_position[gap_index] == 1:
            inserted_gaps.append(gap_alt_base[gap_index])
            gap_index+=1
          if gap_index < len(gap_position):
            inserted_gaps.append(input_base)
            gap_index+=1

        while gap_index < len(gap_position):
          inserted_gaps.append(gap_alt_base[gap_index])
          gap_index+=1

        record.seq = Seq(''.join(inserted_gaps))
        gapped_alignments.append(record)

    output_handle = open(output_fasta_filename, "a")
    AlignIO.write(MultipleSeqAlignment(gapped_alignments), output_handle, "fasta")

    return


    # reparsing a fasta file splits the lines which makes fastml work
  @staticmethod
  def reconvert_fasta_file(input_filename, output_filename):
    input_handle = open(input_filename, "rU")
    output_handle = open(output_filename, "w+")
    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, "fasta")
    output_handle.close()
    input_handle.close()
    return

  @staticmethod
  def pairwise_comparison(filename,base_filename,gubbins_exec,alignment_filename,fastml_exec):
    sequence_names = GubbinsCommon.get_sequence_names_from_alignment(filename)
    GubbinsCommon.create_pairwise_newick_tree(sequence_names, base_filename+".tre")

    subprocess.check_call(GubbinsCommon.generate_fastml_command(fastml_exec, base_filename+".gaps.snp_sites.aln", base_filename+".tre"), shell=True)
    shutil.copyfile(base_filename+'.tre.output_tree',base_filename+".tre")
    shutil.copyfile(base_filename+'.tre.seq.joint.txt', base_filename+".snp_sites.aln")
    subprocess.check_call(gubbins_exec+" -r -v "+base_filename+".vcf -t "+base_filename+".tre -f "+ alignment_filename +" "+ base_filename+".snp_sites.aln", shell=True)

  @staticmethod
  def create_pairwise_newick_tree(sequence_names, output_filename):
    tree = Phylo.read(StringIO('('+sequence_names[0]+','+sequence_names[1]+')'), "newick")
    Phylo.write(tree, output_filename, 'newick')

  @staticmethod
  def delete_files_based_on_list_of_regexes(directory_to_search, regex_for_file_deletions, verbose):
    for dirname, dirnames, filenames in os.walk(directory_to_search):
      for filename in filenames:
        for deletion_regex in regex_for_file_deletions:
          full_path_of_file_for_deletion = os.path.join(directory_to_search, filename)
          if(re.match(str(deletion_regex), filename) != None and os.path.exists(full_path_of_file_for_deletion)):
            if verbose > 0:
              print "Deleting file: "+ os.path.join(directory_to_search, filename) + " regex:"+deletion_regex
            os.remove(full_path_of_file_for_deletion)

  @staticmethod
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

  @staticmethod
  def use_bundled_exec(input_executable, bundled_executable):
    # Pop the first value off the input_executable and replace it with the bundled exec
    executable_and_params = input_executable.split(" ")
    executable_and_params.pop(0)
    executable_and_params.insert(0, bundled_executable)
    return  ' '.join(executable_and_params)

