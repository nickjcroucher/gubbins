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

import os
import sys
import tempfile
import dendropy
import subprocess
import shutil
import time
from random import randint
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

class RAxMLSequenceReconstruction(object):
	def __init__(self, input_alignment_filename, input_tree, output_alignment_filename, output_tree, raxml_internal_sequence_reconstruction_command, verbose = False ):
		self.input_alignment_filename = os.path.abspath(input_alignment_filename)
		self.input_tree = os.path.abspath(input_tree)
		self.output_alignment_filename = os.path.abspath(output_alignment_filename)
		self.output_tree = os.path.abspath(output_tree)
		self.raxml_internal_sequence_reconstruction_command = raxml_internal_sequence_reconstruction_command
		self.verbose = verbose
		
		self.working_dir = tempfile.mkdtemp(dir=os.getcwd())
		self.temp_rooted_tree = self.working_dir +'/' +'rooted_tree.newick'
		self.temp_interal_fasta = self.working_dir +'/' +'internal.fasta'
		self.internal_node_prefix = 'internal_'
	
	def reconstruct_ancestor_sequences(self):
		self.root_tree(self.input_tree, self.temp_rooted_tree)
		
		self.run_raxml_ancestor_command(self.temp_rooted_tree)
		self.convert_raw_ancestral_states_to_fasta(self.raw_internal_sequence_filename(), self.temp_interal_fasta)
		self.combine_fastas(self.input_alignment_filename, self.temp_interal_fasta,self.output_alignment_filename)
		
		if os.path.exists(self.temp_rooted_tree):
			self.transfer_internal_names_to_tree(self.raw_internal_rooted_tree_filename(), self.temp_rooted_tree, self.output_tree)

		shutil.rmtree(self.working_dir)
	
	def run_raxml_ancestor_command(self,rooted_tree):
		current_directory = os.getcwd()
		if self.verbose > 0:
			print(self.raxml_reconstruction_command(rooted_tree))
		try:
			os.chdir(self.working_dir)
			subprocess.check_call(self.raxml_reconstruction_command(rooted_tree), shell=True)
			os.chdir(current_directory)
		except:
			os.chdir(current_directory)
			sys.exit("Something went wrong while creating the ancestor sequences using RAxML")
		if self.verbose > 0:
			print(int(time.time()))
	
	def raw_internal_sequence_filename(self):
		return self.working_dir +'/RAxML_marginalAncestralStates.internal'
	
	def raw_internal_rooted_tree_filename(self):
		return self.working_dir +'/RAxML_nodeLabelledRootedTree.internal'
	
	def raxml_reconstruction_command(self,rooted_tree):
		verbose_suffix = ''
		if not self.verbose:
			verbose_suffix = '> /dev/null 2>&1'
		
		return " ".join([self.raxml_internal_sequence_reconstruction_command, ' -s', self.input_alignment_filename, '-t', rooted_tree, '-n', 'internal' ,verbose_suffix ])
	
	def write_tree(self, tree, output_tree):
		output_tree_string = tree.as_string(
			schema='newick',
			suppress_leaf_taxon_labels=False,
			suppress_leaf_node_labels=True,
			suppress_internal_taxon_labels=False,
			suppress_internal_node_labels=False,
			suppress_rooting=False,
			suppress_edge_lengths=False,
			unquoted_underscores=True,
			preserve_spaces=False,
			store_tree_weights=False,
			suppress_annotations=True,
			annotations_as_nhx=False,
			suppress_item_comments=True,
			node_label_element_separator=' '
			)
		with open(output_tree, 'w+') as output_file:
			output_file.write(output_tree_string.replace('\'', ''))
			output_file.closed
			
			return output_tree
			
	def transfer_internal_names_to_tree(self, source_tree, destination_tree, output_tree):
		source_tree_obj  = dendropy.Tree.get_from_path(source_tree, 'newick', preserve_underscores=True)
		source_internal_node_labels = []
		for source_internal_node in source_tree_obj.internal_nodes():
			if source_internal_node.label:
				source_internal_node_labels.append(source_internal_node.label)
			else:
				source_internal_node_labels.append('')
		
		destination_tree_obj  = dendropy.Tree.get_from_path(destination_tree, 'newick', preserve_underscores=True)
		for index, destination_internal_node in enumerate(destination_tree_obj.internal_nodes()):
			destination_internal_node.label = None
			destination_internal_node.taxon = dendropy.Taxon(self.internal_node_prefix + str(source_internal_node_labels[index]))
		self.write_tree( destination_tree_obj, output_tree)
		
	
	def root_tree(self, input_tree_filename, output_tree):
		# split bi nodes and root tree
		tree  = dendropy.Tree.get_from_path(input_tree_filename, 'newick', preserve_underscores=True)
		self.split_all_non_bi_nodes(tree.seed_node)
		self.write_tree( tree, output_tree)

	
	def convert_raw_ancestral_states_to_fasta(self, input_filename, output_filename):
		with open(input_filename, 'r') as infile:
			with open(output_filename, 'w+') as outfile:
				for sequence_line in infile:
					[sequence_name, sequence_bases] = sequence_line.split(' ')
					sequence_bases = sequence_bases.replace('?', 'N')
					outfile.write('>'+sequence_name+'\n')
					outfile.write(sequence_bases)
    
    # Warning - recursion
	def split_all_non_bi_nodes(self, node):
		if node.is_leaf():
			return None
		elif len(node.child_nodes()) > 2:
			self.split_child_nodes(node)
		for child_node in node.child_nodes():
			self.split_all_non_bi_nodes(child_node)
		return None
	
	def split_child_nodes(self,node):
		all_child_nodes = node.child_nodes()
		first_child = all_child_nodes.pop()
		new_child_node = node.new_child(edge_length=0)
		new_child_node.set_child_nodes(all_child_nodes)
		node.set_child_nodes((first_child,new_child_node))
	
	def combine_fastas(self, leaf_node_filename, internl_node_filename, output_file ):
		with open(output_file, 'w') as output_handle:
			# print out leafnodes as is
			with open(leaf_node_filename, 'r') as input_handle:
				alignments = AlignIO.parse(input_handle, "fasta")
				AlignIO.write(alignments,output_handle, "fasta")
				input_handle.closed
			
			with open(internl_node_filename, 'r') as input_handle:
				alignments = AlignIO.parse(input_handle, "fasta")
				output_alignments = []
				for alignment in alignments:
					for record in alignment:
						record.id = self.internal_node_prefix + str(record.id)
						record.description = ''
						output_alignments.append(record)
				
				AlignIO.write(MultipleSeqAlignment(output_alignments),output_handle, "fasta")
				input_handle.closed
				output_handle.closed
	 