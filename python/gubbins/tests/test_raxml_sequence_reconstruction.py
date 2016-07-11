#! /usr/bin/env python3
# encoding: utf-8
#  python3 setup.py test -s gubbins.tests.test_raxml_sequence_reconstruction.TestRAxMLSequenceReconstruction

"""
Tests for reconstructing internal sequences using raxml
"""

import unittest
import re
import os
import subprocess
import filecmp
import shutil
import dendropy
from gubbins.RAxMLSequenceReconstruction import RAxMLSequenceReconstruction
from gubbins.RAxMLExecutable import RAxMLExecutable

os.environ["PATH"] += os.pathsep + 'gubbins/tests/bin'

class TestRAxMLSequenceReconstruction(unittest.TestCase):
	
	def test_ancestor_raxml_command_no_verbose(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('input_alignment.fasta', 'input_tree',
			'output_alignment_filename', 'output_tree',
			'raxmlHPC -f A -p 1 -m GTRGAMMA',
			verbose = False)
		self.assertEqual(raxml_seq_recon.raxml_reconstruction_command(), 'raxmlHPC -f A -p 1 -m GTRGAMMA  -s '+raxml_seq_recon.input_alignment_filename + ' -t ' + raxml_seq_recon.working_dir+'/rooted_tree.newick -n internal ')
		self.cleanup()
	
	def test_ancestor_raxml_command_verbose(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('input_alignment.fasta', 'input_tree',
			'output_alignment_filename', 'output_tree',
			'raxmlHPC -f A -p 1 -m GTRGAMMA',
			verbose = True)
		self.assertEqual(raxml_seq_recon.raxml_reconstruction_command(), 'raxmlHPC -f A -p 1 -m GTRGAMMA  -s '+raxml_seq_recon.input_alignment_filename+' -t ' + raxml_seq_recon.working_dir+'/rooted_tree.newick -n internal > /dev/null 2>&1')
		self.cleanup()
	
	def test_working_directory_construction(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('', '', '', '', '', False)
		self.assertTrue( os.path.exists(raxml_seq_recon.working_dir) )
		self.cleanup()
	
	def test_root_input_tree(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('abc', 'gubbins/tests/data/raxml_sequence_reconstruction/unrooted_tree.newick', 'abc', 'abc', '', False)
		output_tree = raxml_seq_recon.root_tree_and_label_internal_nodes('gubbins/tests/data/raxml_sequence_reconstruction/unrooted_tree.newick',raxml_seq_recon.temp_rooted_tree)
		self.assertTrue(filecmp.cmp(str(raxml_seq_recon.temp_rooted_tree), 'gubbins/tests/data/raxml_sequence_reconstruction/expected_rooted_tree.newick'))
		self.cleanup()
	
	def test_run_raxml_ancestor_reconstruction(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('gubbins/tests/data/raxml_sequence_reconstruction/input_alignment.fasta',
			'gubbins/tests/data/raxml_sequence_reconstruction/unrooted_tree.newick',
			'outputfile', 'output_tree', RAxMLExecutable(1).internal_sequence_reconstruction_command(), False)
		raxml_seq_recon.reconstruct_ancestor_sequences()

		assert os.path.exists('outputfile')
		self.cleanup()
	
	def test_convert_raw_ancestral_file_to_fasta(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('', '', '', '', '', False)
		raxml_seq_recon.convert_raw_ancestral_states_to_fasta('gubbins/tests/data/raxml_sequence_reconstruction/raw_marginalAncestralStates.phylip','outputfile')
		self.assertTrue(filecmp.cmp('outputfile','gubbins/tests/data/raxml_sequence_reconstruction/expected_marginalAncestralStates.fasta'))
		self.cleanup()
	
	def test_merging_fasta_files(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('', '', '', '', '', False)
		raxml_seq_recon.combine_fastas('gubbins/tests/data/raxml_sequence_reconstruction/1.fasta','gubbins/tests/data/raxml_sequence_reconstruction/2.fasta','combined.fasta')
		self.assertTrue(filecmp.cmp('combined.fasta','gubbins/tests/data/raxml_sequence_reconstruction/expected_combined_1_2.fasta'))
		self.cleanup()
		
	def test_add_labels_to_tree(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('', '', '', '', '', False)
		raxml_seq_recon.root_tree_and_label_internal_nodes('gubbins/tests/data/raxml_sequence_reconstruction/unrooted_tree.newick', raxml_seq_recon.temp_rooted_tree)
		
		tree  = dendropy.Tree.get_from_path(raxml_seq_recon.temp_rooted_tree, 'newick', preserve_underscores=True)
		self.assertEqual("((B:0.1,(C:0.1,(D:0.1,E:0.1)10)9)8,(A:0.1,F:0.1)7:0.0)ROOT;\n",tree.as_string(schema='newick'))
	
	def test_transfer_internal_labels(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('', '', '', 'output_tree', '', False)
		raxml_seq_recon.transfer_internal_names_to_tree('gubbins/tests/data/source_tree.tre', 'gubbins/tests/data/destination_tree.tre', 'renamed_output_tree')
		assert os.path.exists('renamed_output_tree')
		self.assertTrue(filecmp.cmp('renamed_output_tree','gubbins/tests/data/expected_renamed_output_tree'))
		os.remove('renamed_output_tree')
	
	def test_more_complex_tree(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('gubbins/tests/data/multiple_recombinations.aln',
			'gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.output_tree',
			'output_alignment', 'output_tree', RAxMLExecutable(1).internal_sequence_reconstruction_command(), False)
		raxml_seq_recon.reconstruct_ancestor_sequences()

		assert os.path.exists('output_alignment')
		assert os.path.exists('output_tree')
		self.cleanup()
	
	def cleanup(self):
		for file_to_delete in ['combined.fasta','outputfile','RAxML_nodeLabelledRootedTree.internal','RAxML_marginalAncestralProbabilities.internal', 'RAxML_info.internal', 'RAxML_flagCheck', 'output_tree']:
			if os.path.exists(file_to_delete):
				os.remove(file_to_delete)

if __name__ == "__main__":
  unittest.main()
  
  
  
  
  