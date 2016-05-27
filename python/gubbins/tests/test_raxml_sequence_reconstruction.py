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
from gubbins.RAxMLSequenceReconstruction import RAxMLSequenceReconstruction
from gubbins.RAxMLExecutable import RAxMLExecutable

os.environ["PATH"] += os.pathsep + 'gubbins/tests/bin'

class TestRAxMLSequenceReconstruction(unittest.TestCase):
	
	def test_ancestor_raxml_command_no_verbose(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('input_alignment.fasta', 'input_tree',
			'output_alignment_filename', 'output_tree',
			'raxmlHPC -f A -p 1 -m GTRGAMMA',
			verbose = False)
		self.assertEqual(raxml_seq_recon.raxml_reconstruction_command(), 'raxmlHPC -f A -p 1 -m GTRGAMMA  -s '+raxml_seq_recon.input_alignment_filename+' -t '+raxml_seq_recon.working_dir+'/rooted_tree.newick -n internal ')
		self.cleanup()
	
	def test_ancestor_raxml_command_verbose(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('input_alignment.fasta', 'input_tree',
			'output_alignment_filename', 'output_tree',
			'raxmlHPC -f A -p 1 -m GTRGAMMA',
			verbose = True)
		self.assertEqual(raxml_seq_recon.raxml_reconstruction_command(), 'raxmlHPC -f A -p 1 -m GTRGAMMA  -s '+raxml_seq_recon.input_alignment_filename+' -t '+raxml_seq_recon.working_dir+'/rooted_tree.newick -n internal > /dev/null 2>&1')
		self.cleanup()
	
	def test_working_directory_construction(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('', '', '', '', '', False)
		self.assertTrue( os.path.exists(raxml_seq_recon.working_dir) )
		self.cleanup()
	
	def test_root_input_tree(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('abc', 'gubbins/tests/data/raxml_sequence_reconstruction/unrooted_tree.newick', 'abc', 'abc', '', False)
		output_tree = raxml_seq_recon.root_tree('gubbins/tests/data/raxml_sequence_reconstruction/unrooted_tree.newick',raxml_seq_recon.temp_rooted_tree)
		self.assertTrue(filecmp.cmp(str(raxml_seq_recon.temp_rooted_tree), 'gubbins/tests/data/raxml_sequence_reconstruction/expected_rooted_tree.newick'))
		self.cleanup()
	
	def test_run_raxml_ancestor_reconstruction(self):
		raxml_seq_recon = RAxMLSequenceReconstruction('gubbins/tests/data/raxml_sequence_reconstruction/input_alignment.fasta',
			'gubbins/tests/data/raxml_sequence_reconstruction/unrooted_tree.newick',
			'outputfile', 'output_tree', RAxMLExecutable(1).internal_sequence_reconstruction_command(), False)
		output_tree = raxml_seq_recon.root_tree('gubbins/tests/data/raxml_sequence_reconstruction/unrooted_tree.newick',raxml_seq_recon.temp_rooted_tree)
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
	
	def cleanup(self):
		for file_to_delete in ['combined.fasta','outputfile','RAxML_nodeLabelledRootedTree.internal','RAxML_marginalAncestralProbabilities.internal', 'RAxML_info.internal', 'RAxML_flagCheck', 'output_tree']:
			if os.path.exists(file_to_delete):
				os.remove(file_to_delete)

if __name__ == "__main__":
  unittest.main()
  
  
  
  
  