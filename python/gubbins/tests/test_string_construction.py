#! /usr/bin/env python3
# encoding: utf-8

"""
Tests the construction of file names and commands with no dependencies on external applications.
"""

import unittest
import os
import subprocess
from gubbins import common


class TestStringConstruction(unittest.TestCase):
  def test_raxml_base_name(self):
    assert common.GubbinsCommon.raxml_base_name('AAA', 1234) == 'AAA.1234iteration_'
    assert common.GubbinsCommon.raxml_base_name('AAA', '') == 'AAA.iteration_'

  def test_raxml_current_tree_name(self):
    assert common.GubbinsCommon.raxml_current_tree_name('AAA', 1234, 5) == 'RAxML_result.AAA.1234iteration_5'
    assert common.GubbinsCommon.raxml_current_tree_name('AAA', '', 5) == 'RAxML_result.AAA.iteration_5'

  def test_raxml_previous_tree_name(self):
    assert common.GubbinsCommon.raxml_previous_tree_name('AAA', 'BBB', 1234, 5) == 'RAxML_result.AAA.1234iteration_4'
    assert common.GubbinsCommon.raxml_previous_tree_name('AAA', 'BBB', '', 5) == 'RAxML_result.AAA.iteration_4'

  def test_raxml_previous_tree(self):
    assert common.GubbinsCommon.raxml_previous_tree('AAA', 'BBB', 1234, 5, 'ABC') == '-t ABC'
    assert common.GubbinsCommon.raxml_previous_tree('AAA', 'BBB', 1234, 1, 'ABC') == ''

  def test_raxml_tree_building_command(self):
    assert common.GubbinsCommon.raxml_tree_building_command(5, 'AAA', 'BBB', 1234, 'raxml_exec', 'CCC', 0) \
           == 'raxml_exec -s CCC.phylip -n AAA.1234iteration_5 -t CCC> /dev/null 2>&1'
    assert common.GubbinsCommon.raxml_tree_building_command(1, 'AAA', 'BBB', 1234, 'raxml_exec', 'CCC', 0) \
           == 'raxml_exec -s CCC.phylip -n AAA.1234iteration_1 > /dev/null 2>&1'
    assert common.GubbinsCommon.raxml_tree_building_command(2, 'AAA', 'BBB', '', 'raxml_exec', 'CCC', 0) \
           == 'raxml_exec -s CCC.phylip -n AAA.iteration_2 -t CCC> /dev/null 2>&1'
    assert common.GubbinsCommon.raxml_tree_building_command(2, 'AAA', 'BBB', '', 'raxml_exec', 'CCC', 1) \
           == 'raxml_exec -s CCC.phylip -n AAA.iteration_2 -t CCC'

  def test_raxml_gubbins_command(self):
    assert common.GubbinsCommon.raxml_gubbins_command('AAA', 'BBB', '', 5, 'CCC', 'DDD', 3, 'EEE', 10, 100) \
           == 'DDD -r -v BBB.vcf -a 10 -b 100 -f EEE -t RAxML_result.AAA.iteration_5 -m 3 BBB.snp_sites.aln'
    assert common.GubbinsCommon.raxml_gubbins_command('base_filename_without_ext', 'starting_base_filename', 1234, 5,
               'alignment_filename', 'gubbins_exec', 3, 'original_aln', 10, 100) \
           == 'gubbins_exec -r -v starting_base_filename.vcf -a 10 -b 100 -f original_aln -t ' \
              'RAxML_result.base_filename_without_ext.1234iteration_5 -m 3 starting_base_filename.snp_sites.aln'

  def test_raxml_regex_for_file_deletions(self):
    assert common.GubbinsCommon.raxml_regex_for_file_deletions('AAA', 1234, 'BBB', 5) \
           == ['^RAxML_(bestTree|info|log|parsimonyTree).AAA.1234iteration_', '^RAxML_result.AAA.1234iteration_1\\.',
               '^RAxML_result.AAA.1234iteration_1$', '^RAxML_result.AAA.1234iteration_2\\.',
               '^RAxML_result.AAA.1234iteration_2$', '^RAxML_result.AAA.1234iteration_3\\.',
               '^RAxML_result.AAA.1234iteration_3$', '^RAxML_result.AAA.1234iteration_4\\.',
               '^RAxML_result.AAA.1234iteration_4$', '^RAxML_result.AAA.1234iteration_5.ancestor.tre',
               '^RAxML_result.AAA.1234iteration_5.seq.joint.txt', '^RAxML_result.AAA.1234iteration_5.prob.joint.txt']

  def test_raxml_regex_for_file_deletions_10(self):
    assert common.GubbinsCommon.raxml_regex_for_file_deletions('AAA', 1234, 'BBB', 10) \
           == ['^RAxML_(bestTree|info|log|parsimonyTree).AAA.1234iteration_', '^RAxML_result.AAA.1234iteration_1\\.',
               '^RAxML_result.AAA.1234iteration_1$', '^RAxML_result.AAA.1234iteration_2\\.',
               '^RAxML_result.AAA.1234iteration_2$', '^RAxML_result.AAA.1234iteration_3\\.',
               '^RAxML_result.AAA.1234iteration_3$', '^RAxML_result.AAA.1234iteration_4\\.',
               '^RAxML_result.AAA.1234iteration_4$', '^RAxML_result.AAA.1234iteration_5\\.',
               '^RAxML_result.AAA.1234iteration_5$', '^RAxML_result.AAA.1234iteration_6\\.',
               '^RAxML_result.AAA.1234iteration_6$', '^RAxML_result.AAA.1234iteration_7\\.',
               '^RAxML_result.AAA.1234iteration_7$', '^RAxML_result.AAA.1234iteration_8\\.',
               '^RAxML_result.AAA.1234iteration_8$', '^RAxML_result.AAA.1234iteration_9\\.',
               '^RAxML_result.AAA.1234iteration_9$', '^RAxML_result.AAA.1234iteration_10.ancestor.tre',
               '^RAxML_result.AAA.1234iteration_10.seq.joint.txt', '^RAxML_result.AAA.1234iteration_10.prob.joint.txt'
               ]

  def test_fasttree_regex_for_file_deletions(self):
    assert common.GubbinsCommon.fasttree_regex_for_file_deletions('AAA', 5) \
        == ['^AAA.iteration_1[$|\\.]', '^AAA.iteration_2[$|\\.]', '^AAA.iteration_3[$|\\.]', '^AAA.iteration_4[$|\\.]']

  def test_starting_files_regex(self):
    assert common.GubbinsCommon.starting_files_regex('AAA') == 'AAA.(gaps|vcf|snp_sites|phylip|start)'
    assert common.GubbinsCommon.starting_files_regex('^') == '^.(gaps|vcf|snp_sites|phylip|start)'
    assert common.GubbinsCommon.starting_files_regex('') == '.(gaps|vcf|snp_sites|phylip|start)'

  def test_translation_of_fasttree_filenames_to_final_filenames(self):
    assert common.GubbinsCommon.translation_of_fasttree_filenames_to_final_filenames('AAA', 5, 'test') == {
      'AAA.iteration_5.vcf':             'test.summary_of_snp_distribution.vcf',
      'AAA.iteration_5.tab':             'test.recombination_predictions.embl',
      'AAA.iteration_5.branch_snps.tab': 'test.branch_base_reconstruction.embl',
      'AAA.iteration_5.gff':             'test.recombination_predictions.gff',
      'AAA.iteration_5.snp_sites.aln':   'test.filtered_polymorphic_sites.fasta',
      'AAA.iteration_5.phylip':          'test.filtered_polymorphic_sites.phylip',
      'AAA.iteration_5.stats':           'test.per_branch_statistics.csv',
      'AAA.iteration_5.output_tree':     'test.node_labelled.tre',
      'AAA.iteration_5':                 'test.final_tree.tre'
    }

  def test_translation_of_raxml_filenames_to_final_filenames(self):
    assert common.GubbinsCommon.translation_of_raxml_filenames_to_final_filenames('AAA', 1234, 10, 'test') == {
      'RAxML_result.AAA.1234iteration_10.snp_sites.aln':   'test.filtered_polymorphic_sites.fasta',
      'RAxML_result.AAA.1234iteration_10.branch_snps.tab': 'test.branch_base_reconstruction.embl',
      'RAxML_result.AAA.1234iteration_10.tab':             'test.recombination_predictions.embl',
      'RAxML_result.AAA.1234iteration_10.gff':             'test.recombination_predictions.gff',
      'RAxML_result.AAA.1234iteration_10':                 'test.final_tree.tre',
      'RAxML_result.AAA.1234iteration_10.phylip':          'test.filtered_polymorphic_sites.phylip',
      'RAxML_result.AAA.1234iteration_10.stats':           'test.per_branch_statistics.csv',
      'RAxML_result.AAA.1234iteration_10.output_tree':     'test.node_labelled.tre',
      'RAxML_result.AAA.1234iteration_10.vcf':             'test.summary_of_snp_distribution.vcf'
    }

  def test_rename_files(self):
    subprocess.check_call('touch temp_file; touch  another_file', shell=True)
    common.GubbinsCommon.rename_files({'temp_file': 'output_file', 'another_file': 'another_output_file'})
    assert os.path.exists('output_file')
    assert os.path.exists('another_output_file')
    os.remove('output_file')
    os.remove('another_output_file')

  def test_fasttree_current_tree_name(self):
    assert common.GubbinsCommon.fasttree_current_tree_name('AAA', 1)  == 'AAA.iteration_1'
    assert common.GubbinsCommon.fasttree_current_tree_name('AAA', 5)  == 'AAA.iteration_5'
    assert common.GubbinsCommon.fasttree_current_tree_name('AAA', '') == 'AAA.iteration_'

  def test_fasttree_previous_tree_name(self):
    assert common.GubbinsCommon.fasttree_previous_tree_name('AAA', 5) == 'AAA.iteration_4'

  def test_fasttree_tree_building_command(self):
    assert common.GubbinsCommon.fasttree_tree_building_command(5, 'AAA', 'BBB', 'CCC', 'EEE', 'FFF', 'GGG', 'HHH') \
           == 'FFF -intree AAA GGG CCC.snp_sites.aln > HHH.iteration_5'

  def test_fasttree_gubbins_command(self):
    assert common.GubbinsCommon.fasttree_gubbins_command('AAA', 'BBB', 5, 'CCC', 'DDD', 3, 'EEE', 10, 100) \
           == 'DDD -r -v BBB.vcf -a 10 -b 100 -f EEE -t AAA.iteration_5 -m 3 BBB.snp_sites.aln'


if __name__ == "__main__":
  unittest.main()
