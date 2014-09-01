#! /usr/bin/env python
# encoding: utf-8

"""
Tests the construction of file names and commands with no dependancies on external applications.
"""

import unittest
import re
from gubbins import common

class TestStringConstruction(unittest.TestCase):
  def test_raxml_base_name(self):
    assert common.GubbinsCommon.raxml_base_name('AAA',1234) == 'AAA.1234iteration_'
    assert common.GubbinsCommon.raxml_base_name('AAA','') == 'AAA.iteration_'

  def test_raxml_current_tree_name(self):
    assert common.GubbinsCommon.raxml_current_tree_name('AAA',1234, 5) == 'RAxML_result.AAA.1234iteration_5'
    assert common.GubbinsCommon.raxml_current_tree_name('AAA','', 5) == 'RAxML_result.AAA.iteration_5'

  def test_raxml_previous_tree_name(self):
    assert common.GubbinsCommon.raxml_previous_tree_name('AAA','BBB', 1234,5) == 'RAxML_result.AAA.1234iteration_4'
    assert common.GubbinsCommon.raxml_previous_tree_name('AAA','BBB', '',5) == 'RAxML_result.AAA.iteration_4'

  def test_raxml_previous_tree(self):
    assert common.GubbinsCommon.raxml_previous_tree('AAA','BBB', 1234,5,'ABC') == '-t ABC'
    assert common.GubbinsCommon.raxml_previous_tree('AAA','BBB', 1234,1,'ABC') == ''

  def test_raxml_tree_building_command(self):
    assert common.GubbinsCommon.raxml_tree_building_command(5,'AAA','BBB',1234, 'raxml_exec','CCC', 0) == 'raxml_exec -s CCC.phylip -n AAA.1234iteration_5 -t CCC> /dev/null 2>&1'
    assert common.GubbinsCommon.raxml_tree_building_command(1,'AAA','BBB',1234, 'raxml_exec','CCC', 0) == 'raxml_exec -s CCC.phylip -n AAA.1234iteration_1 > /dev/null 2>&1'
    assert common.GubbinsCommon.raxml_tree_building_command(2,'AAA','BBB','', 'raxml_exec','CCC', 0) == 'raxml_exec -s CCC.phylip -n AAA.iteration_2 -t CCC> /dev/null 2>&1'
    assert common.GubbinsCommon.raxml_tree_building_command(2,'AAA','BBB','', 'raxml_exec','CCC', 1) == 'raxml_exec -s CCC.phylip -n AAA.iteration_2 -t CCC'

  def test_raxml_gubbins_command(self):
    assert common.GubbinsCommon.raxml_gubbins_command('AAA','BBB','', 5,'CCC','DDD',3, 'EEE') == 'DDD -r -v BBB.vcf -f EEE -t RAxML_result.AAA.iteration_5 -m 3 BBB.snp_sites.aln'
    assert common.GubbinsCommon.raxml_gubbins_command('base_filename_without_ext','starting_base_filename',1234, 5,'alignment_filename','gubbins_exec',3, 'original_aln') == 'gubbins_exec -r -v starting_base_filename.vcf -f original_aln -t RAxML_result.base_filename_without_ext.1234iteration_5 -m 3 starting_base_filename.snp_sites.aln'

  def test_raxml_regex_for_file_deletions(self):
    assert common.GubbinsCommon.raxml_regex_for_file_deletions('AAA',1234,'BBB', 5) == ['^RAxML_(bestTree|info|log|parsimonyTree).AAA.1234iteration_', '^RAxML_result.AAA.1234iteration_1\\.', '^RAxML_result.AAA.1234iteration_1$', '^RAxML_result.AAA.1234iteration_2\\.', '^RAxML_result.AAA.1234iteration_2$', '^RAxML_result.AAA.1234iteration_3\\.', '^RAxML_result.AAA.1234iteration_3$', '^RAxML_result.AAA.1234iteration_4\\.', '^RAxML_result.AAA.1234iteration_4$', '^RAxML_result.AAA.1234iteration_5.ancestor.tre', '^RAxML_result.AAA.1234iteration_5.seq.joint.txt', '^RAxML_result.AAA.1234iteration_5.prob.joint.txt', '^RAxML_result.AAA.1234iteration_5.output_tree']

  def test_raxml_regex_for_file_deletions_10(self):
    assert common.GubbinsCommon.raxml_regex_for_file_deletions('AAA',1234,'BBB', 10) == ['^RAxML_(bestTree|info|log|parsimonyTree).AAA.1234iteration_', '^RAxML_result.AAA.1234iteration_1\\.', '^RAxML_result.AAA.1234iteration_1$', '^RAxML_result.AAA.1234iteration_2\\.', '^RAxML_result.AAA.1234iteration_2$', '^RAxML_result.AAA.1234iteration_3\\.', '^RAxML_result.AAA.1234iteration_3$', '^RAxML_result.AAA.1234iteration_4\\.', '^RAxML_result.AAA.1234iteration_4$', '^RAxML_result.AAA.1234iteration_5\\.', '^RAxML_result.AAA.1234iteration_5$', '^RAxML_result.AAA.1234iteration_6\\.', '^RAxML_result.AAA.1234iteration_6$', '^RAxML_result.AAA.1234iteration_7\\.', '^RAxML_result.AAA.1234iteration_7$', '^RAxML_result.AAA.1234iteration_8\\.', '^RAxML_result.AAA.1234iteration_8$', '^RAxML_result.AAA.1234iteration_9\\.', '^RAxML_result.AAA.1234iteration_9$', '^RAxML_result.AAA.1234iteration_10.ancestor.tre', '^RAxML_result.AAA.1234iteration_10.seq.joint.txt', '^RAxML_result.AAA.1234iteration_10.prob.joint.txt', '^RAxML_result.AAA.1234iteration_10.output_tree']

  def test_fasttree_regex_for_file_deletions(self):
    assert common.GubbinsCommon.fasttree_regex_for_file_deletions('AAA', 5) == ['^AAA.iteration_1[$|\\.]', '^AAA.iteration_2[$|\\.]', '^AAA.iteration_3[$|\\.]', '^AAA.iteration_4[$|\\.]']

  def test_starting_files_regex(self):
    assert common.GubbinsCommon.starting_files_regex('AAA') == 'AAA.(gaps|vcf|snp_sites|phylip|aln.start)'
    assert common.GubbinsCommon.starting_files_regex('^')   == '^.(gaps|vcf|snp_sites|phylip|aln.start)'
    assert common.GubbinsCommon.starting_files_regex('')    == '.(gaps|vcf|snp_sites|phylip|aln.start)'

  def test_translation_of_fasttree_filenames_to_final_filenames(self):
    assert common.GubbinsCommon.translation_of_fasttree_filenames_to_final_filenames('AAA', 5, 'test') == {
     'AAA.iteration_5.vcf':             'test.summary_of_snp_distribution.vcf',
     'AAA.iteration_5.tab':             'test.recombination_predictions.embl',
     'AAA.iteration_5.branch_snps.tab': 'test.branch_base_reconstruction.embl',
     'AAA.iteration_5.gff':             'test.recombination_predictions.gff',
     'AAA.iteration_5.snp_sites.aln':   'test.filtered_polymorphic_sites.fasta',
     'AAA.iteration_5.phylip':          'test.filtered_polymorphic_sites.phylip',
     'AAA.iteration_5.stats':           'test.per_branch_statistics.csv',
     'AAA.iteration_5':                 'test.final_tree.tre'}

  def test_translation_of_raxml_filenames_to_final_filenames(self):
    assert common.GubbinsCommon.translation_of_raxml_filenames_to_final_filenames('AAA',1234, 10, 'test') == {
    'RAxML_result.AAA.1234iteration_10.snp_sites.aln':   'test.filtered_polymorphic_sites.fasta',
    'RAxML_result.AAA.1234iteration_10.branch_snps.tab': 'test.branch_base_reconstruction.embl',
    'RAxML_result.AAA.1234iteration_10.tab':             'test.recombination_predictions.embl',
    'RAxML_result.AAA.1234iteration_10.gff':             'test.recombination_predictions.gff',
    'RAxML_result.AAA.1234iteration_10':                 'test.final_tree.tre',
    'RAxML_result.AAA.1234iteration_10.phylip':          'test.filtered_polymorphic_sites.phylip',
    'RAxML_result.AAA.1234iteration_10.stats':           'test.per_branch_statistics.csv',
    'RAxML_result.AAA.1234iteration_10.vcf':             'test.summary_of_snp_distribution.vcf'}

  def test_fasttree_current_tree_name(self):
    assert common.GubbinsCommon.fasttree_current_tree_name('AAA', 1)  == 'AAA.iteration_1'
    assert common.GubbinsCommon.fasttree_current_tree_name('AAA', 5)  == 'AAA.iteration_5'
    assert common.GubbinsCommon.fasttree_current_tree_name('AAA', '') == 'AAA.iteration_'

  def test_fasttree_previous_tree_name(self):
    assert common.GubbinsCommon.fasttree_previous_tree_name('AAA', 5) == 'AAA.iteration_4'

  def test_fasttree_tree_building_command(self):
    assert common.GubbinsCommon.fasttree_tree_building_command(5, 'AAA', 'BBB','CCC', 'EEE','FFF', 'GGG','HHH' ) == 'FFF  -intree AAA GGG CCC.snp_sites.aln   > HHH.iteration_5'

  def test_fasttree_gubbins_command(self):
    assert common.GubbinsCommon.fasttree_gubbins_command('AAA','BBB', 5,'CCC','DDD',3,'EEE') == 'DDD -r -v BBB.vcf -f EEE -t AAA.iteration_5 -m 3 BBB.snp_sites.aln'

  def test_fasttree_fastml_command(self):
    assert common.GubbinsCommon.fasttree_fastml_command('AAA', 'BBB', 'CCC',2) == 'AAA -s BBB -t CCC.iteration_2 -x CCC.iteration_2.output_tree -y CCC.iteration_2.ancestor.tre -j CCC.iteration_2.seq.joint.txt -k CCC.iteration_2.seq.marginal.txt -d CCC.iteration_2.prob.joint.txt -e CCC.iteration_2.prob.marginal.txt'

  def test_raxml_fastml_command(self):
    assert common.GubbinsCommon.raxml_fastml_command('AAA', 'BBB', 'CCC',1234, 5) == 'AAA -s BBB -t RAxML_result.CCC.1234iteration_5 -x RAxML_result.CCC.1234iteration_5.output_tree -y RAxML_result.CCC.1234iteration_5.ancestor.tre -j RAxML_result.CCC.1234iteration_5.seq.joint.txt -k RAxML_result.CCC.1234iteration_5.seq.marginal.txt -d RAxML_result.CCC.1234iteration_5.prob.joint.txt -e RAxML_result.CCC.1234iteration_5.prob.marginal.txt'

  def test_generate_fastml_command(self):
    assert common.GubbinsCommon.generate_fastml_command('AAA', 'BBB', 'CCC') == 'AAA -s BBB -t CCC -x CCC.output_tree -y CCC.ancestor.tre -j CCC.seq.joint.txt -k CCC.seq.marginal.txt -d CCC.prob.joint.txt -e CCC.prob.marginal.txt'

if __name__ == "__main__":
  unittest.main()