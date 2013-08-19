#! /usr/bin/env python
# encoding: utf-8

"""
Integration testing of external dependancies. Likely to be the most brittle tests, but the most important.
"""

import unittest
import re
import os
import argparse
from gubbins import common

class TestExternalDependancies(unittest.TestCase):

  def test_which(self):
    # the location of ls varies depending on OS so just check end
    assert re.match('.*/ls$', common.GubbinsCommon.which('ls')) != None
    # Strip parameters
    assert re.match('.*/ls$', common.GubbinsCommon.which('ls -alrt')) != None
    assert common.GubbinsCommon.which('non_existant_program') == None

  def test_parse_and_run(self):

    parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
    parser.add_argument('alignment_filename',       help='Multifasta alignment file')
    parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
    parser.add_argument('--starting_tree',    '-s', help='Starting tree')
    parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
    parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging',default = 1)
    parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
    parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
    parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
    parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
    parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 25)
    
    #  multiple recombinations
    gubbins_runner  = common.GubbinsCommon(parser.parse_args(['gubbins/tests/data/multiple_recombinations.aln']))
    gubbins_runner.parse_and_run()

    actual_file_content1    = open('multiple_recombinations.aln.start',   'U').readlines()
    actual_file_content2    = open('RAxML_result.multiple_recombinations.iteration_5.vcf',   'U').readlines()
    actual_file_content3    = open('RAxML_result.multiple_recombinations.iteration_5.tab',   'U').readlines()
    actual_file_content4    = open('RAxML_result.multiple_recombinations.iteration_5.stats',   'U').readlines()
    actual_file_content5    = open('RAxML_result.multiple_recombinations.iteration_5.snp_sites.aln',   'U').readlines()
    actual_file_content6    = open('RAxML_result.multiple_recombinations.iteration_5.phylip',   'U').readlines()
    actual_file_content7    = open('RAxML_result.multiple_recombinations.iteration_5.output_tree',   'U').readlines()
    actual_file_content8   = open('RAxML_result.multiple_recombinations.iteration_5.gff',   'U').readlines()
    actual_file_content9   = open('RAxML_result.multiple_recombinations.iteration_5.branch_snps.tab',   'U').readlines()
    actual_file_content10   = open('RAxML_result.multiple_recombinations.iteration_5',   'U').readlines()
    
    expected_file_content1  = open('gubbins/tests/data/expected_multiple_recombinations.aln.start',   'U').readlines()
    expected_file_content2  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.vcf',   'U').readlines()
    expected_file_content3  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.tab',   'U').readlines()
    expected_file_content4  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.stats',   'U').readlines()
    expected_file_content5  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.snp_sites.aln',   'U').readlines()
    expected_file_content6  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.phylip',   'U').readlines()
    expected_file_content7  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.output_tree',   'U').readlines()
    expected_file_content8 = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.gff',   'U').readlines()
    expected_file_content9 = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.branch_snps.tab',   'U').readlines()
    expected_file_content10 = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5',   'U').readlines()
    
    assert actual_file_content1 == expected_file_content1
    assert actual_file_content2 == expected_file_content2
    assert actual_file_content3 == expected_file_content3
    assert actual_file_content4 == expected_file_content4
    assert actual_file_content5 == expected_file_content5
    assert actual_file_content6 == expected_file_content6
    assert actual_file_content7 == expected_file_content7
    assert actual_file_content8 == expected_file_content8
    assert actual_file_content9 == expected_file_content9
    assert actual_file_content10 == expected_file_content10
    
    os.remove('multiple_recombinations.aln.start')
    os.remove('log.txt')
    os.remove('latest_tree.multiple_recombinations.tre')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.vcf')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.tab')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.stats')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.snp_sites.aln')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.seq.joint.txt')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.prob.joint.txt')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.phylip')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.output_tree')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.gff')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.branch_snps.tab')
    os.remove('RAxML_result.multiple_recombinations.iteration_5.ancestor.tre')
    os.remove('RAxML_result.multiple_recombinations.iteration_5')
  
  #def test_running_fasttree(self):
  #  assert 0 == 1
  
  def test_pairwise_comparison(self):
    common.GubbinsCommon.pairwise_comparison('gubbins/tests/data/pairwise.aln','gubbins/tests/data/pairwise.aln','../src/gubbins','gubbins/tests/data/pairwise.aln','../external/fastml/programs/fastml/fastml  -mg -qf -b ')
    # Check the tree file is as expected
    actual_file_content   = open('gubbins/tests/data/pairwise.aln.tre',   'U').readlines()
    expected_file_content = open('gubbins/tests/data/pairwise_expected.tre', 'U').readlines()
    assert actual_file_content == expected_file_content
    
    # Check the VCF file is as expected
    actual_file_content   = open('gubbins/tests/data/pairwise.aln.tre.vcf',   'U').readlines()
    expected_file_content = open('gubbins/tests/data/pairwise.aln.tre.vcf_expected', 'U').readlines()
    assert actual_file_content == expected_file_content
    
    # Check the reconstruction of internal nodes
    actual_file_content   = open('gubbins/tests/data/pairwise.aln.snp_sites.aln',   'U').readlines()
    expected_file_content = open('gubbins/tests/data/pairwise.aln.snp_sites.aln_expected', 'U').readlines()
    assert actual_file_content == expected_file_content
    
    os.remove('gubbins/tests/data/pairwise.aln.snp_sites.aln')
    os.remove('gubbins/tests/data/pairwise.aln.tre')
    os.remove('gubbins/tests/data/pairwise.aln.tre.ancestor.tre')
    os.remove('gubbins/tests/data/pairwise.aln.tre.branch_snps.tab')
    os.remove('gubbins/tests/data/pairwise.aln.tre.gff')
    os.remove('gubbins/tests/data/pairwise.aln.tre.output_tree')
    os.remove('gubbins/tests/data/pairwise.aln.tre.phylip')
    os.remove('gubbins/tests/data/pairwise.aln.tre.prob.joint.txt')
    os.remove('gubbins/tests/data/pairwise.aln.tre.seq.joint.txt')
    os.remove('gubbins/tests/data/pairwise.aln.tre.snp_sites.aln')
    os.remove('gubbins/tests/data/pairwise.aln.tre.stats')
    os.remove('gubbins/tests/data/pairwise.aln.tre.tab')
    os.remove('gubbins/tests/data/pairwise.aln.tre.vcf')
    os.remove('log.txt')

  
  def test_delete_files_based_on_list_of_regexes(self):
    open('gubbins/tests/data/AAA', 'w').close()
    open('gubbins/tests/data/BBB', 'w').close() 
    open('gubbins/tests/data/BBBAAA', 'w').close() 
    open('gubbins/tests/data/AAABBB', 'w').close() 
    
    common.GubbinsCommon.delete_files_based_on_list_of_regexes('gubbins/tests/data', ['AAA'], 0)
    assert not os.path.exists('gubbins/tests/data/AAA')
    assert     os.path.exists('gubbins/tests/data/BBB')
    assert     os.path.exists('gubbins/tests/data/BBBAAA')
    assert not os.path.exists('gubbins/tests/data/AAABBB')
    os.remove('gubbins/tests/data/BBB')
    os.remove('gubbins/tests/data/BBBAAA')
  
  
  def test_use_bundled_exec(self):
    assert re.search('../external/standard-RAxML/raxmlHPC -f d -p 1 -m GTRGAMMA',common.GubbinsCommon.use_bundled_exec('raxmlHPC -f d -p 1 -m GTRGAMMA', '../external/standard-RAxML/raxmlHPC')) != None
    assert re.search('../external/fastml/programs/fastml/fastml -mg -qf -b ',common.GubbinsCommon.use_bundled_exec('fastml -mg -qf -b ', '../external/fastml/programs/fastml/fastml')) != None
    assert re.search('../src/gubbins',common.GubbinsCommon.use_bundled_exec('gubbins', '../src/gubbins')) != None

if __name__ == "__main__":
  unittest.main()

