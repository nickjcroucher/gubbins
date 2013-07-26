#! /usr/bin/env python
# encoding: utf-8

"""
Integration testing of external dependancies. Likely to be the most brittle tests, but the most important.
"""

import unittest
import re
import os
from gubbins import common

class TestExternalDependancies(unittest.TestCase):

  def test_which(self):
    # the location of ls varies depending on OS so just check end
    assert re.match('.*/ls$', common.GubbinsCommon.which('ls')) != None
    # Strip parameters
    assert re.match('.*/ls$', common.GubbinsCommon.which('ls -alrt')) != None
    assert common.GubbinsCommon.which('non_existant_program') == None

  #def test_parse_and_run(self):
  #  assert 1 == 0
  #
  
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

