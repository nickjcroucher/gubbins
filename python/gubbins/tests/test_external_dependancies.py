#! /usr/bin/env python
# encoding: utf-8

"""
Integration testing of external dependancies. Likely to be the most brittle tests, but the most important.
"""

import unittest
import re
import shutil
import sys
import os
import glob
import argparse
import pkg_resources
from gubbins import common

class TestExternalDependancies(unittest.TestCase):

  def test_which(self):
    # the location of ls varies depending on OS so just check end
    assert re.match('.*/ls$', common.GubbinsCommon.which('ls')) != None
    # Strip parameters
    assert re.match('.*/ls$', common.GubbinsCommon.which('ls -alrt')) != None
    assert common.GubbinsCommon.which('non_existant_program') == None
    
    
  def test_change_window_size(self):
    parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
    parser.add_argument('alignment_filename',       help='Multifasta alignment file')
    parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
    parser.add_argument('--starting_tree',    '-s', help='Starting tree')
    parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
    parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging',default = 0)
    parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
    parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
    parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
    parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
    parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 15)
    parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')
    parser.add_argument('--threads',          '-c', help='Number of threads to run with RAXML, but only if a PTHREADS version is available', type=int,  default = 1)
    parser.add_argument('--converge_method',  '-z', help='Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]',  default = 'weighted_robinson_foulds')
    parser.add_argument('--version',                action='version', version=str(pkg_resources.get_distribution("gubbins").version))
    parser.add_argument('--min_window_size',  '-a', help='Minimum window size, default 100', type=int,  default = 80)
    parser.add_argument('--max_window_size',  '-b', help='Maximum window size, default 10000', type=int,  default = 1000)
    gubbins_runner  = common.GubbinsCommon(parser.parse_args(["--prefix", "different_prefix",'gubbins/tests/data/multiple_recombinations.aln']))
    gubbins_runner.parse_and_run()
    
    assert os.path.exists('different_prefix.summary_of_snp_distribution.vcf')
    assert os.path.exists('different_prefix.recombination_predictions.embl')
    assert os.path.exists('different_prefix.per_branch_statistics.csv')
    assert os.path.exists('different_prefix.filtered_polymorphic_sites.fasta')
    assert os.path.exists('different_prefix.filtered_polymorphic_sites.phylip')
    assert os.path.exists('different_prefix.recombination_predictions.gff')
    assert os.path.exists('different_prefix.branch_base_reconstruction.embl')
    assert os.path.exists('different_prefix.final_tree.tre')

    os.remove('different_prefix.summary_of_snp_distribution.vcf')
    os.remove('different_prefix.node_labelled.tre')
    os.remove('different_prefix.recombination_predictions.embl')
    os.remove('different_prefix.per_branch_statistics.csv')
    os.remove('different_prefix.filtered_polymorphic_sites.fasta')
    os.remove('different_prefix.filtered_polymorphic_sites.phylip')
    os.remove('different_prefix.recombination_predictions.gff')
    os.remove('different_prefix.branch_base_reconstruction.embl')
    os.remove('different_prefix.final_tree.tre')
  
  def test_rename_final_output(self):
    parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
    parser.add_argument('alignment_filename',       help='Multifasta alignment file')
    parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
    parser.add_argument('--starting_tree',    '-s', help='Starting tree')
    parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
    parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging',default = 0)
    parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
    parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
    parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
    parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
    parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 15)
    parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')
    parser.add_argument('--threads',          '-c', help='Number of threads to run with RAXML, but only if a PTHREADS version is available', type=int,  default = 1)
    parser.add_argument('--converge_method',  '-z', help='Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]',  default = 'weighted_robinson_foulds')
    parser.add_argument('--version',                action='version', version=str(pkg_resources.get_distribution("gubbins").version))
    parser.add_argument('--min_window_size',  '-a', help='Minimum window size, default 100', type=int,  default = 100)
    parser.add_argument('--max_window_size',  '-b', help='Maximum window size, default 10000', type=int,  default = 10000)
    
    gubbins_runner  = common.GubbinsCommon(parser.parse_args(["--prefix", "different_prefix",'gubbins/tests/data/multiple_recombinations.aln']))
    gubbins_runner.parse_and_run()
    
    assert os.path.exists('different_prefix.summary_of_snp_distribution.vcf')
    assert os.path.exists('different_prefix.recombination_predictions.embl')
    assert os.path.exists('different_prefix.per_branch_statistics.csv')
    assert os.path.exists('different_prefix.filtered_polymorphic_sites.fasta')
    assert os.path.exists('different_prefix.filtered_polymorphic_sites.phylip')
    assert os.path.exists('different_prefix.recombination_predictions.gff')
    assert os.path.exists('different_prefix.branch_base_reconstruction.embl')
    assert os.path.exists('different_prefix.final_tree.tre')

    os.remove('different_prefix.summary_of_snp_distribution.vcf')
    os.remove('different_prefix.node_labelled.tre')
    os.remove('different_prefix.recombination_predictions.embl')
    os.remove('different_prefix.per_branch_statistics.csv')
    os.remove('different_prefix.filtered_polymorphic_sites.fasta')
    os.remove('different_prefix.filtered_polymorphic_sites.phylip')
    os.remove('different_prefix.recombination_predictions.gff')
    os.remove('different_prefix.branch_base_reconstruction.embl')
    os.remove('different_prefix.final_tree.tre')
    
  def test_recombination_convergence(self):
    parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
    parser.add_argument('alignment_filename',       help='Multifasta alignment file')
    parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
    parser.add_argument('--starting_tree',    '-s', help='Starting tree')
    parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
    parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging',default = 0)
    parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
    parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
    parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
    parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
    parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 15)
    parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')
    parser.add_argument('--threads',          '-c', help='Number of threads to run with RAXML, but only if a PTHREADS version is available', type=int,  default = 1)
    parser.add_argument('--converge_method',  '-z', help='Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]',  default = 'weighted_robinson_foulds')
    parser.add_argument('--min_window_size',  '-a', help='Minimum window size, default 100', type=int,  default = 100)
    parser.add_argument('--max_window_size',  '-b', help='Maximum window size, default 10000', type=int,  default = 10000)
    
    gubbins_runner  = common.GubbinsCommon(parser.parse_args(["--converge_method", "recombination", "--iterations", '15','--no_cleanup', 'gubbins/tests/data/multiple_recombinations.aln']))
    gubbins_runner.parse_and_run()
    
    assert common.GubbinsCommon.have_recombinations_been_seen_before('multiple_recombinations.recombination_predictions.embl', ['RAxML_result.multiple_recombinations.iteration_1.tab','RAxML_result.multiple_recombinations.iteration_2.tab','RAxML_result.multiple_recombinations.iteration_3.tab','RAxML_result.multiple_recombinations.iteration_4.tab','RAxML_result.multiple_recombinations.iteration_5.tab','RAxML_result.multiple_recombinations.iteration_6.tab','RAxML_result.multiple_recombinations.iteration_7.tab','RAxML_result.multiple_recombinations.iteration_8.tab','RAxML_result.multiple_recombinations.iteration_9.tab','RAxML_result.multiple_recombinations.iteration_10.tab','RAxML_result.multiple_recombinations.iteration_11.tab','RAxML_result.multiple_recombinations.iteration_12.tab','RAxML_result.multiple_recombinations.iteration_13.tab','RAxML_result.multiple_recombinations.iteration_14.tab','RAxML_result.multiple_recombinations.iteration_15.tab']) == 1 
    
    os.remove('log.txt')
    r = glob.glob("RAxML_*")
    for i in r:
       os.remove(i)
    r = glob.glob("multiple_recombinations.aln.*")
    for i in r:
       os.remove(i)
       
       
  def test_robinson_foulds_convergence(self):
    parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
    parser.add_argument('alignment_filename',       help='Multifasta alignment file')
    parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
    parser.add_argument('--starting_tree',    '-s', help='Starting tree')
    parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
    parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging',default = 0)
    parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
    parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
    parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
    parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
    parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 15)
    parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')
    parser.add_argument('--threads',          '-c', help='Number of threads to run with RAXML, but only if a PTHREADS version is available', type=int,  default = 1)
    parser.add_argument('--converge_method',  '-z', help='Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]',  default = 'weighted_robinson_foulds')
    parser.add_argument('--min_window_size',  '-a', help='Minimum window size, default 100', type=int,  default = 100)
    parser.add_argument('--max_window_size',  '-b', help='Maximum window size, default 10000', type=int,  default = 10000)
    
    gubbins_runner  = common.GubbinsCommon(parser.parse_args(["--converge_method", "robinson_foulds", "--iterations", '15','--no_cleanup', 'gubbins/tests/data/multiple_recombinations.aln']))
    gubbins_runner.parse_and_run()

    assert common.GubbinsCommon.has_tree_been_seen_before(['RAxML_result.multiple_recombinations.iteration_1','RAxML_result.multiple_recombinations.iteration_2','RAxML_result.multiple_recombinations.iteration_3','RAxML_result.multiple_recombinations.iteration_4','RAxML_result.multiple_recombinations.iteration_5','RAxML_result.multiple_recombinations.iteration_6','RAxML_result.multiple_recombinations.iteration_7','RAxML_result.multiple_recombinations.iteration_8','RAxML_result.multiple_recombinations.iteration_9','RAxML_result.multiple_recombinations.iteration_10','multiple_recombinations.final_tree.tre'],'robinson_foulds') == 1

    os.remove('log.txt')
    r = glob.glob("RAxML_*")
    for i in r:
       os.remove(i)
    r = glob.glob("multiple_recombinations.aln.*")
    for i in r:
       os.remove(i)

  def test_rename_final_output_fasttree(self):
    parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
    parser.add_argument('alignment_filename',       help='Multifasta alignment file')
    parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
    parser.add_argument('--starting_tree',    '-s', help='Starting tree')
    parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
    parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging',default = 0)
    parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
    parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
    parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
    parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
    parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 15)
    parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')
    parser.add_argument('--threads',          '-c', help='Number of threads to run with RAXML, but only if a PTHREADS version is available', type=int,  default = 1)
    parser.add_argument('--converge_method',  '-z', help='Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]',  default = 'weighted_robinson_foulds')
    parser.add_argument('--min_window_size',  '-a', help='Minimum window size, default 100', type=int,  default = 100)
    parser.add_argument('--max_window_size',  '-b', help='Maximum window size, default 10000', type=int,  default = 10000)
    
    gubbins_runner  = common.GubbinsCommon(parser.parse_args(["--prefix", "ft_prefix","--tree_builder", "fasttree",'gubbins/tests/data/multiple_recombinations.aln']))
    gubbins_runner.parse_and_run()

    assert os.path.exists('ft_prefix.summary_of_snp_distribution.vcf')
    assert os.path.exists('ft_prefix.recombination_predictions.embl')
    assert os.path.exists('ft_prefix.per_branch_statistics.csv')
    assert os.path.exists('ft_prefix.filtered_polymorphic_sites.fasta')
    assert os.path.exists('ft_prefix.filtered_polymorphic_sites.phylip')
    assert os.path.exists('ft_prefix.recombination_predictions.gff')
    assert os.path.exists('ft_prefix.branch_base_reconstruction.embl')
    assert os.path.exists('ft_prefix.final_tree.tre')

    os.remove('ft_prefix.summary_of_snp_distribution.vcf')
    os.remove('ft_prefix.recombination_predictions.embl')
    os.remove('ft_prefix.per_branch_statistics.csv')
    os.remove('ft_prefix.filtered_polymorphic_sites.fasta')
    os.remove('ft_prefix.filtered_polymorphic_sites.phylip')
    os.remove('ft_prefix.recombination_predictions.gff')
    os.remove('ft_prefix.branch_base_reconstruction.embl')
    os.remove('ft_prefix.final_tree.tre')
    os.remove('ft_prefix.node_labelled.tre')
    
  def test_rename_final_output_hybrid(self):
    parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
    parser.add_argument('alignment_filename',       help='Multifasta alignment file')
    parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
    parser.add_argument('--starting_tree',    '-s', help='Starting tree')
    parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
    parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging',default = 0)
    parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
    parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
    parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
    parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
    parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 15)
    parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')
    parser.add_argument('--threads',          '-c', help='Number of threads to run with RAXML, but only if a PTHREADS version is available', type=int,  default = 1)
    parser.add_argument('--converge_method',  '-z', help='Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]',  default = 'weighted_robinson_foulds')
    parser.add_argument('--min_window_size',  '-a', help='Minimum window size, default 100', type=int,  default = 100)
    parser.add_argument('--max_window_size',  '-b', help='Maximum window size, default 10000', type=int,  default = 10000)
    
    gubbins_runner  = common.GubbinsCommon(parser.parse_args(["--prefix", "hybrid_prefix","--tree_builder", "hybrid",'gubbins/tests/data/multiple_recombinations.aln']))
    gubbins_runner.parse_and_run()
 
    assert os.path.exists('hybrid_prefix.summary_of_snp_distribution.vcf')
    assert os.path.exists('hybrid_prefix.recombination_predictions.embl')
    assert os.path.exists('hybrid_prefix.per_branch_statistics.csv')
    assert os.path.exists('hybrid_prefix.filtered_polymorphic_sites.fasta')
    assert os.path.exists('hybrid_prefix.filtered_polymorphic_sites.phylip')
    assert os.path.exists('hybrid_prefix.recombination_predictions.gff')
    assert os.path.exists('hybrid_prefix.branch_base_reconstruction.embl')
    assert os.path.exists('hybrid_prefix.final_tree.tre')
 
    os.remove('hybrid_prefix.summary_of_snp_distribution.vcf')
    os.remove('hybrid_prefix.node_labelled.tre')
    os.remove('hybrid_prefix.recombination_predictions.embl')
    os.remove('hybrid_prefix.per_branch_statistics.csv')
    os.remove('hybrid_prefix.filtered_polymorphic_sites.fasta')
    os.remove('hybrid_prefix.filtered_polymorphic_sites.phylip')
    os.remove('hybrid_prefix.recombination_predictions.gff')
    os.remove('hybrid_prefix.branch_base_reconstruction.embl')
    os.remove('hybrid_prefix.final_tree.tre')
    
  def test_fasttree_default_output_names(self):
    parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
    parser.add_argument('alignment_filename',       help='Multifasta alignment file')
    parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
    parser.add_argument('--starting_tree',    '-s', help='Starting tree')
    parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
    parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging',default = 0)
    parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
    parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
    parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
    parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
    parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 15)
    parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')
    parser.add_argument('--threads',          '-c', help='Number of threads to run with RAXML, but only if a PTHREADS version is available', type=int,  default = 1)
    parser.add_argument('--converge_method',  '-z', help='Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]',  default = 'weighted_robinson_foulds')
    parser.add_argument('--min_window_size',  '-a', help='Minimum window size, default 100', type=int,  default = 100)
    parser.add_argument('--max_window_size',  '-b', help='Maximum window size, default 10000', type=int,  default = 10000)
    
    gubbins_runner  = common.GubbinsCommon(parser.parse_args(["--tree_builder", "fasttree",'gubbins/tests/data/multiple_recombinations.aln']))
    gubbins_runner.parse_and_run()

    assert os.path.exists('multiple_recombinations.summary_of_snp_distribution.vcf')
    assert os.path.exists('multiple_recombinations.recombination_predictions.embl')
    assert os.path.exists('multiple_recombinations.per_branch_statistics.csv')
    assert os.path.exists('multiple_recombinations.filtered_polymorphic_sites.fasta')
    assert os.path.exists('multiple_recombinations.filtered_polymorphic_sites.phylip')
    assert os.path.exists('multiple_recombinations.recombination_predictions.gff')
    assert os.path.exists('multiple_recombinations.branch_base_reconstruction.embl')
    assert os.path.exists('multiple_recombinations.final_tree.tre')
    assert os.path.exists('multiple_recombinations.node_labelled.tre')

    os.remove('multiple_recombinations.node_labelled.tre')
    os.remove('multiple_recombinations.summary_of_snp_distribution.vcf')
    os.remove('multiple_recombinations.recombination_predictions.embl')
    os.remove('multiple_recombinations.per_branch_statistics.csv')
    os.remove('multiple_recombinations.filtered_polymorphic_sites.fasta')
    os.remove('multiple_recombinations.filtered_polymorphic_sites.phylip')
    os.remove('multiple_recombinations.recombination_predictions.gff')
    os.remove('multiple_recombinations.branch_base_reconstruction.embl')
    os.remove('multiple_recombinations.final_tree.tre')

  def test_hybrid_default_output_names(self):
    parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
    parser.add_argument('alignment_filename',       help='Multifasta alignment file')
    parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
    parser.add_argument('--starting_tree',    '-s', help='Starting tree')
    parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
    parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging',default = 0)
    parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
    parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
    parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
    parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
    parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 15)
    parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')
    parser.add_argument('--threads',          '-c', help='Number of threads to run with RAXML, but only if a PTHREADS version is available', type=int,  default = 1)
    parser.add_argument('--converge_method',  '-z', help='Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]',  default = 'weighted_robinson_foulds')
    parser.add_argument('--min_window_size',  '-a', help='Minimum window size, default 100', type=int,  default = 100)
    parser.add_argument('--max_window_size',  '-b', help='Maximum window size, default 10000', type=int,  default = 10000)
    
    gubbins_runner  = common.GubbinsCommon(parser.parse_args(["--tree_builder", "hybrid",'gubbins/tests/data/multiple_recombinations.aln']))
    gubbins_runner.parse_and_run()

    assert os.path.exists('multiple_recombinations.summary_of_snp_distribution.vcf')
    assert os.path.exists('multiple_recombinations.recombination_predictions.embl')
    assert os.path.exists('multiple_recombinations.per_branch_statistics.csv')
    assert os.path.exists('multiple_recombinations.filtered_polymorphic_sites.fasta')
    assert os.path.exists('multiple_recombinations.filtered_polymorphic_sites.phylip')
    assert os.path.exists('multiple_recombinations.recombination_predictions.gff')
    assert os.path.exists('multiple_recombinations.branch_base_reconstruction.embl')
    assert os.path.exists('multiple_recombinations.final_tree.tre')

    os.remove('multiple_recombinations.summary_of_snp_distribution.vcf')
    os.remove('multiple_recombinations.recombination_predictions.embl')
    os.remove('multiple_recombinations.per_branch_statistics.csv')
    os.remove('multiple_recombinations.filtered_polymorphic_sites.fasta')
    os.remove('multiple_recombinations.filtered_polymorphic_sites.phylip')
    os.remove('multiple_recombinations.recombination_predictions.gff')
    os.remove('multiple_recombinations.branch_base_reconstruction.embl')
    os.remove('multiple_recombinations.final_tree.tre')
    os.remove('multiple_recombinations.node_labelled.tre')

  def test_parse_and_run(self):

    parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
    parser.add_argument('alignment_filename',       help='Multifasta alignment file')
    parser.add_argument('--outgroup',         '-o', help='Outgroup name for rerooting')
    parser.add_argument('--starting_tree',    '-s', help='Starting tree')
    parser.add_argument('--use_time_stamp',   '-u', action='count', help='Use a time stamp in file names')
    parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging',default = 0)
    parser.add_argument('--no_cleanup',       '-n', action='count', help='Dont cleanup intermediate files')
    parser.add_argument('--tree_builder',     '-t', help='Application to use for tree building (raxml, fasttree, hybrid), default RAxML', default = "raxml")
    parser.add_argument('--iterations',       '-i', help='Maximum No. of iterations, default is 5', type=int,  default = 5)
    parser.add_argument('--min_snps',         '-m', help='Min SNPs to identify a recombination block, default is 3', type=int,  default = 3)
    parser.add_argument('--filter_percentage','-f', help='Filter out taxa with more than this percentage of gaps, default is 25', type=int,  default = 15)
    parser.add_argument('--prefix',           '-p', help='Add a prefix to the final output filenames')
    parser.add_argument('--threads',          '-c', help='Number of threads to run with RAXML, but only if a PTHREADS version is available', type=int,  default = 1)
    parser.add_argument('--converge_method',  '-z', help='Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]',  default = 'weighted_robinson_foulds')
    parser.add_argument('--min_window_size',  '-a', help='Minimum window size, default 100', type=int,  default = 100)
    parser.add_argument('--max_window_size',  '-b', help='Maximum window size, default 10000', type=int,  default = 10000)
    
    #  multiple recombinations
    gubbins_runner  = common.GubbinsCommon(parser.parse_args(['gubbins/tests/data/multiple_recombinations.aln']))
    gubbins_runner.parse_and_run()

    actual_file_content2    = open('multiple_recombinations.summary_of_snp_distribution.vcf',    encoding='utf-8').readlines()
    actual_file_content3    = open('multiple_recombinations.recombination_predictions.embl',    encoding='utf-8').readlines()
    actual_file_content4    = open('multiple_recombinations.per_branch_statistics.csv',    encoding='utf-8').readlines()
    actual_file_content5    = open('multiple_recombinations.filtered_polymorphic_sites.fasta',    encoding='utf-8').readlines()
    actual_file_content6    = open('multiple_recombinations.filtered_polymorphic_sites.phylip',    encoding='utf-8').readlines()
    actual_file_content8    = open('multiple_recombinations.recombination_predictions.gff',    encoding='utf-8').readlines()
    actual_file_content9    = open('multiple_recombinations.branch_base_reconstruction.embl',    encoding='utf-8').readlines()
    actual_file_content10   = open('multiple_recombinations.final_tree.tre',    encoding='utf-8').readlines()
    
    expected_file_content2  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.vcf',    encoding='utf-8').readlines()
    expected_file_content3  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.tab',    encoding='utf-8').readlines()
    expected_file_content4  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.stats',    encoding='utf-8').readlines()
    expected_file_content5  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.snp_sites.aln',    encoding='utf-8').readlines()
    expected_file_content6  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.phylip',    encoding='utf-8').readlines()
    expected_file_content8  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.gff',    encoding='utf-8').readlines()
    expected_file_content9  = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5.branch_snps.tab',    encoding='utf-8').readlines()
    expected_file_content10 = open('gubbins/tests/data/expected_RAxML_result.multiple_recombinations.iteration_5',    encoding='utf-8').readlines()
    
    # Slightly different values of internal nodes if run on fastml on linux and osx
    #assert actual_file_content2 == expected_file_content2
    #assert actual_file_content3 == expected_file_content3
    #assert actual_file_content4 == expected_file_content4
    #assert actual_file_content5 == expected_file_content5
    #assert actual_file_content6 == expected_file_content6
    #assert actual_file_content8 == expected_file_content8
    #assert actual_file_content9 == expected_file_content9
    #assert actual_file_content10 == expected_file_content10
    
    assert not os.path.exists('multiple_recombinations.aln.start')
    assert not os.path.exists('RAxML_result.multiple_recombinations.iteration_5.ancestor.tre')
    assert not os.path.exists('RAxML_result.multiple_recombinations.iteration_5.seq.joint.txt')
    assert not os.path.exists('RAxML_result.multiple_recombinations.iteration_5.prob.joint.txt')
    assert not os.path.exists('log.txt')
    assert not os.path.exists('latest_tree.multiple_recombinations.tre')
  
    os.remove('multiple_recombinations.summary_of_snp_distribution.vcf')
    os.remove('multiple_recombinations.recombination_predictions.embl')
    os.remove('multiple_recombinations.per_branch_statistics.csv')
    os.remove('multiple_recombinations.filtered_polymorphic_sites.fasta')
    os.remove('multiple_recombinations.filtered_polymorphic_sites.phylip')
    os.remove('multiple_recombinations.recombination_predictions.gff')
    os.remove('multiple_recombinations.branch_base_reconstruction.embl')
    os.remove('multiple_recombinations.final_tree.tre')
    os.remove('multiple_recombinations.node_labelled.tre')
  
  
  
  def test_pairwise_comparison(self):
    shutil.copyfile('gubbins/tests/data/input_pairwise.aln.vcf','gubbins/tests/data/pairwise.aln.vcf' )
    shutil.copyfile('gubbins/tests/data/input_pairwise.aln.phylip','gubbins/tests/data/pairwise.aln.phylip' )
    common.GubbinsCommon.pairwise_comparison('gubbins/tests/data/pairwise.aln','gubbins/tests/data/pairwise.aln','../src/gubbins','gubbins/tests/data/pairwise.aln','fastml  -mg -qf -b ','pairwise')
    # Check the tree file exists
    assert os.path.exists('pairwise.final_tree.tre')
    
    # Check the VCF file is as expected
    actual_file_content   = open('pairwise.summary_of_snp_distribution.vcf',    encoding='utf-8').readlines()
    expected_file_content = open('gubbins/tests/data/pairwise.aln.tre.vcf_expected',  encoding='utf-8').readlines()
    assert actual_file_content == expected_file_content
    
    # Check the reconstruction of internal nodes
    actual_file_content   = open('pairwise.filtered_polymorphic_sites.fasta',    encoding='utf-8').readlines()
    expected_file_content = open('gubbins/tests/data/pairwise.aln.snp_sites.aln_expected',  encoding='utf-8').readlines()
    assert actual_file_content == expected_file_content
    
    
    os.remove('pairwise.summary_of_snp_distribution.vcf')
    os.remove('pairwise.filtered_polymorphic_sites.fasta')
    os.remove('pairwise.filtered_polymorphic_sites.phylip')
    os.remove('pairwise.final_tree.tre')
    
  
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
    assert re.search('raxmlHPC -f d -p 1 -m GTRGAMMA',common.GubbinsCommon.use_bundled_exec('raxmlHPC -f d -p 1 -m GTRGAMMA', 'raxmlHPC')) != None
    assert re.search('fastml -mg -qf -b ',common.GubbinsCommon.use_bundled_exec('fastml -mg -qf -b ', 'fastml')) != None
    assert re.search('../src/gubbins',common.GubbinsCommon.use_bundled_exec('gubbins', '../src/gubbins')) != None

if __name__ == "__main__":
  unittest.main()

