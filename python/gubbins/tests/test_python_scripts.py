#! /usr/bin/env python3
# encoding: utf-8

"""
Testing the python scripts work
"""

import unittest
import os
import subprocess
import hashlib
from gubbins import common

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
working_dir = os.path.join(modules_dir, 'tests')

class TestPythonScripts(unittest.TestCase):

    ## Test the alignment_checker script 

    def test_alignment_checker(self):
        small_aln = os.path.join(data_dir, "valid_alignment.aln")
        output_file = os.path.join(working_dir, "valid_alignment_test")
        output_csv = os.path.join(working_dir, "valid_alignment_test.csv")
        test_csv = os.path.join(data_dir, "test_valid_output.csv")
        aln_cmd = "alignment_checker.py --aln " + small_aln + " --out " + output_file 
        subprocess.check_call(aln_cmd, shell=True)
        assert self.md5_check(output_csv, test_csv)
        os.remove(output_csv)

    def test_clade_extraction(self):
        multiple_aln = os.path.join(data_dir, "multiple_recombinations.aln")
        clade_list = os.path.join(data_dir, "clade_to_extract.txt")
        multiple_gff = os.path.join(data_dir, "multiple_recombinations_gubbins.recombination_predictions.gff")
        multiple_tree = os.path.join(data_dir, "multiple_recombinations_gubbins.node_labelled.final_tree.tre")
        out_aln = os.path.join(data_dir, "multiple_recombinations_extract.aln")
        out_gff = os.path.join(data_dir, "multiple_recombinations_extract.gff")
        out_tree = os.path.join(data_dir, "multiple_recombinations_extract.tree")
        base_path = os.path.join(data_dir, "multiple_recombinations_extract")
        ## Get the test files 
        test_aln = os.path.join(data_dir, "multiple_recombinations_clade_extract.aln")
        test_gff = os.path.join(data_dir, "multiple_recombinations_clade_extract.gff")
        test_tree = os.path.join(data_dir, "multiple_recombinations_clade_extract.tree")
        extract_clade_cmd = "extract_gubbins_clade.py --list " + clade_list + " --aln " + multiple_aln +\
            " --gff " + multiple_gff + " --tree " + multiple_tree + " --out " + base_path +  " --out-fmt fasta"
        subprocess.check_call(extract_clade_cmd, shell=True)
        assert self.md5_check(out_aln, test_aln)
        assert self.md5_check(out_gff, test_gff)
        assert self.md5_check(out_tree, test_tree)
        os.remove(out_aln)
        os.remove(out_gff)
        os.remove(out_tree)


    def test_masking_aln(self):
        multiple_aln = os.path.join(data_dir, "multiple_recombinations.aln")
        multiple_gff = os.path.join(data_dir, "multiple_recombinations_gubbins.recombination_predictions.gff")
        out_aln = os.path.join(data_dir, "multiple_recombinations_mask.aln")
        ## Get the test file 
        test_aln = os.path.join(data_dir, "masking_multiple.aln")
        
        extract_clade_cmd = "mask_gubbins_aln.py --aln " + multiple_aln +\
            " --gff " + multiple_gff + " --out " + out_aln +  " --out-fmt fasta"
        subprocess.check_call(extract_clade_cmd, shell=True)
        assert self.md5_check(out_aln, test_aln)
        os.remove(out_aln)
        
    # def test_clade_stats(self):
    #     multiple_aln = os.path.join(data_dir, "multiple_recombinations.aln")
    #     clade_list = os.path.join(data_dir, "clade_to_extract.txt")
    #     multiple_gff = os.path.join(data_dir, "multiple_recombinations_gubbins.recombination_predictions.gff")
    #     multiple_tree = os.path.join(data_dir, "multiple_recombinations_gubbins.node_labelled.final_tree.tre")
    #     out_aln = os.path.join(data_dir, "multiple_recombinations_extract.aln")
    #     out_gff = os.path.join(data_dir, "multiple_recombinations_extract.gff")
    #     out_tree = os.path.join(data_dir, "multiple_recombinations_extract.tree")
    #     base_path = os.path.join(data_dir, "multiple_recombinations_extract")
    #     extract_clade_cmd = "extract_gubbins_clade.py --list " + clade_list + " --aln " + multiple_aln +\
    #         " --gff " + multiple_gff + " --tree " + multiple_tree + " --out " + base_path +  " --out-fmt fasta"
    #     subprocess.check_call(extract_clade_cmd, shell=True)
    #     assert self.md5_check(out_aln, "f9a1645bc6f668eecd4fe828df5d5d0f")
    #     assert self.md5_check(out_gff, "adc2afaebd2adc7fe6690e8a63b73c8b")
    #     assert self.md5_check(out_tree, "cafd551664494ba6caab287c7360aa67")
    #     os.remove(out_aln)
    #     os.remove(out_gff)
    #     os.remove(out_tree)


    @staticmethod
    def md5_check(file_path, correct_output):
        file_sum = hashlib.md5(open(file_path, "rb").read()).hexdigest()
        correct_sum = hashlib.md5(open(correct_output, "rb").read()).hexdigest()
        return file_sum == correct_sum

if __name__ == "__main__":
    unittest.main(buffer=True)

