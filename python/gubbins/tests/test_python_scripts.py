#! /usr/bin/env python3
# encoding: utf-8

"""
Testing the python scripts work
"""

import unittest
import os
import subprocess
import hashlib
import glob
from gubbins import common, run_gubbins

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
preprocess_dir = os.path.join(data_dir, 'preprocessfasta')
working_dir = os.path.join(modules_dir, 'tests')

class TestPythonScripts(unittest.TestCase):

    ## Test the alignment_checker script 

    def test_alignment_checker(self):
        small_aln = os.path.join(data_dir, "valid_alignment.aln")
        output_file = os.path.join(working_dir, "valid_alignment_test")
        output_csv = os.path.join(working_dir, "valid_alignment_test.csv")
        test_csv = os.path.join(data_dir, "test_valid_output.csv")
        aln_cmd = "gubbins_alignment_checker.py --aln " + small_aln + " --out " + output_file
        subprocess.check_call(aln_cmd, shell=True)
        assert self.md5_check(output_csv, test_csv)
        os.remove(output_csv)

    ## Test the clade extraction script

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
        # Script name
        extract_clade_cmd = "extract_gubbins_clade.py --list " + clade_list + " --aln " + multiple_aln +\
            " --gff " + multiple_gff + " --tree " + multiple_tree + " --out " + base_path +  " --out-fmt fasta"
        subprocess.check_call(extract_clade_cmd, shell=True)
        assert self.md5_check(out_aln, test_aln)
        assert self.md5_check(out_gff, test_gff)
        assert self.md5_check(out_tree, test_tree)
        os.remove(out_aln)
        os.remove(out_gff)
        os.remove(out_tree)

    ## Test the masking aln script 
    def test_masking_aln(self):
        multiple_aln = os.path.join(data_dir, "multiple_recombinations.aln")
        multiple_gff = os.path.join(data_dir, "multiple_recombinations_gubbins.recombination_predictions.gff")
        out_aln = os.path.join(data_dir, "multiple_recombinations_mask.aln")
        ## Get the test file 
        test_aln = os.path.join(data_dir, "masking_multiple.aln")
        # Script name
        extract_clade_cmd = "extract_gubbins_clade.py --aln " + multiple_aln +\
            " --gff " + multiple_gff + " --out " + out_aln +  " --out-fmt fasta"
        subprocess.check_call(extract_clade_cmd, shell=True)
        assert self.md5_check(out_aln, test_aln)
        os.remove(out_aln)

    ## Test the ska alignment generator 
    def test_generate_ska_alignment(self):
        exit_code = 1
        ## Change to the preprocess dir to run the tests
        os.chdir(preprocess_dir) 
        ## Run the generate_ska_alignment script
        fasta_loc = './ska_fasta_list.txt'
        ref_seq = os.path.join(preprocess_dir, 'sequence_t1.fasta')
        aln_out = os.path.join(preprocess_dir, 'ska_test_aln.aln')
        # Script name
        ska_cmd = "generate_ska_alignment.py --fasta " + fasta_loc +\
            " --reference " + ref_seq + " --out " + aln_out +\
                " --k 6"
        subprocess.check_call(ska_cmd, shell=True)
        ## Now run gubbins on the aln and check all the output is produced 
        parser = run_gubbins.parse_input_args()
        common.parse_and_run(parser.parse_args(["--prefix", "ska_test",
                                                    "--verbose", "--mar",
                                                    aln_out]))
        exit_code = self.check_for_output_files('ska_test')
        self.cleanup('ska_test')
        os.remove(aln_out)
        assert exit_code == 0

    @staticmethod
    def check_for_output_files(prefix):
        assert os.path.exists(prefix + '.summary_of_snp_distribution.vcf')
        assert os.path.exists(prefix + '.recombination_predictions.embl')
        assert os.path.exists(prefix + '.per_branch_statistics.csv')
        assert os.path.exists(prefix + '.filtered_polymorphic_sites.fasta')
        assert os.path.exists(prefix + '.filtered_polymorphic_sites.phylip')
        assert os.path.exists(prefix + '.recombination_predictions.gff')
        assert os.path.exists(prefix + '.branch_base_reconstruction.embl')
        assert os.path.exists(prefix + '.final_tree.tre')
        assert os.path.exists(prefix + '.node_labelled.final_tree.tre')
        return 0

    @staticmethod
    def cleanup(prefix):
        #os.chdir(working_dir)
        regex_to_remove = prefix + ".*"
        for file in glob.glob(regex_to_remove):
            os.remove(file)
        tmp_to_remove = "./tmp*/*"
        for file in glob.glob(tmp_to_remove):
            os.remove(file)
        for dir in glob.glob("./tmp*"):
            if os.path.isdir(dir):
                os.rmdir(dir)

    
    @staticmethod
    def md5_check(file_path, correct_output):
        file_sum = hashlib.md5(open(file_path, "rb").read()).hexdigest()
        correct_sum = hashlib.md5(open(correct_output, "rb").read()).hexdigest()
        return file_sum == correct_sum

if __name__ == "__main__":
    unittest.main(buffer=True)

