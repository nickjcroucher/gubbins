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
        aln_cmd = "alignment_checker.py --aln " + small_aln + " --out " + output_file 
        subprocess.check_call(aln_cmd, shell=True)
        assert self.md5_check(output_csv, "b409830c1ee7772d551c13abe5cde7a7")
        os.remove(output_csv)

    @staticmethod
    def md5_check(file_path, correct_md5sum):
        file_sum = hashlib.md5(open(file_path, "rb").read()).hexdigest()
        return file_sum == correct_md5sum

if __name__ == "__main__":
    unittest.main(buffer=True)

