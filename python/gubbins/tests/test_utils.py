#! /usr/bin/env python3
# encoding: utf-8

"""
Tests of utilities.
"""

import unittest
import os
import subprocess
import re
import io
from contextlib import redirect_stdout
from gubbins import common, utils

modules_dir = os.path.dirname(os.path.abspath(utils.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestUtilities(unittest.TestCase):

    def test_gubbins_command(self):
        assert common.create_gubbins_command('AAA', 'BBB', 'CCC', 'DDD', 'EEE', 5, 10, 200, 0.05, 1.0, 0) \
               == 'AAA -r -v CCC -a 10 -b 200 -f EEE -t DDD -m 5 -p 0.05 -i 1.0 BBB'

    def test_translation_of_filenames_to_final_filenames(self):
        assert common.translation_of_filenames_to_final_filenames('AAA', 'test') == {
            'AAA.vcf':             'test.summary_of_snp_distribution.vcf',
            'AAA.branch_snps.tab': 'test.branch_base_reconstruction.embl',
            'AAA.tab':             'test.recombination_predictions.embl',
            'AAA.gff':             'test.recombination_predictions.gff',
            'AAA.stats':           'test.per_branch_statistics.csv',
            'AAA.snp_sites.aln':   'test.filtered_polymorphic_sites.fasta',
            'AAA.phylip':          'test.filtered_polymorphic_sites.phylip',
            'AAA.internal':        'test.node_labelled.final_tree.tre',
            'AAA':                 'test.final_tree.tre',
            'AAA.bootstrapped':    'test.final_bootstrapped_tree.tre',
            'AAA.sh_support':      'test.final_SH_support_tree.tre'
        }

    def test_check_and_fix_window_size(self):
        class ArgsObject:
            pass
        arguments = ArgsObject
        arguments.min_window_size = 1
        arguments.max_window_size = 2000000
        common.check_and_fix_window_size(arguments)
        assert arguments.min_window_size == 3
        assert arguments.max_window_size == 1000000
        arguments.min_window_size = 20
        arguments.max_window_size = 10
        common.check_and_fix_window_size(arguments)
        assert arguments.min_window_size == 10
        assert arguments.max_window_size == 20

    def test_which(self):
        # the location of ls varies depending on OS so just check end
        assert re.match('.*/ls$', utils.which('ls')) is not None
        # Strip parameters
        assert re.match('.*/ls$', utils.which('ls -alrt')) is not None
        assert utils.which('non_existent_program') is None

    def test_is_executable(self):
        program = utils.which('ls')
        assert utils.is_executable(program)
        assert not utils.is_executable('non_existent_program')

    def test_replace_executable(self):
        assert 'raxmlHPC -f d -p 1 -m GTRGAMMA' in utils.replace_executable('raxml -f d -p 1 -m GTRGAMMA', 'raxmlHPC')
        assert '../src/gubbins' in utils.replace_executable('gubbins', '../src/gubbins')

    def test_rename_files(self):
        subprocess.check_call('touch temp_file; touch another_file', shell=True)
        utils.rename_files({'temp_file': 'output_file', 'another_file': 'another_output_file'})
        assert os.path.exists('output_file')
        assert os.path.exists('another_output_file')
        os.remove('output_file')
        os.remove('another_output_file')

    def test_delete_files(self):
        open(os.path.join(data_dir, 'AAA.rex'), 'w').close()
        open(os.path.join(data_dir, 'BBB.rex'), 'w').close()
        open(os.path.join(data_dir, 'AAA.rox'), 'w').close()
        open(os.path.join(data_dir, 'BBB.rox'), 'w').close()
        utils.delete_files(data_dir, ['AAA'], ".r[eo]x")
        assert not os.path.exists(os.path.join(data_dir, 'AAA.rex'))
        assert not os.path.exists(os.path.join(data_dir, 'AAA.rox'))
        assert os.path.exists(os.path.join(data_dir, 'BBB.rex'))
        assert os.path.exists(os.path.join(data_dir, 'BBB.rox'))
        os.remove(os.path.join(data_dir, 'BBB.rex'))
        os.remove(os.path.join(data_dir, 'BBB.rox'))

    def test_do_files_exist(self):
        open(os.path.join(data_dir, 'AAA.rex'), 'w').close()
        open(os.path.join(data_dir, 'BBB.rox'), 'w').close()
        assert utils.do_files_exist(data_dir, ['AAA', 'BBB'], ".r[eo]x")
        assert not utils.do_files_exist(data_dir, ['AAA', 'BBB'], ".rax")
        os.remove(os.path.join(data_dir, 'AAA.rex'))
        os.remove(os.path.join(data_dir, 'BBB.rox'))

    def test_verbose_printer(self):
        printer = utils.VerbosePrinter(True, "-")
        assert printer.is_verbose()
        assert printer.separator() == "-"
        f = io.StringIO()
        with redirect_stdout(f):
            printer.print(["AAA", "BBB"])
            printed = f.getvalue()
        assert printed == "AAA-BBB\n"

    def test_seed(self):
        set_seed_val = utils.set_seed(42)
        assert set_seed_val == "42"
        random_seed_val = utils.set_seed(None)
        assert(int(random_seed_val) < 10001)
