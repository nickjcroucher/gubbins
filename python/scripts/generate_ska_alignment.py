#! python

# encoding: utf-8
# Wellcome Trust Sanger Institute and Imperial College London
# Copyright (C) 2020 Imperial College London and Imperial College London
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

# Generic imports
import sys
import os
import argparse
import subprocess
from multiprocessing import Pool
from functools import partial
from shutil import which

# command line parsing
def get_options():

    parser = argparse.ArgumentParser(description='Generate a ska2 alignment from a list '
                                                    'of assemblies',
                                     prog='generate_ska_alignment')

    # input options
    parser.add_argument('--reference',
                        help = 'Name of reference sequence to use for alignment',
                        required = True)
    parser.add_argument('--input',
                        help = 'List of sequence data; one row per isolate, with first column being the isolate name',
                        default = None,
                        required = False)
    parser.add_argument('--out',
                        help = 'Name of output alignment',
                        required = True)
    parser.add_argument('--k',
                        help = 'Split kmer size',
                        type = int,
                        default = 17)
    parser.add_argument('--threads',
                        help = 'Number of threads to use',
                        type = int,
                        default = 1)
    parser.add_argument('--no-cleanup',
                        help = 'Do not remove intermediate files',
                        action = 'store_true',
                        default = False)

    return parser.parse_args()

def map_fasta_sequence(seq, k = None, names = None):
    subprocess.check_output('ska fasta -o ' + names[seq]  + ' -k ' + str(k) + ' ' + seq,
                            shell = True)

def map_fastq_sequence(read_pair, k = None, names = None):
    subprocess.check_output('ska fastq -o ' + names[read_pair[0]]  + ' -k ' + str(k) + ' ' + \
                            read_pair[0] + ' ' + read_pair[1],
                            shell = True)

def ska_map_sequences(seq, k = None, ref = None):
    subprocess.check_output('ska map -o ' + seq + '.map -k ' + str(k) + ' -r ' + ref + \
                            ' ' + seq + '.skf',
                            shell = True)

# main code
if __name__ == "__main__":

    __spec__ = None

    # Get command line options
    args = get_options()

    # Check if ska is installed
    if which('ska') is None:
        sys.stderr.write('ska2 cannot be found on PATH; install with "conda install ska2"')
        sys.exit(1)

    # Check if k value is acceptable:
    if (args.k % 2) == 0 or args.k < 5 or args.k > 63:
        sys.stderr.write('k must be odd and between 5 and 63\n')
        sys.exit(1)

    # Build ska sketch
    subprocess.check_output('ska build -o ' + args.out + ' -k ' + str(args.k) + \
                            ' -f ' + args.input + ' --threads ' + str(args.threads),
                            shell = True)

    # Run ska mapping
    if os.path.exists(args.out + '.skf'):
        tmp_aln = os.path.join(os.path.dirname(args.out), 'tmp.' + os.path.basename(args.out))
        subprocess.check_output('ska map -o ' + tmp_aln + ' --threads ' + str(args.threads) + ' ' +  \
                                    args.reference + ' ' + args.out + '.skf',
                                shell = True)
    else:
        sys.stderr.write('ska building failed\n')
        sys.exit(1)
                            
    # Clean alignment to prep for Gubbins
    subprocess.check_output('gubbins_alignment_checker.py --aln ' + tmp_aln + ' --out-aln '  + args.out + \
                             ' --out ' + args.out + '.csv',
                            shell = True)
    
    sys.stderr.write("Completed generating alignment with ska2 (https://github.com/bacpop/ska.rust)\n")

    # Clean up
    if not args.no_cleanup:
        subprocess.check_output('rm ' + args.out + '.skf ' + tmp_aln,
                                shell = True)
