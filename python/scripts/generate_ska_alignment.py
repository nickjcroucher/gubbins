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

    parser = argparse.ArgumentParser(description='Generate a ska alignment from a list '
                                                    'of assemblies',
                                     prog='generate_ska_alignment')

    # input options
    parser.add_argument('--reference',
                        help = 'Name of reference sequence to use for alignment',
                        required = True)
    parser.add_argument('--fasta',
                        help = 'Two column list of names and FASTA files to include in alignment',
                        default = None,
                        required = False)
    parser.add_argument('--fastq',
                        help = 'Two/three column list of names and of FASTQ files to include in alignment',
                        default = None,
                        required = False)
    parser.add_argument('--out',
                        help = 'Name of output alignment',
                        required = True)
    parser.add_argument('--k',
                        help = 'Split kmer size',
                        type = int,
                        default = 15)
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
        sys.stderr.write('SKA cannot be found on PATH; install with "conda install ska"')
        sys.exit(1)

    # Dictionary for sequence names
    seq_names = {}
    all_names = []

    # Make split kmers from assemblies
    fasta_names = []
    if args.fasta is not None:
        # Read in FASTA assemblies
        with open(args.fasta,'r') as fasta_list:
            for line in fasta_list.readlines():
                info = line.strip().split()
                if os.path.isfile(info[1]):
                    fasta_names.append(info[1])
                    seq_names[info[1]] = info[0]
                    all_names.append(info[0])
                else:
                    sys.stderr.write('Unable to find file ' + info[1] + '\n')
        # Sketch into split kmers
        with Pool(processes = args.threads) as pool:
            pool.map(partial(map_fasta_sequence,
                                k = args.k,
                                names = seq_names),
                                fasta_names)
    
    # Make split kmers from FASTQs
    fastq_names = []
    if args.fastq is not None:
        # Read in FASTQ reads
        with open(args.fastq,'r') as fastq_list:
            for line in fastq_list.readlines():
                info = line.strip().split()
                if os.path.isfile(fns[1]) and os.path.isfile(fns[2]):
                    fastq_names.append((info[1],info[2]))
                    seq_names[info[1]] = info[0]
                    all_names.append(info[0])
                else:
                    sys.stderr.write('Unable to find files ' + fns[1] + ' and ' + fns[2] + '\n')
        # Sketch into split kmers
        with Pool(processes = args.threads) as pool:
            return_codes = pool.map(partial(map_fastq_sequence,
                                            k = args.k,
                                            names = seq_names),
                                            fastq_names)

    # Map sequences
    with Pool(processes = args.threads) as pool:
        return_codes = pool.map(partial(ska_map_sequences,
                                        k = args.k,
                                        ref = args.reference),
                                        all_names)

    # Generate alignment
    subprocess.check_output('cat ' + ' '.join([seq + '.map.aln' for seq in all_names]) + ' > ' + args.out,
                            shell = True)

    # Clean up
    if not args.no_cleanup:
        subprocess.check_output('rm ' + ' '.join([seq + '.map.aln' for seq in all_names]),
                                shell = True)
        subprocess.check_output('rm ' + ' '.join([seq + '.skf' for seq in all_names]),
                                shell = True)
