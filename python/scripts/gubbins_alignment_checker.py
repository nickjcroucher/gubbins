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

import argparse
from asyncore import write
import re
from collections import Counter

def parse_input_args():

    parser = argparse.ArgumentParser(
        description='Script to output bases in an alignment by isolate',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--aln', '-a', dest="aln",
                        help='Multifasta alignment file', required=True)
    parser.add_argument('--out-aln', '-a', dest="aln",
                        help='Reformatted alignment file', required=False)
    parser.add_argument('--out', '-o',
                        dest="out",
                        help="Out csv name for writing results too", required=True)

    return parser.parse_args()

def main(input_args):

    ## Lets set up the csv
    row_num = 0
    tot_lines = 0
    with open(input_args.aln, "r") as aln_file:
        for line in aln_file:
            if re.search("^>", line):
                row_num += 1
            tot_lines += 1

    # Going to use counter in a first pass to store sequence counts and the range of sequences
    total_base_counts = []
    iso_data = []
    total_headers = []
    tot_length_str = str(tot_lines)
    print("Running through alignment file: %s" % input_args.aln)
    print()
    with open(input_args.aln, "r") as aln_file:
        for index,line in enumerate(aln_file):
            num_zeros = len(tot_length_str) - len(str(index + 1))
            fmt_index = (("0" * num_zeros) + str(index + 1))
            
            if re.search("^>", line):
                iso = re.sub("\n","",(re.sub("^>","",line)))
                iso_data.append(iso)
                
            else:
                current_count = Counter(line.strip())
                iso_headers = sorted(list(current_count))
                add_list = list(set(iso_headers) - set(total_headers))
                if len(add_list) > 0:
                    total_headers.extend(add_list)
                total_base_counts.append(current_count)
            print("Counted %s of %s rows" % (fmt_index, tot_length_str), end="\r", flush=True)

    ## Second pass to line up all the counts across the isolates 
    print("")
    print("Assessing counts...")
    isolate_bases = []
    for count in total_base_counts:
        iso_row = []
        for header in total_headers:
            iso_row.append(count[header])
        isolate_bases.append(iso_row)

    total_headers.insert(0, "isolate")
    write_out_headers = [re.sub("-","gap",i) for i in total_headers]
    print("Writing out results...")
    with open((input_args.out), "w") as output:
        output.write(",".join(write_out_headers) + "\n")
        for i, aln_row in enumerate(isolate_bases):
            aln_row_str = list(map(str, aln_row))
            aln_row_str.insert(0, iso_data[i])
            output.write(",".join(aln_row_str) + "\n")

    print("Finished")

if __name__ == '__main__':
    input_args = parse_input_args()
    main(input_args)


