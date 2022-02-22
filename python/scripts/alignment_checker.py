import argparse
import re
#import pandas
#from tqdm import tqdm

def parse_input_args():

    parser = argparse.ArgumentParser(
        description='Script to output bases in an alignment by isolate',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--aln', '-a', dest="aln",
                        help='Multifasta alignment file', required=True)
    parser.add_argument('--out', '-o',
                        dest="out",
                        help="Out csv name for writing results too", required=True)

    return parser.parse_args()

def main(input_args):

    ## Lets set up the pandas csv
    row_num = 0
    tot_lines = 0
    with open(input_args.aln, "r") as aln_file:
        for line in aln_file:
            if re.search("^>", line):
                row_num += 1
            tot_lines += 1

    # Going to need 11 columns: id, A, C, G, T, a, c, g, t, N, -
    # Lets append to a list and go from there

    row_data = []
    iso_data = []
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
                if index > 0:
                    row_data.append([a_num,A_num, t_num, T_num,
                                     c_num, C_num, g_num, G_num,
                                     N_num, num_dash])
                a_num = 0
                A_num = 0
                c_num = 0
                C_num = 0
                t_num = 0
                T_num = 0
                g_num = 0
                G_num = 0
                N_num = 0
                num_dash = 0
            else:
                a_num += line.count("a")
                A_num += line.count("A")
                t_num += line.count("t")
                T_num += line.count("T")
                c_num += line.count("c")
                C_num += line.count("C")
                g_num += line.count("g")
                G_num += line.count("G")
                N_num += line.count("N")
                num_dash += line.count("-")
                if index == (tot_lines - 1):
                    row_data.append([a_num, A_num, t_num, T_num,
                                     c_num, C_num, g_num, G_num,
                                     N_num, num_dash])
            print("Finished %s of %s rows" % (fmt_index, tot_length_str), end="\r", flush=True)

    
    print("Writing out results...")
    with open((input_args.out + ".csv"), "w") as output:
        output.write(",".join(['isolate','a','A','t','T','c','C','g','G','N','gap']) + "\n")
        for i, aln_row in enumerate(row_data):
            aln_row_str = list(map(str, aln_row))
            aln_row_str.insert(0, iso_data[i])
            output.write(",".join(aln_row_str) + "\n")

    print("Finished")

if __name__ == '__main__':
    input_args = parse_input_args()
    main(input_args)


