# encoding: utf-8
# Wellcome Trust Sanger Institute
# Copyright (C) 2013  Wellcome Trust Sanger Institute
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

import sys
from gubbins import utils


class FastTree:
    """Class for operations with the FastTree executable"""

    def __init__(self, verbose=False):
        """Initialises the object"""
        self.verbose = verbose
        self.tree_prefix = ""
        self.tree_suffix = ".tre"

        self.potential_executables = ["FastTree", "fasttree"]
        self.executable = utils.choose_executable(self.potential_executables)
        if self.executable is None:
            sys.exit("No usable version of FastTree could be found.")
        self.tree_building_parameters = ["-nosupport", "-gtr", "-gamma", "-nt"]

    def tree_building_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the FastTree executable"""
        output_tree = basename + self.tree_suffix
        command = [self.executable]
        command.extend(self.tree_building_parameters)
        if input_tree:
            command.extend(["-intree", input_tree])
        command.extend(["-out", output_tree])
        command.append(alignment_filename)
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)


class IQTree:
    """Class for operations with the IQTree executable"""

    def __init__(self, threads: int, internal_node_prefix="", verbose=False):
        """Initialises the object"""
        self.verbose = verbose
        self.threads = threads
        self.tree_prefix = ""
        self.tree_suffix = ".treefile"
        self.asr_prefix = ""
        self.asr_suffix = ".state"
        self.asr_tree_prefix = ""
        self.asr_tree_suffix = ".treefile"
        self.internal_node_prefix = internal_node_prefix

        self.executable = "iqtree"
        if utils.which(self.executable) is None:
            sys.exit("No usable version of IQTree could be found.")
        self.tree_building_parameters = ["-safe -m GTR+G4"]
        self.internal_sequence_reconstruction_parameters = ["-safe -asr -m GTR+G4"]

    def tree_building_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the IQTree executable"""
        command = [self.executable]
        command.extend(self.tree_building_parameters)
        command.extend(["-s", alignment_filename, "-pre", basename])
        if self.threads:
            command.extend(["-nt", str(self.threads)])
        else:
            command.extend(["-nt", "AUTO"])
        if input_tree:
            command.extend(["-t", input_tree])
        if not self.verbose:
            command.append("-quiet")
        return " ".join(command)

    def internal_sequence_reconstruction_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the IQTree executable for ancestral sequence reconstruction"""
        command = [self.executable]
        command.extend(self.internal_sequence_reconstruction_parameters)
        command.extend(["-s", alignment_filename, "-pre", basename])
        if self.threads:
            command.extend(["-nt", str(self.threads)])
        if input_tree:
            command.extend(["-te", input_tree])
        if not self.verbose:
            command.append("-quiet")
        return " ".join(command)

    def convert_raw_ancestral_states_to_fasta(self, input_filename, output_filename):
        """Converts the file containing ancestral sequences into FASTA format"""
        raw_sequences = {}
        with open(input_filename, 'r') as infile:
            for line in infile:
                if not line.startswith("#") and not line.startswith("Node\tSite\tState"):
                    elements = line.split("\t")
                    name = elements[0]
                    base = elements[2]
                    if name in raw_sequences:
                        raw_sequences[name].append(base)
                    else:
                        raw_sequences[name] = [base]

        with open(output_filename, 'w+') as outfile:
            for sequence_name, sequence_bases in raw_sequences.items():
                outfile.write(">" + self.replace_internal_node_label(sequence_name) + "\n")
                outfile.write("".join(sequence_bases) + "\n")

    def replace_internal_node_label(self, label):
        """Changes the label of internal nodes"""
        return self.internal_node_prefix + label.replace("Node", "")


class RAxML:
    """Class for operations with the RAxML executable"""

    def __init__(self, threads: int, model='GTRCAT', internal_node_prefix="", verbose=False):
        """Initialises the object"""
        self.verbose = verbose
        self.threads = threads
        self.model = model
        self.tree_prefix = "RAxML_result."
        self.tree_suffix = ""
        self.asr_prefix = "RAxML_marginalAncestralStates."
        self.asr_suffix = ""
        self.asr_tree_prefix = "RAxML_nodeLabelledRootedTree."
        self.asr_tree_suffix = ""
        self.internal_node_prefix = internal_node_prefix

        self.single_threaded_executables = ['raxmlHPC-AVX2', 'raxmlHPC-AVX', 'raxmlHPC-SSE3', 'raxmlHPC']
        self.multi_threaded_executables = ['raxmlHPC-PTHREADS-AVX2', 'raxmlHPC-PTHREADS-AVX',
                                           'raxmlHPC-PTHREADS-SSE3', 'raxmlHPC-PTHREADS']
        self.executable = self.select_executable_based_on_threads()
        if self.executable is None:
            sys.exit("No usable version of RAxML could be found.")

        self.tree_building_parameters = ["-f", "d", "-p", str(1)]
        if self.model == "GTRGAMMA":
            self.tree_building_parameters.extend(["-m", "GTRGAMMA"])
        else:
            self.tree_building_parameters.extend(["-m", "GTRCAT", "-V"])
        self.internal_sequence_reconstruction_parameters = ["-f", "A", "-p", str(1), "-m", "GTRGAMMA"]

    def tree_building_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the RAxML executable for tree building"""
        command = [self.executable]
        command.extend(self.tree_building_parameters)
        command.extend(["-s", alignment_filename, "-n", basename])
        if self.threads > 1:
            command.extend(["-T", str(self.threads)])
        if input_tree:
            command.extend(["-t", input_tree])
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)

    def internal_sequence_reconstruction_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the RAxML executable for ancestral sequence reconstruction"""
        command = [self.executable]
        command.extend(self.internal_sequence_reconstruction_parameters)
        command.extend(["-s", alignment_filename, "-n", basename])
        if self.threads > 1:
            command.extend(["-T", str(self.threads)])
        command.extend(["-t", input_tree])
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)

    def select_executable_based_on_threads(self):
        """Chooses an appropriate executable"""
        if self.threads == 1:
            single_threaded_exec = utils.choose_executable_based_on_processor(
                self.single_threaded_executables)
            if single_threaded_exec is not None:
                return single_threaded_exec
            else:
                print("Trying multithreaded version of RAxML because no single threaded version of RAxML could be "
                      "found. Just to warn you, this requires 2 threads.\n")
                self.threads = 2

        if self.threads > 1:
            multi_threaded_exec = utils.choose_executable_based_on_processor(
                self.multi_threaded_executables)
            if multi_threaded_exec is not None:
                return multi_threaded_exec
            else:
                return None

    def convert_raw_ancestral_states_to_fasta(self, input_filename, output_filename):
        """Converts the file containing ancestral sequences into FASTA format"""
        with open(input_filename, 'r') as infile:
            with open(output_filename, 'w+') as outfile:
                for sequence_line in infile:
                    [sequence_name, sequence_bases] = sequence_line.split(' ')
                    sequence_bases = sequence_bases.replace('?', 'N')
                    outfile.write('>' + self.replace_internal_node_label(sequence_name) + '\n')
                    outfile.write(sequence_bases)

    def replace_internal_node_label(self, label):
        """Changes the label of internal nodes"""
        return self.internal_node_prefix + label
