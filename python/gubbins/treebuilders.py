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
import os
import subprocess
from random import randint

from Bio import SeqIO

from gubbins import utils

class Star:
    """Class for constructing star phylogenies"""
    
    def __init__(self):
        self.executable = "star phylogeny"
        self.tree_prefix = ""
        self.tree_suffix = ".tre"
    
    def tree_building_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        # Extract taxon names from alignment
        taxon_names = SeqIO.index(alignment_filename,"fasta")

        # Write tree
        star_tree_string = "("
        star_tree_string = star_tree_string + ':0.9,'.join(taxon_names.keys())
        star_tree_string = star_tree_string + ':1);' # mid point rooting fails with equidistant taxa

        # Print to file
        output_tree = basename + self.tree_suffix
        with open(output_tree,'w') as out_file:
            out_file.write(star_tree_string + '\n')
        return output_tree

class RapidNJ:
    """Class for operations with the rapidNJ executable"""

    def __init__(self, threads = int, model='GTRCAT', bootstrap = 0, verbose=False, additional_args = None):
        """Initialises the object"""
        self.verbose = verbose
        self.threads = threads
        self.tree_prefix = ""
        self.tree_suffix = ".tre"
        self.model = model
        self.additional_args = additional_args
        self.bootstrap = bootstrap

        self.executable = "rapidnj"
        if utils.which(self.executable) is None:
            sys.exit("No usable version of rapidnj could be found.")
        command = [self.executable]
        command.extend(["-i fa", "-t d", "-n"])
        command.extend(["-c", str(self.threads)])
        if self.model == 'JC':
            command.extend(["-a", "jc"])
        elif self.model == 'K2P':
            command.extend(["-a", "kim"])
        else:
            command.extend(["-a", self.model])
        # Additional arguments
        if self.additional_args is not None:
            command.extend([self.additional_args])
        self.base_command = command

    def tree_building_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the rapidNJ executable"""
        command = self.base_command.copy()
        # Alignment file needs to be first argument
        executable = command.pop(0)
        command.insert(0,alignment_filename)
        command.insert(0,executable)
        # Specify output file
        output_tree = basename + self.tree_suffix
        command.extend(["-x", output_tree])
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)
        
    def bootstrapping_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Runs a bootstrapping analysis and annotates the nodes of a summary tree"""
        command = self.base_command.copy()
        # Alignment file needs to be first argument
        executable = command.pop(0)
        command.insert(0,alignment_filename)
        command.insert(0,executable)
        # Specify output file
        output_tree = basename + self.tree_suffix + '.bootstrapped'
        # Number of bootstraps
        command.extend(["-b", str(self.bootstrap)])
        command.extend(["-x", output_tree])
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)

class FastTree:
    """Class for operations with the FastTree executable"""

    def __init__(self, threads: int, bootstrap = 0, model='GTRCAT', verbose=False, additional_args = None):
        """Initialises the object"""
        self.verbose = verbose
        self.threads = threads
        self.model = model
        self.tree_prefix = ""
        self.tree_suffix = ".tre"
        self.bootstrap = bootstrap
        self.additional_args = additional_args

        # Identify executable
        self.potential_executables = ["FastTree", "fasttree"]
        self.executable = utils.choose_executable(self.potential_executables)
        if self.executable is None:
            sys.exit("No usable version of FastTree could be found.")
        
        # Function for returning base command
        command = [self.executable]
        command.extend(["-nt"])
        if self.model == 'JC':
            pass # default model
        elif self.model == 'GTR':
            command.extend(["-gtr","-nocat"])
        elif self.model == 'GTRGAMMA':
            command.extend(["-gtr","-gamma"])
        elif self.model == 'GTRCAT':
            command.extend(["-gtr"])
        else:
            command.extend([self.model])
        # Additional arguments
        if self.additional_args is not None:
            command.extend([self.additional_args])
        self.base_command = command
        
        # Set the number of threads for parallelisation
        omp_threads_command = 'export OMP_NUM_THREADS=' + str(self.threads)
        try:
            subprocess.check_call(omp_threads_command, shell=True)
        except subprocess.SubprocessError:
            sys.exit("Failed to set number of threads for fasttree with command " + omp_threads_command)

    def tree_building_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the FastTree executable"""
        command = self.base_command.copy()
        if input_tree:
            command.extend(["-intree", input_tree])
        output_tree = basename + self.tree_suffix
        command.extend(["-nosupport"])
        command.extend(["-out", output_tree])
        command.extend(["-log", basename + '.log'])
        command.append(alignment_filename)
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)
    
    def model_fitting_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Fits a nucleotide substitution model to a tree and an alignment"""
        command = self.base_command.copy()
        command.extend(["-mllen","-nome"])
        command.extend(["-nosupport"])
        command.extend(["-intree",input_tree])
        command.extend(["-log", basename + ".log"])
        command.extend(["-out", basename + ".treefile"])
        command.extend([alignment_filename])
        return " ".join(command)
        
    def bootstrapping_command(self, alignment_filename: str, input_tree: str, basename: str, tmp: str) -> str:
        """Runs a bootstrapping analysis and annotates the nodes of a summary tree"""
        command = self.base_command.copy()
        output_tree = basename + self.tree_suffix
        command.extend(["-nosupport"])
        command.extend(["-out", tmp + "/" + basename + ".bootstrapped_trees"])
        command.extend(["-log", basename + ".log"])
        command.extend(["-n", str(self.bootstrap)])
        command.append(alignment_filename + ".bootstrapping.aln")
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)
    
    def sh_test(self, alignment_filename: str, input_tree: str, basename: str, tmp: str) -> str:
        """Runs a single branch support test"""
        command = self.base_command.copy()
        command.extend(["-mllen","-nome"])
        command.extend(["-intree",input_tree])
        command.extend(["-out",input_tree + ".sh_support"])
        command.extend([alignment_filename])
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)

class IQTree:
    """Class for operations with the IQTree executable"""

    def __init__(self, threads: 1, model: str, bootstrap = 0, internal_node_prefix="", verbose=False, additional_args = None):
        """Initialises the object"""
        self.verbose = verbose
        self.threads = threads
        self.model = model
        self.tree_prefix = ""
        self.tree_suffix = ".treefile"
        self.asr_prefix = ""
        self.asr_suffix = ".state"
        self.asr_tree_prefix = ""
        self.asr_tree_suffix = ".treefile"
        self.internal_node_prefix = internal_node_prefix
        self.bootstrap = bootstrap
        self.additional_args = additional_args
    
        # Construct base command
        self.executable = "iqtree"
        if utils.which(self.executable) is None:
            sys.exit("No usable version of IQTree could be found.")
        command = [self.executable]
        
        # Set parallelisation
        command.extend(["-nt", str(self.threads)])

        # Add flags
        command.extend(["-safe"])
        if self.model == 'JC':
            command.extend(["-m", "JC"])
        elif self.model == 'K2P':
            command.extend(["-m", "K2P"])
        elif self.model == 'HKY':
            command.extend(["-m", "HKY"])
        elif self.model == 'GTR':
            command.extend(["-m","GTR"])
        elif self.model == 'GTRGAMMA':
            command.extend(["-m","GTR+G4"])
        else:
            command.extend(["-m",self.model])
        # Additional arguments
        if self.additional_args is not None:
            command.extend([self.additional_args])
        self.base_command = command

    def tree_building_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the IQTree executable"""
        command = self.base_command.copy()
        command.extend(["-s", alignment_filename, "-pre", basename])
        if input_tree:
            command.extend(["-t", input_tree])
        if not self.verbose:
            command.append("-quiet")
        return " ".join(command)

    def internal_sequence_reconstruction_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the IQTree executable for ancestral sequence reconstruction"""
        command = self.base_command.copy()
        command.extend(["-asr"])
        command.extend(["-s", alignment_filename, "-pre", basename])
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

    def model_fitting_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Fits a nucleotide substitution model to a tree and an alignment"""
        # Using http://www.iqtree.org/doc/Advanced-Tutorial#user-defined-substitution-models
        command = self.base_command.copy()
        command.extend(["-s", alignment_filename, "-t", input_tree, "--prefix", basename, " -n 0 --mlrate", "-redo"])
        return " ".join(command)
    
    def bootstrapping_command(self, alignment_filename: str, input_tree: str, basename: str, tmp: str) -> str:
        """Runs a bootstrapping analysis"""
        command = self.base_command.copy()
        command.extend(["-s", alignment_filename, "-t", input_tree, "--prefix", tmp + "/" + basename + ".bootstrapped", "-B", str(self.bootstrap), "-wbt"])
        return " ".join(command)

    def sh_test(self, alignment_filename: str, input_tree: str, basename: str, tmp: str) -> str:
        """Runs a single branch support test"""
        command = self.base_command.copy()
        command.extend(["-s", alignment_filename])
        command.extend(["--prefix", tmp + "/" + input_tree + ".sh_support"])
        command.extend(["-te", input_tree])
        command.extend(["-alrt 0"])
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)

class RAxML:
    """Class for operations with the RAxML executable"""

    def __init__(self, threads: 1, model='GTRCAT', bootstrap = 0, internal_node_prefix="", verbose=False, additional_args = None):
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
        self.bootstrap = bootstrap
        self.additional_args = additional_args

        self.single_threaded_executables = ['raxmlHPC-AVX2', 'raxmlHPC-AVX', 'raxmlHPC-SSE3', 'raxmlHPC']
        self.multi_threaded_executables = ['raxmlHPC-PTHREADS-AVX2', 'raxmlHPC-PTHREADS-AVX',
                                           'raxmlHPC-PTHREADS-SSE3', 'raxmlHPC-PTHREADS']
        self.executable = self.select_executable_based_on_threads()
        if self.executable is None:
            sys.exit("No usable version of RAxML could be found.")
        command = [self.executable]
        
        # Set parallelisation
        if self.threads > 1:
            command.extend(["-T", str(self.threads)])

        # Add flags
        command.extend(["-safe"])
        if self.model == 'JC':
            command.extend(["-m", "GTRCAT","--JC69"])
        elif self.model == 'K2P':
            command.extend(["-m", "GTRCAT","--K80"])
        elif self.model == 'HKY':
            command.extend(["-m", "GTRCAT","--HKY85"])
        elif self.model == 'GTRCAT':
            command.extend(["-m","GTRCAT", "-V"])
        elif self.model == 'GTRGAMMA':
            command.extend(["-m","GTRGAMMA"])
        else:
            command.extend(["-m", self.model])
        # Additional arguments
        if self.additional_args is not None:
            command.extend([self.additional_args])
        self.base_command = command

    def tree_building_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the RAxML executable for tree building"""
        command = self.base_command.copy()
        command.extend(["-f", "d", "-p", str(1)])
        command.extend(["-s", alignment_filename, "-n", basename])
        if input_tree:
            command.extend(["-t", input_tree])
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)

    def internal_sequence_reconstruction_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the RAxML executable for ancestral sequence reconstruction"""
        command = self.base_command.copy()
        command.extend(["-f", "A", "-p", str(1)])
        command.extend(["-s", alignment_filename, "-n", basename])
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
    
    def model_fitting_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Fits a nucleotide substitution model to a tree and an alignment"""
        command = self.base_command.copy()
        command.extend(["-s", alignment_filename, "-n", os.path.basename(basename) + '_reconstruction', "-t", input_tree])
        command.extend(["-f e"])
        command.extend(["-w",os.path.dirname(basename)])
        return " ".join(command)
    
    def generate_alignments_for_bootstrapping(self, alignment_filename: str, basename: str, tmp: str) -> str:
        """Generates subsampled alignments for bootstrap analysis with FastTree"""
        # Generate alignments
        command = self.base_command.copy()
        command.extend(["-s", os.path.basename(alignment_filename)])
        command.extend(["-f j"])
        p_seed = str(randint(0, 10000))
        command.extend(["-b",p_seed])
        command.extend(["-#",str(self.bootstrap)])
        command.extend(["-n",basename + ".bootstrapping"])
        command.extend(["-w",tmp])
        # Then concatenate
        command.extend(["; cat", tmp + "/" + os.path.basename(alignment_filename) + ".BS* >", tmp + "/" + basename + ".bootstrapping.aln"])
        return " ".join(command)
    
    def bootstrapping_command(self, alignment_filename: str, input_tree: str, basename: str, tmp: str) -> str:
        """Runs a bootstrapping analysis and annotates the nodes of a summary tree"""
        # Run bootstraps
        command = self.base_command.copy()
        command.extend(["-s", alignment_filename, "-n", basename + ".bootstrapped_trees"])
        command.extend(["-w",tmp])
        p_seed = str(randint(0, 10000))
        command.extend(["-p",p_seed])
        command.extend(["-x",p_seed])
        command.extend(["-#",str(self.bootstrap)])
        # Output
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        command.extend([";"])
        return " ".join(command)
        
    def annotate_tree_using_bootstraps_command(self, alignment_filename: str, input_tree: str, bootstrapped_trees: str, basename: str, tmp: str) -> str:
        # Annotate tree with bootstraps
        command = self.base_command.copy()
        p_seed = str(randint(0, 10000))
        command.extend(["-p",p_seed])
        command.extend(["-f","b"])
        command.extend(["-t",input_tree])
        command.extend(["-z",bootstrapped_trees]) # "RAxML_bootstrap." + basename + ".bootstrapped_trees"
        command.extend(["-n",basename + ".bootstrapped"])
        command.extend(["-w",tmp])
        # Output
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        command.extend([";"])
        # Rename final file
        command.extend(["cp",tmp + "/RAxML_bipartitions." + basename + ".bootstrapped", basename + ".tre.bootstrapped"])
        return " ".join(command)

    def sh_test(self, alignment_filename: str, input_tree: str, basename: str, tmp: str) -> str:
        """Runs a single branch support test"""
        command = self.base_command.copy()
        p_seed = str(randint(0, 10000))
        command.extend(["-p",p_seed])
        command.extend(["-f", "J"])
        command.extend(["-s", alignment_filename, "-n", input_tree + ".sh_support"])
        command.extend(["-t", input_tree])
        command.extend(["-w",tmp])
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)

class RAxMLNG:
    """Class for operations with the RAxML executable"""

    def __init__(self, threads: 1, model: str, bootstrap = 0, internal_node_prefix = "", verbose = False, additional_args = None):
        """Initialises the object"""
        self.verbose = verbose
        self.threads = threads
        self.model = model
        self.tree_prefix = ""
        self.tree_suffix = ".raxml.bestTree"
        self.asr_prefix = "RAxML_marginalAncestralStates."
        self.asr_suffix = ""
        self.asr_tree_prefix = "RAxML_nodeLabelledRootedTree."
        self.asr_tree_suffix = ""
        self.internal_node_prefix = internal_node_prefix
        self.bootstrap = bootstrap
        self.additional_args = additional_args

        self.single_threaded_executables = ['raxml-ng']
        self.multi_threaded_executables = ['raxml-ng']
        self.executable = self.select_executable_based_on_threads()
        if self.executable is None:
            sys.exit("No usable version of RAxML could be found.")
        command = [self.executable]
        
        # Set parallelisation
        if self.threads > 1:
            command.extend(["--threads", str(self.threads)])

        # Add model
        command.extend(["-model"])
        if self.model == 'JC':
            command.extend(["JC"])
        elif self.model == 'K2P':
            command.extend(["K80"])
        elif self.model == 'HKY':
            command.extend(["HKY"])
        elif self.model == 'GTR':
            command.extend(["GTR"])
        elif self.model == 'GTRGAMMA':
            command.extend(["GTR+G"])
        else:
            command.extend([self.model])
        # Additional arguments
        if self.additional_args is not None:
            command.extend([self.additional_args])
        self.base_command = command

    def tree_building_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the RAxML executable for tree building"""
        command = self.base_command.copy()
        command.extend(["--search"])
        command.extend(["--msa", alignment_filename, "--prefix", basename])
        if input_tree:
            command.extend(["--tree", input_tree])
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        return " ".join(command)

    def internal_sequence_reconstruction_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Constructs the command to call the RAxML executable for ancestral sequence reconstruction"""
        command = self.base_command.copy()
        command.extend(["--ancestral"])
        command.extend(["--msa", alignment_filename, "--prefix", basename])
        command.extend(["--tree", input_tree])
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

    def model_fitting_command(self, alignment_filename: str, input_tree: str, basename: str) -> str:
        """Fits a nucleotide substitution model to a tree and an alignment"""
        command = self.base_command.copy()
        command.extend(["--evaluate"])
        command.extend(["--msa", alignment_filename, "--prefix", os.path.basename(basename) + '_reconstruction', "--tree", input_tree])
        command.extend([])
#        command.extend(["-w",os.path.dirname(basename)])
        return " ".join(command)

    def generate_alignments_for_bootstrapping(self, alignment_filename: str, basename: str, tmp: str) -> str:
        """Generates subsampled alignments for bootstrap analysis with FastTree"""
        # Generate alignments
        command = self.base_command.copy()
        command.extend(["--bsmsa"])
        command.extend(["--msa", alignment_filename, "--prefix", tmp + "/" + basename + ".bootstrapping"])
        command.extend(["--bs-trees",str(self.bootstrap)])
        # Then concatenate
        command.extend(["; cat", tmp + "/" + os.path.basename(alignment_filename) + ".BS* >", tmp + "/" + basename + ".bootstrapping.aln"])
        return " ".join(command)

    def bootstrapping_command(self, alignment_filename: str, input_tree: str, basename: str, tmp: str) -> str:
        """Runs a bootstrapping analysis and annotates the nodes of a summary tree"""
        # Run bootstraps
        command = self.base_command.copy()
        command.extend(["--bootstrap"])
        command.extend(["-s", alignment_filename, "-n", tmp + "/" + basename + ".bootstrapped_trees"])
        command.extend(["--bs-trees",str(self.bootstrap)])
        # Output
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        command.extend([";"])
        return " ".join(command)
        
    def annotate_tree_using_bootstraps_command(self, alignment_filename: str, input_tree: str, bootstrapped_trees: str, basename: str, tmp: str) -> str:
        # Annotate tree with bootstraps
        command = self.base_command.copy()
        command.extend(["--support"])
        command.extend(["--bs-trees",bootstrapped_trees])
        command.extend(["--tree",input_tree])
        command.extend(["-n",tmp + "/" + basename + ".bootstrapped"])
        # Output
        if not self.verbose:
            command.extend([">", "/dev/null", "2>&1"])
        command.extend([";"])
        # Rename final file
        command.extend(["cp",tmp + "/RAxML_bipartitions." + basename + ".bootstrapped", basename + ".tre.bootstrapped"])
        return " ".join(command)
