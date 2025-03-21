#!/usr/bin/env python
# encoding: utf-8
#
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

# Generic imports
import os
import re
import shutil
import subprocess
import sys
import tempfile
import gzip
import time
# Phylogenetic imports
import dendropy
# Biopython imports
from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo import Consensus
from Bio.Seq import Seq
# Gubbins imports
from gubbins.PreProcessFasta import PreProcessFasta
from gubbins.ValidateFastaAlignment import ValidateFastaAlignment
from gubbins.treebuilders import FastTree, VeryFastTree, IQTree, RAxML, RAxMLNG, RapidNJ, Star
from gubbins.pyjar import jar, get_base_patterns
from gubbins import utils
from gubbins.__init__ import version
from gubbins.pyjar import jar, get_base_patterns, Pyjar
from gubbins.treebuilders import FastTree, VeryFastTree, IQTree, RAxML, RAxMLNG, RapidNJ, Star

# Phylogenetic models valid for each algorithm
tree_models = {
    'star': ['JC','GTRCAT','GTRGAMMA'],
    'raxml': ['JC','K2P','HKY','GTRCAT','GTRGAMMA'],
    'raxmlng': ['JC','K2P','HKY','GTR','GTRGAMMA'],
    'iqtree': ['JC','K2P','HKY','GTR','GTRGAMMA'],
    'fasttree': ['JC','GTRCAT','GTRGAMMA'],
    'rapidnj': ['JC','K2P'],
    'veryfasttree': ['JC','GTRCAT','GTRGAMMA']
}

def parse_and_run(input_args, program_description=""):
    """Main function of the Gubbins program"""
    start_time = time.time()
    current_directory = os.getcwd()
    printer = utils.VerbosePrinter(True, "\n")

    # Process input options
    input_args = process_input_arguments(input_args)
    # Check if the Gubbins C-program is available. If so, print a welcome message. Otherwise exit.
    os.environ["PATH"] = os.environ["PATH"] + ":/usr/lib/gubbins/"
    gubbins_exec = 'gubbins'
    if utils.which(gubbins_exec) is None:
        # Check if the Gubbins C-program is available in its source directory (for tests/CI)
        gubbins_python_dir = os.path.dirname(os.path.abspath(utils.__file__))
        gubbins_bundled_exec = os.path.abspath(os.path.join(gubbins_python_dir, '../../src/gubbins'))
        if utils.which(gubbins_bundled_exec) is None:
            sys.exit(gubbins_exec + " is not in your path")
        else:
            gubbins_exec = utils.replace_executable(gubbins_exec, gubbins_bundled_exec)
    program_version = version()
    printer.print(["\n--- Gubbins " + program_version + " ---\n", program_description])
    # Log algorithms used
    methods_log = {property:[] for property in ['citation','process','version','algorithm','model']}
    methods_log['algorithm'].append("Gubbins")
    methods_log['citation'].append("https://doi.org/10.1093/nar/gku1196")
    methods_log['process'].append("Overall")
    methods_log['version'].append(program_version)
    methods_log['model'].append("-")

    # Initialize tree builder and check if all required dependencies are available
    printer.print("\nChecking dependencies and input files...")
    current_tree_name = input_args.starting_tree
    tree_file_names = []
    internal_node_label_prefix = "internal_"
    
    # Select the algorithms used for the first iteration
    current_tree_builder, current_model_fitter, current_model, current_recon_model, extra_tree_arguments, extra_model_arguments, custom_model, custom_recon_model = return_algorithm_choices(input_args,1)
    check_model_validity(current_model,current_tree_builder,input_args.mar,current_recon_model,current_model_fitter,custom_model, custom_recon_model)
    # Initialise tree builder
    tree_builder = return_algorithm(current_tree_builder, current_model, input_args, node_labels = internal_node_label_prefix, extra = extra_tree_arguments)
    alignment_suffix = tree_builder.alignment_suffix
    methods_log = update_methods_log(methods_log, method = tree_builder, step = 'Tree constructor (1st iteration)')
    # Initialise model fitter
    model_fitter = return_algorithm(current_model_fitter, current_recon_model, input_args, node_labels = internal_node_label_prefix, extra = extra_model_arguments)
    methods_log = update_methods_log(methods_log, method = model_fitter, step = 'Model fitter (1st iteration)')
    # Initialise sequence reconstruction if MAR
    if input_args.mar:
        sequence_reconstructor = return_algorithm(input_args.seq_recon, current_recon_model, input_args, node_labels = internal_node_label_prefix, extra = input_args.seq_recon_args)
        methods_log = update_methods_log(methods_log, method = sequence_reconstructor, step = 'Sequence reconstructor (1st iteration)')

    # Initialise IQtree
    tree_dater = IQTree(threads = input_args.threads,
                            model = current_model,
                            verbose = input_args.verbose
                        )

    # Check - and potentially correct - further input parameters
    check_and_fix_window_size(input_args)

    # Get the base filename
    (base_directory, base_filename) = os.path.split(input_args.alignment_filename)
    (basename, extension) = os.path.splitext(base_filename)
    if input_args.use_time_stamp:
        time_stamp = str(int(time.time()))
        basename = basename + "." + time_stamp
    snp_alignment_filename = base_filename + ".snp_sites.aln"
    gaps_alignment_filename = base_filename + ".gaps.snp_sites.aln"
    gaps_vcf_filename = base_filename + ".gaps.vcf"
    joint_sequences_filename = base_filename + ".seq.joint.aln"

    # If restarting from a previous run
    starting_iteration = 1
    if input_args.resume is not None:
        search_itr = re.search(r'iteration_(\d+)', input_args.resume)
        if search_itr is None:
            sys.stderr.write('Resuming a Gubbins run requires a tree file name containing the phrase "iteration_X"\n')
            exit(1)
        else:
            starting_iteration = int(search_itr.group(1)) + 1
            if starting_iteration > input_args.iterations:
                sys.stderr.write('Run has already reached the number of specified iterations\n')
                exit(1)
            else:
                sys.stderr.write('Resuming Gubbins analysis at iteration ' + str(starting_iteration) + '\n')
            input_args.starting_tree = input_args.resume
            current_tree_name = input_args.starting_tree

    # Check if intermediate files from a previous run exist
    intermediate_files = [basename + ".iteration_"]
    if not input_args.no_cleanup and input_args.resume is None:
        utils.delete_files(".", intermediate_files, "", input_args.verbose)
    if utils.do_files_exist(".", intermediate_files, "", input_args.verbose) and input_args.resume is None:
        sys.exit("Intermediate files from a previous run exist. Please rerun without the --no_cleanup option "
                 "to automatically delete them or with the --use_time_stamp to add a unique prefix.")

    # Filter the input alignment and save as temporary alignment file
    # Create temporary directory for storing working copies of input files
    temp_working_dir = tempfile.mkdtemp(dir=os.getcwd())

    # Check if the input files exist and have the right format
    printer.print("\nChecking input alignment file...")
    if not os.path.exists(input_args.alignment_filename):
        sys.exit("The input alignment file " + input_args.alignment_filename + " does not exist")
    temp_alignment_filename = temp_working_dir + "/" + base_filename
    shutil.copyfile(input_args.alignment_filename, temp_alignment_filename)
    input_args.alignment_filename = temp_alignment_filename
    if not ValidateFastaAlignment(temp_alignment_filename).is_input_fasta_file_valid():
        sys.exit("The input alignment file " + input_args.alignment_filename + " is invalid")
    # Filter the input alignment and save as temporary alignment file
    printer.print("\nFiltering input alignment...")
    pre_process_fasta = PreProcessFasta(temp_alignment_filename,
                                        input_args.verbose,
                                        input_args.filter_percentage)
    taxa_removed = pre_process_fasta.remove_duplicate_sequences_and_sequences_missing_too_much_data(
        temp_alignment_filename, input_args.remove_identical_sequences)

    # Check on number of sequences remaining in alignment after validation and processing
    if input_args.pairwise:
        if number_of_sequences_in_alignment(input_args.alignment_filename) != 2:
            sys.exit("Pairwise mode should only be used for two sequences.")
    else:
        if number_of_sequences_in_alignment(input_args.alignment_filename) < 3:
            sys.exit("Three or more sequences are required for a meaningful phylogenetic analysis.")

    # If outgroup is specified, check it is still in the alignment
    if input_args.outgroup is not None:
        if input_args.outgroup in taxa_removed:
            sys.stderr.write('Outgroup removed due to proportion of missing bases\n')
            sys.exit(1)

    # Initialise tree dating algorithm if dates supplied
    if input_args.date is not None:
        if os.path.isfile(input_args.date):
            # Get sequence names from alignment
            sequence_names_in_alignment = pre_process_fasta.get_sequence_names()
            # Edit taxon names as in tree
            new_date_file = os.path.join(temp_working_dir,basename + '.dates')
            with open(input_args.date,'r') as in_dates, open(new_date_file,'w') as out_dates:
                for line in in_dates.readlines():
                    info = line.rstrip().split()
                    if len(info) == 2:
                        new_name = utils.process_sequence_names(info[0])
                        if new_name in sequence_names_in_alignment:
                            out_dates.write(new_name + '\t' + info[1] + '\n')
            input_args.date = new_date_file
        else:
            sys.stderr.write('Cannot open dates file ' + input_args.date + '\n')
            sys.exit(1)

    # If a starting tree has been provided check its validity
    # Also make sure that taxa filtered out in the previous step are removed from it
    if input_args.starting_tree is not None and input_args.starting_tree != "":
        if not os.path.exists(input_args.starting_tree):
            sys.exit("The starting tree " + input_args.starting_tree + " does not exist")
        (tree_base_directory, tree_base_filename) = os.path.split(input_args.starting_tree)
        temp_starting_tree = temp_working_dir + '/' + tree_base_filename
        shutil.copyfile(input_args.starting_tree, temp_starting_tree)
        input_args.starting_tree = temp_starting_tree
        if not is_starting_tree_valid(temp_starting_tree):
            sys.exit("The starting tree " + input_args.starting_tree + " is invalid")
        if not do_the_names_match_the_fasta_file(temp_starting_tree, temp_alignment_filename):
            sys.exit("The names in the starting tree do not match the names in the alignment file")
        filter_out_removed_taxa_from_tree(temp_starting_tree,
                                            temp_starting_tree,
                                            taxa_removed)

    printer.print("...done. Run time: {:.2f} s".format(time.time() - start_time))

    # Find all SNP sites with Gubbins
    gubbins_command = f"{gubbins_exec} \"{input_args.alignment_filename}\""
    printer.print(["\nRunning Gubbins to detect SNPs...", gubbins_command])
    try:
        subprocess.check_call(gubbins_command, shell=True)
    except subprocess.SubprocessError:
        sys.exit("Gubbins crashed, please ensure you have enough free memory")
    printer.print("...done. Run time: {:.2f} s".format(time.time() - start_time))
    reconvert_fasta_file(snp_alignment_filename, snp_alignment_filename)
    reconvert_fasta_file(gaps_alignment_filename, base_filename + ".start")
    # Start the main loop
    printer.print("\nEntering the main loop.")
    for i in range(starting_iteration, input_args.iterations+1):
        printer.print("\n*** Iteration " + str(i) + " ***")

        # Define file names
        if i == 1:
            previous_tree_name = input_args.starting_tree
            alignment_filename = base_filename + alignment_suffix
        else:
            previous_tree_name = current_tree_name
            alignment_filename = previous_tree_name + alignment_suffix

        # 1.1. Construct the tree-building command depending on the iteration and employed options
        if i == 2 or input_args.resume is not None:
            # Select the algorithms used for the subsequent iterations
            current_tree_builder, current_model_fitter, current_model, current_recon_model, extra_tree_arguments, extra_model_arguments, custom_model, custom_recon_model = return_algorithm_choices(input_args,i)
            # Pick best model through ML tests
            if input_args.best_model:
                printer.print("\nSelecting best phylogenetic model")
                current_model = select_best_models(alignment_filename,
                                                    basename,
                                                    current_tree_builder,
                                                    input_args)
                input_args.model = current_model
                if current_tree_builder != 'iqtree':
                    check_model_validity(current_model,current_tree_builder,input_args.mar,current_recon_model,current_model_fitter,custom_model, custom_recon_model)
                printer.print("Phylogeny will be constructed with a " + current_model + " model")
            # Initialise tree builder
            tree_builder = return_algorithm(current_tree_builder, current_model, input_args, node_labels = internal_node_label_prefix, extra = extra_tree_arguments)
            alignment_suffix = tree_builder.alignment_suffix
            methods_log = update_methods_log(methods_log, method = tree_builder, step = 'Tree constructor (later iterations)')
            # Update date model (should not make a difference)
            if input_args.date is not None:
                tree_dater.model = current_model

        current_basename = basename + ".iteration_" + str(i)
        current_tree_name = current_basename + ".tre"
        if previous_tree_name and input_args.first_tree_builder != "star":
            tree_building_command = tree_builder.tree_building_command(
                os.path.abspath(alignment_filename), os.path.abspath(previous_tree_name), current_basename)
        else:
            tree_building_command = tree_builder.tree_building_command(
                os.path.abspath(alignment_filename), "", current_basename)
        built_tree = temp_working_dir + "/" + tree_builder.tree_prefix + current_basename + tree_builder.tree_suffix

        # 1.2. Construct the phylogenetic tree
        if input_args.starting_tree is not None and i == 1:
            printer.print("\nCopying the starting tree...")
            shutil.copyfile(input_args.starting_tree, current_tree_name)
        else:

            printer.print(["\nConstructing the phylogenetic tree with " + tree_builder.executable + "...",
                           tree_building_command])
            if current_tree_builder == "star":
                # Move star tree into temp dir
                shutil.move(tree_builder.tree_prefix + current_basename + tree_builder.tree_suffix,
                            built_tree)
            else:
                try:
                    os.chdir(temp_working_dir)
                    subprocess.check_call(tree_building_command, shell=True)
                except subprocess.SubprocessError:
                    sys.exit("Failed while building the tree.")
                os.chdir(current_directory)
            shutil.copyfile(built_tree, current_tree_name)
        printer.print("...done. Run time: {:.2f} s".format(time.time() - start_time))

        # 2. Re-root the tree
        reroot_tree(str(current_tree_name), input_args.outgroup)
        temp_rooted_tree = temp_working_dir + "/" + current_tree_name + ".rooted"
        if input_args.tree_builder == "iqtree":
            shutil.copyfile(current_tree_name, temp_rooted_tree)
        else:
            root_tree(current_tree_name, temp_rooted_tree)

        # 3.1. Construct the command for ancestral state reconstruction depending on the iteration and employed options
        ancestral_sequence_basename = current_basename + ".internal"
        current_tree_name_with_internal_nodes = current_tree_name + ".internal"

        if not input_args.mar:
        
            # 3.2a. Joint ancestral reconstruction
            printer.print(["\nReconstructing ancestral sequences with pyjar..."])
            
            if i == starting_iteration:

                # 3.3a. Read alignment and identify unique base patterns in first iteration only
                alignment_filename = base_filename + ".start"
                alignment_type = 'fasta' # input starting polymorphism alignment file assumed to be fasta format
                ordered_sequence_names, base_pattern_bases_array, base_pattern_positions_array, max_pos = \
                                                            get_base_patterns(base_filename,
                                                                                input_args.verbose,
                                                                                threads = input_args.threads)
                # 3.3b. Record in methods log (just once)
                pyjar_method = Pyjar(current_model)
                methods_log = update_methods_log(methods_log, method = pyjar_method, step = 'Sequence reconstructor')

            # 3.4a. Re-fit full polymorphism alignment to new tree
            model_fitting_command = model_fitter.model_fitting_command(snp_alignment_filename,
                                                                os.path.abspath(temp_rooted_tree),
                                                                temp_working_dir + '/' + current_basename)
            printer.print(["\nFitting substitution model to tree...", model_fitting_command])
            try:
                subprocess.check_call(model_fitting_command, shell = True)
            except:
                sys.exit("Unable to fit model to data")

            # 3.5a. Joint ancestral reconstruction with new tree and info file in each iteration
            info_filename = model_fitter.get_info_filename(temp_working_dir,current_basename)
            recontree_filename = model_fitter.get_recontree_filename(temp_working_dir,current_basename)
            # If requested, use a time-calibrated tree for sequence reconstruction
            if input_args.date is not None and input_args.recon_with_dates:
                dating_command = tree_dater.run_time_tree(snp_alignment_filename,
                                                recontree_filename,
                                                input_args.date,
                                                temp_working_dir,
                                                base_filename,
                                                outgroup = input_args.outgroup)
                try:
                    subprocess.check_call(dating_command, shell=True)
                    recontree_filename = os.path.join(temp_working_dir,base_filename + '.timetree.nwk')
                    # Set root of reconstruction tree to match that of the current tree
                    # Cannot just midpoint root both, because the branch lengths differ between them
                    harmonise_roots(recontree_filename, temp_rooted_tree, algorithm = model_fitter.name)
                except subprocess.SubprocessError:
                    # If this fails, continue to generate rest of output
                    sys.stderr.write("Unable to use time calibrated tree for sequence reconstruction in "
                    " iteration " + str(i))
            else:
                # Set root of reconstruction tree to match that of the current tree
                # Cannot just midpoint root both, because the branch lengths differ between them
                harmonise_roots(recontree_filename, temp_rooted_tree, algorithm = model_fitter.name)
            
            printer.print(["\nRunning joint ancestral reconstruction with pyjar"])
            jar(sequence_names = ordered_sequence_names, # complete polymorphism alignment
                base_patterns = base_pattern_bases_array, # array of unique base patterns in alignment
                base_pattern_positions = base_pattern_positions_array, # nparray of positions of unique base patterns in alignment
                alignment_filename = base_filename + ".start", # gap and SNP alignment file name
                tree_filename = recontree_filename, # tree generated by model fit
                info_filename = info_filename, # file containing evolutionary model parameters
                info_filetype = input_args.model_fitter, # model fitter - format of file containing evolutionary model parameters
                output_prefix = temp_working_dir + "/" + ancestral_sequence_basename, # output prefix
                outgroup_name = input_args.outgroup, # outgroup for rooting and reconstruction
                threads = input_args.threads, # number of cores to use
                verbose = input_args.verbose,
                max_pos = max_pos)
            gaps_alignment_filename = temp_working_dir + "/" + ancestral_sequence_basename + ".joint.aln"
            raw_internal_rooted_tree_filename = temp_working_dir + "/" + ancestral_sequence_basename + ".joint.tre"
            printer.print(["\nTransferring pyjar results onto original recombination-corrected tree"])
            transfer_internal_node_labels_to_tree(raw_internal_rooted_tree_filename,
                                                  temp_rooted_tree,
                                                  current_tree_name_with_internal_nodes,
                                                  "pyjar")
            printer.print(["\nDone transfer"])
            
        else:

            # 3.2b. Marginal ancestral reconstruction with RAxML, RAxML-NG or IQTree
            sequence_reconstruction_command = sequence_reconstructor.internal_sequence_reconstruction_command(
                os.path.abspath(base_filename + alignment_suffix), os.path.abspath(temp_rooted_tree),
                ancestral_sequence_basename)
            raw_internal_sequence_filename \
                = temp_working_dir + "/" + sequence_reconstructor.asr_prefix \
                + ancestral_sequence_basename + sequence_reconstructor.asr_suffix
            processed_internal_sequence_filename = temp_working_dir + "/" + ancestral_sequence_basename + ".aln"
            raw_internal_rooted_tree_filename \
                = temp_working_dir + "/" + sequence_reconstructor.asr_tree_prefix \
                + ancestral_sequence_basename + sequence_reconstructor.asr_tree_suffix

            # 3.3b. Reconstruct the ancestral sequence
            printer.print(["\nReconstructing ancestral sequences with " + sequence_reconstructor.executable + "...",
                           sequence_reconstruction_command])
            os.chdir(temp_working_dir)
            try:
                subprocess.check_call(sequence_reconstruction_command, shell=True)
            except subprocess.SubprocessError:
                sys.exit("Failed while reconstructing the ancestral sequences.")
            os.chdir(current_directory)
            # 3.4b. Join ancestral sequences with given sequences
            current_tree_name_with_internal_nodes = current_tree_name + ".internal"
            sequence_reconstructor.convert_raw_ancestral_states_to_fasta(raw_internal_sequence_filename,
                                                                         processed_internal_sequence_filename)
            concatenate_fasta_files([snp_alignment_filename, processed_internal_sequence_filename],
                                    joint_sequences_filename)

            if input_args.seq_recon == "raxml":
                transfer_internal_node_labels_to_tree(raw_internal_rooted_tree_filename, temp_rooted_tree,
                                                  current_tree_name_with_internal_nodes, sequence_reconstructor)
            elif input_args.seq_recon == "iqtree" or input_args.seq_recon == "raxmlng":
                # IQtree returns an unrooted tree
                harmonise_roots(raw_internal_rooted_tree_filename, temp_rooted_tree)
                transfer_internal_node_labels_to_tree(raw_internal_rooted_tree_filename,
                                                 temp_rooted_tree,
                                                  current_tree_name_with_internal_nodes,
                                                  sequence_reconstructor,
                                                  use_root = False)
            else:
                sys.stderr.write("Unrecognised sequence reconstruction command: " + input_args.seq_recon + '\n')
                sys.exit()
            printer.print("...done. Run time: {:.2f} s".format(time.time() - start_time))
            # 3.5b. Reinsert gaps (cp15 note: something is wonky here, the process is at the very least terribly inefficient)
            printer.print("\nReinserting gaps into the alignment...")
            shutil.copyfile(base_filename + ".start", gaps_alignment_filename)
            reinsert_gaps_into_fasta_file(joint_sequences_filename, gaps_vcf_filename, gaps_alignment_filename)
            if not os.path.exists(gaps_alignment_filename) \
                    or not ValidateFastaAlignment(gaps_alignment_filename).is_input_fasta_file_valid():
                sys.exit("There is a problem with your FASTA file after running internal sequence reconstruction. "
                         "Please check this intermediate file is valid: " + gaps_alignment_filename)

        # Ancestral reconstruction complete
        printer.print("...done. Run time: {:.2f} s".format(time.time() - start_time))
        # 4. Detect recombination sites with Gubbins (cp15 note: copy file with internal nodes back and forth to
        # ensure all created files have the desired name structure and to avoid fiddling with the Gubbins C program)
        shutil.copyfile(current_tree_name_with_internal_nodes, current_tree_name)
        gubbins_command = create_gubbins_command(
            gubbins_exec, gaps_alignment_filename, gaps_vcf_filename, current_tree_name,
            input_args.alignment_filename, input_args.min_snps, input_args.min_window_size, input_args.max_window_size,
            input_args.p_value, input_args.trimming_ratio, input_args.extensive_search, input_args.threads)
        printer.print(["\nRunning Gubbins to detect recombinations...", gubbins_command])
        try:
            subprocess.check_call(gubbins_command, shell=True)
        except subprocess.SubprocessError:
            sys.exit("Failed while running Gubbins. Please ensure you have enough free memory")
        printer.print("...done. Run time: {:.2f} s".format(time.time() - start_time))
        shutil.copyfile(current_tree_name, current_tree_name_with_internal_nodes)
        # 5. Check for convergence
        printer.print("\nChecking for convergence...")
        remove_internal_node_labels_from_tree(current_tree_name_with_internal_nodes, current_tree_name)
        tree_file_names.append(current_tree_name)
        if i > 1:
            if input_args.converge_method == 'recombination':
                current_recomb_file, previous_recomb_files = get_recombination_files(tree_file_names)
                if have_recombinations_been_seen_before(current_recomb_file, previous_recomb_files):
                    printer.print("Convergence after " + str(i) + " iterations: Recombinations observed before.")
                    break
            else:
                if has_tree_been_seen_before(tree_file_names, input_args.converge_method):
                    printer.print("Convergence after " + str(i) + " iterations: Tree observed before.")
                    break
        printer.print("...done. Run time: {:.2f} s".format(time.time() - start_time))
    else:
        printer.print("Maximum number of iterations (" + str(input_args.iterations) + ") reached.")
    printer.print("\nExiting the main loop.")

    # 6. Run bootstrap analysis if requested
    final_aln = current_basename + ".tre" + alignment_suffix # For use with bootstrap and SH tests
    if input_args.bootstrap > 0:
        printer.print(["\nRunning bootstrap analysis..."])
        shutil.copyfile(final_aln, temp_working_dir + "/" + final_aln)
        # NJ bootstraps
        if current_tree_builder == "rapidnj":
            # Bootstraps for NJ tree have to be run in a single command - deterministic algorithm means tree assumed to be the same
            # as the final tree
            bootstrap_command = tree_builder.bootstrapping_command(os.path.abspath(final_aln), os.path.abspath(current_tree_name), temp_working_dir + "/" + current_basename)
            try:
                subprocess.check_call(bootstrap_command, shell=True)
            except subprocess.SubprocessError:
                sys.exit("Failed while running bootstrap analysis.")
            transfer_bootstraps_to_tree(temp_working_dir + "/" + current_basename + ".tre.bootstrapped",
                                                    os.path.abspath(current_tree_name),
                                                    current_basename + ".tre.bootstrapped",
                                                    outgroups = input_args.outgroup)
        # ML bootstraps
        else:
            # Define alignment and a RAxML object for bootstrapping utilities
            bootstrap_aln = final_aln
            if current_tree_builder == "raxmlng":
                bootstrap_utility = tree_builder
            else:
                bootstrap_utility = return_algorithm("raxmlng", current_model, input_args, node_labels = "")
            # Generate alignments for bootstrapping if FastTree being used
            if current_tree_builder == "fasttree":
                bootstrap_aln = generate_bootstrap_alignments(bootstrap_aln,
                                                                input_args.bootstrap,
                                                                temp_working_dir + "/" + current_basename)
            # Generate bootstrap trees
            bootstrap_command = tree_builder.bootstrapping_command(os.path.abspath(bootstrap_aln), os.path.abspath(current_tree_name), current_basename, os.path.abspath(temp_working_dir))
            try:
                subprocess.check_call(bootstrap_command, shell=True)
            except subprocess.SubprocessError:
                sys.exit("Failed while running bootstrap analysis.")
            # Annotate the final tree using the bootstraps
            bootstrapped_trees_file = tree_builder.get_bootstrapped_trees_file(temp_working_dir,current_basename)
            annotation_command = bootstrap_utility.annotate_tree_using_bootstraps_command(os.path.abspath(final_aln),
                                                                                            os.path.abspath(current_tree_name),
                                                                                            bootstrapped_trees_file,
                                                                                            current_basename,
                                                                                            os.path.abspath(temp_working_dir),
                                                                                            transfer = input_args.transfer_bootstrap)
            try:
                subprocess.check_call(annotation_command, shell=True)
            except subprocess.SubprocessError:
                sys.exit("Failed while annotating final tree with bootstrapping results.")
        printer.print("...done. Run time: {:.2f} s".format(time.time() - start_time))

    # 7. Run node branch support analysis if requested
    if input_args.sh_test:
        sh_test_command = tree_builder.sh_test(final_aln,
                                                current_tree_name,
                                                current_basename,
                                                os.path.abspath(temp_working_dir))
        try:
            subprocess.check_call(sh_test_command, shell=True)
        except subprocess.SubprocessError:
            sys.exit("Failed while running SH test.")
        reformat_sh_support(current_tree_name,
                            os.path.abspath(temp_working_dir),
                            current_tree_name,
                            algorithm = current_tree_builder,
                            outgroup = input_args.outgroup)

    # 8. Run time calibration of final tree
    if input_args.date is not None:
        dating_command = tree_dater.run_time_tree(final_aln,
                                                    current_tree_name,
                                                    input_args.date,
                                                    temp_working_dir,
                                                    basename,
                                                    outgroup = input_args.outgroup)
        try:
            subprocess.check_call(dating_command, shell=True)
        except subprocess.SubprocessError:
            # If this fails, continue to generate rest of output
            sys.stderr.write("Failed running tree time calibration with LSD.")
            input_args.date = None

    # Create the final output
    printer.print("\nCreating the final output...")
    if input_args.prefix is None:
        input_args.prefix = basename
    output_filenames_to_final_filenames = translation_of_filenames_to_final_filenames(
        current_tree_name, input_args.prefix)
    utils.rename_files(output_filenames_to_final_filenames)
    if input_args.date is not None:
        # Save final output files
        output_dating_filenames_to_final_dating_filenames = \
            translation_of_dating_filenames_to_final_filenames(temp_working_dir,
                                                                basename,
                                                                input_args.prefix)
        utils.rename_files(output_dating_filenames_to_final_dating_filenames)
        # Add dating method to methods log
        tree_dater.model = "LSD"
        tree_dater.citation = "https://doi.org/10.1093/sysbio/syv068"
        methods_log = update_methods_log(methods_log, method = tree_dater, step = 'Time calibration of tree')
    print_log(methods_log, input_args.prefix)

    # Cleanup intermediate files
    if not input_args.no_cleanup:
        shutil.rmtree(temp_working_dir)
        utils.delete_files(".", tree_file_names[:-1], intermediate_files_regex(), input_args.verbose)
        utils.delete_files(".", [base_filename], starting_files_regex(), input_args.verbose)
    printer.print("...finished. Total run time: {:.2f} s".format(time.time() - start_time))

#############
# Functions #
#############

def process_input_arguments(input_args):
    # Alter settings if pairwise comparison of sequences
    if input_args.pairwise:
        input_args.iterations = 1
        input_args.model_fitter = 'fasttree'
        input_args.tree_builder = 'star'
        input_args.model = 'GTRGAMMA'
        if input_args.mar:
            sys.stderr.write('Cannot use marginal reconstruction for a pairwise comparison; switching to joint reconstruction')
            input_args.mar = False
    else:
        if input_args.tree_builder == 'hybrid':
            input_args.first_tree_builder = "fasttree"
            input_args.tree_builder = "raxml"
        # Parse IQtree in fast mode
        elif input_args.tree_builder == "iqtree-fast":
            input_args.tree_builder = "iqtree"
            input_args.tree_args = utils.extend_args(input_args.tree_args,"--fast")
        if input_args.first_tree_builder == "iqtree-fast":
            input_args.first_tree_builder = "iqtree"
            input_args.first_tree_args = utils.extend_args(input_args.first_tree_args,"--fast")
        # Make model fitting consistent with tree building
        if input_args.model_fitter is None:
            if input_args.best_model:
                input_args.model_fitter = "iqtree"
            elif input_args.tree_builder in ['raxml', 'raxmlng', 'iqtree', 'fasttree']:
                input_args.model_fitter = input_args.tree_builder
            else:
                # Else use IQtree where not possible - raxml failed on some unrealistic test datasets
                input_args.model_fitter = 'iqtree'
        # Make sequence reconstruction consistent with tree building where possible
        if input_args.seq_recon is None:
            if input_args.best_model:
                input_args.seq_recon = "iqtree"
            elif input_args.tree_builder in ['raxml', 'raxmlng', 'iqtree']:
                input_args.seq_recon = input_args.tree_builder
            else:
                # Else use RAxML where not possible
                input_args.seq_recon = 'iqtree'
        elif not input_args.mar:
            sys.stderr.write('Sequence reconstruction uses pyjar unless the '
            '--mar flag is specified\n')
            sys.exit()
        # Only allow time calibration to be used where it makes sense
        if input_args.recon_with_dates:
            if input_args.date is None or input_args.mar:
                sys.stderr.write("Reconstruction using dates is only possible with joint reconstruction and a dates file\n")
                sys.exit()
        # Check on arguments for measures of branch support
        if input_args.bootstrap > 0 and input_args.bootstrap < 1000 and input_args.tree_builder == "iqtree":
            sys.stderr.write("IQtree requires at least 1,000 bootstrap replicates\n")
            sys.exit()
        if input_args.sh_test and input_args.tree_builder not in ["raxml","iqtree","fasttree"]:
            sys.stderr.write("SH test only available for RAxML, IQtree or Fasttree\n")
            sys.exit()
    return input_args

def check_model_validity(current_model,current_tree_builder,mar,recon_model,model_fitter,custom_model,custom_recon_model):
    # Check substitution model consistent with tree building algorithm
    invalid_model = False
    # Determine model to be used for subsequent iterations
    if not custom_model and current_model not in tree_models[current_tree_builder]:
        sys.stderr.write('Evolutionary model ' + current_model +
                        ' and algorithm ' + current_tree_builder +
                         ' are incompatible\n')
        invalid_model = True
    # Determine model to be used for ancestral state reconstruction
    if not mar:
        if not custom_recon_model and recon_model not in tree_models[model_fitter]:
            sys.stderr.write('Evolutionary model ' + recon_model +
                            ' and algorithm ' + model_fitter +
                             ' are incompatible\n')
            invalid_model = True
    # Information for rectifying incompatible combinations
    if invalid_model:
        sys.stderr.write('Available combinations are:\n')
        for algorithm in tree_models:
            models = ', '.join(tree_models[algorithm])
            sys.stderr.write(algorithm + ':\t' + models + '\n')
        sys.exit()

def return_algorithm_choices(args,i):
    # Check that at least 2 iterations will be run if customised options for 1st iteration
    if args.iterations == 1:
        if args.first_tree_builder is not None or args.first_model \
            or args.custom_first_model is not None or args.first_tree_args is not None:
            sys.stderr.write('Please do not use options specific to the first iteration when'
                             ' only one iteration is to be run\n')
            sys.exit()
    # Check on first tree builder
    if args.first_tree_builder is not None:
        # Raise error if first tree builder and starting tree
        if args.starting_tree is not None:
            sys.stderr.write('Initial tree builder is not used if a starting tree is provided\n')
            sys.exit()
    # Pick tree builder
    if args.first_tree_builder is not None and i==1:
        current_tree_builder = args.first_tree_builder
        extra_tree_arguments = args.first_tree_args
    else:
        current_tree_builder = args.tree_builder
        # Get tree builder arguments
        extra_tree_arguments = args.tree_args
    # If RAXML-NG first tree builder use search1 option to decrease runtime
    if i == 1 and current_tree_builder == "raxmlng":
        if extra_tree_arguments is None:
            extra_tree_arguments = "--search1"
        else:
            extra_tree_arguments = [extra_tree_arguments]
            extra_tree_arguments.extend(["--search1"])
            extra_tree_arguments = " ".join(extra_tree_arguments)
    # Current model
    custom_model = False
    if i == 1:
        if args.custom_first_model is not None:
            current_model = args.custom_first_model
            custom_model = True
        elif args.first_model is not None:
            current_model = args.first_model
        elif args.custom_model is not None:
            current_model = args.custom_model
            custom_model = True
        elif args.model is not None:
            current_model = args.model
        elif current_tree_builder == "rapidnj":
            current_model = "JC"
        else:
            current_model = "GTRGAMMA"
    else:
        if args.custom_model is not None:
            current_model = args.custom_model
            custom_model = True
        elif args.model is not None:
            current_model = args.model
        elif current_tree_builder == "rapidnj":
            current_model = "JC"
        else:
            current_model = "GTRGAMMA"
    # Pick model fitter and model
    custom_recon_model = False
    current_model_fitter = args.model_fitter
    if args.custom_recon_model is not None:
        current_recon_model = args.custom_recon_model
        custom_recon_model = True
    else:
        current_recon_model = args.recon_model
    extra_recon_arguments = args.model_fitter_args
    # Return choices
    return current_tree_builder, current_model_fitter, current_model, current_recon_model, extra_tree_arguments, extra_recon_arguments, custom_model, current_recon_model

def return_algorithm(algorithm_choice, model, input_args, node_labels = None, extra = None):
    initialised_algorithm = None
    if algorithm_choice == "fasttree":
        initialised_algorithm = FastTree(threads = input_args.threads, model = model, seed = input_args.seed, bootstrap = input_args.bootstrap, verbose = input_args.verbose, additional_args = extra)
    elif algorithm_choice == "veryfasttree":
        initialised_algorithm = VeryFastTree(threads = input_args.threads, model = model, seed = input_args.seed, bootstrap = input_args.bootstrap, verbose = input_args.verbose, additional_args = extra)
    elif algorithm_choice == "raxml":
        initialised_algorithm = RAxML(threads = input_args.threads, model = model, seed = input_args.seed, bootstrap = input_args.bootstrap, internal_node_prefix = node_labels, verbose = input_args.verbose, additional_args = extra)
    elif algorithm_choice == "raxmlng":
        initialised_algorithm = RAxMLNG(threads = input_args.threads, model = model, seed = input_args.seed, bootstrap = input_args.bootstrap, internal_node_prefix = node_labels, verbose = input_args.verbose, additional_args = extra)
    elif algorithm_choice == "iqtree":
        initialised_algorithm = IQTree(threads = input_args.threads, model = model, seed = input_args.seed, bootstrap = input_args.bootstrap, internal_node_prefix = node_labels, verbose = input_args.verbose, use_best = (model is None and input_args.best_model), additional_args = extra)
    elif algorithm_choice == "rapidnj":
        initialised_algorithm = RapidNJ(threads = input_args.threads, model = model, bootstrap = input_args.bootstrap, verbose = input_args.verbose, additional_args = extra)
    elif algorithm_choice == "star":
        initialised_algorithm = Star()
    else:
        sys.stderr.write("Unrecognised algorithm: " + algorithm_choice + "\n")
        sys.exit()
    return initialised_algorithm

def select_best_models(snp_alignment_filename,basename,current_tree_builder,input_args):
    model_tester = IQTree(threads = input_args.threads,
                            model = 'GTR',
                            verbose = input_args.verbose
                    )
    model_test_command = model_tester.run_model_comparison(snp_alignment_filename,basename)
    try:
        subprocess.check_call(model_test_command, shell=True)
    except subprocess.SubprocessError:
        sys.exit("Unable to identify best-fitting model")
    current_model = None
    iqtree_specific_model = None
    with gzip.open(basename + '.model.gz','rb') as model_file:
        for line in model_file:
            if line.decode().startswith('best_model_list_BIC:'):
                model_list = line.decode().rstrip().split()
                iqtree_specific_model = model_list[1]
                if current_tree_builder == 'iqtree':
                    current_model = iqtree_specific_model
                else:
                    for model in model_list:
                        model_aspects = model.split('+')
                        if model_aspects[0] == 'GTR' and 'G4' in model_aspects:
                            model_name = 'GTRGAMMA'
                        elif model_aspects[0] == 'GTR' and \
                            ('R2' in model_aspects or 'R3' in model_aspects):
                            model_name = 'GTRCAT'
                        else:
                            model_name = model_aspects[0]
                        if model_name in tree_models[current_tree_builder]:
                            current_model = model_name
                            break
    return current_model

def create_gubbins_command(gubbins_exec, alignment_filename, vcf_filename, current_tree_name,
                           original_alignment_filename, min_snps, min_window_size, max_window_size,
                           p_value, trimming_ratio, extensive_search, num_threads):
    command = [gubbins_exec, "-r", "-v", vcf_filename, "-a", str(min_window_size),
               "-b", str(max_window_size), "-f", original_alignment_filename, "-t", current_tree_name,
               "-m", str(min_snps), "-p", str(p_value), "-i", str(trimming_ratio), "-n", str(num_threads)]
    if extensive_search:
            command.append("-x")
    command.append(alignment_filename)
    return " ".join(command)

def number_of_sequences_in_alignment(filename):
    return len(get_sequence_names_from_alignment(filename))


def get_sequence_names_from_alignment(filename):
    sequence_names = []
    with open(filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence_names.append(record.id)
    return sequence_names


def is_starting_tree_valid(starting_tree):
    try:
        Phylo.read(starting_tree, 'newick')
        dendropy.Tree.get_from_path(starting_tree, 'newick', preserve_underscores=True)
    except Exception:
        print("Error with the input starting tree: Is it a valid Newick file?")
        return False
    return True


def do_the_names_match_the_fasta_file(starting_tree, alignment_filename):

    # Extract sequence names from alignment
    sequence_names = set()
    with open(alignment_filename, "r") as input_handle:
        alignments = AlignIO.parse(input_handle, "fasta")
        for alignment in alignments:
            for record in alignment:
                sequence_names.add(record.name)

    # Extract sequence names from tree
    tree = dendropy.Tree.get_from_path(starting_tree, 'newick', preserve_underscores=True)
    leaf_names = set()
    for lf in tree.leaf_nodes():
        lf.taxon.label = utils.process_sequence_names(lf.taxon.label)
        leaf_names.add(lf.taxon.label)

    # Write out modified tree
    output_tree_string = tree_as_string(tree, suppress_internal=False)
    with open(starting_tree, 'w+') as output_file:
        output_file.write(output_tree_string.replace('\'', ''))

    # Check if alignment names are a subset of the tree names
    # Superfluous taxa can be pruned from the tree
    if not sequence_names.issubset(leaf_names):
        missing_seqs = list(sequence_names - leaf_names)
        print("Error: The following sequences are present in your alignment, but not the "\
        "starting tree: ")
        for name in missing_seqs:
            print(name)
        return False

    return True


def check_and_fix_window_size(input_args):
    if input_args.min_window_size < 3:
        input_args.min_window_size = 3
    if input_args.max_window_size > 1000000:
        input_args.max_window_size = 1000000
    if input_args.min_window_size > input_args.max_window_size:
        input_args.max_window_size, input_args.min_window_size = input_args.min_window_size, input_args.max_window_size


def reconvert_fasta_file(input_filename, output_filename):
    with open(input_filename, "r") as input_handle:
        alignment = AlignIO.read(input_handle, "fasta")
    with open(output_filename, "w+") as output_handle:
        AlignIO.write(alignment, output_handle, "fasta")


def concatenate_fasta_files(input_filenames, output_filename):
    alignments = []
    for input_filename in input_filenames:
        with open(input_filename, "r") as input_handle:
            alignments.append(AlignIO.read(input_handle, "fasta"))

    with open(output_filename, "w+") as output_handle:
        for alignment in alignments:
            AlignIO.write(alignment, output_handle, "fasta")


def starting_files_regex():
    return "\\.(gaps|snp_sites|phylip|vcf|start|seq)"


def intermediate_files_regex():
    return "($|\\.(gff|vcf|snp_sites|branch_snps|phylip|stats|tab|internal))"


def root_tree(input_filename, output_filename):
    # split bi nodes and root tree
    tree = dendropy.Tree.get_from_path(input_filename, 'newick', preserve_underscores=True)
    split_all_non_bi_nodes(tree.seed_node)
    output_tree_string = tree_as_string(tree, suppress_internal=False, suppress_rooting=False)
    with open(output_filename, 'w+') as output_file:
        output_file.write(output_tree_string.replace('\'', ''))


def reroot_tree(tree_name, outgroups):
    if outgroups:
        reroot_tree_with_outgroup(tree_name, outgroups.split(','))
    else:
        reroot_tree_at_midpoint(tree_name)


def reroot_tree_with_outgroup(tree_name, outgroups):
    clade_outgroups = get_monophyletic_outgroup(tree_name, outgroups)
    tree = dendropy.Tree.get_from_path(tree_name, 'newick', preserve_underscores=True)
    outgroup_mrca = tree.mrca(taxon_labels=clade_outgroups)
    print('Edge length is: ' + str(outgroup_mrca.edge.length))
    tree.reroot_at_edge(outgroup_mrca.edge,
                        length1 = outgroup_mrca.edge.length/2,
                        length2 = outgroup_mrca.edge.length/2,
                        update_bipartitions=False)
    tree.update_bipartitions()
    output_tree_string = tree_as_string(tree, suppress_internal=False)
    with open(tree_name, 'w+') as output_file:
        output_file.write(output_tree_string.replace('\'', ''))

def reroot_tree_at_midpoint(tree_name):
    tree = dendropy.Tree.get_from_path(tree_name, 'newick', preserve_underscores=True)
    split_all_non_bi_nodes(tree.seed_node)
    tree.update_bipartitions()
    tree.deroot()
    tree.reroot_at_midpoint()
    tree.update_bipartitions()
    output_tree_string = tree_as_string(tree, suppress_internal=False)
    with open(tree_name, 'w+') as output_file:
        output_file.write(output_tree_string.replace('\'', ''))

def unroot_tree(input_filename, output_filename):
    tree = dendropy.Tree.get_from_path(input_filename, 'newick', preserve_underscores=True)
    tree.deroot()
    output_tree_string = tree_as_string(tree, suppress_internal=False)
    with open(output_filename, 'w+') as output_file:
        output_file.write(output_tree_string.replace('\'', ''))

def harmonise_roots(new_tree_fn, tree_for_root_fn, algorithm = None):
    # Read in tree and get nodes adjacent to root
    taxa = dendropy.TaxonNamespace()
    tree_for_root = dendropy.Tree.get_from_path(tree_for_root_fn,
                                                'newick',
                                                preserve_underscores=True,
                                                taxon_namespace=taxa,
                                                rooting="force-rooted")
    tree_for_root.encode_bipartitions()
    root = tree_for_root.seed_node
    root_adjacent_bipartitions = []
    for root_adjacent_node in root.child_nodes():
        root_adjacent_bipartitions.append(root_adjacent_node.edge.bipartition)
    # Now search the new tree for these nodes
    new_tree = dendropy.Tree.get_from_path(new_tree_fn,
                                            'newick',
                                            preserve_underscores=True,
                                            taxon_namespace=taxa,
                                            rooting="force-rooted")
    new_tree.encode_bipartitions()
    for node in new_tree.preorder_node_iter():
        if node.edge.bipartition in root_adjacent_bipartitions:
            half_branch_length = node.edge_length/2
            new_tree.reroot_at_edge(node.edge,
                                    length1 = half_branch_length,
                                    length2 = half_branch_length)
            break
    new_tree_string = tree_as_string(new_tree,
                                        suppress_internal=False,
                                        suppress_rooting=False)

    # Check both trees are topologically identical
    missing_bipartitions = dendropy.calculate.treecompare.find_missing_bipartitions(tree_for_root,
                                                                                    new_tree,
                                                                                    is_bipartitions_updated=False)
    
    if len(missing_bipartitions) > 0:
        sys.stderr.write('Bipartitions missing when harmonising roots between trees ' + new_tree_fn + ' and ' + tree_for_root_fn + '\n')
        sys.stderr.write('The missing bipartitions are: ' + str([str(x) for x in missing_bipartitions]) + '\n')
        if algorithm == 'FastTree':
            sys.stderr.write('This is a known issue when using FastTree to fit a phylogenetic model; use an alternative algorithm\n')
        sys.exit(1)

    # Write output
    with open(new_tree_fn, 'w+') as output_file:
        output_file.write(new_tree_string.replace('\'', ''))

def filter_out_removed_taxa_from_tree(input_filename, output_filename, taxa_removed):
    tree = dendropy.Tree.get_from_path(input_filename, 'newick', preserve_underscores=True)
    tree.prune_taxa_with_labels(taxa_removed, update_bipartitions=True)
    tree.prune_leaves_without_taxa(update_bipartitions=True)
    tree.deroot()
    output_tree_string = tree_as_string(tree)
    with open(output_filename, 'w+') as output_file:
        output_file.write(output_tree_string.replace('\'', ''))


def tree_as_string(tree, suppress_internal=True, suppress_rooting=True):

    return tree.as_string(
        schema='newick',
        suppress_leaf_taxon_labels=False,
        suppress_leaf_node_labels=True,
        suppress_internal_taxon_labels=suppress_internal,
        suppress_internal_node_labels=suppress_internal,
        suppress_rooting=suppress_rooting,
        suppress_edge_lengths=False,
        unquoted_underscores=True,
        preserve_spaces=False,
        store_tree_weights=False,
        suppress_annotations=True,
        annotations_as_nhx=False,
        suppress_item_comments=True,
        node_label_element_separator=' '
    )


def split_all_non_bi_nodes(node):
    if node.is_leaf():
        return None
    elif len(node.child_nodes()) > 2:
        split_child_nodes(node)

    for child_node in node.child_nodes():
        split_all_non_bi_nodes(child_node)

    return None


def split_child_nodes(node):
    all_child_nodes = node.child_nodes()
    # skip over the first node
    first_child = all_child_nodes.pop()
    # create a placeholder node to hang everything else off
    new_child_node = node.new_child(edge_length=0)
    # move the subtree into the placeholder
    new_child_node.set_child_nodes(all_child_nodes)
    # probably not really necessary
    node.set_child_nodes((first_child, new_child_node))


def get_monophyletic_outgroup(tree_name, outgroups):
    if len(outgroups) == 1:
        return outgroups

    tree = dendropy.Tree.get_from_path(tree_name, 'newick', preserve_underscores=True)
    tree.deroot()
    tree.update_bipartitions()

    for leaf_node in tree.mrca(taxon_labels=outgroups).leaf_nodes():
        if leaf_node.taxon.label not in outgroups:
            print("Your outgroups do not form a clade.\n  Using the first taxon " + str(outgroups[0]) +
                  " as the outgroup.\n  Taxon " + str(leaf_node.taxon.label) +
                  " is in the clade but not in your list of outgroups.")
            return [outgroups[0]]

    return outgroups


def transfer_internal_node_labels_to_tree(source_tree_filename, destination_tree_filename, output_tree_filename,
                                          sequence_reconstructor, use_root = True):
    # read source tree and extract node labels, to match with the ancestral sequence reconstruction
    taxa = dendropy.TaxonNamespace()
    source_tree = dendropy.Tree.get_from_path(source_tree_filename,
                                              'newick',
                                              preserve_underscores=True,
                                              taxon_namespace=taxa,
                                              rooting='force-rooted')
    source_tree.encode_bipartitions()
    source_internal_node_dict = {}
    for source_internal_node in source_tree.internal_nodes():
        node_bipartition = source_internal_node.edge.bipartition
        if source_internal_node.label:
            source_internal_node_dict[node_bipartition] = source_internal_node.label
        else:
            source_internal_node_dict[node_bipartition] = ''
    # read original tree and add in the labels from the ancestral sequence reconstruction
    destination_tree = dendropy.Tree.get_from_path(destination_tree_filename,
                                                    'newick',
                                                    preserve_underscores=True,
                                                    taxon_namespace=taxa,
                                                    rooting='force-rooted')
    destination_tree.encode_bipartitions()
    
    # Check both trees are topologically identical
    missing_bipartitions = dendropy.calculate.treecompare.find_missing_bipartitions(source_tree,
                                                                                    destination_tree,
                                                                                    is_bipartitions_updated=False)
    if len(missing_bipartitions) > 0:
        sys.stderr.write('Bipartitions missing when transferring node labels: ' + str([str(x) for x in missing_bipartitions]))
        sys.exit(1)
    
    root_alternative = ''
    for destination_internal_node in destination_tree.internal_nodes():
        if destination_internal_node != destination_tree.seed_node or use_root:
            node_bipartition = destination_internal_node.edge.bipartition
            if sequence_reconstructor == 'pyjar':
                try:
                    new_label = source_internal_node_dict[node_bipartition]
                except:
                    sys.stderr.write('Unable to find bipartition ' + str(node_bipartition) + '\n')
                    sys.exit(1)
            else:
                new_label = sequence_reconstructor.replace_internal_node_label(str(source_internal_node_dict[node_bipartition]))
            destination_internal_node.label = None
            destination_internal_node.taxon = dendropy.Taxon(new_label)

    # output final tree
    if not use_root:
        alt_root_name = ''
        for root_child_node in destination_tree.seed_node.child_nodes():
            alt_root_name = root_child_node.taxon.label
            break
        destination_tree.seed_node.label = None
        destination_tree.seed_node.taxon = dendropy.Taxon(alt_root_name)
    destination_tree_string = destination_tree.as_string(schema="newick",
                                                            suppress_rooting=True,
                                                            unquoted_underscores=True,
                                                            suppress_internal_node_labels=True).replace("'","").replace('\'', '')

    with open(output_tree_filename, 'w+') as output_file:
        print(destination_tree_string,
                 file=output_file,
                 end='')

def remove_internal_node_labels_from_tree(input_filename, output_filename):
    tree = dendropy.Tree.get_from_path(input_filename, 'newick', preserve_underscores=True)
    output_tree_string = tree_as_string(tree)
    with open(output_filename, 'w+') as output_file:
        output_file.write(output_tree_string.replace('\'', ''))


def reinsert_gaps_into_fasta_file(input_fasta_filename, input_vcf_file, output_fasta_filename):
    # find out where the gaps are located
    # PyVCF removed for performance reasons
    with open(input_vcf_file) as vcf_file:

        sample_names = []
        gap_position = []
        gap_alt_base = []

        for vcf_line in vcf_file:
            if re.match('^#CHROM', vcf_line) is not None:
                sample_names = vcf_line.rstrip().split('\t')[9:]
            elif re.match(r'^\d', vcf_line) is not None:
                # If the alternate is only a gap it wont have a base in this column
                if re.match('^([^\t]+\t){3}([ACGTacgt])\t([^ACGTacgt])\t', vcf_line) is not None:
                    m = re.match('^([^\t]+\t){3}([ACGTacgt])\t([^ACGTacgt])\t', vcf_line)
                    gap_position.append(1)
                    gap_alt_base.append(m.group(2))
                elif re.match('^([^\t]+\t){3}([^ACGTacgt])\t([ACGTacgt])\t', vcf_line) is not None:
                    # sometimes the ref can be a gap only
                    m = re.match('^([^\t]+\t){3}([^ACGTacgt])\t([ACGTacgt])\t', vcf_line)
                    gap_position.append(1)
                    gap_alt_base.append(m.group(3))
                else:
                    gap_position.append(0)
                    gap_alt_base.append('-')

        gapped_alignments = []
        # interleave gap only and snp bases
        with open(input_fasta_filename, "r") as input_handle:
            alignments = AlignIO.parse(input_handle, "fasta")
            for alignment in alignments:
                for record in alignment:
                    inserted_gaps = []
                    if record.id in sample_names:
                        # only apply to internal nodes
                        continue
                    gap_index = 0
                    for input_base in record.seq:
                        while gap_index < len(gap_position) and gap_position[gap_index] == 1:
                            inserted_gaps.append(gap_alt_base[gap_index])
                            gap_index += 1
                        if gap_index < len(gap_position):
                            inserted_gaps.append(input_base)
                            gap_index += 1

                    while gap_index < len(gap_position):
                        inserted_gaps.append(gap_alt_base[gap_index])
                        gap_index += 1

                    record.seq = Seq(''.join(inserted_gaps))
                    gapped_alignments.append(record)

        with open(output_fasta_filename, "a") as output_handle:
            AlignIO.write(MultipleSeqAlignment(gapped_alignments), output_handle, "fasta")
            output_handle.close()
    return


def get_recombination_files(basenames):
    previous_files = []
    for name in basenames:
        previous_files.append(name + ".tab")
    current_file = previous_files.pop()
    return current_file, previous_files


def have_recombinations_been_seen_before(current_file, previous_files):
    if not os.path.exists(current_file):
        return False
    current_file_recombinations = extract_recombinations_from_embl(current_file)

    for previous_file in previous_files:
        if not os.path.exists(previous_file):
            continue
        previous_file_recombinations = extract_recombinations_from_embl(previous_file)
        if current_file_recombinations == previous_file_recombinations:
            return True
    return False


def extract_recombinations_from_embl(filename):
    with open(filename, "r") as fh:
        sequences_to_coords = {}
        start_coord = -1
        end_coord = -1
        for line in fh:
            search_obj = re.search(r'misc_feature    ([\d]+)..([\d]+)$', line)
            if search_obj is not None:
                start_coord = int(search_obj.group(1))
                end_coord = int(search_obj.group(2))
                continue

            if start_coord >= 0 and end_coord >= 0:
                search_taxa = re.search(r'taxa\=\"([^"]+)\"', line)
                if search_taxa is not None:
                    taxa_names = search_taxa.group(1).strip().split(' ')
                    for taxa_name in taxa_names:
                        if taxa_name in sequences_to_coords:
                            sequences_to_coords[taxa_name].append([start_coord, end_coord])
                        else:
                            sequences_to_coords[taxa_name] = [[start_coord, end_coord]]

                    start_coord = -1
                    end_coord = -1
                continue
        fh.close()
    return sequences_to_coords


def has_tree_been_seen_before(tree_file_names, converge_method):
    if len(tree_file_names) <= 2:
        return False

    tree_files_which_exist = []
    for tree_file_name in tree_file_names:
        if os.path.exists(tree_file_name):
            tree_files_which_exist.append(tree_file_name)

    for tree_file_name in tree_files_which_exist:
        if tree_file_name is not tree_files_which_exist[-1]:
            if converge_method == 'weighted_robinson_foulds':
                current_rf_distance = robinson_foulds_distance(
                    tree_file_name, tree_files_which_exist[-1])
                if current_rf_distance == 0.0:
                    return True
            else:
                current_rf_distance = symmetric_difference(
                    tree_file_name, tree_files_which_exist[-1])
                if current_rf_distance == 0.0:
                    return True

    return False


def robinson_foulds_distance(input_tree_name, output_tree_name):
    tns = dendropy.TaxonNamespace()
    input_tree = dendropy.Tree.get_from_path(input_tree_name, 'newick', taxon_namespace=tns)
    output_tree = dendropy.Tree.get_from_path(output_tree_name, 'newick', taxon_namespace=tns)
    input_tree.encode_bipartitions()
    output_tree.encode_bipartitions()
    return dendropy.calculate.treecompare.weighted_robinson_foulds_distance(input_tree, output_tree)


def symmetric_difference(input_tree_name, output_tree_name):
    tns = dendropy.TaxonNamespace()
    input_tree = dendropy.Tree.get_from_path(input_tree_name, 'newick', taxon_namespace=tns)
    output_tree = dendropy.Tree.get_from_path(output_tree_name, 'newick', taxon_namespace=tns)
    input_tree.encode_bipartitions()
    output_tree.encode_bipartitions()
    return dendropy.calculate.treecompare.symmetric_difference(input_tree, output_tree)

def generate_bootstrap_alignments(bootstrap_aln, n, output_aln_prefix):
    snp_aln = list(AlignIO.parse(bootstrap_aln, "fasta"))[0]
    with open(output_aln_prefix + '.bootstrapping.aln', "w+") as output_handle:
        for aln in (Phylo.Consensus.bootstrap(snp_aln, n)):
            AlignIO.write(aln, output_handle, "phylip")
    return output_aln_prefix

def transfer_bootstraps_to_tree(source_tree_filename, destination_tree_filename, output_tree_filename, outgroups = None):
    # read source tree and extract bootstraps as node labels, match with bipartition
    reroot_tree(source_tree_filename, outgroups = outgroups)
    source_tree = dendropy.Tree.get_from_path(source_tree_filename, 'newick', preserve_underscores=True)
    source_bootstraps = {}
    for source_internal_node in source_tree.internal_nodes():
        leaves = []
        for leaf in source_internal_node.leaf_iter():
            leaves.append(leaf.taxon.label.replace("'",""))
        descendant_taxa = frozenset(leaves)
        if source_internal_node.label:
            source_bootstraps[descendant_taxa] = source_internal_node.label
        else:
            source_bootstraps[descendant_taxa] = ''
    # read original tree and add in the bootstrap values
    destination_tree = dendropy.Tree.get_from_path(destination_tree_filename, 'newick', preserve_underscores=True)
    for destination_internal_node in destination_tree.internal_nodes():
        leaves = []
        for descendant in destination_internal_node.leaf_iter():
            leaves.append(descendant.taxon.label.replace("'",""))
        descendant_taxa = frozenset(leaves)
        if descendant_taxa in source_bootstraps:
            destination_internal_node.label = source_bootstraps[descendant_taxa]
        else:
            sys.stderr.write('Cannot find the internal node with descendants ' + str(descendant_taxa) + '\n')
            destination_internal_node.label = "NA"
    # output final tree
    output_tree_string = tree_as_string(destination_tree, suppress_internal=False, suppress_rooting=False)
    with open(output_tree_filename, 'w+') as output_file:
        output_file.write(output_tree_string.replace('\'', ''))

def reformat_sh_support(tree_name, tmpdir, final_tree_fn, algorithm = "raxml", outgroup = None):
    # Tree file name
    outtree_fn = tree_name + ".sh_support"
    # Read in tree
    if algorithm == "raxml":
        intree_fn = tmpdir + "/RAxML_fastTreeSH_Support." + tree_name + ".sh_support"
    elif algorithm == "iqtree":
        intree_fn = tmpdir + "/" + tree_name + ".sh_support.treefile"
    # Change SH support value formatting
    if algorithm == "raxml" or algorithm == "iqtree":
        with open(intree_fn,'r') as intree, open(tmpdir + "/" + outtree_fn,'w') as outtree:
            for line in intree.readlines():
                if algorithm == "raxml":
                    new_line = re.sub(r':(\d*[.]?\d*)\[(\d+)\]', '\\2:\\1', line)
                elif algorithm == "iqtree":
                    new_line = re.sub(r'\/', '', line)
                outtree.write(new_line)
    # Transfer SH support to final tree
    transfer_bootstraps_to_tree(tmpdir + "/" + outtree_fn,
                                final_tree_fn,
                                outtree_fn,
                                outgroups = outgroup)
    
    

def update_methods_log(log, method = None, step = ''):
    """Record methods used at each step"""
    log['citation'].append(method.citation)
    log['process'].append(step)
    log['version'].append(method.version)
    log['algorithm'].append(method.executable)
    log['model'].append(method.model)
    return log

def print_log(log, prefix):
    """Print a records of the methods used"""
    log_file_name = prefix + ".log"
    with open(log_file_name,'w') as log_file:
        log_file.write("Process,Algorithm,Version,Model,Citation\n")
        for index,process in enumerate(log['process']):
            log_file.write(process + "," + log['algorithm'][index] + "," + log['version'][index] + "," + log['model'][index] + "," + log['citation'][index] + "\n")

def translation_of_filenames_to_final_filenames(input_prefix, output_prefix):
    input_names_to_output_names = {
        str(input_prefix) + ".vcf":             str(output_prefix) + ".summary_of_snp_distribution.vcf",
        str(input_prefix) + ".branch_snps.tab": str(output_prefix) + ".branch_base_reconstruction.embl",
        str(input_prefix) + ".tab":             str(output_prefix) + ".recombination_predictions.embl",
        str(input_prefix) + ".gff":             str(output_prefix) + ".recombination_predictions.gff",
        str(input_prefix) + ".stats":           str(output_prefix) + ".per_branch_statistics.csv",
        str(input_prefix) + ".snp_sites.aln":   str(output_prefix) + ".filtered_polymorphic_sites.fasta",
        str(input_prefix) + ".phylip":          str(output_prefix) + ".filtered_polymorphic_sites.phylip",
        str(input_prefix) + ".internal":        str(output_prefix) + ".node_labelled.final_tree.tre",
        str(input_prefix) + ".bootstrapped":    str(output_prefix) + ".final_bootstrapped_tree.tre",
        str(input_prefix) + ".sh_support":      str(output_prefix) + ".final_SH_support_tree.tre",
        str(input_prefix):                      str(output_prefix) + ".final_tree.tre"
    }
    return input_names_to_output_names

def translation_of_dating_filenames_to_final_filenames(temp_working_dir,
                                                                basename,
                                                                prefix):
    dating_files = {
        os.path.join(temp_working_dir,basename + '.timetree.lsd'): prefix + '.lsd.out',
        os.path.join(temp_working_dir,basename + '.timetree.nwk'): prefix + '.final_tree.timetree.tre'
    }
    return dating_files
