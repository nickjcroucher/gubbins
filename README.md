# Gubbins
Genealogies Unbiased By recomBinations In Nucleotide Sequences

PLEASE NOTE: we currently do not have the resources to provide support for Gubbins, so please do not expect a reply if you flag any issue.

[![Unmaintained](http://unmaintained.tech/badge.svg)](http://unmaintained.tech/)   
[![Build Status](https://travis-ci.org/sanger-pathogens/gubbins.svg?branch=master)](https://travis-ci.org/sanger-pathogens/gubbins)   
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-brightgreen.svg)](https://github.com/sanger-pathogens/gubbins/blob/master/LICENSE)   
[![status](https://img.shields.io/badge/NAR-10.1093-brightgreen.svg)](https://academic.oup.com/nar/article/43/3/e15/2410982)   
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/gubbins/README.html)  
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/gubbins)  
[![Docker Build Status](https://img.shields.io/docker/build/sangerpathogens/gubbins.svg)](https://hub.docker.com/r/sangerpathogens/gubbins)  
[![Docker Pulls](https://img.shields.io/docker/pulls/sangerpathogens/gubbins.svg)](https://hub.docker.com/r/sangerpathogens/gubbins)  
[![codecov](https://codecov.io/gh/sanger-pathogens/gubbins/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/gubbins)

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Required dependencies](#required-dependencies)
    * [OSX/Linux \- conda](#osxlinux---conda)
    * [Linux \- Ubuntu Xenial (16\.04) &amp; Debian (unstable)](#linux---ubuntu-xenial-1604--debian-unstable)
    * [OSX/Linux/Cloud/Windows \- Docker](#osxlinuxcloudwindows---docker)
    * [OSX/Linux \- from source](#osxlinux---from-source)
    * [OSX/Linux/Windows \- Virtual Machine](#osxlinuxwindows---virtual-machine)
    * [Running the tests](#running-the-tests)
  * [Usage](#usage)
    * [Generating Input files](#generating-input-files)
    * [Output files](#output-files)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)
  * [Further Information](#further-information)
    * [Data from the paper](#data-from-the-paper)
    * [Midpoint rerooting](#midpoint-rerooting)
    * [Ancestral sequence reconstruction](#ancestral-sequence-reconstruction)

## Introduction
Since the introduction of high-throughput, second-generation DNA sequencing technologies, there has been an enormous
increase in the size of datasets being used for estimating bacterial population phylodynamics.
Although many phylogenetic techniques are scalable to hundreds of bacterial genomes, methods which have been used
for mitigating the effect of mechanisms of horizontal sequence transfer on phylogenetic reconstructions cannot cope
with these new datasets. Gubbins (Genealogies Unbiased By recomBinations In Nucleotide Sequences) is an algorithm
that iteratively identifies loci containing elevated densities of base substitutions while concurrently constructing
a phylogeny based on the putative point mutations outside of these regions. Simulations demonstrate the algorithm
generates highly accurate reconstructions under realistic models of short-term bacterial evolution, and can be run
in only a few hours on alignments of hundreds of bacterial genome sequences.

## Installation
Before you do anything, please have a look at the [Gubbins webpage](http://sanger-pathogens.github.io/gubbins/).

### Required dependencies
* [FastTree](http://www.microbesonline.org/fasttree/#Install) ( >=2.1.4 )
* [RAxML](https://github.com/stamatak/standard-RAxML) ( >=8.0 )
* Python modules: Biopython (> 1.59), DendroPy (>=4.0), nose
* Standard build environment tools (e.g. python3, pip3, make, autoconf, libtool, gcc, check, etc...)

There are a number of ways to install Gubbins and details are provided below. If you encounter an issue when installing Gubbins please contact your local system administrator.

### OSX/Linux - conda
Install conda and enable the bioconda channels.

```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install gubbins
```

### Linux - Ubuntu Xenial (16.04) & Debian (unstable)
Gubbins has been packaged by the Debian Med team and is trivial to install using apt.

    sudo apt-get install gubbins

### OSX/Linux/Cloud/Windows - Docker
We have a docker container which gets automatically built from the latest version of Gubbins in Debian Med. To install it:

    docker pull sangerpathogens/gubbins

To use it you would use a command such as this:

    docker run --rm -it -v <local_dir>:<remote_dir> sangerpathogens/gubbins run_gubbins.py -p <remote_dir>/<prefix> [<options>] <remote_dir>/<alignment_file>

The flag `-v` synchronizes a directory on the host machine (here denoted as `<local_dir>`) with a directory in the Docker container (here denoted as `<remote_dir>`).
`<remote_dir>` does not need to exist in the container before the run, a common choice is `/data`.
Note that both `<local_dir>` and `<remote_dir>` must be absolute paths.

The input alignment file must be present in `<local_dir>` (or in one of its subdirectories).
In order to retrieve the files produced by Gubbins, run the program with option `-p`; the argument of this option must consist of `<remote_dir>`,
followed by an arbitrary identifier (here denoted as `<prefix>`).

### OSX/Linux - from source
This is the most difficult method and is only suitable for someone with advanced computing skills. Please consider using Conda instead.

Install the dependencies and include them in your `PATH`.

Clone or download the source code from GitHub and run the following commands to install Gubbins:
```
autoreconf -i
./configure [--prefix=$PREFIX]
make
[sudo] make install
cd python
[sudo] python3 setup.py install
```
Use `sudo` to install Gubbins system-wide. If you don't have the permissions, run `configure` with a prefix to install Gubbins in your home directory.

### OSX/Linux/Windows - Virtual Machine
Gubbins wont run natively on Windows but we have created a virtual machine which has all of the software setup, along with the test datasets from the paper.
It is based on [Bio-Linux 8](http://environmentalomics.org/bio-linux/).  You need to first install [VirtualBox](https://www.virtualbox.org/),
then load the virtual machine, using the 'File -> Import Appliance' menu option. The root password is 'manager'.

* ftp://ftp.sanger.ac.uk/pub/pathogens/pathogens-vm/pathogens-vm.latest.ova

More importantly though, if your trying to do bioinformatics on Windows, your not going to get very far and you should seriously consider upgrading to Linux.

### Running the tests
The test can be run from the top level directory:  

`make check`

## Usage
To run Gubbins with default settings:

    run_gubbins.py [FASTA alignment]

Input options:

    --outgroup, -o

The name of a sequence in the alignment on which to root the tree

    --starting_tree, -s

A Newick-format starting tree on which to perform the first iteration analysis. The default is to compute a starting tree using RAxML

    --filter_percentage -f

Filter out taxa with more than this percentage of missing data. Default is 25%

Processing options:

    --tree_builder, -t

The algorithm to use in the construction of phylogenies in the analysis; can be ‘raxml’, to use RAxML, ‘fasttree’,
to use Fasttree, or ‘hybrid’, to use Fasttree for the first iteration and RAxML in all subsequent iterations. Default is raxml

    --iterations, -i	
    
The maximum number of iterations to perform; the algorithm will stop earlier than this if it converges on the same tree in two successive iterations. Default is 5.

    --min_snps, -m

The minimum number of base substitutions required to identify a recombination. Default is 3.

    --converge_method, -z

Criteria to use to know when to halt iterations [weighted_robinson_foulds|robinson_foulds|recombination]. Default is weighted_robinson_foulds.

Output options:

    --use_time_stamp, -u

Include a time stamp in the name of output files to avoid overwriting previous runs on the same input file. Default is to not include a time stamp.

    --prefix, -p

Specifiy a prefix for output files. If none is provided it defaults to the name of the input FASTA alignment

    --verbose, -v

Print debugging messages. Default is off.

    --no_cleanup, -n

Do not remove files from intermediate iterations. This option will also keep other files created by RAxML and fasttree,
which would otherwise be deleted. Default is to only keep files from the final iteration.

    --raxml_model, -r

Change the model used by RAxML. The default it GTRCAT (with -V). You can set it to GTRGAMMA.

### Generating Input files
Any application which can generate a whole genome multi-FASTA alignment can be used with Gubbins, such as [Snippy](https://github.com/tseemann/snippy).

### Output files

__Prefix__
If a prefix is not defined with the –prefix option, the default prefix of the output files is:
X.Y

where:
X = Prefix taken from the input fasta file
Y = Time stamp. NOTE: This will only be included in the output file prefix if the –u flag has been selected

__Output file suffices:__

    .recombination_predictions.embl

Recombination predictions in EMBL tab file format.

    .recombination_predictions.gff

Recombination predictions in GFF3 format

    .branch_base_reconstruction.embl

Base substitution reconstruction in EMBL tab format.

    .summary_of_snp_distribution.vcf

VCF file summarising the distribution of SNPs

    .per_branch_statistics.csv

Per branch reporting of the base substitutions inside and outside recombinations events.

    .filtered_polymorphic_sites.fasta

FASTA format alignment of filtered polymorphic sites used to generate the phylogeny in the final iteration.

    .filtered_polymorphic_sites.phylip

Phylip format alignment of filtered polymorphic sites used to generate the phylogeny in the final iteration.

    .final_tree.tre

Final phylogenetic tree in newick format.

    .node_labelled.final_tree.tre

Final phylogenetic tree in newick format but with internal node labels.

## License
Gubbins is free software, licensed under [GPLv2](https://github.com/sanger-pathogens/gubbins/blob/master/LICENSE).

## Feedback/Issues
We currently do not have the resources to provide support for Gubbins. However, the community might be able to help you out if you report any issues about usage of the software to the [issues page](https://github.com/sanger-pathogens/gubbins/issues).

## Citation
If you use this software please cite:
[Croucher N. J., Page A. J., Connor T. R., Delaney A. J., Keane J. A., Bentley S. D., Parkhill J., Harris S.R.
"Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins". doi:10.1093/nar/gku1196, Nucleic Acids Research, 2014.]
(http://nar.oxfordjournals.org/content/43/3/e15)

## Further Information
For more information on this software see the [Gubbins webpage](http://sanger-pathogens.github.io/gubbins/).

### Data from the paper
* ftp://ftp.sanger.ac.uk/pub/project/pathogens/gubbins/PMEN1.aln.gz
* ftp://ftp.sanger.ac.uk/pub/project/pathogens/gubbins/ST239.aln.gz

### Midpoint rerooting
From version 1.3.5 (25/6/15) to version 1.4.6 (29/2/16) trees were not midpoint rerooted by default.
This doesnt have any effect on the recombination detection, but the output trees may not look as expected. Users are advised to upgrade to the latest version.

### Ancestral sequence reconstruction
From version 2.0.0 onwards, RAxML is used to reconstruction ancestral sequences instead of fastML.
RAxML doesnt always produce results as you would expect so the results can be lower quaility than fastML.
If you would like to stick with fastML for ancestral sequence reconstruction, please checkout and install v1.4.9.
