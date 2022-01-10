# Gubbins
**G**enealogies **U**nbiased **B**y recom**B**inations **I**n **N**ucleotide **S**equences

![build](https://github.com/nickjcroucher/gubbins/workflows/build/badge.svg)  
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-brightgreen.svg)](https://github.com/nickjcroucher/gubbins/blob/master/LICENSE)   
[![status](https://img.shields.io/badge/NAR-10.1093-brightgreen.svg)](https://academic.oup.com/nar/article/43/3/e15/2410982)   
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/gubbins/README.html)  
[![codecov](https://codecov.io/gh/nickjcroucher/gubbins/branch/master/graph/badge.svg)](https://codecov.io/gh/nickjcroucher/gubbins)

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Required dependencies](#required-dependencies)
    * [OSX/Linux \- conda](#osxlinux---conda)
    * [Linux \- Ubuntu Xenial (16\.04) &amp; Debian (unstable)](#linux---ubuntu-xenial-1604--debian-unstable)
    * [OSX/Linux \- from source](#osxlinux---from-source)
    * [OSX/Linux/Windows \- Virtual Machine](#osxlinuxwindows---virtual-machine)
    * [Running the tests](#running-the-tests)
  * [Usage](#usage)
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
for mitigating the effect of horizontal sequence transfer on phylogenetic reconstructions cannot cope
with these new datasets. Gubbins (Genealogies Unbiased By recomBinations In Nucleotide Sequences) is an algorithm
that iteratively identifies loci containing elevated densities of base substitutions while concurrently constructing
a phylogeny based on the putative point mutations outside of these regions. Simulations demonstrate the algorithm
generates highly accurate reconstructions under realistic models of short-term bacterial evolution, and can be run
in only a few hours on alignments of hundreds of bacterial genome sequences.

## Installation
<<<<<<< HEAD
Before starting your analysis, please have a look at the [Gubbins webpage](http://sanger-pathogens.github.io/gubbins/), [publication](https://academic.oup.com/nar/article/43/3/e15/2410982) and [manual](docs/gubbins_manual.md).

### Required dependencies

=======
Before starting your analysis, please have a look at the [Gubbins webpage](http://nickjcroucher.github.io/gubbins/), [publication](https://academic.oup.com/nar/article/43/3/e15/2410982) and [manual](docs/gubbins_manual.md).

### Required dependencies
>>>>>>> master
Phylogenetic software:
* [RAxML](https://doi.org/10.1093/bioinformatics/btu033)
* [IQTree](https://doi.org/10.1093/molbev/msaa015)
* [RAxML-NG](https://doi.org/10.1093/bioinformatics/btz305)
* [FastTree](https://doi.org/10.1371/journal.pone.0009490)
* [Rapidnj](https://doi.org/10.1007/978-3-540-87361-7_10)

Python modules:
<<<<<<< HEAD
* Biopython (> 1.59),
=======
* Biopython (>1.59),
>>>>>>> master
* DendroPy (>=4.0)
* Scipy
* Numpy
* Multiprocessing
<<<<<<< HEAD
* Nose

See environment.yml for details. These are in addition to standard build environment tools (e.g. python >=3.8, pip3, make, autoconf, libtool, gcc, check, etc...). There are a number of ways to install Gubbins and details are provided below. If you encounter an issue when installing Gubbins please contact your local system administrator.

### OSX/Linux - conda
Install conda and enable the bioconda channels.
=======
* Numba

See environment.yml for details. These are in addition to standard build environment tools (e.g. python >=3.8, pip3, make, autoconf, libtool, gcc, check, etc...). There are a number of ways to install Gubbins and details are provided below. If you encounter an issue when installing Gubbins please contact your local system administrator.

### Recommended installation method - conda
Install conda and enable the bioconda channels. This can be done using the normal command line (Linux), with Terminal (OSX) or the Powershell (Windows versions >=10).
>>>>>>> master

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

### OSX/Linux - from source
Install the dependencies and include them in your `PATH`. Clone or download the source code from GitHub and run the following commands to install Gubbins:
```
autoreconf -i
./configure [--prefix=$PREFIX]
make
[sudo] make install
cd python
[sudo] python3 setup.py install
```
Use `sudo` to install Gubbins system-wide. If you don't have the permissions, run `configure` with a prefix to install Gubbins in your home directory.

### OSX/Linux - installing from the repository
The easiest way to install the latest version of the code from this repository is to set up a conda environment with the packages needed for installation, then remove gubbins:
```
conda create -c bioconda -n gubbins_git gubbins python=3.9
conda activate gubbins_git
conda install -c conda-forge libtool autoconf-archive automake pkg-config
conda remove --force gubbins
```

Then download and install the repository in the same environment:

```
git clone https://github.com/nickjcroucher/gubbins
cd gubbins
autoreconf -i
chmod +x configure 
./configure --prefix=$CONDA_PREFIX
make
sudo make install
cd python
python setup.py install
```

### OSX/Linux/Windows - Virtual Machine
Gubbins can be run through the Powershell in Windows versions >=10. We have also created a virtual machine which has all of the software setup, along with the test datasets from the paper.
It is based on [Bio-Linux 8](http://environmentalomics.org/bio-linux/).  You need to first install [VirtualBox](https://www.virtualbox.org/),
then load the virtual machine, using the 'File -> Import Appliance' menu option. The root password is 'manager'.

* ftp://ftp.sanger.ac.uk/pub/pathogens/pathogens-vm/pathogens-vm.latest.ova

### Running the tests
The test can be run from the top level directory:  

`make check`

## Usage
To run Gubbins with default settings:

    run_gubbins.py [FASTA alignment]

Information on on further options can be found in the [manual](docs/gubbins_manual.md).

## License
Gubbins is free software, licensed under [GPLv2](https://github.com/nickjcroucher/gubbins/blob/master/LICENSE).

## Feedback/Issues
There is no specific support for development or maintenance of Gubbins. However, we will try to help you out if you report any issues about usage of the software to the [issues page](https://github.com/nickjcroucher/gubbins/issues).

## Citation
If you use this software please cite:
[Croucher N. J., Page A. J., Connor T. R., Delaney A. J., Keane J. A., Bentley S. D., Parkhill J., Harris S.R.
"Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins". doi:10.1093/nar/gku1196, Nucleic Acids Research, 2014.]
(http://nar.oxfordjournals.org/content/43/3/e15)

## Further Information
For more information on this software see the [Gubbins webpage](http://nickjcroucher.github.io/gubbins/).

### Data from the paper
* ftp://ftp.sanger.ac.uk/pub/project/pathogens/gubbins/PMEN1.aln.gz
* ftp://ftp.sanger.ac.uk/pub/project/pathogens/gubbins/ST239.aln.gz

### Midpoint rerooting
From version 1.3.5 (25/6/15) to version 1.4.6 (29/2/16) trees were not midpoint rerooted by default.
This doesnt have any effect on the recombination detection, but the output trees may not look as expected. Users are advised to upgrade to the latest version.

### Ancestral sequence reconstruction
<<<<<<< HEAD
From version 3.0.0 onwards, Gubbins will use joint ancestral reconstructions with pyjar. Version 2 used marginal ancestral reconstruction with RAxML. Version 1 used joint ancestral reconstruction with fastML.
=======
From version 3.0.0 onwards, Gubbins will use joint ancestral reconstructions with a modified version of [pyjar](https://github.com/simonrharris/pyjar) by default. Version 2 used marginal ancestral reconstruction with RAxML; this is still available in version 3, using the `--mar` flag (IQtree can also be used for reconstruction in version >3.0.0). This may useful in cases where memory use is limiting. Version 1 used joint ancestral reconstruction with [fastML](http://fastml.tau.ac.il/).
>>>>>>> master
