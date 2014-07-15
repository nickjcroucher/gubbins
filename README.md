Gubbins Install
=======
This piece of software is not the version that was used in Croucher et. al. "Rapid Pneumococcal Evolution in Response to Clinical Interventions", Science 2011, 331 (6016): 430-434.
It is a work in progress and as such is not yet recommended for use by the community. We expect to make an official release of the software soon.

# Install #

== Prep work ==
As listed in the dependancies section below, you need
 * Python 2.7.

== Included external applications ==
The following applications should be installed before building gubbins:

* FastTree (>=2.1.4) http://www.microbesonline.org/fasttree/,
* a modified version of FASTML (3.0.1) http://fastml.tau.ac.il/, and
* RAxML (>=7.2.8, the full version) http://sco.h-its.org/exelixis/software.html,

On Ubuntu Saucy systems these dependencies can be installed as follows:

Install the DendroPy dependancy:

``` bash
$ wget  http://pypi.python.org/packages/source/D/DendroPy/DendroPy-3.12.0.tar.gz
$ tar xzvf DendroPy-3.12.0.tar.gz
$ cd DendroPy-3.12.0
$ sudo python setup.py install
```

Then install gubbins
``` bash
$ sudo apt-get-repository ppa:a-j-delaney/gubbins-ppa
$ sudo apt-get update
$ sudo apt-get install fasttree raxml fastml2
```

Finally, the raxmlHPC binary should be named "raxmlHPC".  If you have a modern CPU you 
can get vastly increased performance by installing RAxML yourself and selecting
the most appropriate makefile.

== Installing from a git checkout ==
``` bash
$ autoreconf -i
$ ./configure
$ su -
$ make
$ make install
$ cd python
$ python setup.py install
```

== Installing from a released source tarball ==

Please ensure that the following package is installed on your system before proceeding with the installation:

python-dev

To use this version of Gubbins you will need the following version of Python installed on your machine:

Python 2.7

You will also need the following python packages installed, if they are not then they will be downloaded and installed automatically during the Gubbins installation process:

Biopython ( >=1.59 )
DendroPy ( >=3.11.1 )
ReportLab ( >= 2.5 )

Gubbins also depends on the following pieces of software. These will be downloaded and installed automatically during the Gubbins installation process. 

FastTree ( >=2.1.4 ) http://www.microbesonline.org/fasttree/
RAxML ( >=7.2.8, the full version ) http://sco.h-its.org/exelixis/software.html
FASTML ( 2.02 ) http://fastml.tau.ac.il/

=== Installation ===

The Gubbins software can be installed directly from a pre-compiled package or from source.

==== Installing from pre-compiled package (x86-64 Ubuntu) ====

Check out a version of the repository from GitHub 

$ git clone https://github.com/sanger-pathogens/gubbins

Run the installation script to install it in your home directory

$ ./install-userspace.sh

Source the .bashrc file
	
$source ~/.bashrc

==== Installing from a git checkout ====

Check out a version of the repository from github and 

$ git clone https://github.com/sanger-pathogens/gubbins

Run the following commands to compile and install the software:

$ autoreconf -i
$ ./configure
$ su -
$ make
$ make install

==== Installing into a specific user directory ====

If you do not have permission to install the software as root and instead want to install it in a local user directory then the following commands can be used instead:

$ ./configure --prefix=~/.local
$ make
$ make install

=== Running Gubbins ===

For bash users ensure you run

$ source ~/.bashrc
 
To run the Gubbins application use:

$ run_gubbins my_alignment.fa

To see full usage of this script run: 

$ run_gubbins -h
