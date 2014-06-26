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

``` bash
$ tar xzvf gubbins-0.1.3.tar.gz
$ cd gubbins-0.1.3
$ ./configure
$ make
$ su -
$ make install
$ cd python
$ python setup.py install
```


To run the script with all of the default bundled executables.
```bash
$ cd python
$ scripts/run_gubbins my_alignment.fa
```

To run the script anywhere you must first install the external dependances must be available in your PATH:

```bash
$ which raxmlHPC
$ which fastml
$ which fastTree
```

Then you can run this from any directory:
$ run_gubbins my_alignment.fa

The main executable should be in ./src/gubbins and in the system-wide directory
for locally compiled applications.  On Linux this is /usr/local/bin.

On Linux, you should ensure that /usr/local/bin is in your $PATH and that
/usr/local/lib (or /usr/local/lib64 on 64bit machines) are in $LD_LIBRARY_PATH.

== Installing into a specific user directory ==

It is common to want to build and install an application so that it is not
system-wide and is only available to your user.  On a Linux system this is
generally achieved by installing into the .local directory in your home
directory.

$ ./configure --prefix=~/.local
$ make
$ make install

On Linux, you should ensure that ~/.local/bin is in your $PATH and that
~/.local/lib (or ~/.local/lib64 on 64bit machines) are in $LD_LIBRARY_PATH.

If you are in the situation where your home directory is mounted over NFS (or
similar) and shared between multiple Linux distros (for example, 32bit Debian on
a cluster and 64bit Ubuntu on a desktop) then installing per-user in this manner
is probably not a good idea.

== Dependancies ==

If these applications are not available in your PATH, less efficient ones will be used:
FastTree (>=2.1.4) http://www.microbesonline.org/fasttree/
RAxML (>=7.2.8, the full version) http://sco.h-its.org/exelixis/software.html
FASTML (2.02) http://fastml.tau.ac.il/

== Python ==
Python 2.7 is required

These Python packages must be in your PYTHONPATH or will be automatically 
installed:
Biopython (>=1.59)
DendroPy  (>=3.11.1)
ReportLab (>= 2.5)

If you get errors from Biopython make sure you have the following packages installed (as named by Ubuntu):
python-dev
python-reportlab
python-numpy
python-support

== Running the Python wrapper ==

There is a Python wrapper script that automates some of the preparatory work
that is necessary to run Gubbins.

To run the application:
run_gubbbins.py my_alignment.aln
