Gubbins Install
=======
This piece of software is not the version that was used in Croucher et. al. "Rapid Pneumococcal Evolution in Response to Clinical Interventions", Science 2011, 331 (6016): 430-434.
It is a work in progress and as such is not yet recommended for use by the community. We expect to make an official release of the software soon.

## Prep work ##
As listed in the dependancies section below, you need
 * Python 2.7.
 * Python development headers (the python-dev package on Debian/Ubuntu).

You will also need the following python packages installed, if they are not then they will be downloaded and installed automatically during the Gubbins installation process:

 * Biopython ( >=1.59 )
 * DendroPy ( >=3.11.1 )
 * ReportLab ( >= 2.5 )

Gubbins also depends on the following pieces of software.  Depending on your installation choice these may be downloaded and installed automatically during the Gubbins installation process. 

* [FastTree](http://www.microbesonline.org/fasttree/) ( >=2.1.4 )
* [RAxML](http://sco.h-its.org/exelixis/software.html) ( >=7.2.8, the full version )
* [FASTML](http://fastml.tau.ac.il/) ( 2.02 )

If you choose to install from source, you should follow the FastTree and RAxML instructions to do so.  In the case of installing FASTML from source, we recommend using our slight modifications to the FASTML build system available at https://github.com/AidanDelaney/fastml2.

## Install ##

There are multiple ways to install gubbins depending on your requirements

1. install system-wide from binaries,
2. install per-user from binaries,
3. install system-wied from source, and
4. install per-user from source.

Each of the system-wide cases assumes you have permissions to _sudo_.  The per-user cases do not make that assumption.

### System-wide from binaries ###

We currently only support Ubuntu 14.04 x86_64 as a system-wide binary install.  Other architectures will be added on request.

Install the DendroPy dependancy:

``` bash
$ wget  http://pypi.python.org/packages/source/D/DendroPy/DendroPy-3.12.0.tar.gz
$ tar xzvf DendroPy-3.12.0.tar.gz
$ cd DendroPy-3.12.0
$ sudo python setup.py install
```

Then install gubbins

``` bash
$ sudo apt-get-repository ppa:ap13/gubbins
$ sudo apt-get update
$ sudo apt-get install fasttree raxml fastml2 gubbins
```

If you have your own version of the raxml binary, then you can omit it from the list of packages to install.  Many users have reported vastly increased performance by installing RAxML from source by selecting the most appropriate makefile.


### System-wide from binaries for other debian based systems ###
This might work on other Debian based systems and other versions of Ubuntu, but is untested.

Install the DendroPy dependancy:

``` bash
$ wget  http://pypi.python.org/packages/source/D/DendroPy/DendroPy-3.12.0.tar.gz
$ tar xzvf DendroPy-3.12.0.tar.gz
$ cd DendroPy-3.12.0
$ sudo python setup.py install
```

Then install gubbins


```bash
echo "deb http://ppa.launchpad.net/ap13/gubbins/ubuntu trusty main" >> /etc/apt/sources.list
echo "deb-src http://ppa.launchpad.net/ap13/gubbins/ubuntu trusty main" >> /etc/apt/sources.list

sudo apt-get update
sudo apt-get install fasttree raxml fastml2 gubbins
```

### Per-user from binaries  ###

Again, we currently only support Ubuntu 14.04 x86_64 as a binary install option

Check out a version of the repository from GitHub

> $ git clone https://github.com/sanger-pathogens/gubbins

Run the installation script to install it in your home directory, spectifically into ~/.local.  The installation script ensures that the Python dependencies are also installed locally.

> $ ./install-userspace.sh

Source the .bashrc file

> $ source ~/.bashrc

### System-wide from source ###

You must first enssure that the dependencies are installed.

On a Debian/Ubuntu system
``` bash
$ sudo apt-get install python-biopython python-setuptools
$ easy_install -U dendropy
```

Alternatively, if you need to install the dependencies from source:
``` bash
$ wget https://bootstrap.pypa.io/ez_setup.py -O - | sudo python -
$ sudo easy_install -U biopython
$ sudo easy_install -U dendropy
```

Check out a version of the repository from github and

> $ git clone https://github.com/sanger-pathogens/gubbins

Run the following commands to compile and install the software:

``` bash
$ autoreconf -i
$ ./configure
$ make
$ sudo make install
```

### Per-user from source ###

If you do not have permission to install the software as root and instead want to install it in a local user directory then the following commands can be used instead:

``` bash
$ wget https://bootstrap.pypa.io/ez_setup.py -O - | python - --user
$ ~/.local/bin/easy_install -U --user biopython
$ ~/.local/bin/easy_install -U --user dendropy
$ ./configure --prefix=~/.local
$ make
$ make install
```

Ensure that

> LD_LIBRARY_PATH=~/.local/lib

and

> PATH=~/.local/bin:$PATH

This is best achieved by adding them to ~/.bashrc in the usual manner.

## Running Gubbins ##

For bash users ensure you run

> $ source ~/.bashrc

To run the Gubbins application use:

> $ run_gubbins my_alignment.fa

To see full usage of this script run:

> $ run_gubbins -h
