Quick start
=======
Before you do anything, please have a look at the [Gubbins webpage](http://sanger-pathogens.github.io/gubbins/). It contains links to the latest precompiled binaries.


Gubbins Install
=======

There are a few ways to Install/Use Gubbins, with detailed instructions below:

1. OSX - using homebrew
2. Ubuntu Trusty - using apt-get
3. Linux - using precompiled binaries or compiling from source
4. Windows - using our virtual machine



## Prep work ##
As listed in the dependancies section below, you need
 * Python 3
 * Python development headers (the python-dev package on Debian/Ubuntu).

You will also need the following python packages installed, if they are not then they will be downloaded and installed automatically during the Gubbins installation process:

 * Biopython ( >=1.59 )
 * DendroPy ( >=4.0.0 )
 * ReportLab ( >= 2.5 )

Gubbins also depends on the following pieces of software.  Depending on your installation choice these may be downloaded and installed automatically during the Gubbins installation process. 

* [FastTree](http://www.microbesonline.org/fasttree/) ( >=2.1.4 )
* [RAxML](http://sco.h-its.org/exelixis/software.html) ( >=7.2.8, the full version )
* [FASTML](http://fastml.tau.ac.il/) ( 2.02 )

If you choose to install from source, you should follow the FastTree and RAxML instructions to do so.  In the case of installing FASTML from source, we recommend using our slight modifications to the FASTML build system available at https://github.com/sanger-pathogens/fastml.

## Install ##

There are multiple ways to install gubbins depending on your requirements

1. install system-wide from binaries,
2. install per-user from binaries,
3. install system-wied from source, and
4. install per-user from source.

Each of the system-wide cases assumes you have permissions to _sudo_.  The per-user cases do not make that assumption.  Please note this software does not work on Windows (and never will), only on Linux and OSX (*nix).



### System-wide from source ###

You must first enssure that the dependencies are installed.

On a Debian/Ubuntu system
``` bash
$ sudo apt-get install python-biopython python-setuptools
$ easy_install -U dendropy
```

Alternatively, if you need to install the dependencies from source:
``` bash
$ wget https://bootstrap.pypa.io/ez_setup.py -O - | sudo python3 -
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

