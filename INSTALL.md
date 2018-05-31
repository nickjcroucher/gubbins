Before you do anything, please have a look at the [Gubbins webpage](http://sanger-pathogens.github.io/gubbins/).

# Installation
There are a few ways to install Gubbins and its dependancies. The simpliest way is using conda.

* OSX/Linux - Bioconda
* Linux - Ubuntu Xenial (16.04) & Debian (unstable)
* OSX/Linux/Cloud/Windows - Docker
* OSX/Linux - from source
* OSX/Linux/Windows - Virtual Machine

## OSX/Linux - conda
Install conda and enable the bioconda channels.

```
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install gubbins
```

## Linux - Ubuntu Xenial (16.04) & Debian (unstable)
Gubbins has been packaged by the Debian Med team and is trivial to install using apt.
```
sudo apt-get install gubbins
```

## OSX/Linux/Cloud/Windows - Docker
We have a docker container which gets automatically built from the latest version of Gubbins in Debian Med. To install it:

```
docker pull sangerpathogens/gubbins
```

To use it you would use a command such as this (substituting in your directories), where your input files are  assumed to be stored in /home/ubuntu/data:
```
docker run --rm -it -v /home/ubuntu/data:/data sangerpathogens/gubbins run_gubbins --help
```

## OSX/Linux - from source
This is the most difficult method and is only suitable for someone with advanced computing skills. Please consider using Conda instead.

Install the dependances and include them in your PATH:
* [FastTree](http://www.microbesonline.org/fasttree/#Install) ( >=2.1.4 )
* [RAxML](https://github.com/stamatak/standard-RAxML) ( >=8.0 )
* Python modules: Biopython (> 1.59), DendroPy (>=4.0), Reportlab, nose, pillow
* Standard build environment tools (e.g. python3, pip3, make, autoconf, libtool, gcc, check, etc...)

```
autoreconf -i
./configure
make
sudo make install
cd python
sudo python3 setup.py install
```

## OSX/Linux/Windows - Virtual Machine
Roary wont run natively on Windows but we have created virtual machine which has all of the software setup, along with the test datasets from the paper. 
It is based on [Bio-Linux 8](http://environmentalomics.org/bio-linux/).  You need to first install [VirtualBox](https://www.virtualbox.org/), 
then load the virtual machine, using the 'File -> Import Appliance' menu option. The root password is 'manager'.

* ftp://ftp.sanger.ac.uk/pub/pathogens/pathogens-vm/pathogens-vm.latest.ova

More importantly though, if your trying to do bioinformatics on Windows, your not going to get very far and you should seriously consider upgrading to Linux.

