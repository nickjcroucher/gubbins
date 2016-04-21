Before you do anything, please have a look at the [Gubbins webpage](http://sanger-pathogens.github.io/gubbins/).

# Installation
There are a few ways to install Gubbins and its dependancies. The simpliest way is using HomeBrew (OSX) or LinuxBrew.

* OSX - Mavericks (10.9) & Yosemite (10.10) & El Capitan (10.11)
* OSX - Mountain Lion (10.8)
* Linux - Ubuntu Trusty (14.04) & Precise (12.04)
* Linux - Ubuntu Xenial (16.04) & Debian (unstable)
* Linux - CentOS 7
* Linux - CentOS 6
* OSX/Linux - from source
* OSX/Linux/Windows - Virtual Machine


## OSX - Mavericks (10.9) & Yosemite (10.10) & El Capitan (10.11)
Install [HomeBrew](http://brew.sh/). It requires a minimum of Xcode 5.1.1 (xcodebuild -version). Then run:
```
brew tap homebrew/science
brew install gubbins
```

## OSX - Mountain Lion (10.8)
Install [HomeBrew](http://brew.sh/). It requires a minimum of Xcode 5.1.1 (xcodebuild -version).

Manually install [FastML](http://fastml.tau.ac.il/source.php) and include the binary in your PATH. For example:
```
wget 'http://fastml.tau.ac.il/source/FastML.v3.1.tgz'
tar -xzf FastML.v3.1.tgz
cd FastML.v3.1
export PATH=${HOME}/FastML.v3.1/bin:$PATH
```
Then run:
```
brew tap homebrew/science
brew install python3
brew install gubbins --without-fastml
```

## OSX - It failed to install
* Run 'brew doctor' and correct any errors with your homebrew setup.
* Make sure Xcode is 5.1.1 or greater (xcodebuild -version). 
* Install it in /usr/local. Homebrew warn 'Pick another prefix at your peril!'.
* Run 'brew install -vd gubbins' and try and correct any errors.

## Linux - Ubuntu Trusty (14.04) & Precise (12.04)
Tested on Ubuntu Trusty (14.04) and Precise (12.04). Install [LinuxBrew](http://brew.sh/linuxbrew/). Then run:

```
sudo apt-get install gfortran
brew tap homebrew/science
brew install python3
brew install gubbins
```

## Linux - Ubuntu Xenial (16.04) & Debian (unstable)
Gubbins has been packaged by the Debian Med team and is trivial to install using apt.
```
sudo apt-get install gubbins
```

## Linux - CentOS/RHEL 7
Enable EPEL and make sure compilers are installed.
```
sudo yum install epel-release gcc gcc-c++ automake
```
Install [LinuxBrew](http://brew.sh/linuxbrew/).
```
brew tap homebrew/science
brew tap homebrew/dupes	
ln -s $(which gcc) ~/.linuxbrew/bin/gcc-4.8
ln -s $(which g++) ~/.linuxbrew/bin/g++-4.8
ln -s $(which gfortran) ~/.linuxbrew/bin/gfortran-4.8
brew install ruby gpatch python3
brew install gubbins
```

## Linux - CentOS/RHEL 6.6
Enable EPEL and make sure compilers are installed.
```
sudo yum install epel-release gcc gcc-c++ automake ruby-irb
```
Install [LinuxBrew](http://brew.sh/linuxbrew/).
```
brew tap homebrew/science
ln -s $(which gcc) ~/.linuxbrew/bin/gcc-4.4
ln -s $(which g++) ~/.linuxbrew/bin/g++-4.4
ln -s $(which gfortran) ~/.linuxbrew/bin/gfortran-4.4
brew install ruby python3
brew install gubbins
```

## OSX/Linux - from source
This is the most difficult method and is only suitable for someone with advanced computing skills. Please consider using HomeBrew/LinuxBrew instead.

Install the dependances and include them in your PATH:
* [FastTree](http://www.microbesonline.org/fasttree/#Install) ( >=2.1.4 )
* [RAxML](https://github.com/stamatak/standard-RAxML) ( >=8.0 )
* [FASTML](http://fastml.tau.ac.il/source.php) ( >=2.02 )
* Python modules: Biopython (> 1.59), DendroPy (>=4.0), Reportlab, nose, pillow
* Standard build environment tools (e.g. python3, pip3, make, autoconf, libtool, gcc, check, etc...)

```
autoreconf -i
./configure
make
sudo make install
```

## OSX/Linux/Windows - Virtual Machine
Roary wont run natively on Windows but we have created virtual machine which has all of the software setup, along with the test datasets from the paper. 
It is based on [Bio-Linux 8](http://environmentalomics.org/bio-linux/).  You need to first install [VirtualBox](https://www.virtualbox.org/), 
then load the virtual machine, using the 'File -> Import Appliance' menu option. The root password is 'manager'.

* ftp://ftp.sanger.ac.uk/pub/pathogens/pathogens-vm/pathogens-vm.latest.ova

More importantly though, if your trying to do bioinformatics on Windows, your not going to get very far and you should seriously consider upgrading to Linux.

