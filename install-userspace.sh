#!/bin/bash

# This setup script installs gubbins and dependencies on Ubuntu 14.04 x86_64.
# In particular, we set things up in the users' home directory in ~/.local (as
# per the XDG spec).  We also add ~/.local/bin to PATH and ~/.local/lib to
# LD_LIBRARY_PATH.
#
# If you'd prefer to install this system wide, have a look at the documentation
# at https://github.com/AidanDelaney/gubbins
#
# Contact Aidan Delaney <aidan@ontologyengineering.org> for issues with this
# script.
#

py_pkgs=( "biopython" "dendropy" )
deb_urls=( "http://uk.archive.ubuntu.com/ubuntu/pool/universe/r/raxml/raxml_7.2.8-2_amd64.deb" "https://launchpad.net/~ap13/+archive/ubuntu/gubbins/+files/fastml2_2.3~trusty1_amd64.deb" "https://launchpad.net/~ap13/+archive/ubuntu/gubbins/+files/gubbins_1.1.3~trusty1_amd64.deb" )

function check_platform {
    # Ubuntu 14.04
    echo -n "Checking Platform "
    source /etc/lsb-release
    case $DISTRIB_DESCRIPTION in
        "Ubuntu 14.04")
            echo "pass"
            ;;
        *)
            echo "This setup script has only been tested on Ubuntu 14.04 but will continue with the installation process."
            ;;
    esac

    case "$(uname -m)" in
        "x86_64")
            echo "pass"
            ;;
        *)
            echo "This setup script only works on Ubuntu x86_64"
            exit
            ;;
    esac
}

function check_dependencies {
    # Python 2.7
    echo -n "Checking Python 2.7 "
    case "$(python --version 2>&1)" in
    *" 2.7"*)
        echo "pass"
        ;;
    *)
        echo "Wrong Python version! Please ensure you are using Python 2.7"
        exit 1
        ;;
    esac

    echo -n "Checking Python development headers "
    # ensure python-dev is installed
    dpkg -s python-dev > /dev/null 2>&1
    if [ $? -eq 0 ] ; then
        echo "pass"
    else
        echo "No Python.h on this system. Please ensure python-dev is installed."
        exit 1
    fi

    # wget
    echo -n "Checking wget "
    if hash wget 2>/dev/null; then
        echo "pass"
    else
        echo "No wget found"
        exit 1
    fi
}

function install {
    # install Python setuptools
    wget https://bootstrap.pypa.io/ez_setup.py -O - | python - --user

    for pkg in ${py_pkgs[@]}
    do
        ${HOME}/.local/bin/easy_install -U --user $pkg
    done

    for deb in ${deb_urls[@]}
    do
        wget $deb
    done

    # extract all debs to a subfolder
    mkdir expanded

    # expand each deb
    for deb in *.deb
    do
        dpkg -x ${deb} expanded/.
    done

    # copy binaries to ~/.local/bin
    cp expanded/usr/bin/* ~/.local/bin
    cp expanded/usr/lib/x86_64-linux-gnu/libgubbins.* ~/.local/lib
    cp -R expanded/usr/lib/python2.7/ ~/.local/lib

    rm -rf expanded
}

function verify_path {
    # tell the user session to pick up the installed tools
    cat >> ~/.bashrc <<EOF
LD_LIBRARY_PATH=~/.local/lib${LD_LIBRARY_PATH:+:}${LD_LIBRARY_PATH:-}
export LD_LIBRARY_PATH
PATH=~/.local/bin${PATH:+:}${PATH:-}
export PATH
EOF
    echo "" >> ~/.bashrc
    echo "run 'source ~/.bashrc' to ensure you can find gubbins"
}

check_platform
check_dependencies
install
verify_path
