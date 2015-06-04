#!/bin/bash

set -x
set -e

start_dir=$(pwd)

RAXML_VERSION="8.1.20"
FASTML_VERSION="2.3"
FASTTREE_VERSION="2.1.8"

RAXML_DOWNLOAD_URL="https://github.com/stamatak/standard-RAxML/archive/v${RAXML_VERSION}.tar.gz"
FASTML_DOWNLOAD_URL="https://github.com/sanger-pathogens/fastml/archive/v${FASTML_VERSION}.tar.gz"
FASTTREE_DOWNLOAD_URL="http://www.microbesonline.org/fasttree/FastTree-${FASTTREE_VERSION}.c"

# Make an install location
if [ ! -d 'build' ]; then
  mkdir build
fi
cd build
build_dir=$(pwd)

# DOWNLOAD ALL THE THINGS
download () {
  url=$1
  download_location=$2

  if [ -e $download_location ]; then
    echo "Skipping download of $url, $download_location already exists"
  else
    echo "Downloading $url to $download_location"
    wget $url -O $download_location
  fi
}

download $RAXML_DOWNLOAD_URL "raxml-${RAXML_VERSION}.tgz"
download $FASTML_DOWNLOAD_URL "fastml-${FASTML_VERSION}.tgz"
download $FASTTREE_DOWNLOAD_URL "fasttree-${FASTTREE_VERSION}.c"

# Update dependencies
if [ "$TRAVIS" = 'true' ]; then
  echo "Using Travis's apt plugin"
else
  sudo apt-get update -q
  sudo apt-get install -y -q autoconf \
                             check \
                             g++ \
                             libtool \
                             pkg-config \
                             python-dev
fi

# Build all the things
cd $build_dir

## RAxML
raxml_dir=$(pwd)/"standard-RAxML-${RAXML_VERSION}"
if [ ! -d $raxml_dir ]; then
  tar xzf raxml-${RAXML_VERSION}.tgz
fi
cd $raxml_dir
if [ -e "${raxml_dir}/raxmlHPC" ]; then
  echo "Already build RAxML; skipping build"
else
  make -f Makefile.gcc
fi

cd $build_dir

## FASTML
fastml_dir=$(pwd)/"fastml-${FASTML_VERSION}"
if [ ! -d $fastml_dir ]; then
  tar xzf fastml-${FASTML_VERSION}.tgz
fi
cd $fastml_dir
if [ -e "${fastml_dir}/programs/fastml/fastml" ]; then
  echo "Already build FASTML; skipping build"
else
  make
fi

cd $build_dir

## FastTree
fasttree_dir=${build_dir}/fasttree-${FASTTREE_VERSION}
if [ ! -d $fasttree_dir ]; then
  mkdir $fasttree_dir
fi
cd $fasttree_dir
if [ -e "${fasttree_dir}/FastTree" ]; then
  echo "Skipping, FastTree already exists"
else
  gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree ${build_dir}/fasttree-${FASTTREE_VERSION}.c  -lm
fi

# Setup environment variables
update_path () {
  new_dir=$1
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path ${raxml_dir}
update_path ${fastml_dir}/programs/fastml
update_path ${fasttree_dir}

cd $start_dir

set +x
set +e
