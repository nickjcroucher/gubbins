#!/bin/bash

set -x
set -e

start_dir=$(pwd)

RAXML_VERSION="8.2.12"
FASTTREE_VERSION="2.1.10"
IQTREE_VERSION="1.6.6"

RAXML_DIR="standard-RAxML-$RAXML_VERSION"
RAXML_ZIP_FILE="$RAXML_DIR.tar.gz"

IQTREE_DIR="iqtree-$IQTREE_VERSION-Linux"
IQTREE_ZIP_FILE="$IQTREE_DIR.tar.gz"

FASTTREE_DIR="FastTree-$FASTTREE_VERSION"
FASTTREE_SOURCE="$FASTTREE_DIR.c"

RAXML_DOWNLOAD_URL="https://github.com/stamatak/standard-RAxML/archive/v$RAXML_VERSION.tar.gz"
FASTTREE_DOWNLOAD_URL="http://www.microbesonline.org/fasttree/$FASTTREE_SOURCE"
IQTREE_DOWNLOAD_URL="https://github.com/Cibiv/IQ-TREE/releases/download/v$IQTREE_VERSION/$IQTREE_ZIP_FILE"

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

download $RAXML_DOWNLOAD_URL $RAXML_ZIP_FILE
download $FASTTREE_DOWNLOAD_URL $FASTTREE_SOURCE
download $IQTREE_DOWNLOAD_URL $IQTREE_ZIP_FILE

# Update dependencies
if [ "$TRAVIS" = 'true' ]; then
  echo "Using Travis's apt plugin"
else
  sudo apt-get update -q
  sudo apt-get install -y -q autoconf \
                             check \
                             g++ \
                             libtool \
                             libsubunit-dev \
                             pkg-config \
                             python-dev
fi

# Build all the things

## RAxML
cd $build_dir
if [ ! -d $RAXML_DIR ]; then
  tar xzf $RAXML_ZIP_FILE
fi
cd $RAXML_DIR
if [ -e "raxmlHPC" ]; then
  echo "Already built single-processor RAxML; skipping build"
else
  make -f Makefile.gcc
fi
if [ -e "raxmlHPC-PTHREADS" ]; then
  echo "Already built multi-processor RAxML; skipping build"
else
  make -f Makefile.PTHREADS.gcc
fi

## FastTree
cd $build_dir
if [ ! -d $FASTTREE_DIR ]; then
  mkdir $FASTTREE_DIR
fi
cd $FASTTREE_DIR
if [ -e "$FASTTREE_DIR/FastTree" ]; then
  echo "Skipping, FastTree already exists"
else
  gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree $build_dir/$FASTTREE_SOURCE  -lm
fi

## IQTree
cd $build_dir
if [ ! -d $IQTREE_DIR ]; then
  tar xzf $IQTREE_ZIP_FILE
fi
cd $IQTREE_DIR
if [ -e "bin/iqtree" ]; then
  cp bin/iqtree iqtree
fi

# Setup environment variables
update_path () {
  new_dir=$1
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path $build_dir/$RAXML_DIR
update_path $build_dir/$FASTTREE_DIR
update_path $build_dir/$IQTREE_DIR

cd $start_dir

set +x
set +e
