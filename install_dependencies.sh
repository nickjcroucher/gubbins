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

RAXML_DOWNLOAD_URL="https://github.com/stamatak/standard-RAxML/archive/v$RAXML_VERSION.tar.gz"
FASTTREE_DOWNLOAD_URL="http://www.microbesonline.org/fasttree/FastTree-$FASTTREE_VERSION.c"
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
download $FASTTREE_DOWNLOAD_URL "fasttree-${FASTTREE_VERSION}.c"
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
if [ -e "$RAXML_DIR/raxmlHPC" ]; then
  echo "Already build RAxML; skipping build"
else
  make -f Makefile.gcc
fi


## FastTree
cd $build_dir
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

## IQTree
cd $build_dir
if [ ! -d $IQTREE_DIR ]; then
  tar xzf $IQTREE_ZIP_FILE
fi
cd $IQTREE_DIR
if [ -e "$IQTREE_DIR/bin/iqtree" ]; then
  cp bin/iqtree iqtree
fi

# Setup environment variables
update_path () {
  new_dir=$1
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path $RAXML_DIR
update_path ${fasttree_dir}
update_path $IQTREE_DIR

cd $start_dir

set +x
set +e
