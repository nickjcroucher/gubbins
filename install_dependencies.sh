#!/bin/bash

set -x
set -e

start_dir=$(pwd)

# Determine OS type
OS="Linux"
if [[ $OSTYPE == "darwin"* ]]; then
    OS="OSX"
fi

# Alias wget on OSX
if [[ $OS == "OSX" ]]; then
    function _wget() { curl "${1}" -o $(basename "${1}") ; };
    alias wget='_wget'
fi

# Get tree builder options
RAXML_VERSION="8.2.12"
FASTTREE_VERSION="2.1.11"
IQTREE_VERSION="2.0.3"
RAXMLNG_VERSION="1.0.1"
RAPIDNJ_VERSION="2.3.2"

RAXML_DIR="standard-RAxML-$RAXML_VERSION"
RAXML_ZIP_FILE="$RAXML_DIR.tar.gz"

IQTREE_DIR=""
IQTREE_ZIP_FILE=""
if [[ $OS == "Linux" ]]; then
    IQTREE_DIR="iqtree-$IQTREE_VERSION-Linux"
    IQTREE_ZIP_FILE="$IQTREE_DIR.tar.gz"
elif [[ $OS == "OSX" ]]; then
    IQTREE_DIR="iqtree-$IQTREE_VERSION-MacOSX"
    IQTREE_ZIP_FILE="$IQTREE_DIR.zip"
fi

FASTTREE_DIR="FastTree-$FASTTREE_VERSION"
FASTTREE_SOURCE="$FASTTREE_DIR.c"

RAXMLNG_DIR="raxmlng_dir"
RAXMLNG_ZIP_FILE=""
if [[ $OS == "Linux" ]]; then
    RAXMLNG_ZIP_FILE="raxml-ng_v${RAXMLNG_VERSION}_linux_x86_64.zip"
elif [[ $OS == "OSX" ]]; then
    RAXMLNG_ZIP_FILE="raxml-ng_v${RAXMLNG_VERSION}_macos_x86_64.zip"
fi

RAPIDNJ_DIR="rapidnj-$RAPIDNJ_VERSION"
RAPIDNJ_ZIP_FILE="$RAPIDNJ_VERSION.zip"

RAXML_DOWNLOAD_URL="https://github.com/stamatak/standard-RAxML/archive/v$RAXML_VERSION.tar.gz"
FASTTREE_DOWNLOAD_URL="http://www.microbesonline.org/fasttree/$FASTTREE_SOURCE"
IQTREE_DOWNLOAD_URL="https://github.com/Cibiv/IQ-TREE/releases/download/v$IQTREE_VERSION/$IQTREE_ZIP_FILE"
RAXMLNG_DOWNLOAD_URL="https://github.com/amkozlov/raxml-ng/releases/download/$RAXMLNG_VERSION/$RAXMLNG_ZIP_FILE"
RAPIDNJ_DOWNLOAD_URL="https://github.com/johnlees/rapidnj/archive/$RAPIDNJ_ZIP_FILE"

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
download $RAXMLNG_DOWNLOAD_URL $RAXMLNG_ZIP_FILE
download $RAPIDNJ_DOWNLOAD_URL $RAPIDNJ_ZIP_FILE

# Update dependencies
if [[ "$OS" == 'Linux' ]]; then
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
   if [[ $OS == "Linux" ]]; then
      tar xzf $IQTREE_ZIP_FILE
   elif [[ $OS == "OSX" ]]; then
      unzip $IQTREE_ZIP_FILE
   fi
fi
cd $IQTREE_DIR
if [ -e "bin/iqtree" ]; then
  cp bin/iqtree iqtree
fi

## RAxML-NG
cd $build_dir
if [ ! -d $RAXMLNG_DIR ]; then
  mkdir $RAXMLNG_DIR
fi
unzip $RAXMLNG_ZIP_FILE
mv raxml-ng $RAXMLNG_DIR

## RapidNJ
cd $build_dir
if [ ! -d $RAPIDNJ_DIR ]; then
  unzip $RAPIDNJ_ZIP_FILE
fi
cd $RAPIDNJ_DIR
make
if [ -e "bin/rapidnj" ]; then
  cp bin/rapidnj rapidnj
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
update_path $build_dir/$RAXMLNG_DIR
update_path $build_dir/$RAPIDNJ_DIR

cd $start_dir

set +x
set +e
