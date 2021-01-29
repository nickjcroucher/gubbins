# Use a LTS Ubuntu version as parent image
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
ARG PYTHONWARNINGS="ignore:Unverified HTTPS request"
ENV LD_LIBRARY_PATH='/usr/local/lib'
WORKDIR /opt

# Install general dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    curl \
    build-essential \
    automake \
    pkg-config \
    zlib1g-dev \
    unzip \
    autoconf \
    check \
    libtool \
    libsubunit-dev \
    python3.8 \
    python3.8-dev \
    python3-setuptools \
    python3-pip \
    gdb \
    valgrind \
    locales

# Set locales
RUN   sed -i -e 's/# \(en_GB\.UTF-8 .*\)/\1/' /etc/locale.gen && \
      touch /usr/share/locale/locale.alias && \
      locale-gen

ENV   LANG     en_GB.UTF-8
ENV   LANGUAGE en_GB:en
ENV   LC_ALL   en_GB.UTF-8

# Get tree builder versions
ARG raxml_version='8.2.12'
ARG fasttree_version='2.1.11'
ARG iqtree_version='2.0.3'
ARG raxmlng_version='1.0.1'
ARG rapidnj_version='2.3.2'

# Install RAxML
RUN curl -L https://github.com/stamatak/standard-RAxML/archive/v${raxml_version}.tar.gz -o standard-RAxML-${raxml_version}.tar.gz \
  && tar xzf standard-RAxML-${raxml_version}.tar.gz \
  && cd standard-RAxML-${raxml_version} \
  && make -f Makefile.gcc \
  && make -f Makefile.PTHREADS.gcc \
  && cp raxml* /usr/local/bin/ \
  && cd .. \
  && rm -rf standard-RAxML-${raxml_version} standard-RAxML-${raxml_version}.tar.gz

# Install FastTree
RUN curl http://www.microbesonline.org/fasttree/FastTree-${fasttree_version}.c -o FastTree.c \
  && gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm \
  && mv FastTree /usr/local/bin/ \
  && rm FastTree.c
  
# Install IQTree
RUN curl -L https://github.com/Cibiv/IQ-TREE/releases/download/v${iqtree_version}/iqtree-${iqtree_version}-Linux.tar.gz -o iqtree-${iqtree_version}-Linux.tar.gz \
  && tar xzf iqtree-${iqtree_version}-Linux.tar.gz \
  && cp iqtree-${iqtree_version}-Linux/bin/iqtree /usr/local/bin \
  && rm -rf iqtree-${iqtree_version}-Linux iqtree-${iqtree_version}-Linux.tar.gz

# Install RAxML-NG
RUN curl -LO https://github.com/amkozlov/raxml-ng/releases/download/${raxmlng_version}/raxml-ng_v${raxmlng_version}_linux_x86_64.zip \
  && unzip raxml-ng_v${raxmlng_version}_linux_x86_64.zip -d raxml-ng_v${raxmlng_version} \
  && cp raxml-ng_v${raxmlng_version}/raxml-ng /usr/local/bin \
  && rm -rf raxml-ng_v${raxmlng_version} raxml-ng_v${raxmlng_version}_linux_x86_64.zip

# Install RapidNJ
RUN curl -LO https://github.com/johnlees/rapidnj/archive/${rapidnj_version}.zip \
  && unzip ${rapidnj_version}.zip \
  && cd rapidnj-${rapidnj_version} \
  && make \
  && cp bin/rapidnj /usr/local/bin \
  && cd .. \
  && rm -rf rapidnj-${rapidnj_version} ${rapidnj_version}.zip

# Install Python dependencies
RUN pip3 install \
  pytest \
  pytest-cov \
  biopython==1.78 \
  multiprocess==0.70.11 \
  scipy==1.6.0 \
  numpy==1.19.5 \
  dendropy==4.5.1

# Install Gubbins
ENV BUILD_DIR /opt/gubbins
RUN mkdir -p ${BUILD_DIR}
COPY . ${BUILD_DIR}
RUN cd ${BUILD_DIR} \
  && autoreconf -i \
  && ./configure \
  && make \
  && make install \
  && make check
