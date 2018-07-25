# Use a LTS Ubuntu version as parent image
FROM ubuntu:16.04
ARG DEBIAN_FRONTEND=noninteractive
ARG PYTHONWARNINGS="ignore:Unverified HTTPS request"
ENV LD_LIBRARY_PATH='/usr/local/lib'
WORKDIR /opt

# Install general dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
	curl \
	make \
	automake \
	gcc \
	pkg-config \
	zlib1g-dev \
	autoconf \
	check \
	libtool \
	libsubunit-dev \
	git \
    python3 \
    python3-dev \
    python3-setuptools \
    python3-pip

# Install python dependencies 
RUN pip3 install --trusted-host pypi.python.org --upgrade pip
RUN pip3 install certifi \
  && pip3 install wheel \
  && pip3 install nose \
  && pip3 install pillow \
  && pip3 install dendropy \
  && pip3 install biopython
  
# Install RAxML
ARG raxml_version='8.2.12'
RUN curl -L https://github.com/stamatak/standard-RAxML/archive/v${raxml_version}.tar.gz -o standard-RAxML-${raxml_version}.tar.gz \
  && tar xzf standard-RAxML-${raxml_version}.tar.gz \
  && cd standard-RAxML-${raxml_version} \
  && make -f Makefile.gcc \
  && make -f Makefile.PTHREADS.gcc \
  && cp raxmlHPC /usr/local/bin/ \
  && cp raxmlHPC-PTHREADS /usr/local/bin/ \
  && cd .. \
  && rm -rf standard-RAxML-${raxml_version}

# Install FastTree
ARG fasttree_version='2.1.10'
RUN curl http://www.microbesonline.org/fasttree/FastTree-${fasttree_version}.c -o FastTree.c \
  && gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm \
  && mv FastTree /usr/local/bin/ \
  && rm FastTree.c
  
# Install IQTree
ARG iqtree_version='1.6.6'
RUN curl -L https://github.com/Cibiv/IQ-TREE/releases/download/v${iqtree_version}/iqtree-${iqtree_version}-Linux.tar.gz -o iqtree-${iqtree_version}-Linux.tar.gz \
  && tar xzf iqtree-${iqtree_version}-Linux.tar.gz \
  && cp iqtree-${iqtree_version}-Linux/bin/iqtree /usr/local/bin \
  && rm -rf iqtree-${iqtree_version}-Linux

# Install Gubbins
RUN git clone https://github.com/sanger-pathogens/gubbins.git \
  && cd gubbins \
  && autoreconf -i \
  && ./configure \
  && make \
  && make check \
  && make install \
  && cd python \
  && python3 setup.py install
  
# Rename executable for backwards compatibility
RUN mv /usr/local/bin/run_gubbins.py /usr/local/bin/run_gubbins
