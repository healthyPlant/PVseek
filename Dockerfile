FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update && apt-get -qq -y install \
    automake \
    build-essential \
    bzip2 \
    curl \
    wget \
    unzip \
    make \
    git \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libz-dev \
    libfreetype6-dev \
    libpng-dev \
    libcurl4-gnutls-dev \
    libdeflate-dev \
    libncurses-dev \
    autogen \
    libtool \
    shtool \
    libblas-dev \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-setuptools \
    python3-versioneer \
    python3-matplotlib \
    python3-distutils \
    rsync \
    texlive-latex-base \
    tzdata \
    apt-utils

## set up tool config and deployment area:
ENV SRC /opt
ENV BIN /usr/local/bin
ENV LD_LIBRARY_PATH=/usr/local/lib

#rename python3
RUN ln -sf /bin/python3 /bin/python

#apt-get install fastqc
WORKDIR $SRC
RUN wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp && \
    mv ./fastp $BIN/

#install snakemake
WORKDIR $SRC
RUN pip install Cython
RUN wget https://github.com/snakemake/snakemake/archive/refs/tags/v7.18.0.tar.gz && \
    tar -xvzf v7.18.0.tar.gz && \
    cd snakemake-7.18.0 && \
    python setup.py install

#install htslib
WORKDIR $SRC
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar -xvjf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    autoheader  && \
    autoconf -Wno-syntax  && \
    ./configure && \
    make && make install

#install samtools
WORKDIR $SRC
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar -xvjf samtools-1.16.1.tar.bz2 && \
    cd samtools-1.16.1 && \
    autoheader  && \
    autoconf -Wno-syntax  && \
    ./configure --with-htslib=/usr/local && \
    make  && \
    make install

WORKDIR $SRC
#install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 && \
    tar -xvjf bcftools-1.16.tar.bz2 && \
    cd bcftools-1.16 && \
    autoheader  && \
    autoconf -Wno-syntax  && \
    ./configure --with-htslib=/usr/local && \
    make  && \
    make install

WORKDIR $SRC
#install Bowtie
RUN wget https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/bowtie-1.3.1-linux-x86_64.zip && \
    unzip bowtie-1.3.1-linux-x86_64.zip && \
    mv bowtie-1.3.1-linux-x86_64 bowtie-1.3.1 && \
    ln -s $SRC/bowtie-1.3.1/bowtie* $BIN/. 

WORKDIR $SRC
#install Bowtie2
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip && \
    unzip bowtie2-2.4.2-sra-linux-x86_64.zip  && \
    mv bowtie2-2.4.2-sra-linux-x86_64 bowtie2-2.4.2  && \
    ln -s $SRC/bowtie2-2.4.2/bowtie2* $BIN/.

WORKDIR $SRC
#install bedtools
RUN mkdir bedtools2 && \
    cd bedtools2  && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary  && \
    mv bedtools.static.binary bedtools  && \
    chmod a+x bedtools  && \
    cd .. && \
    ln -s $SRC/bedtools2/bedtools $BIN/bedtools

WORKDIR $SRC
#install seqtk
#RUN apt-get -qq -y install seqtk
RUN git clone https://github.com/lh3/seqtk.git  && \
    cd seqtk  && \
    make  && \
    ln -s $SRC/seqtk/seqtk $BIN/seqtk 

WORKDIR $SRC
#Install BWA
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa; make  && \
    cd ..  && \
    ln -s $SRC/bwa/bwa $BIN/bwa

WORKDIR $SRC
#instal NCBI EDirect
RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz  && \
    gunzip -c edirect.tar.gz | tar xf -  && \
    cd edirect && \
    wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz &&\
    gunzip -f xtract.Linux.gz && \
    chmod +x xtract.Linux && \
    cd .. && \
    find $SRC/edirect -type f -executable -exec ln -s {} $BIN/ \;

WORKDIR $SRC
#instal minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 && \
    tar -xvjf minimap2-2.26_x64-linux.tar.bz2 && \
    mv minimap2-2.26_x64-linux minimap2-2.26 && \
    ln -s $SRC/minimap2-2.26/minimap2 $BIN/.

#install biopython
RUN pip install biopython

#install pandas
RUN pip install pandas

RUN apt-get -qq -y remove git && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN rm -rf $SRC/htslib-1.16
RUN rm -rf $SRC/samtools-1.16.1
RUN rm -rf $SRC/bcftools-1.16
RUN rm -rf $SRC/snakemake-7.18.0

# some cleanup
WORKDIR $SRC
RUN rm -r *.tar.gz *.zip *.bz2

# copy PVseek code
# pvseek folder contains the application
ADD PVseek $SRC/PVseek

################################################
## be sure this is last!
COPY Dockerfile $SRC/PVseek/.
