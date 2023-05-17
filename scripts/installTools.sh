#/usr/bin/env bash
# This script is for installing PVseek required tools on a Ubuntu 20.04 machine

SRC=$1 #/opt

BIN=/usr/local/bin
sudo apt-get update  
sudo apt-get -qq -y install git curl wget unzip make apt-utils build-essential libssl-dev zlib1g-dev libbz2-dev liblzma-dev libz-dev libfreetype6-dev libpng-dev autogen libtool shtool 

#sudo apt-get -qq -y install pkg-config 
#Configuring tzdata
#Geographic area: 2
#Time zone: 106

sudo apt-get -qq -y install python3-dev python3-pip python3-versioneer
#rename python3
#sudo apt-get update && sudo apt-get -qq -y install python3
#sudo ln -sf /bin/python3 /bin/python

#install biopython
sudo pip install biopython

#install fastp
cd $SRC
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
sudo mv ./fastp $BIN/

#install snakemake
wget https://github.com/snakemake/snakemake/archive/refs/tags/v7.18.0.tar.gz 
tar -xvzf v7.18.0.tar.gz
cd snakemake-7.18.0
sudo python setup.py install
cd  ..

#install htslib
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xvjf htslib-1.16.tar.bz2
cd htslib-1.16 
autoheader 
autoconf -Wno-syntax 
./configure 
make
sudo make install
cd ..

#install samtools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xvjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
autoheader 
autoconf -Wno-syntax 
./configure --with-htslib=/usr/local
make 
sudo make install
cd ..

#install bcftools
wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
tar -xvjf bcftools-1.16.tar.bz2
cd bcftools-1.16
autoheader 
autoconf -Wno-syntax 
./configure --with-htslib=/usr/local
make 
sudo make install
cd ..

#install bedtools
mkdir bedtools2
cd bedtools2
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
mv bedtools.static.binary bedtools
chmod a+x bedtools
cd ..
sudo ln -s $SRC/bedtools2/bedtools $BIN/bedtools

#install seqtk
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
sudo cp $SRC/seqtk/seqtk $BIN/seqtk
cd ..

#instal NCBI EDirect
wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
gunzip -c edirect.tar.gz | tar xf -
cd edirect
./nquire -dwn ftp.ncbi.nlm.nih.gov entrez/entrezdirect xtract.Linux.gz
gunzip -f xtract.Linux.gz
chmod +x xtract.Linux
cd ..
sudo find $SRC/edirect -type f -executable -exec ln -s {} $BIN/ \;

#install Bowtie
wget https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/bowtie-1.3.1-linux-x86_64.zip
unzip bowtie-1.3.1-linux-x86_64.zip
mv bowtie-1.3.1-linux-x86_64 bowtie-1.3.1
sudo ln -s $SRC/bowtie-1.3.1/bowtie* $BIN/.

#install Bowtie2
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-sra-linux-x86_64.zip
unzip bowtie2-2.4.2-sra-linux-x86_64.zip
mv bowtie2-2.4.2-sra-linux-x86_64 bowtie2-2.4.2
sudo ln -s $SRC/bowtie2-2.4.2/bowtie2* $BIN/.

#Install BWA
git clone https://github.com/lh3/bwa.git
cd bwa; make
cd ..
sudo ln -s $SRC/bwa/bwa $BIN/bwa

#Install minimap2
wget https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2
tar -xvjf minimap2-2.26_x64-linux.tar.bz2
mv minimap2-2.26_x64-linux minimap2-2.26
sudo ln -s $SRC/minimap2-2.26/minimap2 $BIN/.

#install pandas
sudo pip install pandas
#install pysam
sudo pip install pysam
#install matplotlib
sudo pip install matplotlib

# some cleanup
rm -r *.tar.gz *.zip *.bz2
rm -rf $SRC/htslib-1.13
rm -rf $SRC/snakemake-7.18.0
rm -rf $SRC/seqtk
cd ..

sudo apt-get -qq -y autoremove
sudo apt-get clean 

#Check software versions
echo "Python3 version"
python --version
echo "FastQC version"
fastqc -v 
echo "snakemake version"
snakemake -v
echo "BWA version"
bwa
echo "Bowtie2 version"
bowtie2 --version
echo "Bcftools version"
bcftools --version
echo "samtools version"
samtools version
echo "bedtools version"
bedtools --version
echo "Seqtk version"
seqtk
echo "EDirect version"
esearch -h
