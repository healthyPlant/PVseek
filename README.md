# PVseek: a lightweight cross-platform plant viral diagnostic pipeline for high throughput sequencing data

PVseek is an open-source bioinformatics pipeline for plant virus detection using virous high throughput sequencing data, including Illumina RNA-seq, small RNA-seq, AmpliSeq/HiPlex, and Oxford Nanopore data. It's a fast analysis tool for virus diagnosis using a regular computer/laptop. The pipeline is written in [Snakemake](https://snakemake.readthedocs.io) (Köster and Rahmann 2018), a workflow management system for the development of data analysis workflows. PVseek consists of three major stages: (1). read QC; (2). mapping reads to reference sequences; (3). counting mapped reads for viruses and calculating their genome coverages. PVseek is developed by the USDA APHIS Plant Germplasm Quarantine Program (PGQP).  

# Workflow

![scheme of workflow](doc/*VDflow_scheme*.png?raw=true)

## Tools in pipeline

1. Read QC by [`Fastp`](https://github.com/OpenGene/fastp)
2. Mapping tools
   1. Illumina RNA-seq: [`BWA`](https://github.com/lh3/bwa) *||* [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
   2. Illumina small RNA-seq: [`Bowtie`](https://bowtie-bio.sourceforge.net/index.shtml)
   3. Oxford Nanopore: [`minimap2`](https://github.com/lh3/minimap2)
   4. AmpliSeq/HiPlex: [`BWA`](https://github.com/lh3/bwa)
3. Alignment maniapulation tools
   1. Sort and index alignments by [`SAMtools`](https://github.com/samtools/samtools)
   2. Consensus callers by [`BCFTools`](http://samtools.github.io/bcftools/bcftools.html) 
   3. Genome coverage by [`BEDTools`](https://github.com/arq5x/bedtools2)
4. Sequence processing tool: [`seqtk`](https://github.com/lh3/seqtk)

## Quick start
### Installation

For a Windows system, please use [PVseek docker image](https://hub.docker.com/r/healthyplant/pvseek).

For a Linux or Mac system, you can download our code or clone the repository using Git:
```
cd /path/to/software
git clone https://github.com/healthyPlant/PVseek.git
```

Then install dependencies. For an Ubuntu system, simply run
```
sudo bash /path/to/PVseek/scripts/installTools.sh /path/to/software
```
For other Linux or Mac systems, please follow [PVseek wiki](https://github.com/healthyPlant/PVseek/wiki#dependencies) to install dependencies.

***Since some tools conflict in a conda environment and are hard to update, we do not recommend installing all dependencies using conda.*** For a conda user, please check [PVseek wiki](https://github.com/healthyPlant/PVseek/wiki#dependencies) to see how to create a PVseek conda environment.

### Set up configuration
First, pick up the right config file for your sequencing platform.
1. Illumina paired-end read RNA-seq: `config_pe.yaml`
2. Illumina single read RNA-seq: `config_se.yaml`
3. Illumina small RNA-seq: `config_sRNA.yaml`
4. Illumina AmpliSeq/HiPlex: `config_hiplex.yaml`
5. Oxford Nanopore: `config_Nanopore.yaml`

Second, customize the workflow based on your needs in the config file, such as your sequence file extension and database paths. Please see the details in [wiki](https://github.com/healthyPlant/PVseek/wiki).

### Run PVseek
Please check [dependencies requirements](https://github.com/healthyPlant/PVseek/wiki) first using a dry-run (-n flag). 

For fastq.gz reads input dry-run:
```shell
$ snakemake  --configfile /path/to/PVseek/config_*.yaml -s /path/to/PVseek/Snakefile --config workDir=/path/to/output/folder fastqDir=/path/to/input/fastq/folder --cores [number of cores ex. 8] -n 
```

If the dry-run succeeds, please remove '-n' parameter to run the pipeline. If you'd like to run it in the background, please use 'nohup'. For example:
```shell
$ nohup snakemake  --configfile /path/to/PVseek/config.yaml -s /path/to/PVseek/Snakefile --config workDir=/path/to/output/folder fastqDir=/path/to/input/fastq/folder --cores [number of cores ex. 8] &
```
**Important:** For workDir and fastqDir paths, full paths must be used.

You can view progress or errors in the file 'nohup.out' using the command

`more nohup.out`

### PVseek quick test
After the software are ready, you can run a quick test using [the VIROMOCKchallenge Dataset8](https://gitlab.com/ilvo/VIROMOCKchallenge/-/blob/master/Datasets/Dataset8.md), which is in the test/data folder. The steps can be found in [test notebook](test/quick_test.ipynb).  

# PVseek docker image
[The PVseek docker image](https://hub.docker.com/r/healthyplant/pvseek) can be pulled 
```
docker pull healthyplant/pvseek
```
Docker can help you avoid manually installing the software. You can use the docker image on many systems (Linux, Mac, Windows). PVseek docker image usage is in its [docker README](https://hub.docker.com/r/healthyplant/pvseek). PVseek docker tutorial can be found in [test docker](test/pvseek_docker_tutorial.ipynb).

# Documentation

More information on input/output, dependencies, and databases are found in the [wiki](https://github.com/healthyPlant/PVseek/wiki)

# Citation

> The manuscript is under review.
> Will update later. 

Contact
------------
Alex Hu (xiaojun.hu (at) usda.gov)
