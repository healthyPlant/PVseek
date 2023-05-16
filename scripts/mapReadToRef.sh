#!/usr/bin/env bash
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 03/31/2023
#This program maps reads to the reference and draw coverage graph
#Usage: mapReadToRef.sh [script directory] [sample name] [reference name] [viral reference fasta file] [mapping tool] [read1] [read2]
######################################################################

scripts_dir=$1 #PVseek/scripts
sampleName=$2 #sample name
refName=$3  #reference name
viral_db=$4 #viral sequence db
mappingTool=$5  #bwa, bowtie, bowtie2, minimap2
read1=$6 #read1 (R1) file
read2=$7 #read1 (R2) file

#echo $viral_db
#echo $mappingTool
#echo $# arguments 
if [ "$#" -le 5 ]; then
    echo "Usage: bash mapReadToRef.sh [PVseek script directory] [sample name] [reference name] [viral reference database fasta file] [mapping tool] [read1] [read2]"
    echo "Example: bash mapReadToRef.sh /my/PVseek/scripts sampleA MT701608.1 /my/db/virusDetect_db/vrl_Plants_248_U95_annotated.fasta bwa /my/reads/sampleA_R1.fastq.gz /my/reads/sampleA_R2.fastq.gz"
    exit
fi

#1. get the reference sequence
echo $refName >name.list
seqtk subseq $viral_db name.list > $refName.fasta

#2. mapping read to a reference
#Please change mapping tools and their peremeters here
#mappingTool="bwa" # or "bowtie2", bowtie, minimap2

#mapping tool parameters
#bwa for Illumina single-read or paired-end
bwaOption=" -B 4 " #default bwa options. it's good for read length >=150, change mismatch panelty from 1 to 3 (-B)
#bwaOption=" -k 12 -A 1 -B 3 -O 1 -E 1 " # good for read length <=100, Viral-NGS use it.
#bwaOption=" -B 9 -O 16 " #BWA strict mapping
#Bowtie2 for Illumina single-read or paired-end
#bowtieOption=" --mp 20  --score-min L,-0.1,-0.1 " #read almost exactly match to the reference. It's good for mapping very similar virus genomes
bowtie2Option=" --very-sensitive-local -k 100 --score-min L,20,1.0 " #mapping reads with mismatches. Pathoscope2 used it. 
#minimap2 for Nanopore
minimap2Option="-L --secondary=no -ax map-ont"
#bowtie for small RNA-seq
bowtieOption="-n 1 -k 1 -M 100 --best --strata"

#readName=`basename $reads .fastq.gz`
outputName=$sampleName"."$refName

if [ $mappingTool == "bowtie2" ]; then
	bowtie2-build $refName.fasta $refName.fasta
	if [ -z "$read2" ];then
		bowtie2 -p 8 $bowtie2Option -x $refName.fasta -U  $read1 | samtools view -Sb -F 4 - | samtools sort - > $outputName.sorted.bam
	else
		bowtie2 -p 8 $bowtie2Option -x $refName.fasta -1 $read1 -2 $read2 | samtools view -Sb -F 4 - | samtools sort - > $outputName.sorted.bam
	fi
elif [ $mappingTool == "bwa" ]; then
	bwa index $refName.fasta
	#bwa mem -t {threads} {params.param} $ref {input.reads} | samtools sort -@{threads} -o {params.outputName}/$refName.bam - 2>> {log} 1>&2
	if [ -z "$read2" ];then
		bwa mem -t 8 $bwaOption  $refName.fasta $read1 | samtools view -Sb -F 4 - | samtools sort - > $outputName.sorted.bam  #filter out unmapped reads -F 4
	else
		bwa mem -t 8 $bwaOption $refName.fasta $read1 $read2 | samtools view -Sb -F 4 - | samtools sort - > $outputName.sorted.bam  #filter out unmapped reads -F 4
	fi
elif [ $mappingTool == "bowtie" ]; then
    bowtie-build $refName.fasta $refName.fasta
    bowtie -p 8 -q $bowtieOption -x $refName.fasta -S  $read1 | samtools view -Sb -F 4 | samtools sort - -o $outputName.sorted.bam 
elif [ $mappingTool == "minimap2" ]; then
	minimap2 -t 8 $minimap2Option $refName.fasta $read1 | samtools view -Sb -F 4 | samtools sort - -o $outputName.sorted.bam 
else
	echo "Error: invalid mapping tool. Must be bwa, bowtie, bowtie2, or minimap2"
	exit 1
fi

#3. calcualte coverage
#samtools sort -o $outputName.sorted.bam $outputName.bam
samtools depth -aa $outputName.sorted.bam > $outputName.coverage.txt
#get mapped read number
echo -n "mapped read number:"
samtools view -c -F 260 $outputName.sorted.bam
count0=$(grep -w 0$ $outputName.coverage.txt | wc -l)
total=$(cat $outputName.coverage.txt | wc -l)
echo -n "Genome covered %:"
#echo "($total-$count0)/$total*100" | bc -l
printf '%.2f%%\n' "$((($total-$count0)*100/$total))e-2"  #"e-2" get a float number with two decimals, 
echo -n "Reference lenght: "
wc -l  $outputName.coverage.txt
echo -n "Mean coverage: "
awk '{ total += $3; count++ } END { print total/count }' $outputName.coverage.txt

#4. plot coverage graph
echo "Plot coverage graph"
python $scripts_dir/plotViralCoverage.py -i $outputName.coverage.txt -o $outputName.coverage.png

#5. male consensus
echo "Make consensus"
bcftools mpileup -Ou -f $refName.fasta $outputName.sorted.bam | bcftools call -Ou -mv | bcftools norm -f $refName.fasta -Oz -o $outputName.vcf.gz
bcftools index $outputName.vcf.gz
bcftools consensus -f $refName.fasta $outputName.vcf.gz -o $outputName.consensus.fasta
bedtools genomecov -bga -ibam $outputName.sorted.bam -g $outputName.consensus.fasta > $outputName.bed
awk '{if($4 == 0) print}' $outputName.bed > $outputName.0.bed  
if [[ -s $outputName.0.bed  ]]; then
  bcftools consensus -m $outputName.0.bed -f $refName.fasta $outputName.vcf.gz -o $outputName.consensus.N.fasta
fi
