#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 04/18/2023
##This program gets primer positions in the reference
#
#Required Parameters:
#   -n, --reference X......................reference name from the PVseek report
#   -d, --database X.......................plant virus database for PVseek, a fasta file
#   -f, --forward X........................forward primer fasta file
#   -r, --reverse X........................reverse primer fasta file
#   -o, --output X.........................output folder
######################################################################

import argparse
import os


def getRefseq(outFolder, refName, db):
    """
    extract a reference sequence from the database fasta file using its id 
    """
    output = outFolder + "/" + refName
    #if not os.path.exists(output + ".fasta"):
    #use seqtk subseq to get reference sequence
    command = "echo " + refName + " > " + output
    #print(command)
    os.system(command)
    command = "seqtk subseq " + db + " " + output + " > " + output + ".fasta"
    #print(command)
    os.system(command)
        
    if os.stat(output + ".fasta").st_size == 0:
        print('Exits because', refName, 'does not exist in the database', db)
        exit()

def mapPrimer2Ref(outFolder, refName, fPrimerFile, rPrimerFile):
    """
    map primer sequences to the reference, then output their positions in a bed file using bedtools bamtobed
    """
    refSeq = outFolder + "/" + refName + ".fasta"
    #if not os.path.exists(refSeq + ".sa"):
    command = "bwa index " + refSeq
    os.system(command)
    bwa_param="-l 10 -n 8 -k 4 -M 2 -O 5 -E 1"

    fName = os.path.basename(fPrimerFile) #get basename of the file without path
    rName = os.path.basename(rPrimerFile) #get basename of the file without path

    fName = os.path.splitext(fName)[0] #get name without extension
    rName = os.path.splitext(rName)[0]

    fPrimerBamFile =  outFolder + "/" + refName + "." + fName + ".bam" #forward primer bam file
    rPrimerBamFile =  outFolder + "/" + refName + "." + rName + ".bam" #reverse primer bam file

    fPrimerBedFile =  outFolder + "/" + refName + "." + fName + ".bed" #forward primer bed file
    rPrimerBedFile =  outFolder + "/" + refName + "." + rName + ".bed" #reverse primer bed file

    #if not (os.path.exists(fPrimerBedFile) or os.path.exists(rPrimerBedFile)) or os.stat(fPrimerBedFile).st_size == 0 or os.stat(rPrimerBedFile).st_size == 0:
    #map forward primers to the reference
    command = "bwa aln " + bwa_param + " " + refSeq + " " + fPrimerFile + " | bwa samse " + refSeq + " - " + fPrimerFile + " | samtools view -b - > " + fPrimerBamFile
    #print(command)
    os.system(command)
    #convert bam to bed
    command = "bedtools bamtobed -bed12 -i " + fPrimerBamFile + " > " +  fPrimerBedFile
    #print(command)
    os.system(command)

    #map reverse primers to the reference
    command = "bwa aln " + bwa_param + " " + refSeq + " " + rPrimerFile + " | bwa samse " + refSeq + " - " + rPrimerFile + " | samtools view -b - > " + rPrimerBamFile
    #print(command)
    os.system(command)
    #convert bam to bed
    command = "bedtools bamtobed -bed12 -i " + rPrimerBamFile + " > " + rPrimerBedFile
    #print(command)
    os.system(command)

    #return fPrimerBedFile, rPrimerBedFile

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Get primers positions in a reference sequence")
    parser.add_argument('-n','--reference',dest='reference_name', required=True, default="", help='reference name from the PVseek report')
    parser.add_argument('-d','--database',dest='database', required=True, default="", help='plant virus database for PVseek, a fasta file')
    parser.add_argument('-f','--forward',dest='forward_primer_file', required=True, default="", help='forward primer fasta file')
    parser.add_argument('-r','--reverse',dest='reverse_primer_file', required=True, default="", help='reverse primer fasta file')
    parser.add_argument('-o','--output',dest='output_folder', required=True, default="", help='path to amplicon graphs')
    parser.set_defaults(append=False)
    return parser.parse_args()

def main():
    ### Input arguments
    options = parseArguments()
    refName = options.reference_name
    db = options.database

    fPrimerFile = options.forward_primer_file
    rPrimerFile = options.reverse_primer_file    
    outFolder = options.output_folder

    #if outout folder is not exists, create one
    if not os.path.isdir( outFolder ):
        os.makedirs( outFolder )

    #extract reference sequences
    getRefseq(outFolder, refName, db)

    #get primer positions in the reference
    mapPrimer2Ref(outFolder, refName, fPrimerFile, rPrimerFile)

if __name__ == '__main__':
	main()
