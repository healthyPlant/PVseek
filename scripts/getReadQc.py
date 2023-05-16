#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 03/31/2023
#This program get raw and trimmed read QC: 

#Required Parameters:
#   -w, --workdir X........................working directory
#   -m, --format X.........................read file format, ex. fastq.gz
#   -o, --output X.........................Report file name 
#Optional Parameters:
#   -f, --forward..........................Forward strand surffix in the fastq file, ex. R1
#   -f, --reverse..........................Reverse strand surffix in the fastq file, ex. R2
#######################################################################

import argparse
import glob
import os
import re
import multiprocessing 
import gzip

#################################### Get a raw fastq information function ####################################  
def getQualFromFastq(fastq_file):
    """
    Get raw quality information from a fastq.gz file 
    """

    #@ is the read start 
    # replace ASCII quality with integer
    #ord()	Converts a character to an integer
    #ord(str(quality)) - 33
    totalBases = 0
    totalQual = 0
    q30_count = 0
    totalLength = 0
    lineno = 0
    #print(fastq_file)
    #fin = open(fastq_file, "rt")
    #if fastq_file.endswith(".gz"):
    fin = gzip.open(fastq_file, "rt")
    for line in fin:
        lineno += 1
        #if lineno%4 == 1: flag = (line in ids)
        if lineno%4 == 0:
            line = line.strip() 
            qstr = list(line)
            for q in qstr:
                qual = ord(q) - 33
                if qual >= 30:
                    q30_count += 1
                totalQual += qual
            totalBases += len(line)
            totalLength += len(line)
    totalRead = int(lineno/4)
    avgQ30 = 100 * float(q30_count)/float(totalBases)
    avgQual = float(totalQual)/float(totalBases)
    avgLength = float(totalLength)/float(totalRead)
    sample=os.path.basename(fastq_file).split(".")[0]
    
    avgQualSet = (sample, str(totalRead), "{:.0f}".format(totalBases/1000000), "{:.2f}".format(avgQ30), "{:.2f}".format(avgQual), "{:.2f}".format(avgLength))
    return(avgQualSet)

#################################### Get all raw fastq information function ####################################  
def sumQual(fastqFolder, input_format, strand1, strand2):
    """
    Run parallel to get each fastq file quality
    """
    fastqFiles = [x for x in glob.glob(fastqFolder + '/*.' +  input_format)]
    #print(fastqFiles)
    #Parallel Processing multiple fastq.gz file 
    pool = multiprocessing.Pool() 
    outputs_async = pool.map_async(getQualFromFastq, fastqFiles) 
    outputs = outputs_async.get() #return a list of sets
    #print(outputs)
    #convert set to list
    for i in range(len(outputs)):
        outputs[i] = list(outputs[i])
    #remove _001 in sample name    
    #if strand1 and strand1.endswith("_R1_001"):
    if strand1:
        for output in outputs:
            if output[0].endswith("_001"):
                output[0] = output[0].replace("_001", "") 
               
    #check paired-end or single-end read
    paired = False
    for output in outputs:
        if strand2 and output[0].endswith(strand2):
            paired = True
            break

    #change a list to a dict
    avgQualDict = {}    
    #for paired-end reads
    avgQualDictR1 = {}
    avgQualDictR2 = {}

    if not strand1: #for single-end read without any extensions on sample names
         for output in outputs:
             avgQualDict[output[0]] = output[1:len(output)]
    else:
        for output in outputs:
            if output[0].endswith(strand1):
                sname = output[0].replace("_"+strand1, "")
                avgQualDictR1[sname] = output[1:len(output)]
            elif paired and output[0].endswith(strand2): #paired-end read
                sname = output[0].replace("_"+strand2, "")
                avgQualDictR2[sname] = output[1:len(output)]
            else: #single-end read
                avgQualDict[output[0]] = output[1:len(output)]

    if paired: #combine R1 and R2
        for sname in avgQualDictR2:
            avgQualDict[sname] = (avgQualDictR1[sname][0] + "(F)|" + avgQualDictR2[sname][0] + "(R)",  avgQualDictR1[sname][1] + "(F)|" + avgQualDictR2[sname][1] + "(R)", avgQualDictR1[sname][2] + "(F)|" + avgQualDictR2[sname][2] + "(R)", avgQualDictR1[sname][3] + "(F)|" + avgQualDictR2[sname][3] + "(R)", avgQualDictR1[sname][4] + "(F)|" + avgQualDictR2[sname][4] + "(R)"  )
    elif not avgQualDict:
        avgQualDict = avgQualDictR1

    return(avgQualDict)

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description='Generate text report')
    parser.add_argument("-w", "--workdir", dest='workdir', required=True, help="Working directory, ex. /ppq/data2/pgqp_pipeline/Run45")
    parser.add_argument('-f','--forward',dest='strand1', required=False, default="R1", help='Forward strand, ex. R1')
    parser.add_argument('-r','--reverse',dest='strand2', required=False, default="R2", help='Reverse strand, ex. R2')
    parser.add_argument("-m", "--format", dest='format', required=True, help="read file format, ex. fastq.gz")
    parser.add_argument("-o", "--output", dest='output_file', required=True, help="output file name, ex. qcReadNumber.txt")
    parser.set_defaults(append=False)
    return parser.parse_args()


def main():
    ### Input arguments
    options = parseArguments()
    workdir = options.workdir  #'/ppq/data2/pgqp_pipeline/Run45'
    strand1 = options.strand1  #R1
    strand2 = options.strand2  #R2
    input_format = options.format #"fastq.gz"
    outFile = options.output_file #workdir + '/report/qcReadNumber.txt'

    trimDir = workdir + '/trimmed'
    rawDir = workdir + '/raw'

    #samples = [os.path.basename(x) for x in glob.glob(trimDir + '/*.log')]
    #samples = [x.replace(".log","") for x in samples]
    #print(samples)

    header="Sample\tRawReads\tRawYield(Mbases)\tRawPercent>=Q30Bases\tRawMeanQualityScore\tRawAverageLength\tTrimmedReads\tTrimmedYield(Mbases)\tTrimmedPercent>=Q30Bases\tTrimmedMeanQualityScore\tTrimmedAverageLength\tPercentReadAfterTrim"

    fout = open(outFile, 'w')
    fout.write(header + "\n")

    rawQualDict = sumQual(rawDir, input_format, strand1, strand2)
    trimQualDict = sumQual(trimDir, input_format, strand1, strand2)
    #print(rawQualDict)
    #print(trimQualDict)
    for sample in rawQualDict:
        fout.write(sample + "\t")
        fout.write('\t'.join(rawQualDict[sample]))
        fout.write('\t')
        fout.write('\t'.join(trimQualDict[sample]))
        #print(trimQualDict[sample][0], rawQualDict[sample][0])
        percentTrim = float(trimQualDict[sample][0])/float(rawQualDict[sample][0]) * 100 
        if "(R)" in rawQualDict[sample][0]:
            numbers0 = re.search(r"(\d+)\(F\)\|(\d+)\(R\)", rawQualDict[sample][0])
            numbers1 = re.search(r"(\d+)\(F\)\|(\d+)\(R\)", trimQualDict[sample][0])
            percentTrim = (float(numbers0.group(1)) + float(numbers0.group(2))) / (float(numbers1.group(1)) + float(numbers1.group(2))) * 100
        fout.write('\t' + "{:.2f}".format(percentTrim)+'\n')   
    fout.close()  

if __name__ == "__main__":
    main()

