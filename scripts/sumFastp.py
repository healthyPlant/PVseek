#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 03/31/2023
#This program summarizes fastp json outputs of all samples: 

#Required Parameters:
#   -w, --workdir X........................working directory
#   -o, --output X.........................Report file name 
#######################################################################

import argparse
import glob
import json
import os

#################################### Get raw read information function ####################################  
def parseJson(jsonFile):
    """
    Get raw quality information from a jason file which generate by fastp
    """
    with open(jsonFile) as f:
        data = json.load(f) #json.load() output a dict

    sample = os.path.basename(jsonFile).replace(".fastp.json","")
    rawReads = "{:.0f}".format(data["summary"]["before_filtering"]["total_reads"])
    rawYield = "{:.2f}".format(data["summary"]["before_filtering"]["total_bases"] / 1000000) 
    rawQ20Rate = "{:.2f}".format(data["summary"]["before_filtering"]["q20_rate"])
    rawQ30Rate = "{:.2f}".format(data["summary"]["before_filtering"]["q30_rate"])
    rawAverageLength = "{:.2f}".format(data["summary"]["before_filtering"]["read1_mean_length"])
    trimReads = "{:.0f}".format(data["summary"]["after_filtering"]["total_reads"])
    trimYield = "{:.2f}".format(data["summary"]["after_filtering"]["total_bases"] / 1000000)
    trimQ20Rate = "{:.2f}".format(data["summary"]["after_filtering"]["q20_rate"])
    trimQ30Rate = "{:.2f}".format(data["summary"]["after_filtering"]["q30_rate"])
    trimAverageLength = "{:.2f}".format(data["summary"]["after_filtering"]["read1_mean_length"])   
    percentReadAfterTrim = "{:.2f}".format(float(trimReads) / float(rawReads))

    return [sample, rawReads, rawYield, rawQ20Rate, rawQ30Rate, rawAverageLength, trimReads, trimYield, trimQ20Rate, trimQ30Rate, trimAverageLength, percentReadAfterTrim]

#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description='Generate text report')
    parser.add_argument("-w", "--workdir", dest='workdir', required=True, help="Working directory, ex. /ppq/data2/pgqp_pipeline/Run45")
    parser.add_argument("-o", "--output", dest='output_file', required=True, help="output file name, ex. qcReadNumber.txt")
    parser.set_defaults(append=False)
    return parser.parse_args()


def main():
    ### Input arguments
    options = parseArguments()
    workdir = options.workdir  #'/ppq/data2/pgqp_pipeline/Run45/mapping'
    outFile = options.output_file #workdir + '/report/qcReadNumber.txt'

    jsonFiles = glob.glob(workdir + '/*.fastp.json')
    #sum up QC summary from json files
    summary = [] 
    sum1 = parseJson(jsonFiles[0])
    summary.append(sum1)
    for jf in jsonFiles[1:]:
        sum1 = parseJson(jf)
        #print(sum1)
        summary.append(sum1)

    header="Sample\tRawReads\tRawYield(Mbases)\tRawQ20Rate\tRawQ30Rate\tRawAverageLength\tTrimmedReads\tTrimmedYield(Mbases)\tTrimmedQ20Rate\tTrimmedQ30Rate\tTrimmedAverageLength\tRateReadAfterTrim"
    fout = open(outFile, 'w')
    fout.write(header + "\n")
    for sum in summary:
        fout.write("\t".join(sum) + "\n")
    fout.close

if __name__ == "__main__":
    main()

