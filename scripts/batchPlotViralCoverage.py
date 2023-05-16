#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 04/18/2023
#This program draws all amplicon coverage graphs for detected viruses in all samples and calculate median coverage of the amplicons
#
#Required Parameters:
#   -w, --workdir X........................full path to the PVseek work directory, which has folders: report and trimmed
#   -d, --database X.......................plant virus database for PVseek, a fasta file
#   -s, --script X.........................PVseek script folder
#   -m, --mapping X........................mapping tool: bwa for RNA-seq, bowtie for smallRNA, bowtie2 for RNA-seq, minimap2 for Nanorpore
#   -o, --output X.........................path to output folder 
######################################################################

import argparse
import subprocess
import os
import glob


#Sample  TaxonId ReadInTaxon     Reference       RefCoverage     Virus   TaxonPath       RefDescription
#25_C1   42882   117078  KY510860.1      0.11    Cherry virus A  Viruses; Riboviria; Orthornavirae; Kitrinoviricota; Alsuviricetes; Tymovirales; Betaflexiviridae; Trivirinae; Capillovirus; Cherry virus A;     KY510860.1 Cherry virus A isolate 1286 1A2/13C222_N9, complete genome
def getReport(reportFile):
    """
    get detected virus reference id in samples
    """
    sampleVirusDict = {} #viral taxon id : [[reference, primer_f, primer_R]]
    with open(reportFile, 'r') as fin:
        next(fin) #skip header
        for line in fin:
            cells = line.strip().split("\t")
            sample = cells[0]
            refId = cells[3] 
            if sample in sampleVirusDict:
                sampleVirusDict[sample].append(refId)
            else:
                sampleVirusDict[sample] = [refId]
    return sampleVirusDict

def batchDrawPlot(reportFile, virusDbFasta, scriptFolder, workFolder, outputFolder, mappingTool):
    """
    draw all coverage graphs from a report
    """
    sampleVirusDict = getReport(reportFile)
    #print(sampleVirusDict)

    #subprocess.Popen() is non-blocking, and returns immediately. The shell command is run in the background.
    #extract reference sequences from the database fasta file and map primers to them
    print("\nMapping reads to references. Please wait ...")    
    for sample in sampleVirusDict:
        readFiles = glob.glob(workFolder + "/trimmed/" + sample + "*.trimmed.fastq.gz")
        reads = " ".join(readFiles)
        print("Read files:", readFiles)
        for refName in sampleVirusDict[sample]:
            print("Plotting ", refName, " in ", sample)
            command = ['bash', scriptFolder + '/mapReadToRef.sh', scriptFolder, sample, refName, virusDbFasta, mappingTool, reads ] 
            print(" ".join(command))
            pipe = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) #, shell=True
            output, error = pipe.communicate()
            print(output)
            #output = output.decode("utf-8").strip()
            #print(pipe.stderr.readline())
            #print(pipe.stdout.readline())
            if pipe.returncode != 0: 
                print("plotAmplicons.py %d %s %s" % (pipe.returncode, output, error))
            pipe.wait()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Draw viral coverage graph.")
    parser.add_argument('-w','--workdir',dest='workdir', required=True, default="", help='full path to the PVseek work directory, which has folders: report and trimmed')
    parser.add_argument('-d','--database',dest='database', required=True, default="", help='plant virus database for PVseek, a fasta file')
    parser.add_argument('-s','--script',dest='script', required=True, default="", help='PVseek script folder')
    parser.add_argument('-m','--mapping',dest='mapping', required=True, default="", help='mapping tool: bwa for RNA-seq, bowtie for smallRNA, bowtie2 for RNA-seq, minimap2 for Nanorpore')
    parser.add_argument('-o','--output',dest='outputFolder', required=True, default="", help='path to folder having amplicon graphs')

    args = parser.parse_args()
    workFolder = args.workdir #"/ppq/data0/test_PVseek/1_A1"
    virusDbFasta = args.database #"/ppq/data0/software/PVseek/db/plantvirus.fa"
    scriptFolder = args.script #"/ppq/data0/software/PVseek/scripts"
    mappingTool = args.mapping #bwa, bowtie, bowtie2, minimap2
    outputFolder = args.outputFolder

    if not (os.path.exists(workFolder + "/report") or os.path.exists(workFolder + "/trimmed")):
        print(workFolder, " does not have folders report or/and trimmed. Exit!")
        exit()

    if not os.path.isfile(virusDbFasta):
        print(virusDbFasta, "does not exist. Exit!")
        exit()
    if not os.path.exists(scriptFolder):
        print(scriptFolder, "does not exist. Exit!")
        exit()        

    reportFile = workFolder + "/report/report.txt"
    if not os.path.isfile(reportFile):
        print(reportFile, "does not exist. Exit!")
        exit()

    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder, exist_ok=True)

    #change working directory to outputFolder
    os.chdir(outputFolder)
    #print("working directory", os.getcwd())

    #draw all amplicon plots
    batchDrawPlot(reportFile, virusDbFasta, scriptFolder, workFolder, outputFolder, mappingTool)
    print("\n****************************************************")
    print("All coverage grpahs are plotted")
    print("Please check output files: *.consensus.N.fasta, *.coverage.txt and *.coverage.png in the output foler ", outputFolder)

    #clean files
    command = "rm -rf " + outputFolder + "/*.bam*"
    os.system(command)
    command = "rm -rf " + outputFolder + "/*.fasta.*"
    os.system(command)
