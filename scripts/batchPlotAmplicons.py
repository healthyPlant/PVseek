#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 04/18/2023
#This program draws all amplicon coverage graphs for detected viruses in all samples and calculate median coverage of the amplicons
#
#Required Parameters:
#   -w, --workdir X........................full path to the PVseek work directory, which has folders: report and trimmed
#   -v, --infor X..........................path to the virus information file, which has reference name, primer file names
#   -d, --database X.......................plant virus database for PVseek, a fasta file
#   -s, --script X.........................PVseek script folder
#   -p, --primer X ........................path to primer folder having all primer files
#   -j, --job X............................numbre of jobs of parallel running, ex. 16
#   -o, --output X.........................path to output folder 
######################################################################

import argparse
import subprocess
import os

#Name	Acronym	TaxonId	Reference	ForwardPrimerFile	ReversePrimerFile
#Apple chlorotic leaf spot virus	ACLSV	12175	NC_001409.1	ACLSV_primer_F.fa	ACLSV_primer_R.fa
#Tomato ringspot virus RNA1	TomRSV	12280	NC_003840.1	TomRSV-RNA1_primer_F.fa	TomRSV-RNA1_primer_R.fa
#Tomato ringspot virus RNA2	TomRSV	12280	NC_003839.2	TomRSV-RNA2_primer_F.fa	TomRSV-RNA2_primer_R.fa
def getviralInfor(viralInforFile):
    """
    get viral information: taxon id, reference, forward primer, reverse primer
    """
    virusDict = {} #viral taxon id : [[reference, primer_f, primer_r]]
    with open(viralInforFile, 'r') as fin:
        next(fin)
        for line in fin:
            cells = line.strip().split("\t")
            if len(cells) > 5: 
                taxonId = cells[2]
                refName = cells[3]
                primerF = cells[4]
                primerR = cells[5]
                infor = [refName, primerF, primerR]
                if taxonId in virusDict:
                    virusDict[taxonId].append(infor) #for multiple segment viruses, save viral information into a list of lists
                else:
                    virusDict[taxonId] = [infor]
            else:
                print("Error: reference or/and primer information are missed for the following line")
                print(line)
    return virusDict

#Sample  TaxonId ReadInTaxon     Reference       RefCoverage     Virus   TaxonPath       RefDescription
#25_C1   42882   117078  KY510860.1      0.11    Cherry virus A  Viruses; Riboviria; Orthornavirae; Kitrinoviricota; Alsuviricetes; Tymovirales; Betaflexiviridae; Trivirinae; Capillovirus; Cherry virus A;     KY510860.1 Cherry virus A isolate 1286 1A2/13C222_N9, complete genome
def getReport(reportFile):
    """
    get detected virus taxon id in samples
    """
    sampleVirusDict = {} #viral taxon id : [[reference, primer_f, primer_R]]
    with open(reportFile, 'r') as fin:
        next(fin)
        for line in fin:
            cells = line.strip().split("\t")
            sample = cells[0]
            taxonId = cells[1]
            if sample in sampleVirusDict:
                sampleVirusDict[sample].append(taxonId)
            else:
                sampleVirusDict[sample] = [taxonId]
    return sampleVirusDict

def callSubprocess(parallelJobs, jobs):
    """
    call subprocess.Popen to run parallel jobs, commands in the parameter parallelJobs
    """
    #run parallel jobs
    ps = [] #for processes
    threadNum = jobs  #the number of parrallel jobs
    for i in range(len(parallelJobs)):
        command = parallelJobs[i]
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True) #, shell=True
        ps.append(pipe)
        threadNum -= 1

        #to solve "OSError: [Errno 24] Too many open files", limit the number of parallel jobs
        #wait for number of jobs to finish, then submit another number of jobs
        if threadNum == 0:
            threadNum = jobs
            #wait for parallel processes to finish
            for p in ps:
                output, error = p.communicate()
                #print(p.returncode, output, error)
                if p.returncode != 0: 
                    print("plotAmplicons.py %d %s %s" % (p.returncode, output, error))
                p.wait()
            ps = []

def batchDrawPlot(viralInforFile, reportFile, virusDbFasta, scriptFolder, primerFolder, workFolder, outputFolder, jobs):
    """
    draw all amplicon coverage graph from a report
    """
    virusDict = getviralInfor(viralInforFile)
    sampleVirusDict = getReport(reportFile)
    #print(virusDict)
    #print(sampleVirusDict)

    #subprocess.Popen() is non-blocking, and returns immediately. The shell command is run in the background.
    #extract reference sequences from the database fasta file and map primers to them
    print("\nMapping primers to references. Please wait ...")
    parallelJobs = []
    for sample in sampleVirusDict:
        for taxon in sampleVirusDict[sample]:
            if taxon in virusDict:
                for vinfor in virusDict[taxon]:
                    refName = vinfor[0]
                    primerF = primerFolder + "/" + vinfor[1]
                    primerR = primerFolder + "/" + vinfor[2]
                    print("mapping", vinfor[1], "and", vinfor[2], "to", refName)
                    if not (os.path.isfile(primerF) or os.path.isfile(primerR)):
                        print("Primer files do not exist. please check their names in", viralInforFile, "Exit!")
                        exit()
                    command = ['python', scriptFolder + '/mapPrimer2Ref.py',  '-n', refName, '-d', virusDbFasta, '-f', primerF, '-r', primerR, '-o', outputFolder ] 
                    #print(" ".join(command))
                    parallelJobs.append(command)
            else:
                print("Taxon id", taxon, "is not in the file", viralInforFile)
    callSubprocess(parallelJobs, jobs)


    print("\nMapping reads to references. Please wait ...")    
    parallelJobs = []
    for sample in sampleVirusDict:
        readFile = workFolder + "/trimmed/" + sample + ".trimmed.fastq.gz"
        for taxon in sampleVirusDict[sample]:
            if taxon in virusDict:            
                for vinfor in virusDict[taxon]:
                    refName = vinfor[0]
                    primerF = outputFolder + "/" + refName + "." + vinfor[1].replace(".fa",".bed")
                    primerR = outputFolder + "/" + refName + "." + vinfor[2].replace(".fa",".bed")
                    print("Plotting ", refName, " in ", sample)
                    if not (os.path.isfile(primerF) or os.path.isfile(primerR)):
                        print("Primer files do not exist. please check their names in", viralInforFile, "Exit!")
                        exit()
                    command = ['python', scriptFolder + '/plotAmplicons.py',  '-a',  sample, '-n', refName, '-d', virusDbFasta, '-s', readFile, '-f', primerF, '-r', primerR, '-o', outputFolder ] 
                    #print(" ".join(command))
                    parallelJobs.append(command)
            else:
                print("Taxon id", taxon, "is not in the file", viralInforFile)
    
    print("Running parallel jobs now")
    callSubprocess(parallelJobs, jobs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Draw amplicon graph.")
    parser.add_argument('-w','--workdir',dest='workdir', required=True, default="", help='full path to the PVseek work directory, which has folders: report and trimmed')
    parser.add_argument('-v','--infor',dest='viralInforFile', required=True, default="", help='path to the virus information file, which has reference name, primer file names')
    parser.add_argument('-d','--database',dest='database', required=True, default="", help='plant virus database for PVseek, a fasta file')
    parser.add_argument('-s','--script',dest='script', required=True, default="", help='PVseek script folder')
    parser.add_argument('-p','--primer',dest='primerFolder', required=True, default="", help='path to primer folder having all primer files')
    parser.add_argument('-j','--job',dest='job', required=True, default="", help='numbre of jobs of parallel running, ex. 16')
    parser.add_argument('-o','--output',dest='outputFolder', required=True, default="", help='path to folder having amplicon graphs')
 
    args = parser.parse_args()
    viralInforFile = args.viralInforFile #"/ppq/data0/software/PVseek/db/Pomes_HiPlex_viruses.txt"
    virusDbFasta = args.database #"/ppq/data0/software/PVseek/db/plantvirus.fa"
    scriptFolder = args.script #"/ppq/data0/software/PVseek/scripts"
    primerFolder = args.primerFolder #"/ppq/data0/software/PVseek/db/HiPlex_primers"
    outputFolder = args.outputFolder
    jobs = int(args.job)

    workFolder = args.workdir #"/ppq/data0/test_PVseek/1_A1"
    if not (os.path.exists(workFolder + "/report") or os.path.exists(workFolder + "/trimmed")):
        print(workFolder, " does not have folders report or/and trimmed. Exit!")
        exit()
    if not os.path.isfile(viralInforFile):
        print(viralInforFile, "does not exist. Exit!")
        exit()
    if not os.path.isfile(virusDbFasta):
        print(virusDbFasta, "does not exist. Exit!")
        exit()
    if not os.path.exists(scriptFolder):
        print(scriptFolder, "does not exist. Exit!")
        exit()        
    if not os.path.exists(primerFolder):
        print(primerFolder, "does not exist. Exit!")
        exit()

    reportFile = workFolder + "/report/report.txt"
    if not os.path.isfile(reportFile):
        print(reportFile, "does not exist. Exit!")
        exit()

    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder, exist_ok=True)
        
    #draw all amplicon plots
    batchDrawPlot(viralInforFile, reportFile, virusDbFasta, scriptFolder, primerFolder, workFolder, outputFolder, jobs)
    print("\n****************************************************")
    print("All amplicon coverage grpahs are plotted")
    print("Please check files: *.coverage.ampDepth.txt and *.amplicon.png in the output foler ", outputFolder)

    #clean files
    command = "rm -rf " + outputFolder + "/*.fasta*"
    os.system(command)
    command = "rm -rf " + outputFolder + "/*.bam*"
    os.system(command)
