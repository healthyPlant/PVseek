#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 04/18/2023
##This program draws amplicon coverage graphs for detected viruses in a sample and calculate median coverage of the amplicons
#
#Required Parameters:
#   -a, --sample X.........................sample name
#   -n, --reference X......................reference name from the PVseek report
#   -d, --database X.......................plant virus database for PVseek, a fasta file
#   -s, --reads X..........................trimmed reads file
#   -f, --forward X........................forward primer fasta file
#   -r, --reverse X........................reverse primer fasta file
#   -o, --output X.........................coverage grpah 
######################################################################

import argparse
import pandas as pd
import os
import numpy as np

#to solve error: couldn't connect to display "localhost"
import matplotlib
try:
    matplotlib.use("Agg")
except ValueError:
    pass
import matplotlib.pyplot as plt
from matplotlib import collections  as mc

def mapReads2Ref(outFolder, sampleName, refName, db, reads):
    """
    map reads to the reference, then output the coverage file from the bam file 
    """
    output = outFolder + "/" + sampleName + "." + refName

    #if refseq doesn't exit, use seqtk subseq to get reference sequence
    refFile = outFolder + "/" + refName 
    if not os.path.exists(output + ".fasta"):
        command = "echo " + refName + " > " + refFile 
        #print(command)
        os.system(command)
        command = "seqtk subseq " + db + " " + refFile + " > " + refFile + ".fasta"
        #print(command)
        os.system(command)

    refFile = refFile + ".fasta"
    if os.stat(refFile).st_size == 0:
        print(refName, 'does not exist in the file', db)
        exit()
    
    #if bwa index file doesn't exist, build one
    if not os.path.isfile(refFile + ".sa"):
        command = "bwa index " + refFile
        #print(command)
        os.system(command)

    #map reads to the reference using bwa with default parameters    
    command = "bwa mem " + refFile + " " + reads + " | samtools view -Sb -F 4 - | samtools sort -o " + output + ".bam"
    #print(command)
    os.system(command) 
    command = "samtools index " + output + ".bam"
    #print(command)
    os.system(command)
    #get the genome coverage data using bedtools genomecov 
    coverageFile = output + ".coverage.bed"
    #os.system("samtools depth -a " + output + ".bam" + ">" + coverageFile)
    command = "bedtools genomecov -bga -ibam " + output + ".bam" + " -g " + refFile + " > " + coverageFile
    #print(command)
    os.system(command)

    return coverageFile

def getPrimerPosition(primerBedFile, ypPos):
    """
    prepare start and end points for primers
    """
    #lines = [[(0, 1), (1, 1)], [(2, 3), (3, 3)], [(1, 2), (1, 3)], [(x1, y1), (x2, y2)]]
    primerLines = []
    primerNames = []
    cnt = 0
    
    with open(primerBedFile, 'r') as fin:
        for line in fin:
            cells = line.strip().split("\t")

            #adjst y axis position, every three segments a cycle
            ypPos1 = ypPos
            if cnt % 3 == 1:
                ypPos1 -= ypPos*0.05 
            elif cnt % 3 == 2:
                ypPos1 -= ypPos*0.1 

            startPoint = (cells[1], ypPos1)
            endPoint = (cells[2], ypPos1)
            name = cells[3]
            primer = [startPoint, endPoint]
            primerLines.append(primer)
            primerNames.append(name)
            cnt += 1
    return primerLines, primerNames

def getAmpAvgDepth(depth, forwardPrimer, reversePrimer):
    """
    Calculate mean or median for a amplicon 
    """
    amplicons = []
    ampAvgDep = []
    for i in range(0, len(forwardPrimer)):
        start = forwardPrimer[i][0][0]
        end = reversePrimer[i][1][0]
        amplicons.append([start,end])

    for amp in amplicons:
        startIndex = int(amp[0])
        endIndex = int(amp[1]) + 1
        #print(startIndex, endIndex)
        ampDepth = depth[startIndex : endIndex]
        if ampDepth:
            median_dep = np.median(ampDepth)
        else:
            median_dep = 0 #in case no depth in this range, ampDepth is empty
        ampAvgDep.append(median_dep)
        #mean_dep = np.mean(ampDepth)
        #ampAvgDep.append(mean_dep)        
    return ampAvgDep

#KY510860.1      0       54      0
#KY510860.1      54      69      1
def getCoverageFromBed(covBedFile):
    """
    get depth of coverage from a coverage file by bedtools genomecov, the start position is included, the end position is excluded
    """
    #title = covBedFile.replace(".bed","")
    depth = []
    refName = ""
    with open(covBedFile, 'r') as fin:
        for line in fin:
            cells = line.strip().split() #"\t"
            refName = cells[0]
            start = int(cells[1])
            end = int(cells[2])
            for i in range(start, end):
                depth.append(int(cells[3]))
    
    return refName, depth 


def getCoverage(covFile):
    """
    get depth of coverage from a coverage file by samtools depth/pileup
    """
    #base = os.path.basename(covFile) #get basename of the file without path
    #title = base.replace(".txt","")
    #title = covFile.replace(".txt","")   
    #title = tabFile.replace(".pileup.tab","")
    #read coverage data
    table = pd.read_csv(covFile, sep='\t', header=None)
    #print(table)
    alignment = pd.DataFrame(data=table)
    refName = alignment.iloc[0,0]
    position = list(alignment.iloc[:,1]) #.to_dict()
    depth = alignment.iloc[:,2].to_dict()
    #print(position)
    #print(depth)

    #some positions with 0 depth are missing, fill them up
    newDepth = []
    for i in range(1, max(position)):
        if i in position:
            idx = position.index(i)
            newDepth.append(depth[idx])
        else:
            newDepth.append(0)

    return refName, newDepth 


#################################### Plot coverage function ####################################
def plot_pileup_coverage(covFile, title, fPrimerFile, rPrimerFile, graphFile, ampDepthFile):
    """
    plot coverage graph 
    """
    #get virus coverage from read mapping
    #refName, depth = getCoverage(covFile)
    refName, depth = getCoverageFromBed(covFile)
    #print(depth)
    max_dep = max(depth)
    min_dep = min(depth)
    #mean_dep = sum(depth)/len(depth)
    #mean_dep = np.mean(depth)
    #median_dep = np.median(depth)

    #make primer lines using positions from a bed file
    ypPos = max_dep*1.10 #primer position on y axis
    if max_dep < 100 and max_dep > 60:
        ypPos = 120
    elif max_dep < 60:
        ypPos = 70
    #lines = [[(0, 1), (1, 1)], [(2, 3), (3, 3)], [(1, 2), (1, 3)]]
    forwardPrimerPos, fwdPrimName = getPrimerPosition(fPrimerFile, ypPos)
    reversePrimerPos, rvsPrimName = getPrimerPosition(rPrimerFile, ypPos)
    amplicon = []
    for i in range(0, len(forwardPrimerPos)):
        #print(forwardPrimer[i], reversePrimer[i])
        startPoint = forwardPrimerPos[i][1]
        endPoint = reversePrimerPos[i][0]
        amplicon.append([startPoint,endPoint])
    #print(forwardPrimer)
    #print(reversePrimer)
    #print(amplicon)

    #output amplicon depth
    outputDepth(ampDepthFile, refName, depth, fwdPrimName, rvsPrimName,forwardPrimerPos, reversePrimerPos)

    #use LineCollection to draw primers and amplicons 
    flc = mc.LineCollection(forwardPrimerPos, colors='green', linewidths=8, linestyles='solid') #forward primer
    rlc = mc.LineCollection(reversePrimerPos, colors='red', linewidths=8, linestyles='solid') #reverse primer
    alc = mc.LineCollection(amplicon, colors='orange', linewidths=5, linestyles='solid') #amplicon

    fig, ax = plt.subplots()
    ax.add_collection(flc)
    ax.add_collection(rlc)
    ax.add_collection(alc)
    #ax.autoscale()
    #ax.margins(0.1)

    #plt.plot(position, depth, 'bo', markersize=1)
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                    box.width, box.height * 0.9])    
    
    plt.stackplot(range(0,len(depth)),depth,colors='lightblue',edgecolor='black')
    plt.xlabel('Position in Genome ' + refName)
    plt.ylabel('Depth of Coverage')
    ymax = max_dep*1.25
    if max_dep < 100 and max_dep > 60:
        ymax = 125
    elif max_dep < 60:
        ymax = 75

    plt.ylim(min_dep+0.5, ymax)
    #add mean value
    #plt.axhline(median_dep, color='r', alpha=0.5, linestyle='--', linewidth=1)
    #plt.legend(('median: {:5.0f}'.format(median_dep),'depth'),loc='upper right')
    #draw a horizontal line
    #plt.hlines(y=0.6, xmin=0.0, xmax=1.0, color='b')
    # Put a legend below current axis
    plt.legend(('forward primer', 'reverse primer', 'amplicon', 'depth'), loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=4) 

    plt.title(title)
    #plt.show()

    plt.gcf().set_size_inches(10,6)
    plt.savefig("%s" % graphFile)

    print("Please check amplicon depth figure in the file ", graphFile)


def outputDepth(ampDepthFile, refName, depth, fwdPrimName, rvsPrimName, forwardPrimer, reversePrimer):
    """
    output amplicon depth to a file
    """
    #calcualte average depth for each amplicon
    ampAvgDep = getAmpAvgDepth(depth, forwardPrimer, reversePrimer)
    #print(ampAvgDep)

    sample = os.path.basename(ampDepthFile).split(".")[0]
    #output amplicon depth
    fout = open(ampDepthFile, 'w')
    fout.write("Sample\tReference\t")
    for i in range(0, len(ampAvgDep)):
        fout.write(fwdPrimName[i]+" | " + rvsPrimName[i]+"\t")
    fout.write("\n")
    avgDepth = ["{:.2f}".format(i) for i in ampAvgDep] #format number, map(str,ampAvgDep)
    fout.write(sample + "\t" + refName + "\t" + "\t".join(avgDepth) + "\n")
    fout.close()

    print("Please check amplicon average depth in the file ", ampDepthFile)



#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description = "Draw amplicon graph.")
    parser.add_argument('-a','--sample',dest='sample_name', required=True, default="", help='sample name')
    parser.add_argument('-n','--reference',dest='reference_name', required=True, default="", help='reference name from the PVseek report')
    parser.add_argument('-d','--database',dest='database', required=True, default="", help='plant virus database for PVseek, a fasta file')
    parser.add_argument('-s','--reads',dest='reads', required=True, default="", help='path to sample trimmed reads file')
    parser.add_argument('-f','--forward',dest='forward_primer_file', required=True, default="", help='forward primer bed file')
    parser.add_argument('-r','--reverse',dest='reverse_primer_file', required=True, default="", help='reverse primer bed file')
    parser.add_argument('-o','--output',dest='output_folder', required=True, default="", help='path to amplicon graphs')
    parser.set_defaults(append=False)
    return parser.parse_args()

def main():
    ### Input arguments
    options = parseArguments()
    #covFile = options.coverage
    sampleName = options.sample_name
    refName = options.reference_name
    db = options.database
    reads = options.reads
    fPrimerBedFile = options.forward_primer_file  #forward primers positions in the reference
    rPrimerBedFile = options.reverse_primer_file  #reverse primers positions in the reference  
    #fPrimerBedFile = "outFolder/virus_primer_F.bed"
    #rPrimerBedFile = "outFolder/virus_primer_R.bed"
    outFolder = options.output_folder
    if not os.path.exists(fPrimerBedFile) or not os.path.exists(rPrimerBedFile):
        print(fPrimerBedFile, "and/or", rPrimerBedFile, "do not exist. Exit!")
        exit()
    elif os.stat(fPrimerBedFile).st_size == 0 or os.stat(rPrimerBedFile).st_size == 0:
        print(fPrimerBedFile, "and/or", rPrimerBedFile, "are empty. Exit!")
        exit()        

    #get virus name from primer bed file, which format is {refName}.{virusName}_primer_F.bed
    virusName = fPrimerBedFile.split(".")[-2].replace("_primer_F","")
    title = sampleName + "." + virusName + "." + refName 
    graphFile = outFolder + "/" + title + ".amplicon.png"
    ampDepthFile = outFolder + "/" + title + ".ampDepth.txt"
    #if outout folder is not exists, create one
    if not os.path.isdir( outFolder ):
        os.makedirs( outFolder )

    #map reads to the reference
    covFile = mapReads2Ref(outFolder, sampleName, refName, db, reads)
    #covFile = "outFolder/sample.refName.coverage.bed"

    #draw amplicon graph and output amplicon depth (median)
    plot_pileup_coverage(covFile, title, fPrimerBedFile, rPrimerBedFile, graphFile, ampDepthFile)

if __name__ == '__main__':
	main()
