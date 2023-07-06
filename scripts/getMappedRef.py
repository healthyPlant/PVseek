#!/usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 03/312023
#This program counts mapped reads by taxons, then selects a representive reference with highest mapped reads in a taxon and cauculate its viral genome coverage
#Required Parameters:
#   -m, --map X........................mapping file having read number for references from a bam file
#   -d, --depth X......................depth file having depth for a reference at each position, from samtools depth"
#   -r, --refInfor X...................reference name and sequence length 
#   -n, --readNum X....................minimum read number for slecting a taxon
#   -g, --genomeCov X..................minmum virus genome coverage for selecting a reference in a taxon
#   -t, --taxon X......................taxon db file having taxon id for each reference, ex. plantvirus.gb_taxon.txt
#   -o, --output X.....................output file name 
#######################################################################

import os.path
import argparse   #take arguments

#AB000282.1      85516   Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes; Picornavirales; Secoviridae; Sadwavirus; Satsumavirus; Satsuma dwarf virus; Navel orange infectious mottling virus;
def getTax(taxonFile):
    """
    Get a taxon id for a reference
    """
    taxonDict = {}
    with open(taxonFile) as fin:
        for line in fin:
            cells = line.split("\t") 
            taxonDict[cells[0]] =[cells[1], cells[2]]  #reference name : [ taxon id, taxon path]
    return taxonDict

#502 MH058008.1
def getMappedReadNum(mapFile):
    """
    Get mapped read number for a reference
    """
    refDict = {}
    with open(mapFile) as fin:
        for line in fin:
            cells = line.strip().split()
            refDict[cells[1]] = int(cells[0])  #reference name : mapped read number
    return refDict 

def countTaxon(taxonDict, refDict):
    """
    Count read number for a taxon and the reference with a highest number in a taxon
    """
    taxonRef = {} #keep a reference with the highest number for a taxon 
    taxonTotalNum = {} #total read number for a taxon
    for ref in refDict:
        if ref in taxonDict:
            taxon = taxonDict[ref][0]
            if taxon in taxonTotalNum:
                taxonTotalNum[taxon] += refDict[ref]
            else:
                taxonTotalNum[taxon] = refDict[ref]
            
            if taxon in taxonRef:
                if taxonRef[taxon][1] < refDict[ref]:
                    taxonRef[taxon] = [ref, refDict[ref]]
            else:
                taxonRef[taxon] = [ref, refDict[ref]]
    #print(taxonTotalNum)
    #print(taxonRef)
    return (taxonTotalNum, taxonRef)

#AB000048        2007    protoparvovirus Feline panleukopenia virus gene for nonstructural protein 1, complete cds, isolate: 483.
def getRefInfor(refInforFile):
    """
    Get reference sequence legnth and description
    """
    refInforDict = {} 
    with open(refInforFile) as fin:
        for line in fin:
            cells = line.strip().split("\t")
            refInforDict[cells[0]] = [int(cells[1]), cells[2]]  #reference name : reference length, description
    return refInforDict 

#MN882023        5946    2
def getCoverage(depthFile, refInforDict):
    """
    Get genome coverage for each reference
    """
    covDict = {} 
    with open(depthFile) as fin:
        for line in fin:
            cells = line.split()
            ref = cells[0]
            if ref in covDict:
                covDict[ref] += 1  #count covered positions for each reference
            else:
                covDict[ref] = 1
    
    #calculate coverage: covered positions/reflength
    for ref in covDict:
        if ref in refInforDict:
            coverage = covDict[ref] / refInforDict[ref][0]
            covDict[ref] = f"{coverage:.2f}"
        else:
            print(ref, "does not in the database")
            covDict[ref] = 'nan'
        
    return covDict 



def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description='Generate text report')
    parser.add_argument("-m", "--map", dest='mapFile', required=True, help="mapping file having read number for references from a bam file")
    parser.add_argument("-d", "--depth", dest='depthFile', required=True, help="depth file having depth for a reference at each position, from samtools depth")
    parser.add_argument("-r", "--refinfor", dest='refInforFile', required=True, help="reference sequence length and description")
    parser.add_argument("-t", "--taxon", dest='taxonFile', required=True, help="taxon db file having taxon id for each reference, ex. plantvirus.gb_taxon.txt")
    parser.add_argument("-n", "--readNum", dest='readNumThreshold', required=True, help="minimum read number in a taxon")
    parser.add_argument("-g", "--genomeCov", dest='GenomeCoverageThreshold', required=True, help="minmum virus genome coverage")
    parser.add_argument("-o", "--output", dest='outputFile', required=True, help=".output file name")
    return parser.parse_args()   

def main():
    ### Input arguments
    options = parseArguments()
    mapFile = options.mapFile #"sample.mappedRead.txt"
    taxonFile = options.taxonFile #"/db/plantvirus.gb_taxon.txt"
    depthFile = options.depthFile #"sample.coverage.txt"
    refInforFile = options.refInforFile #db/viralRef.info.txt
    outFile = options.outputFile #"sample.mappedRef"

    coverage_threshold = float(options.GenomeCoverageThreshold) #minimum genome coverage, 0.05
    readNum_threshold = float(options.readNumThreshold) #minimum read in a taxon, 10

    taxonDict = getTax(taxonFile)
    refDict = getMappedReadNum(mapFile)
    refInforDict = getRefInfor(refInforFile)
    taxonTotalRead, taxonRef = countTaxon(taxonDict, refDict)
    covDict = getCoverage(depthFile, refInforDict)
    sample = os.path.basename(mapFile).replace(".mappedRead.txt","")

    #output result
    fout1 = open(outFile+".full.txt", 'w')
    fout2 = open(outFile+".filtered.txt", 'w')
    header = "Sample\tTaxonId\tReadInTaxon\tReference\tRefCoverage\tVirus\tTaxonPath\tRefDescription\n"
    if taxonTotalRead:
        fout1.write(header)
        fout2.write(header)
    #sort dict by value
    #print(sorted(taxonTotalRead, key=taxonTotalRead.get, reverse=True))
    for taxon in sorted(taxonTotalRead, key=taxonTotalRead.get, reverse=True): 
        ref = taxonRef[taxon][0]
        taxonPath = taxonDict[ref][1].strip()
        readNum = taxonTotalRead[taxon]
        if ref in covDict and taxonPath: 
            cov = float(covDict[ref])
            virus = taxonPath.split(";")[-2].strip()
            if ref not in refInforDict:
                fout1.write(sample + "\t" + taxon + "\t" + str(readNum) + "\t" + ref + "\t" + covDict[ref] + "\t" + virus + "\t" + taxonPath + "\tNone\n")
            else:
                fout1.write(sample + "\t" + taxon + "\t" + str(readNum) + "\t" + ref + "\t" + covDict[ref] + "\t" + virus + "\t" + taxonPath + "\t" + refInforDict[ref][1] + "\n")

        #filter viruses: 
        #1. For Amplicon or HiPlex data, the genome coverage is ignored, the only criteria is the mapped read number greater than the thresholds
        #2. if the genome coverage is greater than 75%, keep it; 
        #3. if both the coverage and the mapped read number are greater than the thresholds, keep it
        if ref in covDict and taxonPath:
            if (coverage_threshold == 0.0 and readNum >= readNum_threshold) or (coverage_threshold != 0.0 and cov >= 0.75) or (cov >= coverage_threshold and readNum >= readNum_threshold):
                cov = float(covDict[ref]) 
                virus = taxonPath.split(";")[-2].strip()
                fout2.write(sample + "\t" + taxon + "\t" + str(readNum) + "\t" + ref + "\t" + covDict[ref] + "\t" + virus + "\t" + taxonPath + "\t" + refInforDict[ref][1] + "\n")

    fout1.close()
    fout2.close()

if __name__ == "__main__":
    main()

