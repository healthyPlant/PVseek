#!/usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 04/27/2023
#since some viral refseq missed after CD-Hit cluster, reomve refseq from the cluster result, add all refseq to the database later

#Required Parameters:
#  -i, --input X.......................input file with duplcate viruses
#  -o, --output X......................output file without duplicate viruses
#  -d, --refseq X......................path of viral refseq id, one id per line
#  -c, --cluster X.....................path of CD-Hit cluster file (.fa.clstr)
#  -s, --sequence X....................path of CD-Hit clustered fasta file (.fa)
#  -o, --output X......................path of ouput fasta sequence
#######################################################################

import argparse
import re
from Bio import SeqIO

def getIds(refSeqIdFile):
    """
    Save refseq ids in a dict
    """
    idDict = { line.strip() : 1 for line in open(refSeqIdFile) }

    return idDict


def list2Dict(list, dict):
    """
    Add a list to a dict, make each item in the list as a key, the others as a value 
    """
    list_len = len(list)
    for i in range(list_len):
        key = list[i]
        value = list[0:i] + list[i+1:list_len]
        dict[key] = value


#>Cluster 0
#0       1141nt, >NC_043236.1... at +/98.77%
#1       368683nt, >NC_009898.1... *
def getAccInCluster(clusterFile):
    """
    save accesions in a cluster into a dict, each accession has members in the same cluster
    """
    clusterDict = {} 
    accList = []
    with open(clusterFile) as fin:
        for line in fin:
            if line.startswith(">"): #new cluster start
                if len(accList) > 0:
                    list2Dict(accList, clusterDict) #save accessesions to a dict
                accList = []
                continue
            else:
                #get accessions
                #0       1141nt, >NC_043236.1...
                matched = re.search(r"\d+nt,\s>(\w+)", line) #find accession, (\w+\.\d)
                acc = matched.group(1)
                #print(acc)
                accList.append(acc)
    return clusterDict


def getRemoveId(refSeqIdDict, clusterDict):
    """
    get seq accessions for removing
    """
    removeId = {}
    #removeId includes a refseq or accessions in its cluster 
    for acc in refSeqIdDict:  #refseq accessions
        #print(acc)
        #print(segmentDict[acc])
        #print(clusterDict["NC_010314"])
        if acc in clusterDict: #if refseq in a cluster
            removeId[acc] = 1
            if len(clusterDict[acc]) > 0: #other accessions in the refseq cluster
                for id in clusterDict[acc]:
                    removeId[id] = 1
        #else:
        #    print(acc)
    return removeId

def removeSeq(seqFile, newSeqFile, removeId):
    """
    remove all refseq sequences from a seqFile
    """
    cnt = 0
    total = 0
    with open(newSeqFile, "w") as f:
        for seq in SeqIO.parse(seqFile, 'fasta'):
            id = seq.id.split(".")[0]
            total += 1
            if id not in removeId:
                SeqIO.write(seq, f, "fasta")
            else:
                cnt += 1
            #    print(seq.id)
    print(cnt, " out of ", total, " sequences are removed.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--refseq', dest='refSeqIdFile', required=True,
                        help='path of viral refseq id, one id per line')
    parser.add_argument('-c', '--cluster', dest='clusterFile', required=True,
                        help='path of CD-Hit cluster file (.fa.clstr)')
    parser.add_argument('-s', '--sequence', dest='seqFile', required=True,
                        help='path of CD-Hit clustered fasta file (.fa)')
    parser.add_argument('-o', '--output', dest='output', required=True,
                        help='path of ouput fasta sequence')

    args = parser.parse_args()
    refSeqIdFile = args.refSeqIdFile  #"../plantvirusRefseqId.txt"
    clusterFile = args.clusterFile  #"plantvirus.cdhit.c95.fa.clstr"
    seqFile = args.seqFile  #"plantvirus.cdhit.c95.fa" 
    outFile = args.output #"plantvirus.cdhit.c95.noRefseq.fa"

    refSeqIdDict = getIds(refSeqIdFile)
    clusterDict = getAccInCluster(clusterFile)
    removeId = getRemoveId(refSeqIdDict, clusterDict)
    removeSeq(seqFile, outFile, removeId)




