#!/usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 03/31/2023
#select database sequences, keep complete genome and filter gene, protein, cds

#Required Parameters:
#  -i, --input X.......................fasta fromat sequence file
#  -o, --output X......................output file name 
#######################################################################

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


#plantvirus.gb_taxon.txt
#AB000282.1      85516   Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes; Picornavirales; Secoviridae; Sadwavirus; Satsumavirus; Satsuma dwarf virus; Navel orange infectious mottling virus;
def getSeqTaxon(gb_taxon_file):
    """
    Get sequence taxon infor
    """
    seqTaxonDict = {} 
    with open(gb_taxon_file) as fin:
        for line in fin:
            cells = line.split("\t")
            seqTaxonDict[cells[0]] = [cells[1]]
    return seqTaxonDict

#plantvirus.fa
#>X16637.1 Phaseolus vulgaris endornavirus dsRNA fragment
#TCTGATCGTGATGTTCAAGTCATGAAATACGATTCTGGATGTGACAGCTTTGCATGCATGGTGTTAGTTA
def getSeqDesc(seqFile):
    """
    Get sequence description
    """
    seqObjs = SeqIO.parse(open(args.input),'fasta')
    seqDict = {}
    for seq in seqObjs:
        name = seq.id
        #nuc = str(seq.seq)
        desc = seq.description
        seqDict[name] = desc

    return seqDict

def getRelation(gb_taxon_file, seqFile):
    """
    build a relationship between taxon, sequence and annotation using information from gb_taxon_file and seqFile
    """
    seqTaxonDict = getSeqTaxon(gb_taxon_file)
    seqDict = getSeqDesc(seqFile)

    #add seq description to seqTaxonDict
    for sname in seqTaxonDict:
        if sname in seqDict:
            seqTaxonDict[sname].append(seqDict[sname])
        else:
            print(sname, " is not available in ", gb_taxon_file)
    return seqTaxonDict

def getFilteredSeqId(seqTaxonDict):
    """
    Count sequences in a taxon
    """
    taxonDict = {}
    filteredTaxonDict = {} #filter discription
    filteredSeqId = {}

    for sname in seqTaxonDict:
        taxon = seqTaxonDict[sname][0]
        desc = seqTaxonDict[sname][1]
        if taxon in taxonDict:
            taxonDict[taxon].append(sname)
        else:
            taxonDict[taxon] = [sname]
        
        #select sequences: (1). if it's a complete genome; (2). it's a complete sequence, but not gene, protein, cds
        if "complete genome" in desc:
            if taxon in filteredTaxonDict:
                filteredTaxonDict[taxon].append(sname)
            else:
                filteredTaxonDict[taxon] = [sname]
        elif "complete" in desc and not ("gene" in desc or "protein" in desc or "cds" in desc):
            if taxon in filteredTaxonDict:
                filteredTaxonDict[taxon].append(sname)
            else:
                filteredTaxonDict[taxon] = [sname]
    #print(len(filteredTaxonDict))
    #print(seqTaxonDict)

    #check missed taxon after filtering
    #if the taxon is missed, put all filtered sequences back, leave them to cluster method to filter
    for taxon in taxonDict:
        #print(taxon, len(taxonDict[taxon]), sep="\t", end="")
        if not taxon in filteredTaxonDict:
            #print("\t",'0', end="\t")
            #print(",".join(taxonDict[taxon]))
            filteredTaxonDict[taxon] = taxonDict[taxon] 
        #else:
            #print("\t", len(filteredTaxonDict[taxon]))
    #print(filteredTaxonDict)

    for taxon in filteredTaxonDict:
        for seq in filteredTaxonDict[taxon]:
            filteredSeqId[seq] = 1
    
    return filteredSeqId                

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gbtaxon', dest='gbTaxonFile', required=True,
                        help='path of gb_taxon file, ex. plantvirus.gb_taxon.txt')
    parser.add_argument('-i', '--input', dest='input', required=True,
                        help='path of input fasta sequence')
    parser.add_argument('-o', '--output', dest='output', required=True,
                        help='path of ouput fasta sequence')

    args = parser.parse_args()
    gb_taxon_file = args.gbTaxonFile
    seqFile = args.input
    seqTaxonDict = getRelation(gb_taxon_file, seqFile)
    selectSeqId = getFilteredSeqId(seqTaxonDict)

    #get selected sequences
    totalSequences = SeqIO.parse(open(seqFile),'fasta')
    #print(keyWord)
    cnt = cnt0 = 0
    with open(args.output, "w") as fout:
        for seq in totalSequences:
            cnt0 += 1
            name = seq.id
            nuc = str(seq.seq)
            #desc = seq.description
            #print(desc)
            if name in selectSeqId and len(nuc) >= 150: #minimum sequence length 150nt
                SeqIO.write([seq], fout, "fasta")
                cnt += 1
    print(cnt, " out of ", cnt0, " sequences are selected for clustering")
