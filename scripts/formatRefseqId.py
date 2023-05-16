#!/usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 04/24/2023
#Refseq Ids from Virus-Host-DB and ICTV are needed to format to one id per line and remove duplicates

#Required Parameters:
#  -i, --input X.......................input file for combined refseq ids from Virus-Host-DB and ICTV
#  -o, --output X......................output file for formatted refseq ids
#######################################################################

import argparse

#beet curly top Iran virus	
#Exomis microphylla associated virus	NC_037065
#Potato virus B	2	NC_043447,NC_043448
#horsegram yellow mosaic virus	DNA-A: NC_005635; DNA-B: NC_005636
#olive latent ringspot virus	RNA2: NC_038863
def getSegmentAcc(seqIdFile):
    """
    Get virus refseq accesions 
    """
    idDict = {}
    with open(seqIdFile) as fin:
        next(fin) #skip header
        for line in fin:
            cells = line.strip().split("\t")
            if len(cells) > 1:
                if "," in cells[1]: #host-db format
                    accs = cells[1].split(", ")
                elif ";" in cells[1]: #ICTV format
                    accs = cells[1].split("; ")
                    if ":" in accs[0]: 
                        accs = [x.split(":")[-1].strip() for x in accs]
                elif ":"  in cells[1]:
                    accs = [cells[1].split(":")[-1].strip()]
                else:
                    accs = [cells[1]]

                for acc in accs:
                    idDict[acc] = 1
    return idDict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', required=True,
                        help='path of input refseq id file')
    parser.add_argument('-o', '--output', dest='output', required=True,
                        help='path of ouput file')

    args = parser.parse_args()
    seqIdFile = args.input
    outFile = args.output

    idDict = getSegmentAcc(seqIdFile)
    fout = open(outFile, 'w')
    for id in idDict:
        fout.write(id + "\n")
    fout.close()

