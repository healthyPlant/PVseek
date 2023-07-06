#!/usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 04/21/2023
#get viruses with multiple segments

#Required Parameters:
#  -i, --input X.......................fasta fromat sequence file
#  -o, --output X......................output file name 
#######################################################################

import subprocess
import multiprocessing 
import argparse
import sys

def upcase_first_letter(s):
    """
    Upcase first letter of the sentence
    """
    return s[0].upper() + s[1:]

#Abaca bunchy top virus	NC_010314, NC_010315, NC_010316, NC_010317, NC_010318, NC_010319
#abaca bunchy top virus	DNA-C: NC_010318; DNA-M: NC_010317; DNA-N: NC_010314; DNA-R: NC_010319; DNA-S: NC_010316; DNA-U3: NC_010315
def getVirus(viralInforFile):
    """
    read virus name and its segments, remove duplicates
    """
    viralDict = {} 
    with open(viralInforFile) as fin:
        for line in fin:
            cells = line.strip().split("\t")
            name = upcase_first_letter(cells[0])
            segments = ""
            if ", " in cells[1]:
                segments =  cells[1].split(", ")
            elif "; " in cells[1]:
                segments = cells[1].split("; ")
            if name not in viralDict:    
                viralDict[name] = segments
    return viralDict

def getSegments(virus): #, viralSeqFile=viralSeqFile
    virus1 = virus.replace("[", "\[") #add escape to special character
    virus1 = virus1.replace("]", "\]")
    cmd = "grep \"" + virus1 + "\" " + viralSeqFile
    #print(virus)
    try:
        viralDesc = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
        viralDesc = viralDesc.decode("utf-8").strip() # using decode() function to convert byte string to string
        listDesc = viralDesc.split("\n")
        accessions = []
        descs = []
        for desc in listDesc:
            desc = desc.replace(">","")  #remove ">"
            acc,des = desc.split(' ', 1) #split it to two parts: access and description 
            accessions.append(acc)
            descs.append(des)
        return virus, ";".join(accessions), ";".join(descs)
    except subprocess.CalledProcessError as e:
        print(e.cmd)

def multipleRunGetSegments(viralDict):
    """
    use multiprocessing to call function: getSegments
    """
    viralName = list(viralDict.keys())
    #print(viralName)
    #Parallel Processing multiple fastq.gz file 
    pool = multiprocessing.Pool() 
    outputs_async = pool.map_async(getSegments, viralName)  #map_async for a single arguments
    outputs = outputs_async.get() #return a list of sets:[(virus name, aceesions, descriptions)]
    #print(outputs)
    return outputs

if __name__ == '__main__':
    global viralSeqFile  #set it as a globle viriable, the getSegments function can use it
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infor', dest='viralInforFile', required=True,
                        help='path of viral information file having virus name and its segments')
    parser.add_argument('-s', '--sequence', dest='sequenceFile', required=True,
                        help='path of viral fasta sequence')
    parser.add_argument('-o', '--output', dest='output', required=True,
                        help='path of ouput fasta sequence')

    args = parser.parse_args()
    viralInforFile = args.viralInforFile

    viralSeqFile = args.sequenceFile
    outFile = args.output

    #get virus reference information: {virus name : [segments]}
    viralDict = getVirus(viralInforFile)
    #print(viralDict)

    #get virus sequence description: a list of sets: (virus, accessions, descs)
    viralDesc = multipleRunGetSegments(viralDict)
    #print(viralDesc)

    newViralDict = {}
    fout = open(outFile, 'w')
    header = "Virus\tSegmentNumber\tSequenceNumber\tDescription\n"
    fout.write(header)
    for virus in viralDesc:
        print(virus)
        if virus is not None:
            name = virus[0]
            newViralDict[name] = 1
            segNum = len(viralDict[name])
            accs = virus[1].split(";")
            descs = virus[2].split(";")
            #output a set of segments per line
            i = 0
            while i < len(accs):
                segAcc = accs[i : i+segNum]
                segDesc = descs[i : i+segNum]
                i += segNum 
                fout.write(name+"\t"+str(segNum)+"\t"+str(len(accs))+"\t"+";".join(segAcc)+"\t"+";".join(segDesc)+"\n")
    fout.close()

    #check missed viruses
    print("total virus: ", len(viralDict))
    print("available virus: ", len(newViralDict))
    print("Missed viruses:")
    for virus in viralDict:
        if virus not in newViralDict:
            print(virus)
