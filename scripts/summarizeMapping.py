#! /usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 03/31/2023
#This program summarizes all filtered mapping result: 

#Required Parameters:
#   -w, --workdir X........................working directory
#   -f, --filter X.........................File having unwanted sequence accesions 
#   -o, --output X.........................Report file name 
#######################################################################

import argparse
import glob
import os
import pandas as pd

def catRawContent(mapFiles, outFile):
    """
    Concatenate all mapping unfiltered results
    """
    # create an Empty DataFrame object
    df = pd.DataFrame()
    empty_df = {'Sample': 'NA', 'TaxonId': 'NA', 'ReadInTaxon': 'NA', 'Reference': 'NA', 'RefCoverage': 'NA', 'Virus': 'NA', 'TaxonPath': 'NA', 'RefDescription': 'NA'}
    for mf in mapFiles:
        sample = os.path.basename(mf).replace(".mappedRef.full.txt","")
        #print(sample)
        if os.stat(mf).st_size == 0: #for empty files, file size is 0
            df1 = empty_df
            df1['Sample'] = sample
        else:
            df1 = pd.read_csv(mf, sep='\t', header = 0)
            if len(df1.index) == 0:  #for empty rows
                df1 = empty_df
                df1['Sample'] = sample
        if isinstance(df1, dict): #add a dict to pd
            #df = df.append(df1, ignore_index = True) #As of pandas 2.0, append (previously deprecated) was removed.
            #print(df1)
            df2 = pd.DataFrame.from_dict([df1]) #convert a dict to a dataframe
            df = pd.concat([df, df2], ignore_index = True)
        else:
            df = pd.concat([df, df1], ignore_index = True)

    #df = df[~df.Reference.isin(badIds)] #filter bad reference sequence 

    df.to_csv(outFile, index=False, sep='\t')


def catContent(mapFiles, badIds, outFile):
    """
    Concatenate all mapping results
    """
    # create an Empty DataFrame object
    df = pd.DataFrame()
    empty_df = {'Sample': 'NA', 'TaxonId': 'NA', 'ReadInTaxon': 'NA', 'Reference': 'NA', 'RefCoverage': 'NA', 'Virus': 'NA', 'TaxonPath': 'NA', 'RefDescription': 'NA'}
    for mf in mapFiles:
        sample = os.path.basename(mf).replace(".mappedRef.filtered.txt","")
        #print(sample)
        if os.stat(mf).st_size == 0: #for empty files, file size is 0
            df1 = empty_df
            df1['Sample'] = sample
        else:
            df1 = pd.read_csv(mf, sep='\t', header = 0)
            df1 = df1[~df1.Reference.isin(badIds)] #filter bad reference sequence
            if len(df1.index) == 0:  #for empty rows
                df1 = empty_df
                df1['Sample'] = sample
        if isinstance(df1, dict): #add a dict to pd
            #df = df.append(df1, ignore_index = True) #As of pandas 2.0, append (previously deprecated) was removed.
            #print(df1)
            df2 = pd.DataFrame.from_dict([df1]) #convert a dict to a dataframe
            df = pd.concat([df, df2], ignore_index = True)
        else:
            df = pd.concat([df, df1], ignore_index = True)

    #df = df[~df.Reference.isin(badIds)] #filter bad reference sequence 

    df.to_csv(outFile, index=False, sep='\t')


#################################### Get arguments function ####################################
def parseArguments():
    """
    Get arguments
    """
    parser = argparse.ArgumentParser(description='Generate text report')
    parser.add_argument("-w", "--workdir", dest='workdir', required=True, help="full path working directory, ex. /my/run/45")
    parser.add_argument("-o", "--output", dest='output_file', required=True, help="output file name, ex. report.txt")
    parser.add_argument("-f", "--filter", dest='badId_file', required=True, help="bad sequence accession file, ex. unwantedSeqId.txt")
    parser.set_defaults(append=False)
    return parser.parse_args()


def main():
    ### Input arguments
    options = parseArguments()
    workdir = options.workdir  #'/my/run/45'
    badIdFile = options.badId_file #/db/unwantedSeqId.txt
    outFile = options.output_file #workdir + '/report/report.txt'

    mapDir = workdir + '/mapping'
    #get bad sequence Ids
    badIds = [ line.strip() for line in open(badIdFile) ]

    #cat filtered result
    mapFiles = glob.glob(mapDir + '/*.mappedRef.filtered.txt')
    catContent(mapFiles, badIds, outFile)

    #cal all unfiltered mapping result
    mapFiles = glob.glob(mapDir + '/*.mappedRef.full.txt')
    outFile1 = outFile.replace(".txt", ".unfiltered.txt")
    catRawContent(mapFiles, outFile1)

if __name__ == "__main__":
    main()

