#!/usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 03/31/2023
#remove or retrieve sequences from a fasta file by ids

#Required Parameters:
#  -k, --key X.........................key word for selecting sequences, such as "complete genome", "partial"
#  -p, --operation X...................0 for removing sequences; 1 for extracting sequences                   
#  -i, --input X.......................fasta fromat sequence file
#  -o, --output X......................output file name 
#######################################################################

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--key', dest='key_word', required=True,
                        help='key word for selecting sequences, such as "complete genome", partial')
    parser.add_argument('-p', '--operation', dest='operation', required=True,
                        help='0 for removing sequences; 1 for extracting sequences')                        
    parser.add_argument('-i', '--input', dest='input', required=True,
                        help='path of input fasta sequence')
    parser.add_argument('-o', '--output', dest='output', required=True,
                        help='path of ouput fasta sequence')

    args = parser.parse_args()
    keyWord = args.key_word  #use " " for argument with space
    fasta_sequences = SeqIO.parse(open(args.input),'fasta')
    #print(keyWord)
    with open(args.output, "w") as fout:
        for seq in fasta_sequences:
            #name = seq.id
            #nuc = str(seq.seq)
            desc = seq.description
            #print(desc)
            if args.operation == '0':
                if keyWord not in desc.lower():
                    SeqIO.write([seq], fout, "fasta")
            elif args.operation == '1':
                if keyWord in desc.lower():
                    SeqIO.write([seq], fout, "fasta")                
            else:
                print("Wrong operation! please use --operation 0 for removing; --operation 1 for retrieving")
                sys.exit()
