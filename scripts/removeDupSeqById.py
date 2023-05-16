#!/usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 05/18/2023
#This program remove duplicate sequences by sequence ID
#Usage: python removeDupSeqById.py [input fasta file] [output fasta file]
######################################################################

import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    print("Usage: python removeDupSeqById.py [input fasta file] [output fasta file]")
    print("Argument missed, exit!")
    sys.exit()

inFile = sys.argv[1] #input fasta file
outFile = sys.argv[2] #output fasta file

with open(outFile, 'w') as out:
    record_ids = list()
    for record in SeqIO.parse(inFile, 'fasta'):
        if record.id not in record_ids:
            record_ids.append( record.id )
            SeqIO.write(record, out, 'fasta')

            