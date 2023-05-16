#!/usr/bin/env python
######################################################################
#Alex Hu <xiaojun.hu@usda.gov>
#Updated: 03/31/2023
#This program extracts sequences information: accession, length, description

#Required Parameters:
#   -i, --input X......................fasta fromat sequence file
#   -o, --output X.....................output file name 
#######################################################################

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', required=True,
                        help='path of input fasta sequence')
    parser.add_argument('-o', '--output', dest='output', required=True,
                        help='path of ouput sequence infor')

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.input),'fasta')

    with open(args.output, "w") as fout:
        for seq in fasta_sequences:
            name = seq.id
            #nuc = str(seq.seq)
            length = len(seq)
            desc = seq.description
            desc = desc.replace(name+" ", "")
            fout.write(name + "\t" + str(length) + "\t" + desc + "\n") 

