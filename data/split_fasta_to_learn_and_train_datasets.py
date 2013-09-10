#!/usr/bin/python

import sys
import random
from Bio import SeqIO


def main(argv):
    input_file = argv[1]
    split = float(argv[2])
    train = open('train.fasta', 'w')
    test = open('test.fasta', 'w')
    input = open(input_file, "rU")
    random.seed()
    for record in SeqIO.parse(input, "fasta"):
        if random.randint(0, 100) < split:
            SeqIO.write(record, train, "fasta")
        else:
            SeqIO.write(record, test, "fasta")
    input.close()
    test.close()
    train.close()

if __name__ == "__main__":
   main(sys.argv)