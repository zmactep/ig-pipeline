#!/usr/bin/python

import sys
import random
from Bio import SeqIO


def main(argv):
    input_file = argv[1]
    output_folder = argv[2]
    split = float(argv[3])
    train = open(output_folder + 'train.fasta', 'w')
    test = open(output_folder + 'test.fasta', 'w')
    random.seed()
    for record in SeqIO.parse(input_file, "fasta"):
        if random.uniform(0, 100) < split:
            SeqIO.write(record, train, "fasta")
        else:
            SeqIO.write(record, test, "fasta")
    test.close()
    train.close()

if __name__ == "__main__":
   main(sys.argv)
