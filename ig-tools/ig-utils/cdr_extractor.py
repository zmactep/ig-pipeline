""" Converts the output of region_annotator to two files, one of all CDR and one of just CDR3 regions """

import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse(input):
    while True:
        line = ""
        while line.strip() == "":
            line = input.readline()
            if line == "":
                return

        id = line[1:]
        seq = input.readline()
        nom = input.readline()
        yield id, seq, nom


def get_region_seq(n, seq, nom):
    m = re.search("{0}+".format(n), nom)
    return seq[m.start():m.end()] if m else ""

def get_cdrs(seq, nom):
    nom = nom.strip()
    seq = seq.strip()
    match1 = re.search("1+", nom)
    match3 = re.search("3+", nom)
    match5 = re.search("5+", nom)
    return Seq("{0}${1}${2}".format(get_region_seq(1, seq, nom), get_region_seq(3, seq, nom), get_region_seq(5, seq, nom)))


def get_cdr3(seq, nom):
    nom = nom.strip()
    seq = seq.strip()
    return Seq(get_region_seq(5, seq, nom))


def main():
    with open(sys.argv[1], 'r') as input, open(sys.argv[2], 'w') as cdrs, open(sys.argv[3], 'w') as cdr3s:
        parsed = parse(input)
        for id, seq, nom in parsed:
            SeqIO.write(SeqRecord(get_cdrs(seq, nom), id, description=""), cdrs, 'fasta')
            SeqIO.write(SeqRecord(get_cdr3(seq, nom), id, description=""), cdr3s, 'fasta')



if __name__ == '__main__':
    main()