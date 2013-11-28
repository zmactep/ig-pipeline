import argparse
import os
from Bio import SeqIO
from simple_splitter import split_dataset


def parse_args():
    parser = argparse.ArgumentParser(description='Split dataset sequences according to mids')
    parser.add_argument('filename', help='input filename')
    parser.add_argument('type', help='input file type, a string recognized by SeqIO.parse')
    parser.add_argument('outdir', help='output directory')
    return parser.parse_args()


def main():
    args = parse_args()

    buckets = split_dataset(SeqIO.parse(args.filename, args.type))

    print(buckets)

    for k in buckets.keys():
        with open(os.path.join(args.outdir, k + '.fasta'), 'w') as file:
            SeqIO.write(buckets[k], file, 'fasta')


if __name__ == '__main__':
    main()