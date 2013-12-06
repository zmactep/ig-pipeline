import argparse
import os
from random import random
from Bio import SeqIO
from parsers import kabat

SPLIT_RATIO = .7

def split_files(in_fasta, in_kabat, train_fasta, train_kabat, test_fasta, test_kabat, split_ratio):
    with open(in_kabat) as in_kabat_file, \
         open(train_fasta, 'w') as train_fasta_file, \
         open(train_kabat, 'w') as train_kabat_file, \
         open(test_fasta, 'w') as test_fasta_file, \
         open(test_kabat, 'w') as test_kabat_file:
        nomenclature = kabat.parse(in_kabat_file)
        for record in SeqIO.parse(in_fasta, 'fasta'):
            target_fasta, target_kabat = (test_fasta_file, test_kabat_file) \
                if random() > split_ratio else (train_fasta_file, train_kabat_file)
            SeqIO.write(record, target_fasta, 'fasta')
            target_kabat.write('{0}\t{1}\n'.format(record.id, '\t'.join(
                str(d + 1) for t in nomenclature[record.id] for d in t)))


def split(fasta_path, kabat_path, outdir, split_ratio=SPLIT_RATIO):
    def get_split_paths(path, target_dir):
        filename = os.path.basename(path)
        name_part = filename[:filename.rfind('.')]
        ext = filename[filename.rfind('.'):]
        template = os.path.join(target_dir, name_part + '{0}' + ext)
        return template.format('_train'), template.format('_test')

    train_fasta_path, test_fasta_path = get_split_paths(fasta_path, outdir)
    train_kabat_path, test_kabat_path = get_split_paths(kabat_path, outdir)

    split_files(fasta_path, kabat_path,
               train_fasta_path, train_kabat_path,
               test_fasta_path, test_kabat_path,
               split_ratio)

    return {'train': train_fasta_path, 'test': test_fasta_path}, {'train': train_kabat_path, 'test': test_kabat_path}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('kabat', help='input kabat file')
    parser.add_argument('--out_dir', help='output directory')
    parser.add_argument('--split_ratio', type=float, help='share of training data in resulting split (the rest goes to test)')
    return parser.parse_args()


def main():
    args = parse_args()

    out_dir = args.out_dir if args.out_dir else os.path.dirname(args.fasta)
    split_ratio = args.split_ratio if args.split_ratio else SPLIT_RATIO

    split(args.fasta, args.kabat, out_dir, split_ratio)


if __name__ == '__main__':
    main()