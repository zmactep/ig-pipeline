from IPython.parallel import Client

import argparse
import os
from Bio import SeqIO
from simple_splitter import split_dataset
from itertools import chain

view = Client()[:] 


def parse_args():
    parser = argparse.ArgumentParser(description='Split dataset sequences according to mids')
    parser.add_argument('filename', help='input filename')
    parser.add_argument('type', help='input file type, a string recognized by SeqIO.parse')
    parser.add_argument('outdir', help='output directory')
    return parser.parse_args()


@view.parallel()
def split_dataset_parallel(recs):
    from simple_splitter import split_dataset
    return split_dataset(recs)


def merge_dicts(ds):
    m = {}

    for d in ds:
        for k, v in d.items():
            if k in m:
                m[k].append(iter(v))
            else:
                m[k] = [iter(v)]
    
    for k in m:
        m[k] = chain(*m[k])
    return m


def run_split(recs):
    result = split_dataset_parallel(recs)

    result.wait_interactive()

    return merge_dicts(result.result)


def main():
    args = parse_args()

    buckets = run_split(SeqIO.parse(args.filename, args.type))

    for k in buckets.keys():
        with open(os.path.join(args.outdir, k + '.fasta'), 'w') as file:
            SeqIO.write(buckets[k], file, 'fasta')


if __name__ == '__main__':
    main()