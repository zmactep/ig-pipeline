from IPython.parallel import Client

import argparse
import os
from Bio import SeqIO
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


def reduce_dicts(ds):
    r = {}

    for d in ds:
        for k, v in d.items():
            if k in r:
                r[k].append(iter(v))
            else:
                r[k] = [iter(v)]

    for k in r:
        r[k] = chain(*r[k])

    return r


def run_split(recs):
    result = split_dataset_parallel(recs)

    result.wait_interactive()

    return reduce_dicts(result.result)


def main():
    args = parse_args()

    buckets = run_split(SeqIO.parse(args.filename, args.type))

    for k in buckets.keys():
        with open(os.path.join(args.outdir, k + '.fasta'), 'w') as file:
            SeqIO.write(buckets[k], file, 'fasta')


if __name__ == '__main__':
    main()