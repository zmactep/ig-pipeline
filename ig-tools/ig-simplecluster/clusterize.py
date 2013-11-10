__author__ = 'mactep'

import os
import argparse

from cluster import clustalo
from similarity import alignedsim


def main():
    parser = argparse.ArgumentParser(description="Small clustering utility")
    parser.add_argument('--src', action='store', help='source FASTA')
    parser.add_argument('--sim', action='store_const', metavar='ctype', dest='ctype',
                        default='cluster', const='similatiry', help='cluster by equality')
    parser.add_argument('--skip-first', action='store', metavar='m', dest='m', type=int,
                        help="skip first m letters in all sequences for similarity check (only for similarity mode)")
    parser.add_argument('--use-prct', action='store', metavar='useprct', dest='useprct', type=int,
                        help="use just n%% length of aligned sequences (only for similarity mode)")
    parser.add_argument('--shortest-cons', action='store_const', metavar='cons', dest='cons',
                        default=False, const=True, help='make shortest consensus for similarity mode')
    parser.add_argument('--min-len', action='store', metavar='minlen',
                        dest='minlen', type=int, help="minimal sequence length")
    parser.add_argument('--outdir', action='store', help='output directory')

    args, unknown = parser.parse_known_args()

    if not os.path.isdir(os.path.abspath(args.outdir)):
        os.mkdir(os.path.abspath(args.outdir))

    if args.ctype == 'cluster':
        clustalo.run(args.src, args.outdir, args.minlen)
    else:
        alignedsim.run(args.src, args.outdir, args.minlen, args.m, args.useprct, args.cons)

if __name__ == "__main__":
    main()