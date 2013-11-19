import argparse
from contextlib import ExitStack
import os
from Bio import SeqIO
import itertools
from parsers import kabat
from utils import kmer_generator

CAP_CHAR = '*'
TRAIN_OUTPUT = "train.libsvm"
PREDICT_OUTPUT = "predict.libsvm"
COMMENTS_OUTPUT = "read_names.txt"


def parse_arguments():
    parser = argparse.ArgumentParser(description="""Generates training or prediction datasets in libsvm format.
        For every read in input fasta file, produces a set of k-mers, labeled with region numbers from the corresponding
        kabat file, or zeros if processing new data for prediction, where no kabat is available.""")

    fasta_parent = argparse.ArgumentParser(add_help=False)
    fasta_parent.add_argument('fasta_file', help='input FASTA filename')

    kabat_parent = argparse.ArgumentParser(add_help=False)
    kabat_parent.add_argument('kabat_file', help='input KABAT filename')

    common_parent = argparse.ArgumentParser(add_help=False)
    common_parent.add_argument('k', type=int, help='the resulting k-mer length')
    common_parent.add_argument('kmer_type', choices=["left_aligned", "center_aligned"],
                               help="""specifies how generated k-mers correspond to position in input sequence:
                               k characters starting at given position or surrounding it with k/2 characters
                               left and right""")
    common_parent.add_argument('working_dir', help='working directory')

    subparsers = parser.add_subparsers(dest='mode')
    subparsers.add_parser('train', help='process training data', parents=[fasta_parent, kabat_parent, common_parent])
    subparsers.add_parser('predict', help='process data for prediction', parents=[fasta_parent, common_parent])

    # workaround for argparse bug in Python 3.3: presence of subparser arguments is not checked:
    # http://bugs.python.org/issue9253#msg186387
    subparsers.required = True


    return parser.parse_args()


def get_labels(region_bounds):
    """
    Generates sequence of k-mer labels: region numbers according to "region_bounds" mapping from parsed kabat,
    or zeros if no kabat was provided and "region_bounds" is empty.
    """
    if not region_bounds:
        return itertools.repeat(0)  # for predict data, all labels default to 0

    iterators = []
    last_index = -1
    for region, (start, end) in enumerate(region_bounds):
        if start != last_index + 1:
            raise ValueError("Gaps or overlaps in region coverage")
        iterators.append(itertools.repeat(region, end - start + 1))
        last_index = end
    return itertools.chain(*iterators)


def get_kmers(sequence, k, padded):
    """ Generates k-mers taking care of padding the sequence if requested """
    if padded:
        sequence = '%(cap)s%(seq)s%(cap)s' % {'cap': CAP_CHAR*int(k/2), 'seq': sequence}
    return kmer_generator.get_string_kmers(sequence, k)


def get_dataset(sequence, region_bounds, k, padded):
    """ Generates the dataset as sequence of tuples of region numbers and k-mers. """
    return zip(get_labels(region_bounds), get_kmers(sequence, k, padded))


def write_dataset(dataset, libsvm_file):
    """ Writes the output of get_dataset to file """
    for label, kmer in dataset:
        libsvm_file.write('\t'.join([str(label)] + ["%d:%c" % (index + 1, char) for index, char in enumerate(kmer)]) + '\n')


def process(k, padded, working_dir, libsvm_filename, fasta_filename, kabat_filename=None, comment_filename=None):
    """
    Produces a libsvm format dataset of k-mers labeled with their region numbers. K-mers are generated from input
    fasta file, and labels are taken from the kabat. If kabat is missing, zeros are used for labels. "Padded" specifies
    how k-mers are generated: from the sequence as it is, or from the sequence first padded with k/2 "empty space"
    characters on both ends.
    """
    if k % 2 == 0 or k < 1:
        raise ValueError("k must be an odd positive")

    with ExitStack() as stack:
        fasta_file = stack.enter_context(open(os.path.join(working_dir, fasta_filename), 'rU'))
        libsvm_file = stack.enter_context(open(os.path.join(working_dir, libsvm_filename), 'w'))
        comment_file = stack.enter_context(open(os.path.join(working_dir, COMMENTS_OUTPUT), 'w')) \
            if comment_filename else None
        kabat_file = stack.enter_context(open(os.path.join(working_dir, kabat_filename), 'rU')) \
            if kabat_filename else None

        region_maps = {}
        if kabat_file:
            try:
                region_maps = kabat.parse(kabat_file)
            except kabat.ParseException as e:
                raise RuntimeError from e

        for record in SeqIO.parse(fasta_file, "fasta"):
            if len(record.seq) < k:
                print("Sequence %s shorter than k=%d, skipped\n" % (record.id, k))
                continue
            region_bounds = region_maps[record.id] if record.id in region_maps else []
            dataset = get_dataset(record.seq, region_bounds, k, padded)
            write_dataset(dataset, libsvm_file)
            if comment_file:
                comment_file.write(record.id + '\n')


def main():
    args = parse_arguments()

    padded = True if args.kmer_type == "center_aligned" else False  # argparse ensures there are only two options
    try:
        if args.mode == "train":  # similarly, there are just two options:
            process(args.k, padded, args.working_dir, TRAIN_OUTPUT, args.fasta_file, args.kabat_file)
        else:
            process(args.k, padded, args.working_dir, PREDICT_OUTPUT, args.fasta_file, None, COMMENTS_OUTPUT)
    except Exception as e:
        print(e)


if __name__ == "__main__":
    main()