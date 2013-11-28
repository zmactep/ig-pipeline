import argparse
import os
from itertools import zip_longest
from parsers import kabat


def main():
    parser = argparse.ArgumentParser(description='Calculate prediction score for predicted kabat.')
    parser.add_argument('--ref', nargs=1, help='reference kabat')
    parser.add_argument('--input', nargs=1, help='test kabat')
    parser.add_argument('--output', nargs=1, help='output_dir')
    args = parser.parse_args()

    compare(args.ref.pop(0), args.input.pop(0), args.output.pop(0))


def compare(ref, input, output):
    with open(ref, 'rU') as reference_file:
        ref_dict = kabat.parse(reference_file)

    with open(input, 'rU') as input_file:
        input_dict = kabat.parse(input_file)

    common_keys = set(ref_dict.keys()) & set(input_dict.keys())

    with open(os.path.join(output, 'comparison.kabat'), 'w') as comparison_file:
        for key in common_keys:
            comparison_file.write('%s\t%s\n' % (key, '\t'.join([str(pos) for pair in ref_dict[key] for pos in pair])))
            comparison_file.write('%s\t%s\n' % (key, '\t'.join([str(pos) for pair in input_dict[key] for pos in pair])))

    reference_file.close()
    input_file.close()
    comparison_file.close()

    abs_score = sum([get_abs_score(ref_dict[key], input_dict[key]) for key in common_keys])
    mismatches_per_read = abs_score / len(common_keys) if len(common_keys) != 0 else 1
    print('Absolute score: %s; Mismatches per read: %s' % (abs_score, mismatches_per_read))


def get_abs_score(first, second):
    first_regions = [start - end + 1 for start, end in first]
    second_regions = [start - end + 1 for start, end in second]
    return sum([abs(fst - snd) for fst, snd in zip_longest(first_regions, second_regions, fillvalue=0)])


if __name__ == "__main__":
    main()
