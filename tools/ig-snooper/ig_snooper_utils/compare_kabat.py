import argparse
import os
from itertools import zip_longest

def main():
    parser = argparse.ArgumentParser(description='Calculateprediction score for predicted kabat.')
    parser.add_argument('--ref', nargs=1, help='reference kabat')
    parser.add_argument('--input', nargs=1, help='test kabat')
    parser.add_argument('--output', nargs=1, help='output_dir')
    args = parser.parse_args()

    compare(args.ref.pop(0), args.input.pop(0), args.output.pop(0))


def compare(ref, input, output):
    with open(ref, 'rU') as reference_file:
        ref_dict = dict([parse_kabat_line(line) for line in reference_file])

    with open(input, 'rU') as input_file:
        input_dict = dict([parse_kabat_line(line) for line in input_file])

    common_keys = set(ref_dict.keys()) & set(input_dict.keys())

    with open(os.path.join(output, 'comparison.kabat'), 'w') as comparison_file:
        for key in common_keys:
            comparison_file.write('%s\t%s\n' % (key, '\t'.join([str(i) for i in ref_dict[key]])))
            comparison_file.write('%s\t%s\n' % (key, '\t'.join([str(i) for i in input_dict[key]])))

    reference_file.close()
    input_file.close()
    comparison_file.close()

    abs_score = sum([get_abs_score(ref_dict[key], input_dict[key]) for key in common_keys])
    mismatches_per_read = abs_score / len(common_keys) if len(common_keys) != 0 else 1
    print('Absolute score: %s; Mismatches per read: %s' % (abs_score, mismatches_per_read))


# return: key = sequence_name, value = list of regions
def parse_kabat_line(line):
    tokens = line.split()
    key = tokens.pop(0).strip()
    return key, tokens

def get_abs_score(first, second):
    first_regions = [int(first[i + 1]) - int(first[i]) + 1 for i in range(0, len(first), 2)]
    second_regions = [int(second[i + 1]) - int(second[i]) + 1 for i in range(0, len(second), 2)]
    return sum([abs(fst - snd) for fst, snd in zip_longest(first_regions, second_regions, fillvalue=0)])

if __name__ == "__main__":
   main()
